import sys
import os
import argparse
import json
from dataclasses import dataclass
from collections import defaultdict
from typing import List, Dict, Tuple, Optional, Union, DefaultDict, Mapping, Container
import functools
import warnings
import random

from Bio import PDB
from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue

import rna_bits.utils as utils
import rna_bits.utils.pdb
import rna_bits.utils.bgsu
import rna_bits.utils.ss
from rna_bits.utils.bgsu import UnitId
from rna_bits.utils.data_path import DATA_PATH
from rna_bits.utils.misc import make_warnings_nicer


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description="Convert the output of BayesPairing2 into the input for the 3d builder, including downloading the module PDB files."
    )

    parser.add_argument(
        "--dataset",
        "-d",
        required=True,
        help="Path to the bp2 dataset used, in json form, ie models/ALL.json",
    )
    parser.add_argument(
        "--bp2_result",
        "--result",
        "-r",
        required=True,
        help="Path to the bp2 result file, in json form, ie output.json",
    )
    parser.add_argument(
        "--chefs_choice",
        action="store_true",
        help="Generate the structure from bp2's 'chef's choice' and depicted in the output svg",
    )
    parser.add_argument(
        "--secondary_structures",
        "--secondary_structure",
        "--ss",
        "-ss",
        nargs="+",
        help="Secondary structures to use.",
    )
    parser.add_argument(
        "--output_file",
        "-o",
        default="./bp2b_out/out",
        help="Output .rass file. 'out' will be turned into 'out_0_0.rass when sampling multiple structures",
    )
    parser.add_argument(
        "--motif_directory",
        "-m",
        default="./bp2b_out/motifs",
        help="Directory to store module .pdb's",
    )
    parser.add_argument(
        "--use_bp2_instance_data",
        action="store_true",
        help=(
            "Currently, by default, we fetch instance data (PDB code and "
            "residue ids) from BGSU's webset. However, if this flag is set, "
            "it will take that information from the bp2 dataset file "
            "instead."
        ),
    )
    parser.add_argument(
        "--dont_fill_gaps",
        action="store_true",
        help="Output only the core of the bgsu motif",
    )
    parser.add_argument(
        "--num_outputs",
        "--samples",
        type=int,
        default=10,
        help="How many outputs to generate per provided secondary structure. "
        "10 by default",
    )
    parser.add_argument("--random_seed", type=int)

    args = parser.parse_args(argv[1:])

    return args


def get_bp2_insertion_positions(insertion_data, expand: bool) -> List[int]:
    _seq, pos, _prob, _node_pos, _score, _mapping = insertion_data
    # pos looks like: [[15, 17, 18, 19], [2, 4]]

    flat_pos = []  # Will look like [15,16,17,18,19,2,3,4] or [15,17,18,19,2,4]
    if expand:
        for l in pos:
            flat_pos.extend(range(l[0], l[-1] + 1))
    else:
        for l in pos:
            flat_pos.extend(l)

    return flat_pos


def bp2_numbers_to_expansion(a):
    """
    [0,1,3,100,102,105] -> [[True,True,False,True], [True,False,True,False,False,True]]
    """
    ret = [[]]
    ret[-1].append(True)
    prev = a[0]
    for b in a[1:]:
        pm = prev // 100
        m = b // 100
        if m > pm:
            ret.append([])
            prev = b - 1
        for _ in range(b - prev - 1):
            ret[-1].append(False)
        ret[-1].append(True)
        prev = b
    return ret


def partition_to_sizes(a, s):
    """
    partition_by([0,1,2,3,4,5], [2,3,1]) -> [[0,1],[2,3,4],[5]]
    """
    ret = []
    it = iter(a)
    for ss in s:
        rr = []
        for _ in range(ss):
            rr.append(next(it))
        ret.append(rr)

    try:  # Check that the iterator is empty
        next(it)
        assert False, "Remaining elements after partitioning by sizes"
    except StopIteration:
        pass
    return ret


class MotifInstanceUnsupportedError(Exception):
    pass


def fill_gaps(
    original: List[Residue], chain: Chain
) -> Tuple[List[Residue], List[bool]]:
    """
    Fills in the residues in between the original residues.

    Assumes that the rediues all have the same Chain object as parent.
    Also assumes that the residues are in order.

    Returns a list of the original residues with their gaps filled in, as well
    as a map indicating which of the returned residues are original.
    """
    original_dict = {a.get_id(): a for a in original}

    in_range = False
    ret = []
    expansion: List[bool] = []
    seen_cnt = 0
    for residue in chain:
        residue_id = residue.get_id()
        if residue_id in original_dict:
            in_range = True
        if in_range:
            if residue_id in original_dict:
                ret.append(original_dict[residue_id])
                expansion.append(True)
                seen_cnt += 1
            else:
                ret.append(residue)
                expansion.append(False)
            if seen_cnt == len(original):
                break
    assert expansion[0] and expansion[-1]
    return (ret, expansion)


@functools.lru_cache(maxsize=None)
def fetch_instances_units_from_bgsu(group_id: str) -> List[Tuple[str, List[UnitId]]]:
    data = utils.bgsu.download_motif_group_info(group_id)
    ret = []
    for instance_name, units_ in data["alignment"].items():
        units: list[UnitId] = [UnitId.parse(x) for x in units_]
        ret.append((instance_name, units))
    return ret


# @functools.lru_cache(maxsize=None)
def get_instances_units_from_bp2_PDBs(data) -> List[Tuple[str, List[UnitId]]]:
    """
    data is given in the "PDBs" fields of bp2 motif info and looks like:
    {'3LQX.B': [139, 140, 141, 142, 144, 167, 168, 169],
     '4V9F.0': [548, 549, 550, 551, 553, 606, 607, 608]
    }
    """
    ret = []
    for k, v in data.items():
        pdb_code, chain_id = k.split(".")
        units: List[UnitId] = []
        for residue_id in v:
            units.append(
                UnitId(
                    pdb_code=pdb_code,
                    model_id=1,
                    chain_id=chain_id,
                    residue_id=residue_id,
                )
            )
        ret.append((k, units))

    return ret


class CannotCreateMotifError(Exception):
    pass


def get_motif_including_gap_content(
    structure: Structure, units_by_component: List[List[UnitId]]
) -> Tuple[List[List[Residue]], List[List[bool]]]:
    model_id = units_by_component[0][0].model_id
    model = structure[model_id - 1]

    ret_residues = []
    ret_expansion = []
    for comp_units in units_by_component:
        comp_residues = []
        chain_id = comp_units[0].chain_id
        chain = model[chain_id]
        for unit in comp_units:
            assert unit.model_id == model_id
            assert unit.chain_id == chain_id
            try:
                residue = utils.pdb.get_residue_from_chain(chain, unit)
            except KeyError:
                warnings.warn(
                    f"Residue {str(unit)} does not exist, hence instance cannot be generated."
                )
                raise CannotCreateMotifError

            comp_residues.append(residue)
        comp_residues, comp_expansion = fill_gaps(comp_residues, chain)
        ret_residues.append(comp_residues)
        ret_expansion.append(comp_expansion)
    return (ret_residues, ret_expansion)


def get_motif_directory(bp2_id, args):
    if not args.dont_fill_gaps:
        directory = os.path.join(args.motif_directory, "with_gap_contents", bp2_id)
    else:
        directory = os.path.join(args.motif_directory, "with_no_gap_contents", bp2_id)
    return directory


@dataclass(frozen=True)
class ReorderedDatasetInfo:
    bp2_expansions: Dict[str, List[List[bool]]]
    # instance_id_mapping[pdb_code][unit_id_tuple] = instance_id
    instance_id_mapping: Mapping[str, Mapping[Tuple[UnitId, ...], str]]
    # bp2_ids_mapping[pdb_code][unit_id_tuple] = list_of_bp2_motif_ids
    bp2_ids_mapping: Mapping[str, Mapping[Tuple[UnitId, ...], List[str]]]


def make_reordered_dataset_info(
    dataset, considered_motifs, args
) -> ReorderedDatasetInfo:
    """
    In order extract and process motifs from the PDBs in an efficient way, it
    should be processes in a "PDB-first order".
    I can't fit all the PDBs in memory at the same time, and I don't want to
    load the same PDB multiple times. Therefore, I will load each PDB once
    and process all the motifs that were extracted from it.
    In order to do so, I "reorder" the information given in `dataset` to be
    "PDB-first order".

    """
    bp2_expansions: Dict[str, List[List[bool]]] = {}

    # instance_id_mapping[pdb_code][unit_id_tuple] = instance_id
    # (eg "IL_5TBW_102" or "5TBW.A")
    instance_id_mapping: DefaultDict[str, Dict[Tuple[UnitId, ...], str]] = defaultdict(
        dict
    )
    # bp2_ids_mapping[pdb_code][unit_id_tuple] = list_of_bp2_motif_ids
    # (eg ["12","13"])
    bp2_ids_mapping: DefaultDict[
        str, DefaultDict[Tuple[UnitId, ...], List[str]]
    ] = defaultdict(lambda: defaultdict(list))

    for bp2_motif_id, motif_data in dataset.items():
        if considered_motifs is not None:
            if bp2_motif_id not in considered_motifs:
                continue

        # TODO: Check within bp2 source code to see if using the aln numbers are
        # indeed a correct way to determine the extended shape of the motif.
        bp2_numbers = sorted(map(int, motif_data["aln"].keys()))  # [0,1,3,100,102,105]
        expansion = bp2_numbers_to_expansion(bp2_numbers)  # [[T,T,F,T], [T,F,T,F,F,T]]
        bp2_expansions[bp2_motif_id] = expansion

        if args.use_bp2_instance_data:
            instances = get_instances_units_from_bp2_PDBs(motif_data["PDBs"])
        else:
            group_id = motif_data["atlas_name"]
            instances = fetch_instances_units_from_bgsu(group_id)

        for instance_id, units in instances:
            has_symmetry = any(x.symmetry is not None for x in units)
            if has_symmetry:
                # We don't support symmetry
                continue
            pdb_code = units[0].pdb_code
            model_id = units[0].model_id
            for unit in units:
                assert unit.pdb_code == pdb_code
                assert unit.model_id == model_id
            units_tuple = tuple(units)
            instance_id_mapping[pdb_code][units_tuple] = instance_id
            bp2_ids_mapping[pdb_code][units_tuple].append(bp2_motif_id)

    return ReorderedDatasetInfo(
        bp2_expansions=bp2_expansions,
        instance_id_mapping=instance_id_mapping,
        bp2_ids_mapping=bp2_ids_mapping,
    )


def generate_models(dataset, considered_motifs, args):
    # Check which motif models have already been generated
    models_by_bp2_id = defaultdict(list)
    for bp2_id in considered_motifs:
        directory = get_motif_directory(bp2_id, args)
        if os.path.isfile(os.path.join(directory, "done")):
            _ = models_by_bp2_id[
                bp2_id
            ]  # Force defaultdict to create the entry if it doesn't exist
            for f in os.listdir(directory):
                if f.endswith(".pdb"):
                    model_path = os.path.join(directory, f)
                    # The intention of the None is that perhaps it'll be
                    # replaced by an actual Model object.
                    models_by_bp2_id[bp2_id].append((model_path, None))

    needs_building = set(considered_motifs) - set(models_by_bp2_id)

    # Go through the PDBs and generate the motifs
    a = make_reordered_dataset_info(dataset, needs_building, args)
    tot_stuff = sum(len(v) for v in a.instance_id_mapping.values())
    did = 0
    for pdb_code, units_bp2_ids in a.bp2_ids_mapping.items():
        structure = utils.pdb.fetch_pdb(pdb_code)
        for units, bp2_ids in units_bp2_ids.items():
            instance_id = a.instance_id_mapping[pdb_code][units]
            atlas_names = [dataset[id]["atlas_name"] for id in bp2_ids]
            atlas_name = atlas_names[0]
            some_expansion = a.bp2_expansions[bp2_ids[0]]
            component_sizes = [sum(comp) for comp in some_expansion]

            if sum(component_sizes) != len(units):
                warnings.warn(
                    f"The number of units provided for {atlas_name} {instance_id} "
                    f"does not match the number of units in the dataset "
                    f"json file {repr(bp2_ids[0])} ({bp2_ids})"
                )
                continue

            units_by_component = partition_to_sizes(units, component_sizes)
            try:
                residues, pdb_expansion = get_motif_including_gap_content(
                    structure, units_by_component
                )
            except CannotCreateMotifError:
                warnings.warn(
                    f"Instance {instance_id} of motifs {bp2_ids} could not be generated."
                )
                continue

            num_hits = 0
            model = utils.pdb.build_model_from_lists_of_residues(residues)
            for bp2_id in bp2_ids:
                bp2_expansion = a.bp2_expansions[bp2_id]
                if bp2_expansion == pdb_expansion:
                    num_hits += 1
                    directory = get_motif_directory(bp2_id, args)

                    os.makedirs(directory, exist_ok=True)
                    save_path = os.path.join(directory, instance_id + ".pdb")
                    utils.pdb.save_model_as_pdb(model, save_path)
                    # models_by_bp2_id[bp2_id].append((save_path, model))
                    models_by_bp2_id[bp2_id].append((save_path, None))
            did += 1
            print("Generated", did, "/", tot_stuff, "motif instances")

    # Add a "done" file to indicate which motifs have been processed
    for bp2_id in considered_motifs:
        directory = get_motif_directory(bp2_id, args)
        os.makedirs(directory, exist_ok=True)
        with open(os.path.join(directory, "done"), "w"):
            # touch
            pass

    return models_by_bp2_id


def make_correct_model_path(path, args):
    """
    Generates the model path as it will be included in the rass file.
    """
    return os.path.join(
        ".", os.path.relpath(path, start=os.path.dirname(args.output_file))
    )


def ts(a):
    return tuple(sorted(a))


@dataclass(frozen=True)
class MotifInsertion:
    model: Union[str, PDB.Model.Model]
    positions: List[int]


@dataclass
class Assembly:
    seq: str
    ss: str
    motifs: List[MotifInsertion]


def process_bp2(result, dataset, args) -> Mapping[str, List[Assembly]]:
    if args.chefs_choice:
        hits = result["svg_hits"]["input_seq"]
    else:
        hits = result["all_hits"]["input_seq"]

    # Get the models
    models_by_bp2_id = generate_models(dataset, hits, args)

    # Make the model paths be relative to the output rass file
    for l in models_by_bp2_id.values():
        for k, (path, model) in enumerate(l):
            # path, model = l[k]
            path = make_correct_model_path(path, args)
            l[k] = (path, model)

    # Get the sequence
    _seqs = result["input"]
    assert len(_seqs) <= 1, "Multiple sequences not supported."
    seq = _seqs[0]

    # Get the secondary structures
    if args.secondary_structures is not None:
        sss = args.secondary_structures
    elif "chefs_choice_struct" in result:
        sss = result["chefs_choice_struct"]
    else:
        assert (
            False
        ), "No secondary structure given and no chefs_choice_struct in the bp2 output."

    # For the hits given by bp2, determine the areas where they are inserted
    competing = defaultdict(list)
    for motif_id, insertion_datas in hits.items():
        if len(insertion_datas) == 0:
            continue

        models: List[Tuple[str, Model]] = models_by_bp2_id[motif_id]
        if len(models) == 0:
            warnings.warn(f"There's no models for motif {motif_id}")
            continue

        for ins_d in insertion_datas:
            ins_pos = get_bp2_insertion_positions(
                ins_d, expand=(not args.dont_fill_gaps)
            )
            exp_ins_pos = get_bp2_insertion_positions(ins_d, expand=True)
            ts_eip = ts(exp_ins_pos)

            for model_filename, _ in models:
                competing[ts_eip].append((model_filename, ins_pos))

    # Sort to have deterministic outputs
    # The exact way it's sorted isn't important
    for l in competing.values():
        l.sort()

    # For each of the provided secondary structures, insert the motifs that fit into it
    assemblies = {}
    for ss in sss:
        loop_areas = set()

        loop_infos = utils.ss.Segmenter.from_parens(sss[0], start=0).segment_loops()
        for loop_ss, loop_pos in loop_infos:
            loop_areas.add(ts(loop_pos))

        num_outputs = args.num_outputs

        ss_assemblies = []
        for _ in range(num_outputs):
            ass = Assembly(seq=seq, ss=ss, motifs=[])
            ss_assemblies.append(ass)

        for la in loop_areas:
            possible_models = competing[la]
            if len(possible_models) == 0:
                warnings.warn(f"Loop {la} does not have any models")
                continue
            for ass in ss_assemblies:
                model, positions = random.choice(possible_models)
                ass.motifs.append(MotifInsertion(model, positions))

        assemblies[ss] = ss_assemblies

    return assemblies


def generate_filename(ss_inc, sample_inc, args):
    return args.output_file + "_" + str(ss_inc) + "_" + str(sample_inc) + ".rass"


def write_assembly(assembly, f):
    f.write(assembly.seq)
    f.write("\n")
    f.write(assembly.ss)
    f.write("\n")

    for motif in assembly.motifs:
        shifted = [a + 1 for a in motif.positions]
        assert type(motif.model) is str
        f.write(f"motif:{motif.model}: {','.join(map(str,shifted))}\n")


def main(argv):
    make_warnings_nicer()
    args = parse_args(argv)

    if args.random_seed is not None:
        random.seed(args.random_seed)

    with open(args.bp2_result) as f:
        result = json.load(f)
    with open(args.dataset) as f:
        dataset = json.load(f)

    ss_inc = 0

    def output_ss_assemblies(ss_assemblies):
        nonlocal ss_inc
        for i, ass in enumerate(ss_assemblies):
            fn = generate_filename(ss_inc, i, args)
            with open(fn, "w") as f:
                write_assembly(ass, f)
        ss_inc += 1

    all_assemblies = process_bp2(result, dataset, args)
    if args.secondary_structures is not None:
        for ss in args.secondary_structures:
            output_ss_assemblies(all_assemblies[ss])
    else:
        assert len(all_assemblies) == 1
        # Here I assert so that if it ever happens that more than one ss is
        # used in this case, I remember to change the code to output them in an
        # order that isn't arbitrary.
        for ss_assemblies in all_assemblies.values():
            output_ss_assemblies(ss_assemblies)


def main_wrapper():
    main(sys.argv)


if __name__ == "__main__":
    main_wrapper()
