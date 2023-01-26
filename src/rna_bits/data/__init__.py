import sys, os
import subprocess

_SCRIPT_DIR = os.path.dirname(__file__)


def execute_script(d):
    subprocess.call([sys.executable, d])


def generate_loops():
    old_path = os.getcwd()
    os.chdir(_SCRIPT_DIR + "/loops")

    execute_script("00_download_representative.py")
    execute_script("01_normalize.py")
    execute_script("02_mc_annotate.py")
    execute_script("03_ss_annotate.py")
    execute_script("04_segment.py")
    execute_script("05_extract.py")

    os.chdir(old_path)


def generate_stacks():
    old_path = os.getcwd()
    os.chdir(_SCRIPT_DIR + "/stacks")

    execute_script("00_download.py")
    execute_script("01_mc_annotate.py")
    execute_script("02_ss.py")
    execute_script("03_segment.py")
    execute_script("04_extract.py")
    execute_script("05_copy.py")

    os.chdir(old_path)
