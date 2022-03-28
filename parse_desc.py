import os
import random
from collections import Counter

# class link:
#     def __init__(self, ):
#         self.letter_a
#         self.letter_b
#         self.

MAX_COMP_GAP = 4;

# "1004_A" -> (1004, "A")
def parse_nuc(s):
    n,l = s.split("_")
    return (int(n), l)

all_files = set()
files_with_broken_entrances = set()
all_types = Counter()

CANONICAL = {
    "C +/+ c G",
    "G +/+ c C",
    "U -/- c A",
    "G +/+ c U",
    "U +/+ c G",
    "A -/- c U"
}


def parse_desc(path):
    all_files.add(path)
    file = open(path)
    file.readline()
    a = file.readline()
    b,c = a.split(":")
    assert(b.strip() == "Bases")
    c = [parse_nuc(s) for s in c.strip().split()]

    # get first and last nucleotide position of each component
    chain_ends = []
    current_start = c[0][0]
    prev = c[0][0]
    for (n, l) in c[1:]:
        if (n-prev > MAX_COMP_GAP+1):
            chain_ends.append((current_start, prev))
            current_start = n
        prev = n
    chain_ends.append((current_start, prev))

    entrances = []
    for i in range(len(chain_ends)):
        ii = (i-1)%len(chain_ends)
        entrances.append((chain_ends[i][0], chain_ends[ii][1]))

    links = {}
    for l in file:
        l = l.strip()
        a,m,b = (a.strip() for a in l.split("---"))
        m = m.strip()
        a = parse_nuc(a[1:-1].strip())
        b = parse_nuc(b[1:-1].strip())
        links[(a[0],b[0])] = a[1]+" "+m+" "+b[1]

    #print(links)
    for a,b in entrances:
        if not(a<b):
            a,b=b,a
        if (a,b) in links:
            #print((a,b),links[(a,b)])
            typ = links[(a,b)]
        else:
            typ = "Empty"
        all_types[typ]+=1
        if (typ not in CANONICAL):
            files_with_broken_entrances.add(fileName);




folder = "RNAMoIP/No_Redondance_DESC/";
for motifDesc in (x for x in os.listdir(folder) if x[-5:] == '.desc'):
    fileName = os.path.join(folder, motifDesc)
    parse_desc(fileName)
for n in files_with_broken_entrances:
    print(n)
tot = sum(all_types.values())
tot_can = 0
for (n,i) in all_types.items():
    print(n,i,i/tot)
    if(n in CANONICAL):
        tot_can+=i
print(len(files_with_broken_entrances),'/',len(all_files), len(files_with_broken_entrances)/len(all_files))
print(tot-tot_can, "/", tot, (tot-tot_can)/tot)
