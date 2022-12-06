import pickle
from collections import Counter
f = open("3drins_by_id.pickle", "rb")
p = pickle.load(f)

inteacts = Counter()
for k, v in p.items():
    for gk, gv, in  v["graphs"].adj.items():
        for ggk, ggv in gv.items():
            l = ggv["label"]
            inteacts[l]+=1

for k,v in inteacts.items():
    print(k,v)
