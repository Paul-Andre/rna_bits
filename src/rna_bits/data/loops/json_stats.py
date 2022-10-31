import json
import os
from collections import Counter
from collections import defaultdict

JSON_DIR = "json/"

cnt = defaultdict(list)
fns = [a for a in os.listdir(JSON_DIR) if a.endswith(".json") and "_" not in a]
for i, fn in enumerate(fns):
    print(fn, str(i + 1) + "/" + str(len(fns)))
    dfn = os.path.join(JSON_DIR, fn)
    with open(dfn) as f:
        info = json.load(f)
    for entry in info:
        seqs = tuple(entry["seqs"])
        lfn = entry["file"]
        cnt[seqs].append(lfn)

for n, k, v in sorted((len(v), k, v) for (k, v) in cnt.items()):
    print(k, n, v)
