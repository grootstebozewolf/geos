#!/usr/bin/env python3
"""Re-evaluate the cases behind a findings.jsonl (e.g. against a patched GEOS build)
and report, per check, how many of the original findings still fire.

  python3 replay.py results/findings.jsonl [check ...]
"""
import collections, json, random, sys
import hunt

def main():
    path = sys.argv[1]; want = set(sys.argv[2:])
    cases = collections.defaultdict(set)
    for l in open(path):
        f = json.loads(l)
        if f['cls'] == 'ill_conditioned' or (want and f['check'] not in want):
            continue
        cases[f['check']].add(f['case'])
    allc = sorted(set().union(*cases.values())) if cases else []
    rng = random.Random(1)
    still = collections.Counter()
    for i in range(0, len(allc), 500):
        chunk = allc[i:i + 500]
        _, fs = hunt.evaluate(chunk, rng, use_oracle=False)
        for c, f in zip(chunk, fs):
            for k in set(x.check for x in f if x.cls != 'ill_conditioned'):
                if c in cases.get(k, ()):
                    still[k] += 1
    for k in sorted(cases):
        print('%-32s original %5d  still firing %5d' % (k, len(cases[k]), still[k]))

if __name__ == '__main__':
    main()
