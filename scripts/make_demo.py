#!/usr/bin/env python3
"""Deterministic toy coordinates for trying the CLI; not a physical MD trajectory."""
from pathlib import Path
import random
ROOT=Path(__file__).resolve().parents[1]
rng=random.Random(2026)
target=[(1.1+i*.18,2+.12*(-1)**i,2) for i in range(11)]
water=[]
while len(water)<2000:
    p=tuple(4*rng.random() for _ in range(3))
    ds=[sum((min(abs(a-b),4-abs(a-b)))**2 for a,b in zip(p,q)) for q in target]
    if min(ds)>.23**2: water.append(p)
out=ROOT/'examples';out.mkdir(exist_ok=True)
with (out/'demo.gro').open('w') as f:
    f.write('Synthetic polymer and oxygen sites (not an MD trajectory)\n%d\n'%(len(target)+len(water)))
    for i,p in enumerate(target+water):
        is_target=i<len(target)
        res=1 if is_target else i-len(target)+2
        rn='POL' if is_target else 'SOL'; an='C' if is_target else 'OW'
        f.write(f'{res:5d}{rn:<5}{an:>5}{i+1:5d}{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}\n')
    f.write('4.00000 4.00000 4.00000\n')
with (out/'demo.ndx').open('w') as f:
    for name,indices in [('Polymer',range(1,len(target)+1)),('Water',range(len(target)+1,len(target)+len(water)+1))]:
        f.write(f'[ {name} ]\n')
        ids=list(indices)
        for k in range(0,len(ids),15): f.write(' '.join(map(str,ids[k:k+15]))+'\n')
print('Saved examples/demo.gro and demo.ndx')
