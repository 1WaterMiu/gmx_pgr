#!/usr/bin/env python3
"""Render the measured 596-frame example, without smoothing or interpolation."""
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
g=np.loadtxt(ROOT/'docs/data/example_pgr.xvg',comments=('#','@'))
d=np.loadtxt(ROOT/'docs/data/example_pgr_density.xvg',comments=('#','@'))
rho_bulk=33.9642034
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
fig,axes=plt.subplots(1,2,figsize=(11.4,4.7),layout='constrained')
for ax in axes:
    ax.axvspan(1.6,2,color='#efc46d',alpha=.24,label='Far-field normalization interval')
    ax.set(xlim=(0,2),xlabel='Minimum distance to any polymer atom (nm)')
    ax.grid(alpha=.18)
axes[0].plot(g[:,0],g[:,1],color='#186e8a',lw=1.8,label='Volume-corrected pRDF')
axes[0].axhline(1,color='#585858',lw=1,ls='--')
axes[0].set(ylabel=r'$p_g(r)$',title='A  Proximal radial distribution function',ylim=(0,None))
axes[0].legend(loc='lower right',fontsize=8.5,frameon=False)
axes[1].plot(d[:,0],d[:,1],color='#ad5635',lw=1.8,label=r'Water oxygen density $\rho(r)$')
axes[1].axhline(rho_bulk,color='#585858',lw=1,ls='--',label=r'$\rho_{\rm far}=33.9642$ O nm$^{-3}$')
axes[1].set(ylabel=r'$\rho(r)$ (O nm$^{-3}$)',title='B  Density before far-field normalization',ylim=(0,None))
axes[1].legend(loc='lower right',fontsize=8.5,frameon=False)
fig.suptitle('Polymer–water example | 596 frames · 0.01 nm bins · 200,000 volume samples/frame',fontsize=12)
for ext in ['png','svg']:
    fig.savefig(ROOT/f'docs/images/example-pgr.{ext}',dpi=200)
svg=ROOT/'docs/images/example-pgr.svg'
svg.write_text('\n'.join(line.rstrip() for line in svg.read_text().splitlines())+'\n')
plt.close(fig)
print('Saved docs/images/example-pgr.png and .svg')
