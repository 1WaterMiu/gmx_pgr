#!/usr/bin/env python3
"""Integration and independent geometry/normalization checks (requires numpy)."""
import itertools
import pathlib
import subprocess
import numpy as np
ROOT=pathlib.Path(__file__).resolve().parents[1]
OUT=ROOT/'tests'/'output'
OUT.mkdir(exist_ok=True)
GMX=str(ROOT/'bin'/'gmx')

def run(args, name, stdin=None, ok=True):
    p=subprocess.run([GMX, *args],cwd=ROOT,input=stdin,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT)
    (OUT/(name+'.log')).write_text(p.stdout)
    assert (p.returncode==0)==ok, p.stdout
    return p.stdout

def gro(path, target, water, box):
    with open(path,'w') as f:
        f.write('pgr geometry fixture\n%d\n'%(len(target)+len(water)))
        for i,x in enumerate(np.concatenate([target,water])):
            res=1 if i<len(target) else i-len(target)+2
            rn='POL' if i<len(target) else 'SOL'
            an='C' if i<len(target) else 'OW'
            f.write(f'{res%100000:5d}{rn:<5}{an:>5}{(i+1)%100000:5d}{x[0]:8.3f}{x[1]:8.3f}{x[2]:8.3f}\n')
        f.write(' '.join(f'{v:.6f}' for v in [box[0,0],box[1,1],box[2,2],box[0,1],box[0,2],box[1,0],box[1,2],box[2,0],box[2,1]])+'\n')

def data(path):
    return np.loadtxt(path,comments=('#','@'),ndmin=2)

# Independent 27-image enumeration: multiple target atoms, boundary wrapping,
# triclinic geometry, each water counted once, and exact bulk normalization.
rng=np.random.default_rng(4)
box=np.array([[4.,0,0],[1.,4.,0],[.4,.6,4.]])
target=np.array([[.05,.05,.05],[2.,2.,2.]])
water=np.round(rng.random((400,3))@box,3)
water[:2]=[[3.99,.05,.05],[2.02,2.,2.]]
fixture=OUT/'triclinic.gro'
gro(fixture,target,water,box)
args=['pgr','-f',str(fixture),'-s',str(fixture),'-sel','resname POL; resname SOL','-bin','.1','-rmax','1.8','-bulk-start','1.3','-bulk-end','1.8','-samples','100000','-o',str(OUT/'geometry.xvg')]
run(args,'geometry')
shifts=np.array(list(itertools.product([-1,0,1],repeat=3)))@box
d=np.linalg.norm(water[:,None,None,:]-target[None,:,None,:]-shifts[None,None,:,:],axis=3).min(axis=(1,2))
stats=data(OUT/'geometry_density.xvg')
expected=np.histogram(d[d<1.8],bins=np.arange(19)*.1)[0]
np.testing.assert_array_equal(stats[:,3],expected)
assert stats[:,3].sum()==(d<1.8).sum()
# Bulk boundaries align with bins here; check pooled density, not arithmetic average.
bulk=stats[13:18,3].sum()/stats[13:18,2].sum()
np.testing.assert_allclose(data(OUT/'geometry.xvg')[:,1],stats[:,1]/bulk,rtol=1e-9)
run(args,'reproducibility')
np.testing.assert_array_equal(data(OUT/'geometry_density.xvg'),stats)

# Analytic volume and independent uniform-water test for one point in a cubic box.
box=np.eye(3)*4
water=np.round(rng.random((24000,3))*4,3)
fixture=OUT/'uniform.gro'
gro(fixture,np.array([[2.,2.,2.]]),water,box)
run(['pgr','-f',str(fixture),'-s',str(fixture),'-sel','resname POL; resname SOL','-bin','.25','-rmax','1.5','-bulk-start','1','-bulk-end','1.5','-samples','1000000','-o',str(OUT/'uniform.xvg')],'uniform')
stats=data(OUT/'uniform_density.xvg')
edges=np.arange(7)*.25
analytic=4*np.pi/3*np.diff(edges**3)
np.testing.assert_allclose(stats[:,2],analytic,rtol=.07)
np.testing.assert_allclose(data(OUT/'uniform.xvg')[1:,1],1,atol=.13)
# Bin without sampled volume must be NaN, not an artificial zero.
run(['pgr','-f',str(fixture),'-s',str(fixture),'-sel','resname POL; resname SOL','-bin','.01','-rmax','1.5','-samples','1000','-o',str(OUT/'sparse.xvg')],'sparse')
a=data(OUT/'sparse_density.xvg'); g=data(OUT/'sparse.xvg')
assert np.isnan(g[a[:,2]==0,1]).all() and (a[:,2]==0).any()
# Invalid controls and unsupported PBC fail explicitly.
run(args+['-nopbc'],'reject_nopbc',ok=False)
bad=args.copy(); bad[bad.index('-bin')+1]='0'
run(bad,'reject_bin',ok=False)
bad=args.copy(); bad[bad.index('-bulk-start')+1]='1.9'
run(bad,'reject_bulk',ok=False)
if all((ROOT/'example'/name).exists() for name in ['traj_comp.xtc','topol.tpr','conf.gro','index.ndx']):
    # Real input: TPR + XTC, interactive index selection, GRO topology, TRR input.
    base=['pgr','-f','example/traj_comp.xtc','-s','example/topol.tpr','-n','example/index.ndx','-bin','.02','-samples','20000','-e','0']
    run(base+['-sel','-o',str(OUT/'interactive.xvg')],'interactive',stdin='2\n4\n')
    run(base+['-sel','group 2; group 4','-o',str(OUT/'explicit.xvg')],'explicit')
    np.testing.assert_array_equal(data(OUT/'interactive.xvg'),data(OUT/'explicit.xvg'))
    base[base.index('-s')+1]='example/conf.gro'
    run(base+['-sel','group 2; group 4','-o',str(OUT/'gro_topology.xvg')],'gro_topology')
    np.testing.assert_array_equal(data(OUT/'explicit_density.xvg'),data(OUT/'gro_topology_density.xvg'))
    run(['trjconv','-f','example/traj_comp.xtc','-s','example/topol.tpr','-dump','0','-o',str(OUT/'frame.trr')],'convert_trr',stdin='0\n')
    base[base.index('-f')+1]=str(OUT/'frame.trr')
    run(base+['-sel','group 2; group 4','-o',str(OUT/'trr.xvg')],'trr')
    np.testing.assert_array_equal(data(OUT/'trr_density.xvg'),data(OUT/'gro_topology_density.xvg'))
    print('PASS: optional research data TPR/GRO and XTC/TRR checks.')
else:
    print('SKIP: optional research input files (not distributed).')
# Additional changing-box check: sum of two spherical shell volumes is known.
# Both box sizes exceed twice rmax; a correct NPT implementation gives 2*V_shell.
small=OUT/'npt_small.gro'; large=OUT/'npt_large.gro'
gro(small,np.array([[2.,2.,2.]]),np.round(rng.random((3000,3))*4,3),np.eye(3)*4)
gro(large,np.array([[3.,3.,3.]]),np.round(rng.random((3000,3))*6,3),np.eye(3)*6)
small_text=small.read_text().replace('pgr geometry fixture','pgr geometry fixture t=0')
large_text=large.read_text().replace('pgr geometry fixture','pgr geometry fixture t=1')
multi=OUT/'npt.gro'; multi.write_text(small_text+large_text)
run(['pgr','-f',str(multi),'-s',str(small),'-sel','resname POL; resname SOL','-bin','.25','-rmax','1.5','-bulk-start','1','-bulk-end','1.5','-samples','1000000','-o',str(OUT/'npt.xvg')],'npt')
a=data(OUT/'npt_density.xvg')
np.testing.assert_allclose(a[:,2],2*analytic,rtol=.08)
bulk=a[4:,3].sum()/a[4:,2].sum()
np.testing.assert_allclose(data(OUT/'npt.xvg')[:,1],a[:,3]/a[:,2]/bulk,rtol=1e-9)
assert '# frames 2 ' in (OUT/'npt.xvg').read_text()
print('PASS: changing box volumes and pooled multi-frame density normalization.')

# Public runnable fixture: interactive index selection, XTC/TRR round trip.
base=['pgr','-f','examples/demo.gro','-s','examples/demo.gro','-n','examples/demo.ndx','-bin','.05','-rmax','1.5','-bulk-start','1','-bulk-end','1.5','-samples','20000']
run(base+['-sel','-o',str(OUT/'demo_interactive.xvg')],'demo_interactive',stdin='Polymer\nWater\n')
run(base+['-sel','group "Polymer"; group "Water"','-o',str(OUT/'demo_explicit.xvg')],'demo_explicit')
np.testing.assert_array_equal(data(OUT/'demo_interactive_density.xvg'),data(OUT/'demo_explicit_density.xvg'))
for ext in ['xtc','trr']:
    trajectory=OUT/('demo.'+ext)
    run(['trjconv','-f','examples/demo.gro','-s','examples/demo.gro','-o',str(trajectory)],'demo_convert_'+ext,stdin='0\n')
    args=base.copy();args[args.index('-f')+1]=str(trajectory)
    run(args+['-sel','group "Polymer"; group "Water"','-o',str(OUT/('demo_'+ext+'.xvg'))],'demo_read_'+ext)
    np.testing.assert_array_equal(data(OUT/('demo_'+ext+'_density.xvg'))[:,3],data(OUT/'demo_explicit_density.xvg')[:,3])
print('PASS: public GRO/XTC/TRR fixtures, interactive index selection, and all geometry/normalization checks.')
