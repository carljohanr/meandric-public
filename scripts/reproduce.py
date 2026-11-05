#!/usr/bin/env python3
"""Small from-scratch generation and independent checking. No downloaded certificates."""
import argparse
from fractions import Fraction
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
R = ROOT / 'research'
OUT = ROOT / 'outputs'
for name in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','VECLIB_MAXIMUM_THREADS'):
    os.environ[name] = '1'

def run(args, cwd=ROOT):
    print('+ ' + ' '.join(map(str,args)), flush=True)
    subprocess.run(list(map(str,args)),cwd=cwd,check=True)

def compiler(cpp=False):
    name = os.environ.get('CXX' if cpp else 'CC', 'c++' if cpp else 'cc')
    found = shutil.which(name)
    if not found:raise SystemExit(f'{name} not found; install a C/C++ compiler.')
    return found

def build_upper_b():
    folder=R/'1005_head_edge_target'
    dest=folder/('evaluator.dylib' if sys.platform=='darwin' else 'evaluator.so')
    flags=['-dynamiclib'] if sys.platform=='darwin' else ['-shared','-fPIC']
    run([compiler(True),'-std=c++17','-O2',*flags,folder/'evaluator.cpp','-o',dest])

def b_upper():
    import numpy as np
    build_upper_b()
    source=R/'1201_upper_structures'
    sys.path.insert(0,str(source))
    from core3 import Model3
    m=Model3(4,3,4,4,3)
    start=OUT/'b_upper_start.npy'
    np.save(start,np.zeros(m.nv))
    prefix=OUT/'b_upper'
    run([sys.executable,source/'run_big.py','--a',4,'--cfg','(3,4,4,3)',
         '--start',start,'--betas','100,1000','--iters',100,'--seconds',30,'--out',prefix])
    report=OUT/'b_upper_exact.json'
    run([sys.executable,R/'1251_referee_upper_12297/run_formal.py',
         str(prefix)+'_cand.json',report])
    d=json.loads(report.read_text())
    assert d['rows']==512 and d['violations_vs_claimed']==0
    upper=Fraction(d['Rmax'])**2
    print('CERTIFIED upper bound on mu:',float(upper),'(exact rational in outputs/b_upper_bound.json)')
    (OUT/'b_upper_bound.json').write_text(json.dumps({'mu_upper':str(upper)},indent=2)+'\n')

def b_lower():
    cert=OUT/'b_lower_r3.pkl'
    run([sys.executable,R/'805_entropy_certificates/exact_cert.py',3,23,32,7,16,cert])
    run([sys.executable,R/'805_entropy_certificates/verify_cert.py',cert])

def e_upper():
    solver=OUT/'bic';checker=OUT/'check_e_upper'
    run([compiler(True),'-std=c++17','-O2','-pthread',R/'602_bicubic_upper_bound/bic.cpp','-o',solver])
    run([compiler(),'-O2','-pthread',R/'1390_cross_model_audit/upperE/hp_upperE.c','-o',checker])
    table=OUT/'e_upper_table.bin';cert=OUT/'e_upper_r3.bin'
    run([solver,'solve',3,'0.8915','0.685','0.776',f'out={table}','tol=1e-8','maxit=2000','maxsec=30','threads=1'])
    run([solver,'cert',3,1783,2000,137,200,97,125,1000000000,f'in={table}',f'out={cert}','threads=1'])
    run([checker,cert,1])

def e_lower():
    cert=OUT/'e_lower_r4.bin'
    run([sys.executable,R/'1367_bicubic_note/make_cert.py',4,41,50,400,cert])
    # A modest test target: the checker must establish it, not trust a floating-point score.
    run([sys.executable,R/'1390_cross_model_audit/lowerE/hp_lowerE.py',cert,'10'])

def counts():
    # Keep all generated data in outputs rather than changing source files.
    code=(R/'review_2026_11_04/spot_checks.py').read_text()
    code=code.replace("Path(__file__).with_name('spot_checks.json')", "Path('outputs/counts.json')")
    exec(compile(code,'counting-example','exec'),{'__file__':str(R/'review_2026_11_04/spot_checks.py')})

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('example',choices=['counts','b-upper','b-lower','e-upper','e-lower','all'])
    args=ap.parse_args();OUT.mkdir(exist_ok=True)
    os.chdir(ROOT)
    jobs={'counts':counts,'b-upper':b_upper,'b-lower':b_lower,'e-upper':e_upper,'e-lower':e_lower}
    for name,fn in jobs.items():
        if args.example in (name,'all'):
            print('\n=== '+name+' ===',flush=True);fn()

if __name__=='__main__':main()
