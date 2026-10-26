from pathlib import Path
import sys,json,hashlib,importlib.util,time
from fractions import Fraction
import numpy as np
import argparse
parser=argparse.ArgumentParser(description='Exhaustive compensated-curvature integer model; no circuits.')
parser.add_argument('--data',type=Path,required=True)
args=parser.parse_args()
out=Path(__file__).resolve().parent;p=out.parent
(out/'local').mkdir(exist_ok=True)
sys.path[:0]=[str(p),str(args.data/'verifiers')]
import optimized_quadratic as opt
from invsqrt import verify_submission_invsqrt as v
raw=opt.load_coefficients();h=Fraction(1,64)
orig=[[Fraction(abs(z),2**36) for z in r] for r in raw]
xs=np.arange(1,3133453,dtype=np.uint64)
target,tolerance=v.reference_values(v.ProblemSpecInvSqrt.for_n_p(10),xs)
records=[];tables={}
for fc,compensate in [(8,False),(8,True),(9,True),(10,True),(11,True),(12,True)]:
 rows=[]
 for a,b,c in orig:
  cq=Fraction(round(c*2**fc),2**fc);d=cq-c
  rows.append([round((a+(d*h*h/8 if compensate else 0))*2**30),round((b+(d*h if compensate else 0))*2**23),round(cq*2**fc)])
 rows=np.array(rows,dtype=np.int64);widths=(30,23,fc);fail=0;worst=-1.;worstx=0;capacity=0;overflow=0
 physical=[int(max(rows[:,j])).bit_length() for j in range(3)]
 for start in range(1,3133453,65536):
  x=np.arange(start,min(start+65536,3133453),dtype=np.int64)
  shifts=(22-np.searchsorted(1<<np.arange(23),x,side='right'))//2
  normalized=x<<(2*shifts);address=(normalized>>14)-64
  values=[rows[address,j].copy() for j in range(3)]
  for j in (1,0):
   for m in range(7,21):
    sh=widths[j]-widths[j+1]-m
    term=values[j+1]<<sh if sh>=0 else values[j+1]>>(-sh)
    values[j]-=((normalized>>(20-m))&1)*term
   capacity+=int(np.count_nonzero(values[j]<=0))
   overflow+=int(np.count_nonzero(values[j] >= 2**physical[j]))
  phase=(values[0]<<shifts).astype(np.uint64)<<np.uint64(24)
  sl=slice(start-1,start-1+len(x));used=np.abs((phase-target[sl]).view(np.int64)).astype(np.uint64)
  ratio=used/tolerance[sl];k=int(ratio.argmax())
  if ratio[k]>worst:worst=float(ratio[k]);worstx=int(x[k])
  fail+=int(np.count_nonzero(used>tolerance[sl]))
 name=f'curvature-f{fc}-'+('compensated' if compensate else 'naive')
 records.append(dict(name=name,fractional_bits=list(widths),physical_widths=[int(max(rows[:,j])).bit_length() for j in range(3)],failed_inputs=fail,nonpositive_accumulators=capacity,overflow_accumulators=overflow,worst_error_over_tolerance=worst,worst_x=worstx,status='MATHEMATICAL_MODEL_PASS' if not(fail or capacity or overflow) else 'REJECT',circuit_built=False))
 tables[name]=rows.tolist()
 print(records[-1],flush=True)
np.savez(out/'local/curvature-candidates.npz',**{k:np.array(z,dtype=np.int64) for k,z in tables.items()})
report=dict(scope='EXHAUSTIVE_INTEGER_MODEL_ONLY_NO_CIRCUIT',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),inputs=3133452,coefficient_source_sha256=hashlib.sha256((p/'feasibility/B-quadratic-192.json').read_bytes()).hexdigest(),official_verifier_sha256=hashlib.sha256(Path(v.__file__).read_bytes()).hexdigest(),h='1/64',candidates=records,raw_tables_ignored=True)
(out/'curvature-screen.json').write_text(json.dumps(report,indent=2)+'\n')
