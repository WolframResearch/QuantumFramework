from pathlib import Path
import sys,numpy as np,json,hashlib,argparse
parser=argparse.ArgumentParser(description='Exact interval feasibility for rounding-aware affine coefficients; no circuits.')
parser.add_argument('--data',type=Path,required=True)
args=parser.parse_args()
out=Path(__file__).resolve().parent;p=out.parent
sys.path.insert(0,str(args.data/'verifiers'))
from invsqrt import verify_submission_invsqrt as v
x=np.arange(1,3133453,dtype=np.int64);shift=(22-np.searchsorted(1<<np.arange(23),x,side='right'))//2
n=x<<(2*shift);addr=(n>>14)-64
ref,tol=v.reference_values(v.ProblemSpecInvSqrt.for_n_p(10),x.astype(np.uint64))
k=24+shift;unit=1<<k;q=(ref>>k.astype(np.uint64)).astype(np.int64);q[0]=1<<(64-k[0]);rem=(ref&(unit.astype(np.uint64)-1)).astype(np.int64);tol=tol.astype(np.int64)
lo=q-((tol-rem)//unit);hi=q+(tol+rem)//unit
rows=np.load(p/'representation_search/local/curvature-candidates.npz')['curvature-f9-compensated'].copy()
def offset(ix,b,c):
 factor=np.full(len(ix),b,dtype=np.int64)
 for m in range(7,21):
  factor-=((n[ix]>>(20-m))&1)*((c<<(14-m)) if 14-m>=0 else (c>>(m-14)))
 off=np.zeros(len(ix),dtype=np.int64)
 for m in range(7,21):
  off+=((n[ix]>>(20-m))&1)*(factor>>(m-7))
 return off
fixed=[];failed=[]
for a in range(192):
 ix=np.flatnonzero(addr==a);ra,rb,rc=map(int,rows[a]);chosen=None
 for db in [0]+[z for j in range(1,33) for z in (j,-j)]:
  off=offset(ix,rb+db,rc);l=int((lo[ix]+off).max());h=int((hi[ix]+off).min())
  if l<=h:
   chosen=(min(max(ra,l),h),rb+db,rc);break
 if chosen is None: failed.append(a)
 else:
  if chosen!=(ra,rb,rc):fixed.append((a,chosen[0]-ra,chosen[1]-rb))
  rows[a]=chosen
print('modified',fixed,'unresolved',failed,flush=True)
np.save(out/'local/curvature-rounding-rows.npy',rows)
report=dict(scope='FINITE_AFFINE_ROUNDING_SEARCH_NO_CIRCUIT',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),candidate_source_sha256=hashlib.sha256((out/'local/curvature-candidates.npz').read_bytes()).hexdigest(),inputs=len(x),fractional_bits=[30,23,9],B_offsets=list(range(-32,33)),modified_rows=fixed,unresolved_intervals=failed,status='REJECT' if failed else 'INTEGER_MODEL_FEASIBLE',note='Empty A-interval at all searched B offsets rejects only this finite search, not every possible B or another C.')
(out/'affine-rounding.json').write_text(json.dumps(report,indent=2)+'\n')
