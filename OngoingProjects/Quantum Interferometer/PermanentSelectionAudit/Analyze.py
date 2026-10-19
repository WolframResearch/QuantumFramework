"""Analyze the fixed proposal; do not tune coefficients on validation data."""
import csv, json, math
from pathlib import Path
ROOT = Path(__file__).resolve().parent

def proposed(row):
    n, m, k, r = (row[x] for x in ('Photons','Modes','Terms','Rows'))
    factor = 16 if row['Precision']=='Machine' else .1
    return 'Repeated' if k < 2**n and factor*(m*m+k*m*(r+1)) < n*2**n else 'Direct'

def evaluate(rows):
    complete=[r for r in rows if r['DirectStatus']=='OK' and r['RepeatedStatus']=='OK']
    old=lambda r:'Repeated' if r['Terms'] < 2**r['Photons'] else 'Direct'
    metric={}
    for name, policy in [('Old',old),('Proposed',proposed)]:
        selected=[r for r in rows if policy(r)=='Repeated']
        failures=[r for r in rows if r[policy(r)+'Status']!='OK']
        regret=[r[policy(r)+'Seconds']/min(r['DirectSeconds'],r['RepeatedSeconds']) for r in complete]
        wrong=[r for r in complete if r[policy(r)+'Seconds']>1.1*min(r['DirectSeconds'],r['RepeatedSeconds'])]
        total=sum(r[policy(r)+'Seconds'] for r in complete)
        metric[name]={'RepeatedChoices':len(selected),'SelectedCensoredRoutes':len(failures),
          'CompleteCaseSumSeconds':total,'WorstCompleteCaseRegret':max(regret,default=None),
          'CasesOver10PercentSlowerThanFastest':len(wrong)}
    metric['Rows']=len(rows);metric['BothRoutesCompleted']=len(complete)
    metric['MaxAbsoluteDifference']=max((r['AbsoluteError'] or 0 for r in rows),default=0)
    metric['MaxRelativeDifference']=max((r['RelativeError'] or 0 for r in rows),default=0)
    metric['ProposedMisclassifications']=[{k:r[k] for k in ['Seed','Photons','Modes','Precision','Layout','DirectSeconds','RepeatedSeconds']} for r in complete if r[proposed(r)+'Seconds']>1.1*min(r['DirectSeconds'],r['RepeatedSeconds'])]
    return metric

rows=json.loads((ROOT/'benchmark.json').read_text())
for row in rows:
    row['ProposedRoute']=proposed(row)
    row['OldRoute']='Repeated' if row['Terms']<2**row['Photons'] else 'Direct'
columns=['Seed','Split','Layout','Photons','Modes','Precision','Terms','Rows','InputOccupied','OutputOccupied','Direction','OldRoute','ProposedRoute','DirectStatus','RepeatedStatus','DirectSeconds','RepeatedSeconds','AbsoluteError','RelativeError']
with (ROOT/'benchmark.csv').open('w',newline='') as stream:
    writer=csv.DictWriter(stream,fieldnames=columns,extrasaction='ignore');writer.writeheader();writer.writerows(rows)
report={'Coefficients':{'Machine':16,'OtherApproximate':.1},
 'Formula':'Repeated iff K < 2^n and factor*(m^2+K*m*(r+1)) < n*2^n; exact inputs remain Direct.',
 'Training':evaluate([r for r in rows if r['Split']=='Training']),
 'IndependentSeedMachineValidation':evaluate([r for r in rows if r['Split']=='Validation']),
 'All':evaluate(rows),
 'Interpretation':'Times are normalized-amplitude calls including repeated-row/column matrix extraction and normalization and original whole-W precision conversion. Censored rows are not assigned fabricated runtimes. Summed times are only over cases where both routes completed; they are descriptive, not representative workload weights. Constants are empirical for WL15 ARM, not portable guarantees.'}
(ROOT/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')
print(json.dumps(report,indent=2))
