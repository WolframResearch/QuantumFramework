import re
base=__import__("os").path.normpath(__import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), "..", "..")) + "/"
w=base+'exact/B-work/'
t=open(w+'template.md').read()
outs={}
cur=None
for line in open(w+'cells-out-all.txt'):
    m=re.match(r'=== (\d+) ',line)
    if m: cur=m.group(1); outs[cur]=''; continue
    if cur and line.strip(): outs[cur]+=line.rstrip('\n')
short={
'28':"{{{5/(3*Sqrt[2]), ((-5*I)/3)/Sqrt[2], (2*Sqrt[2])/3, ((2*I)/3)*Sqrt[2]}, {(-I)/Sqrt[2], 1/Sqrt[2], 0, 0}, <the two conjugate rows>}, {<six zero arrays>}}",
'31':"{E^(-((e*r*SuperDagger[\\[FormalA]1]**SuperDagger[\\[FormalA]2]*Tan[Sqrt[-(e*eb*r^2)]])/Sqrt[-(e*eb*r^2)]))**Sec[Sqrt[-(e*eb*r^2)]]^(SuperDagger[\\[FormalA]1]**\\[FormalA]1 + SuperDagger[\\[FormalA]2]**\\[FormalA]2)**E^(<conjugate pair term>)*Sec[Sqrt[-(e*eb*r^2)]], 0}",
'32':"{<S(xi) normal form, as in B-051>, E^(-(f*SuperDagger[\\[FormalA]1]))**E^(\\[FormalA]1*Conjugate[f])/E^((f*Conjugate[f])/2), E^((f^2*Tanh[r])/(2*E^(I*phi)))*E^(-(f*SuperDagger[\\[FormalA]1]))**E^(((-2*\\[FormalA]1*f + \\[FormalA]1**\\[FormalA]1)*Tanh[r])/(2*E^(I*phi))), E^(-1/2*(f*fb) + (f^2*Tanh[r])/(2*E^(I*phi)))*Sqrt[Sech[r]], 0}",
'35':"{E^(-1/2*(E^(I*Arg[xi])*SuperDagger[\\[FormalA]1]**SuperDagger[\\[FormalA]1]*Tanh[Abs[xi]]))**Sech[Abs[xi]]^(1/2 + SuperDagger[\\[FormalA]1]**\\[FormalA]1)**E^((\\[FormalA]1**\\[FormalA]1*Tanh[Abs[xi]])/(2*E^(I*Arg[xi]))), \\[FormalA]1*Cosh[Sqrt[x*xb]] + (x*Sinh[Sqrt[x*xb]]*SuperDagger[\\[FormalA]1])/Sqrt[x*xb], <5 x 5 array of 0>}",
'38':"{NormalOrdered[E^(-(lam*SuperDagger[\\[FormalA]1]**\\[FormalA]1))], <S(-xi) normal form>, E^(y*GeneralizedPower[NonCommutativeMultiply, SuperDagger[\\[FormalA]1], 2]*(Cosh[2*Sqrt[l^2]] + (l*Sinh[2*Sqrt[l^2]])/Sqrt[l^2])), 1/Sqrt[1 - 4*x*yy], 0}",
}
def cell(m):
    return open(w+'cells/'+m.group(1)+'.wl').read().rstrip('\n')
def out(m):
    k=m.group(1)
    return short.get(k,outs[k])
t=re.sub(r'\{\{CELL:(\d+)\}\}',cell,t)
t=re.sub(r'\{\{OUT:(\d+)\}\}',out,t)
rep={
'{{B023CHECK}}':"exact WL computation (one mode, complex A11, A22, A12 = A21 = 1 and complex sources, measure d^2f/pi: the exact `Integrate` and the formula agree; `Simplify` stalls on the nested radicals, and the difference evaluates to 0 at 99 digits)",
'{{B027INT}}':"0",
'{{B048NORM}}':"`FullSimplify` of the difference is 0",
'{{B071ONE}}':"residual 0",
'{{B071TWO}}':"residual 0 to 58 digits",
'{{B071PRINTED}}':"a constant residual of about 0.03 in both cases",
}
for k,v in rep.items(): t=t.replace(k,v)
src=open(base+'B-section.md').read()
nf=src[src.index('## Non-formula records'):]
t=t.replace('{{NONFORMULA}}',nf.rstrip('\n'))
assert not re.search(r'\{\{[A-Z0-9:]+\}\}',t), re.findall(r'\{\{[A-Z0-9:]+\}\}',t)
open(base+'exact/B-section.md','w').write(t+'\n')
