# Inactive | [SpanFromLeft]

> [Inactive](https://reference.wolfram.com/language/ref/Inactive.html)[*f*]  — is an inactive form of `*f*`.

## Details

[Inactive](https://reference.wolfram.com/language/ref/Inactive.html)[*f*][*args*] is effectively a purely symbolic form of `*f*[*args*]`, in which no evaluation associated with `*f*` is done.

[Inactive](https://reference.wolfram.com/language/ref/Inactive.html) is conveniently inserted into expressions using [Inactivate](https://reference.wolfram.com/language/ref/Inactivate.html).

[Inactive](https://reference.wolfram.com/language/ref/Inactive.html)[*f*] displays in [StandardForm](https://reference.wolfram.com/language/ref/StandardForm.html), with `*f*` or any special output form associated with `*f*` shown in gray.

[Inactive](https://reference.wolfram.com/language/ref/Inactive.html) does not affect [](https://reference.wolfram.com/language/ref/.html).

[Inactive](https://reference.wolfram.com/language/ref/Inactive.html) has the attribute [HoldFirst](https://reference.wolfram.com/language/ref/HoldFirst.html) and does not evaluate its first argument.

[Inactive](https://reference.wolfram.com/language/ref/Inactive.html)[*atom*] gives `*atom*` for atoms other than symbols.

## Examples

### Basic Examples

Inactive [Length](https://reference.wolfram.com/language/ref/Length.html):

```wolfram
Inactive[Length][{a,b,c}]
(* Output *)
Length[{a,b,c}]
```

Evaluate the expression:

```wolfram
Activate[%]
(* Output *)
3
```

Inactivate [Plus](https://reference.wolfram.com/language/ref/Plus.html):

```wolfram
Inactivate[2+2]
(* Output *)
2+2
```

Display equality of activated and inactivated forms:

```wolfram
%==Activate[%]
(* Output *)
2+2==4
```

Inactive objects are grayed out in [StandardForm](https://reference.wolfram.com/language/ref/StandardForm.html):

```wolfram
Inactive[Sin][0]
(* Output *)
Sin[0]
```

But not in [](https://reference.wolfram.com/language/ref/.html):

```wolfram
[%]
(* Output *)
sin(0)
```

### Scope

#### Basic Uses

Define an inactive expression:

```wolfram
expr=Inactive[Sin][π/2]
(* Output *)
Sin[(π)/(2)]
```

Evaluate the expression using [Activate](https://reference.wolfram.com/language/ref/Activate.html):

```wolfram
Activate[expr]
(* Output *)
1
```

Create an inactive expression using [Inactivate](https://reference.wolfram.com/language/ref/Inactivate.html):

```wolfram
expr=Inactivate[Sin[Pi/2],Sin]
(* Output *)
Sin[(π)/(2)]
```

Evaluate the expression:

```wolfram
Activate[expr]
(* Output *)
1
```

Expression with two inactive terms:

```wolfram
expr =Inactive[Cos][(π)/(3)]+Inactive[Sin][(π)/(3)]
(* Output *)
Cos[(π)/(3)]+Sin[(π)/(3)]
```

Activate different parts of the expression:

```wolfram
Activate[expr,Sin]
(* Output *)
(Sqrt[3])/(2)+Cos[(π)/(3)]
```

```wolfram
Activate[expr,Cos]
(* Output *)
(1)/(2)+Sin[(π)/(3)]
```

```wolfram
Activate[expr]
(* Output *)
(1)/(2)+(Sqrt[3])/(2)
```

Inactivate the symbol `g` only:

```wolfram
Inactivate[f[g[h[x]],y,g[z]],g]
(* Output *)
f[g[h[x]],y,g[z]]
```

Activate `g`:

```wolfram
Activate[%]
(* Output *)
f[g[h[x]],y,g[z]]
```

Inactivate `g` and `h`:

```wolfram
Inactivate[f[g[h[x]],y,g[z]],g|h]
(* Output *)
f[g[h[x]],y,g[z]]
```

Activate `h`:

```wolfram
Activate[%,h]
(* Output *)
f[g[h[x]],y,g[z]]
```

#### Formal Operations

An inactive form of an integral:

```wolfram
Inactive[Integrate][1/(1+x),x]
(* Output *)
(1)/(1+x)
```

Differentiate the inactive form:

```wolfram
D[%,x]
(* Output *)
(1)/(1+x)
```

Formally differentiating a Laplace transform:

```wolfram
Inactive[LaplaceTransform][a t^2,t,s]
(* Output *)
LaplaceTransform[a t^2,t,s]
```

```wolfram
D[%,s]
(* Output *)
LaplaceTransform[-a t^3,t,s]
```

```wolfram
Activate[%]
(* Output *)
-(6 a)/(s^4)
```

```wolfram
D[LaplaceTransform[a t^2,t,s],s]
(* Output *)
-(6 a)/(s^4)
```

Similarly differentiating wrt `t` and `a`:

```wolfram
e=Inactive[LaplaceTransform][a t^2,t,s]
(* Output *)
LaplaceTransform[a t^2,t,s]
```

```wolfram
D[e,t]
(* Output *)
0
```

```wolfram
D[e,a]
(* Output *)
LaplaceTransform[t^2,t,s]
```

Differentiate operators including [Integrate](https://reference.wolfram.com/language/ref/Integrate.html):

```wolfram
Inactive[Integrate][f[x],{x,a[x],b[x]}]
(* Output *)
f[x]
```

```wolfram
D[%,x]
(* Output *)
-f[a[x]] a^′[x]+f[b[x]] b^′[x]
```

[FourierTransform](https://reference.wolfram.com/language/ref/FourierTransform.html):

```wolfram
Inactive[FourierTransform][f[t],t,ω]
(* Output *)
FourierTransform[f[t],t,ω]
```

```wolfram
D[%,ω]
(* Output *)
FourierTransform[ⅈ t f[t],t,ω]
```

[FourierSinTransform](https://reference.wolfram.com/language/ref/FourierSinTransform.html):

```wolfram
Inactive[FourierSinTransform][f[t],t,ω]
(* Output *)
FourierSinTransform[f[t],t,ω]
```

```wolfram
D[%,ω]
(* Output *)
FourierCosTransform[t f[t],t,ω]
```

[Convolve](https://reference.wolfram.com/language/ref/Convolve.html):

```wolfram
Inactive[Convolve][f[t],g[t],t,s]
(* Output *)
Convolve[f[t],g[t],t,s]
```

```wolfram
D[%,s]
(* Output *)
Convolve[f[t],g^′[t],t,s]
```

[Sum](https://reference.wolfram.com/language/ref/Sum.html):

```wolfram
Inactive[Sum][f[k,y],{k,a,b}]
(* Output *)
f[k,y]
```

```wolfram
D[%,y]
(* Output *)
f^(0,1)[k,y]
```

[ZTransform](https://reference.wolfram.com/language/ref/ZTransform.html):

```wolfram
Inactive[ZTransform][f[k],k,z]
(* Output *)
ZTransform[f[k],k,z]
```

```wolfram
D[%,z]
(* Output *)
-(ZTransform[k f[k],k,z])/(z)
```

Difference operators including [Sum](https://reference.wolfram.com/language/ref/Sum.html):

```wolfram
Inactive[Sum][f[k],k]
(* Output *)
f[k]
```

```wolfram
DifferenceDelta[%,k]
(* Output *)
f[k]
```

[DiscreteConvolve](https://reference.wolfram.com/language/ref/DiscreteConvolve.html):

```wolfram
Inactive[DiscreteConvolve][f[m],g[m],m,n]
(* Output *)
DiscreteConvolve[f[m],g[m],m,n]
```

```wolfram
DifferenceDelta[%,n]
(* Output *)
DiscreteConvolve[f[m],Inactive,m,n]
```

[Integrate](https://reference.wolfram.com/language/ref/Integrate.html):

```wolfram
Inactive[Integrate][f[x],{x,a,b}]
(* Output *)
f[x]
```

```wolfram
DifferenceDelta[%,b]
(* Output *)
f[x]
```

[LaplaceTransform](https://reference.wolfram.com/language/ref/LaplaceTransform.html):

```wolfram
Inactive[LaplaceTransform][f[t],t,s]
(* Output *)
LaplaceTransform[f[t],t,s]
```

```wolfram
DifferenceDelta[%,s]
(* Output *)
LaplaceTransform[(-1+ℯ^-t) f[t],t,s]
```

Other formal operations including [Product](https://reference.wolfram.com/language/ref/Product.html):

```wolfram
Inactive[DiscreteRatio][f[k],k]
(* Output *)
Inactive
```

```wolfram
Product[%,k]
(* Output *)
f[k]
```

[ZTransform](https://reference.wolfram.com/language/ref/ZTransform.html):

```wolfram
Inactive[InverseZTransform][F[z],z,k+1]
(* Output *)
InverseZTransform[F[z],z,1+k]
```

```wolfram
ZTransform[%,k,z]
(* Output *)
z F[z]-z InverseZTransform[F[z],z,0]
```

```wolfram
Inactive[InverseZTransform][F[z],z,k]
(* Output *)
InverseZTransform[F[z],z,k]
```

```wolfram
ZTransform[k %,k,z]
(* Output *)
-z F^′[z]
```

#### Code Transformations

[Inactivate](https://reference.wolfram.com/language/ref/Inactivate.html) a function definition:

```wolfram
Inactivate[forSquares[f_,max_]:=For[x=1,x<max,x++,Print[f[x^2]]]]
(* Output *)
forSquares[f_,max_]:=For[x=1,x<max,Increment[x],Print[f[x^2]]]
```

Transform a [For](https://reference.wolfram.com/language/ref/For.html) statement into a [Do](https://reference.wolfram.com/language/ref/Do.html) statement:

```wolfram
%//.Inactivate[For[i_=init_,i_<max_,i_++,body_]:>Do[body,{i,init,max}]]
(* Output *)
forSquares[f_,max_]:=Do[Print[f[x^2]],{x,1,max}]
```

[Activate](https://reference.wolfram.com/language/ref/Activate.html) and use the definition:

```wolfram
Activate[%]
```

```wolfram
forSquares[Range,3]
(* Output *)
{1}
(* Output *)
{1,2,3,4}
(* Output *)
{1,2,3,4,5,6,7,8,9}
```

Replace [With](https://reference.wolfram.com/language/ref/With.html) by an iterated version, so that later variables can refer to earlier ones:

```wolfram
SetAttributes[IterateWith, HoldAll]
IterateWith[expr_] := Activate[Inactivate[expr]//.Inactive[With][{first_,rest__},body_]:>Inactive[With][{first},Inactive[With][{rest},body]]]
```

```wolfram
With[{a=3,b=a+2,c=a+b},{a,b,c}]
(* Output *)
{3,2+a,a+b}
```

```wolfram
IterateWith[With[{a=3,b=a+2,c=a+b},{a,b,c}]]
(* Output *)
{3,5,8}
```

Inactivate [Length](https://reference.wolfram.com/language/ref/Length.html):

```wolfram
Inactivate[{Length[{a,b,c}], Length[{b,c}], Length[{Length[{a}],Length[{b}],Length[{c}],7}]},Length]
(* Output *)
{Length[{a,b,c}],Length[{b,c}],Length[{Length[{a}],Length[{b}],Length[{c}],7}]}
```

Change the measure to be the square of the length:

```wolfram
%//.Inactive[Length][x_]:>Length[x]^2
(* Output *)
{9,4,16}
```

### Applications

#### Basic Identities

Show a sum identity:

```wolfram
Block[{e=Inactivate[23+11]},e==Activate[e]]
(* Output *)
23+11==34
```

Make a whole table of identities:

```wolfram
Table[Block[{e=Inactivate[n+m]},e==Activate[e]],{n,0,3},{m,0,3}]//Grid//
(* Output *)
| 0+0==0 | 0+1==1 | 0+2==2 | 0+3==3 |
| --- | --- | --- | --- |
| 1+0==1 | 1+1==2 | 1+2==3 | 1+3==4 |
| 2+0==2 | 2+1==3 | 2+2==4 | 2+3==5 |
| 3+0==3 | 3+1==4 | 3+2==5 | 3+3==6 |
```

Show a product identity:

```wolfram
Block[{e=Inactivate[3 5]},e==Activate[e]]
(* Output *)
3*5==15
```

Make a whole table of identities:

```wolfram
Table[Block[{e=Inactivate[n m]},e==Activate[e]],{n,0,3},{m,0,3}]//Grid//
(* Output *)
| 0*0==0 | 0*1==0 | 0*2==0 | 0*3==0 |
| --- | --- | --- | --- |
| 1*0==0 | 1*1==1 | 1*2==2 | 1*3==3 |
| 2*0==0 | 2*1==2 | 2*2==4 | 2*3==6 |
| 3*0==0 | 3*1==3 | 3*2==6 | 3*3==9 |
```

Show algebraic identities:

```wolfram
Block[{e=Inactivate[x+1+5x]},e==Activate[e]]//
(* Output *)
x+1+5*x==6 x+1
```

```wolfram
Block[{e=Inactivate[x x x^3]},e==Activate[e]]
(* Output *)
x*x*x^3==x^5
```

Common trigonometric values:

```wolfram
Table[Inactive[Sin][k π/6]==Sin[k π/6],{k,0,3}]//Column//
(* Output *)
sin(0)==0
sin((π)/(6))==(1)/(2)
sin((π)/(3))==(Sqrt(3))/(2)
sin((π)/(2))==1
```

#### Function Identities

[Sin](https://reference.wolfram.com/language/ref/Sin.html) is an odd function of its argument:

```wolfram
Block[{e=Inactivate[Sin[-x],Sin]},e==Activate[e]]//
(* Output *)
sin(-x)==-sin(x)
```

[Cos](https://reference.wolfram.com/language/ref/Cos.html) is an even function of its argument:

```wolfram
Block[{e=Inactivate[Cos[-x],Cos]},e==Activate[e]]//
(* Output *)
cos(-x)==cos(x)
```

Hyperbolic functions with imaginary arguments are equivalent to trigonometric functions:

```wolfram
Block[{e=Inactivate[Sinh[I x],Sinh]},e==Activate[e]]//
(* Output *)
sinh(ⅈ x)==ⅈ sin(x)
```

```wolfram
Block[{e=Inactivate[Cosh[I x],Cosh]},e==Activate[e]]//
(* Output *)
cosh(ⅈ x)==cos(x)
```

[BesselJ](https://reference.wolfram.com/language/ref/BesselJ.html)[1,*x*] is an odd function of `*x*`:

```wolfram
Block[{e=Inactivate[BesselJ[1,-x],BesselJ]},[e]==Activate[e]]
(* Output *)
BesselJ(1,-x)==-BesselJ[1,x]
```

[BesselJ](https://reference.wolfram.com/language/ref/BesselJ.html)[2,*x*] is an even function of `*x*`:

```wolfram
Block[{e=Inactivate[BesselJ[2,-x],BesselJ]},[e]==Activate[e]]
(* Output *)
BesselJ(2,-x)==BesselJ[2,x]
```

#### Calculus Identities

Show identities, including Leibniz's rule for differentiating integrals:

```wolfram
Inactivate[∂_x∫_a[x]^b[x]f[x]ⅆx,D|Integrate]==∂_x∫_a[x]^b[x]f[x]ⅆx
(* Output *)
Inactive==-f[a[x]] a^′[x]+f[b[x]] b^′[x]
```

Product rule:

```wolfram
Inactivate[D[f[x] g[x],x],D]==D[f[x]g[x],x]
(* Output *)
Inactive==g[x] f^′[x]+f[x] g^′[x]
```

Chain rule:

```wolfram
Inactivate[D[f[g[x]],x],D]==D[f[g[x]],x]
(* Output *)
Inactive==f^′[g[x]] g^′[x]
```

Indefinite integrals:

```wolfram
Inactivate[Integrate[x^n,x],Integrate]==Integrate[x^n,x]
(* Output *)
x^n==(x^(1+n))/(1+n)
```

```wolfram
Inactivate[Integrate[1/Sqrt[1-x^2],x],Integrate]==Integrate[1/Sqrt[1-x^2],x]
(* Output *)
(1)/(Sqrt[1-x^2])==ArcSin[x]
```

Infinite sums and products:

```wolfram
Inactivate[Sum[1/n^2,{n,1,Infinity}],Sum]==Sum[1/n^2,{n,1,Infinity}]//
(* Output *)
(1)/(n^(2))==(π^(2))/(6)
```

```wolfram
Inactivate[Product[1-1/n^4,{n,2,Infinity}],Product]==Product[1-1/n^4,{n,2,Infinity}]//
(* Output *)
(1-(1)/(n^(4)))==(sinh(π))/(4 π)
```

Finite and infinite continued fractions:

```wolfram
Inactive[ContinuedFractionK][1,2,{n,1,m}]==ContinuedFractionK[1,2,{n,1,m}]//
(* Output *)
1==(2)/((2 Sqrt(2) ((-3-2 Sqrt(2))^(m)+1))/((-3-2 Sqrt(2))^(m)-1)+2)
```

```wolfram
Inactive[ContinuedFractionK][n,{n,1,Infinity}]==ContinuedFractionK[n,{n,1,Infinity}]//
(* Output *)
1AutoDelete -> TrueEditable -> True]==(1)/(0)
```

Apply [DifferenceDelta](https://reference.wolfram.com/language/ref/DifferenceDelta.html) to an inactive sum:

```wolfram
Timing[DifferenceDelta[Inactive[Sum][(1 + x)^150, x], x]]
(* Output *)
{0.,(1+x)^150}
```

This is significantly faster than the evaluated version:

```wolfram
Timing[DifferenceDelta[Sum[(1 + x)^150, x], x]]
(* Output *)
{1.70041090000000005844071893079672008753,(1+x)^150}
```

Formula for summation by parts:

```wolfram
sumparts=Inactivate[Sum[f[k] g[k],k]==g[k]Sum[f[k],k]-Sum[DifferenceDelta[g[k],k]DiscreteShift[Sum[f[k],k],k],k],Sum|DifferenceDelta|DiscreteShift]
(* Output *)
f[k] g[k]==g[k] f[k]-Inactive Inactive
```

Verify the formula in a special case:

```wolfram
f[k_]:=k
```

```wolfram
g[k_]:=PolyGamma[k]
```

```wolfram
Activate[sumparts]
(* Output *)
True
```

Evaluate the sum:

```wolfram
Inactive[Sum][k PolyGamma[k],k]==Sum[k PolyGamma[k],k]
(* Output *)
k PolyGamma[0,k]==-(1)/(4) k (1+k)+(1)/(2) (-1+k) k PolyGamma[0,k]
```

Interchange the order of summation and integration:

```wolfram
(id=Inactivate[Sum[Integrate[t^(n-1)E^(-t)/(2n+1)!!,
     {t,0,Infinity}],{n,1,Infinity}] ==
 Integrate[
   Sum[E^(-t)t^(n-1)/(2n+1)!!,{n,1,Infinity}],{t,0,Infinity}],Sum|Integrate])//
(* Output *)
(ℯ^(-t) t^(n-1))/((2 n+1)!!)==(ℯ^(-t) t^(n-1))/((2 n+1)!!)
```

Evaluate both sides of the identity:

```wolfram
{Activate[id[[1]]],Activate[id[[2]]]}//Together
(* Output *)
{(4-π)/(2),(4-π)/(2)}
```

Obtain the same result using the corresponding sum:

```wolfram
∑_{n=1}^{∞}(Gamma[n])/((2 n+1)!!)
(* Output *)
(4-π)/(2)
```

The product rule for the [Laplacian](https://reference.wolfram.com/language/ref/Laplacian.html):

```wolfram
Inactivate[Laplacian[f[x,y] g[x,y],{x,y}]==f[x,y] Laplacian[g[x,y],{x,y}]+g[x,y] Laplacian[f[x,y],{x,y}]+2 Grad[f[x,y],{x,y}].Grad[g[x,y],{x,y}],Laplacian|Grad]//
(* Output *)
(f(x,y) g(x,y))==2 f(x,y).g(x,y)+g(x,y) f(x,y)+f(x,y) g(x,y)
```

```wolfram
Activate[%]//Simplify
(* Output *)
True
```

Vector identities for three-vectors `u`, `v`, and `w`:

```wolfram
u={ux,uy,uz};
v={vx,vy,vz};
w={wx,wy,wz};
```

Antisymmetry of the cross product:

```wolfram
Inactivate[Cross[v,v],Cross]==Cross[v,v]
(* Output *)
{vx,vy,vz}⨯{vx,vy,vz}=={0,0,0}
```

```wolfram
%//Activate
(* Output *)
True
```

Orthogonality of the cross product:

```wolfram
Inactive[Dot][u,Inactive[Cross][u,v]]==Inactive[Dot][v,Inactive[Cross][u,v]]==0
(* Output *)
{ux,uy,uz}.{ux,uy,uz}⨯{vx,vy,vz}=={vx,vy,vz}.{ux,uy,uz}⨯{vx,vy,vz}==0
```

```wolfram
%//Activate//Simplify
(* Output *)
True
```

Scalar triple product:

```wolfram
Inactive[Dot][u,Inactive[Cross][v,w]]==Inactive[Det][{u,v,w}]
(* Output *)
{ux,uy,uz}.{vx,vy,vz}⨯{wx,wy,wz}==Det[{{ux,uy,uz},{vx,vy,vz},{wx,wy,wz}}]
```

```wolfram
%//Activate//Simplify
(* Output *)
True
```

#### Derive Identities

The basic commutation trick for differentiating under the integral or summation sign:

```wolfram
Inactivate[Integrate[D[f[x,a],a],x]==((D[Integrate[f[x,a],x],a]))]
(* Output *)
Inactive==Inactive
```

```wolfram
Inactivate[Sum[D[f[x,a],a],x]==((D[Sum[f[x,a],x],a]))]
(* Output *)
Inactive==Inactive
```

Derive a closed form for $\int x e^{x}\mathrm{d}x$ by differentiating $\int e^{a x}\mathrm{d}x$ with respect to $a$ at $a=1$:

```wolfram
D[Inactive[Integrate][E^(a x),x],a]
(* Output *)
ℯ^(a x) x
```

```wolfram
lhs=%/.a->1
(* Output *)
ℯ^x x
```

Now integrate $\int e^{a x}\mathrm{d}x$ and then differentiate with respect to $a$ at $a=1$:

```wolfram
Integrate[E^(a x),x]
(* Output *)
(ℯ^(a x))/(a)
```

```wolfram
D[%,a]
(* Output *)
-(ℯ^(a x))/(a^2)+(ℯ^(a x) x)/(a)
```

```wolfram
rhs=%/.a->1//FullSimplify
(* Output *)
ℯ^x (-1+x)
```

The final result:

```wolfram
lhs==rhs
(* Output *)
ℯ^x x==ℯ^x (-1+x)
```

Verify the result:

```wolfram
Activate[%]
(* Output *)
True
```

Derive a closed form for $\int_{1}^{\infty}e^{-x} log(x)\mathrm{d}x$ by differentiating $\int_{1}^{\infty}e^{-x} x^{a}\mathrm{d}x$ with respect to $a$ at zero:

```wolfram
D[Inactive[Integrate][E^(-x)x^a,{x,1,∞}],a]
(* Output *)
ℯ^-x x^a Log[x]
```

```wolfram
lhs=%/.a->0
(* Output *)
ℯ^-x Log[x]
```

$\int_{1}^{\infty}e^{-x} x^{a}\mathrm{d}x$ is first integrated and then differentiated with respect to $a$ at zero:

```wolfram
Integrate[E^(-x)x^a,{x,1,∞}]
(* Output *)
Gamma[1+a,1]
```

```wolfram
D[%,a]
(* Output *)
MeijerG[{{},{1,1}},{{0,0,1+a},{}},1]
```

```wolfram
rhs=%/.a->0//FullSimplify
(* Output *)
-ExpIntegralEi[-1]
```

The final result:

```wolfram
lhs==rhs
(* Output *)
ℯ^-x Log[x]==-ExpIntegralEi[-1]
```

Verify the result:

```wolfram
Activate[%]
(* Output *)
True
```

Derive a closed form for $\sum_{k}k e^{-k}$ by differentiating $\sum_{k}e^{a k}$ with respect to $a$ at $a=-1$:

```wolfram
D[Inactive[Sum][E^(a k), k], a]
(* Output *)
ℯ^(a k) k
```

```wolfram
lhs=%/.a->-1
(* Output *)
ℯ^-k k
```

Compute $\sum_{k}e^{a k}$ and then differentiate:

```wolfram
Sum[E^(a k),k]
(* Output *)
(ℯ^(a k))/(-1+ℯ^a)
```

```wolfram
D[%,a]
(* Output *)
-(ℯ^(a+a k))/((-1+ℯ^a)^2)+(ℯ^(a k) k)/(-1+ℯ^a)
```

```wolfram
rhs=%/.a->-1//Simplify
(* Output *)
-(ℯ^(1-k) (1+(-1+ℯ) k))/((-1+ℯ)^2)
```

The result:

```wolfram
lhs==rhs
(* Output *)
ℯ^-k k==-(ℯ^(1-k) (1+(-1+ℯ) k))/((-1+ℯ)^2)
```

Verify the result:

```wolfram
Activate[%]//Simplify
(* Output *)
True
```

Derive a closed form for $\sum_{k=-\infty}^{\infty}k^{2} exp(-k^{2})$ by differentiating $\sum_{k=-\infty}^{\infty}exp(a k^{2})$ wrt $a$ at $a=-1$:

```wolfram
D[Inactive[Sum][Exp[a k^2],{k,-∞,∞}],a]
(* Output *)
ℯ^(a k^2) k^2
```

```wolfram
lhs=%/.a->-1
(* Output *)
ℯ^(-k^2) k^2
```

Compute $\sum_{k=-\infty}^{\infty}exp(a k^{2})$ and then differentiate:

```wolfram
Sum[Exp[a k^2],{k,-∞,∞}]
(* Output *)
EllipticTheta[3,0,ℯ^a]
```

```wolfram
D[%,a]
(* Output *)
ℯ^a EllipticTheta^(0,0,1)[3,0,ℯ^a]
```

```wolfram
rhs=%/.a->-1
(* Output *)
(EllipticTheta^(0,0,1)[3,0,(1)/(ℯ)])/(ℯ)
```

The result:

```wolfram
lhs==rhs
(* Output *)
ℯ^(-k^2) k^2==(EllipticTheta^(0,0,1)[3,0,(1)/(ℯ)])/(ℯ)
```

Generalize to $\sum_{k=-\infty}^{\infty}exp(-k^{2}) k^{2 n}$:

```wolfram
Table[D[Inactive[Sum][Exp[a k^2] ,{k,-∞,∞}],{a,n}]==D[Sum[Exp[a k^2] ,{k,-∞,∞}],{a,n}]/.a->-1,{n,1,3}]//Column
(* Output *)
{{ℯ^(-k^2) k^2==(EllipticTheta^(0,0,1)[3,0,(1)/(ℯ)])/(ℯ)}, {ℯ^(-k^2) k^4==(EllipticTheta^(0,0,1)[3,0,(1)/(ℯ)])/(ℯ)+(EllipticTheta^(0,0,2)[3,0,(1)/(ℯ)])/(ℯ^2)}, {ℯ^(-k^2) k^6==(EllipticTheta^(0,0,1)[3,0,(1)/(ℯ)])/(ℯ)+(3 EllipticTheta^(0,0,2)[3,0,(1)/(ℯ)])/(ℯ^2)+(EllipticTheta^(0,0,3)[3,0,(1)/(ℯ)])/(ℯ^3)}}
```

#### Solve Differential Equations

Solution for the three-dimensional Laplace equation in inactive integral form:

```wolfram
V[x_,y_,z_]:=f[z+ⅈ x Cos[u]+ⅈ y Sin[u],u]
```

Obtain a particular solution by specifying the function `*f*`:

```wolfram
f[a_,b_]:=2a^5+b^4
```

```wolfram
V[x,y,z]
(* Output *)
(u^4+2 (z+ⅈ x Cos[u]+ⅈ y Sin[u])^5)
```

```wolfram
sol=Activate[%]//Simplify
(* Output *)
(2 π^5)/(5)+(1)/(2) π z (15 (x^2+y^2)^2-40 (x^2+y^2) z^2+8 z^4)
```

Visualize the solution:

```wolfram
Row[Table[Plot3D[sol/.{z-> j},{x,-3,3},{y,-3,3},Ticks->{Automatic,Automatic,None}],{j,-2,2}]]
(* Output *)
![image](img/image_001.png)
```

Verify the solution:

```wolfram
Laplacian[sol,{x,y,z}]//Simplify
(* Output *)
0
```

Maxwell's equations in natural Lorentz-Heaviside units:

```wolfram
Inactivate[{gaussE=Div[ℰ[x,y,z,t],{x,y,z}]==ρ[x,y,z,t],gaussB=Div[ℬ[x,y,z,t],{x,y,z}]==0,faraday=Curl[ℰ[x,y,z,t],{x,y,z}]==-D[ℬ[x,y,z,t],t],ampere=Curl[ℬ[x,y,z,t],{x,y,z}]== j[x,y,z,t]+D[ℰ[x,y,z,t],t]},Div|Curl|D];
%//TableForm//
(* Output *)
ℰ(x,y,z,t)==ρ(x,y,z,t)
ℬ(x,y,z,t)==0
ℰ(x,y,z,t)==-(∂ℬ(x,y,z,t))/(∂t)
ℬ(x,y,z,t)==(∂ℰ(x,y,z,t))/(∂t)+j(x,y,z,t)
```

Take the curl of Ampere's law in a vacuum ($j=0$ and $\rho=0$):

```wolfram
Inactive[Curl][#, {x,y,z}]& /@ ampere /. j->(0&)
(* Output *)
(ℬ[x,y,z,t])==Inactive
```

Interchange order of differentiation:

```wolfram
% /. Inactive[Curl][Inactive[D][v_,t_],x_] :> Inactive[D][Inactive[Curl][v,x],t]
(* Output *)
(ℬ[x,y,z,t])==Inactive
```

Substitute in Faraday's law:

```wolfram
% /. Solve[faraday,Inactive[Curl][ℰ[x,y,z,t],{x,y,z}]][[1]]
(* Output *)
(ℬ[x,y,z,t])==Inactive
```

Activate the equation, resulting in the wave equation for the magnetic field:

```wolfram
wave=Activate[%]
(* Output *)
ℬ^(0,0,2,0)[x,y,z,t]+ℬ^(0,2,0,0)[x,y,z,t]+ℬ^(2,0,0,0)[x,y,z,t]==-ℬ^(0,0,0,2)[x,y,z,t]
```

Verify plane-wave solutions of the equation:

```wolfram
wave/.ℬ->(A Exp[{kx,ky,kz}.{#1,#2,#3}-ISqrt[kx^2+ky^2+kz^2]#4]&)//Simplify
(* Output *)
True
```

#### Derive a Least-Squares Solution

Define the sum of squares of the vertical deviations for a given set of data:

```wolfram
squareDeviations=Inactivate[Sum[(y[[x]]-(a x+b))^2,{x,1,n}],Sum|Part]
(* Output *)
(-b-a x+Inactive)^2
```

Set up the least-squares equations:

```wolfram
eqna=D[squareDeviations,a]==0
(* Output *)
-2 x (-b-a x+Inactive)==0
```

```wolfram
eqnb=D[squareDeviations,b]==0
(* Output *)
-2 (-b-a x+Inactive)==0
```

Generate some data:

```wolfram
n=200;y=Table[3x+20+RandomReal[{-15,15}],{x,1,n}];
```

```wolfram
ListPlot[y]
```

*([Graphics])*

Solve the least-squares problem for this data:

```wolfram
a x+ b /.Solve[Activate[{eqna,eqnb}],{a,b}][[1]]
(* Output *)
20.5207717227377+2.991390641741944 x
```

#### Code Transformation

Replace uses of [Module](https://reference.wolfram.com/language/ref/Module.html) using [Block](https://reference.wolfram.com/language/ref/Block.html) with unique variables:

```wolfram
Blockhead[icode_] :=
icode //. Inactive[Module][vars_, body_]:>
(Inactive[Block][vars, body] /. (rule = Thread[#->Unique[#]])&[Replace[vars,(Inactive[Set][var_,val_]| var_):> var,{1}]])
```

```wolfram
code = Inactivate[f[s_] := Module[{x, y = Module[{x = 5 s},x+7]}, x= y ^2; Sin[x]]]
(* Output *)
f[s_]:=Module[{x,y=Module[{x=5*s},x+7]},x=y^2;Sin[x]]
```

Applying the function replaces the [Module](https://reference.wolfram.com/language/ref/Module.html) local variables with unique variables:

```wolfram
bcode = Blockhead[code]
(* Output *)
f[s_]:=Block[{x$11002,y$11002=Block[{x$11002$11003=5*s},x$11002$11003+7]},x$11002=y$11002^2;Sin[x$11002]]
```

Activate the code and the transformed code to make definitions for `f` and `fb`:

```wolfram
Activate[code];
Activate[bcode /. f->fb];
```

Compare values for random test values:

```wolfram
testvalues = RandomReal[1, 5];
{f[#],fb[#]}& /@ testvalues
(* Output *)
{{0.872080176111264,0.872080176111264},{-0.9986741169611179,-0.9986741169611179},{0.5174326011384289,0.5174326011384289},{0.8919841211875692,0.8919841211875692},{0.24778347184546407,0.24778347184546407}}
```

Compare timings for a large set of test values:

```wolfram
testvalues = RandomReal[1,10^5];
{AbsoluteTiming[Total[f /@ testvalues]],AbsoluteTiming[Total[fb /@ testvalues]]}
(* Output *)
{{1.91893530000000001045634689944563433528,-62.055825033666395},{0.88926269999999996151984760217601433396,-62.055825033666395}}
```

### Properties & Relations

Inactive expressions can be created using [Inactivate](https://reference.wolfram.com/language/ref/Inactivate.html):

```wolfram
Inactivate[f[x],f]
(* Output *)
f[x]
```

Inactive expressions can be evaluated using [Activate](https://reference.wolfram.com/language/ref/Activate.html):

```wolfram
f[x_]:=x^2
```

```wolfram
expr=Inactive[f][3]
(* Output *)
f[3]
```

```wolfram
Activate[expr]
(* Output *)
9
```

[Inactive](https://reference.wolfram.com/language/ref/Inactive.html) creates inactive forms of symbols and allows parts of expressions to be inactive:

```wolfram
isin=Inactive[Sin][ArcTan[1]]
(* Output *)
Sin[(π)/(4)]
```

```wolfram
Activate[isin]
(* Output *)
(1)/(Sqrt[2])
```

[Hold](https://reference.wolfram.com/language/ref/Hold.html) maintains the expression in unevaluated form, and all parts are inactive:

```wolfram
esin=Hold[Sin[ArcTan[1]]]
(* Output *)
Hold[Sin[ArcTan[1]]]
```

```wolfram
ReleaseHold[%]
(* Output *)
(1)/(Sqrt[2])
```

Compare an inactive expression with the corresponding [FullForm](https://reference.wolfram.com/language/ref/FullForm.html):

```wolfram
Inactivate[Cos[x  Sin[y]]]
(* Output *)
Cos[x*Sin[y]]
```

```wolfram
FullForm[Cos[x  Sin[y]]]
(* Output *)
Cos[Times[x,Sin[y]]]
```

[Inactive](https://reference.wolfram.com/language/ref/Inactive.html) prevents the attributes of its argument from having any effect:

```wolfram
SetAttributes[f,Listable]
Inactive[f][{1,2,3}]
(* Output *)
f[{1,2,3}]
```

Normally, the [Listable](https://reference.wolfram.com/language/ref/Listable.html) attribute would cause `f` to thread over its arguments:

```wolfram
f[{1,2,3}]
(* Output *)
{f[1],f[2],f[3]}
```

### Possible Issues

[Inactive](https://reference.wolfram.com/language/ref/Inactive.html)[*h*] does not have the attributes of `*h*`, which may lead to evaluation leaks:

```wolfram
k=1;
Sum[f[k],{k,1,n}]
(* Output *)
∑_{k=1}^{n}f[k]
```

The value of `k` leaks through in the following, because [Inactive](https://reference.wolfram.com/language/ref/Inactive.html)[[Sum](https://reference.wolfram.com/language/ref/Sum.html)] lacks the [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html) attribute:

```wolfram
Inactive[Sum][f[k],{k,1,n}]
(* Output *)
f[1]
```

### Neat Examples

Create a gallery of multivariate sums:

```wolfram
multisums=Inactivate[{∑_{m=2}^{∞}∑_{n=1}^{∞}(1)/((2 n)^m),∑_{m=1}^{∞}∑_{n=1}^{∞}(1)/((4 n-1)^(2 m)),∑_{i=1}^{∞}∑_{j=1}^{∞}(Zeta[i+j])/(2^(i+j)),  ∑_{i=1}^{∞}∑_{j=1}^{∞}(1)/((i+j)!),∑_{i=1}^{∞}∑_{j=1}^{∞}(1)/(Max[i,j]!),∑_{i=1}^{∞}∑_{j=1}^{∞}(1)/(Max[i,j]^3),∑_{i=1}^{∞}∑_{j=1}^{∞}∑_{k=1}^{∞}(1)/((i j+j k)^s),∑_{m=1}^{∞}∑_{n=1}^{∞}(Zeta[m+2 n])/(4^((m)/(2)+n)),∑_{i=0}^{∞}∑_{j=0}^{∞}∑_{k=1}^{∞}(1)/((i+2 j+k)^4)},Sum];
```

```wolfram
FormulaGallery[forms_List]:=Grid[ParallelMap[{#==Simplify@Activate[#]}&,forms],Dividers->All,Alignment->{"==",Baseline},Background->{{},{{LightDarkSwitched[<|color -> RGBColor[0.87, 0.94, 1]|>],None}}},Spacings->{4,1}];
```

```wolfram
FormulaGallery[multisums]//
(* Output *)
2^(-m) n^(-m)==log(2)
(4 n-1)^(-2 m)==(log(2))/(4)
2^(-i-j) i+j==(π^(2))/(8)
(1)/((i+j)!)==1
(1)/(max(i,j)!)==1+ℯ
(1)/(max(i,j)^(3))==(1)/(3) (π^(2)-3 3)
(i j+j k)^(-s)==(s-1-s) s
4^(-(m)/(2)-n) m+2 n==(1)/(16) (π^(2)-4)
(1)/((i+2 j+k)^(4))==(1)/(384) (192 3+16 π^(2)+π^(4))
```

## Related Guides ▪Partial Differential Equations ▪Evaluation Control ▪Function Composition & Operator Forms ▪Tuning & Debugging ▪Expressions

[Graphics] | Related Workflows
▪Handle Code Symbolically

## History Introduced in 2014 (10.0)
