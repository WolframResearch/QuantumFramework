# GroebnerBasis | [SpanFromLeft]

> [GroebnerBasis](https://reference.wolfram.com/language/ref/GroebnerBasis.html)[{*poly*_1,*poly*_2,…},{*x*_1,*x*_2,…}] — gives a list of polynomials that form a Gröbner basis for the set of polynomials `*poly*_*i*`.
> [GroebnerBasis](https://reference.wolfram.com/language/ref/GroebnerBasis.html)[{*poly*_1,*poly*_2,…},{*x*_1,*x*_2,…},{*y*_1,*y*_2,…}] — finds a Gröbner basis in which the `*y*_*i*` have been eliminated.

## Details and Options

The set of polynomials in a Gröbner basis have the same collection of roots as the original polynomials.

For polynomials in one variable, [GroebnerBasis](https://reference.wolfram.com/language/ref/GroebnerBasis.html) reduces to [PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html).

For linear functions in any number of variables, [GroebnerBasis](https://reference.wolfram.com/language/ref/GroebnerBasis.html) is equivalent to Gaussian elimination.

The Gröbner basis in general depends on the ordering assigned to monomials. This ordering is affected by the ordering of the `*x*_*i*`.

The following options can be given:

| MonomialOrder | Lexicographic | the criterion used for ordering monomials |
| --- | --- | --- |
| CoefficientDomain | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | the type of objects assumed to be coefficients |
| [Method](https://reference.wolfram.com/language/ref/Method.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | the method to use |
| [Modulus](https://reference.wolfram.com/language/ref/Modulus.html) | 0 | the modulus for numerical coefficients |

Possible settings for `MonomialOrder` are `Lexicographic`, `DegreeLexicographic`, `DegreeReverseLexicographic`, `EliminationOrder`, or an explicit weight matrix. Monomials are specified for the purpose of `MonomialOrder` by lists of the exponents with which the `*x*_*i*` appear in them.

The ordering of the `*x*_*i*` and the setting for `MonomialOrder` can substantially affect the efficiency of [GroebnerBasis](https://reference.wolfram.com/language/ref/GroebnerBasis.html).

Possible settings for `CoefficientDomain` are `InexactNumbers`, [Rationals](https://reference.wolfram.com/language/ref/Rationals.html), `RationalFunctions`, and `Polynomials[x]`.

Possible settings for the [Method](https://reference.wolfram.com/language/ref/Method.html) option include `"Buchberger"` and `"GroebnerWalk"`.

## Examples

### Basic Examples

Compute a Gröbner basis:

```wolfram
GroebnerBasis[{x^2-2y^2,x y-3},{x,y}]
(* Output *)
{-9+2 y^4,3 x-2 y^3}
```

Prove that polynomials have no common roots:

```wolfram
GroebnerBasis[{x+y,x^2-1,y^2-2x},{x,y}]
(* Output *)
{1}
```

### Scope

Polynomials with a finite number of common roots:

```wolfram
GroebnerBasis[{x^2+y^2+z^2-1,x y-z+2,z^2-2x+3y},{x,y,z}]
(* Output *)
{1024-832 z-215 z^2+156 z^3-25 z^4+24 z^5+13 z^6+z^8,-11552+2560 y+2197 z+2764 z^2+443 z^3+728 z^4+169 z^5+32 z^6+13 z^7,-34656+5120 x+6591 z+5732 z^2+1329 z^3+2184 z^4+507 z^5+96 z^6+39 z^7}
```

Polynomials with an infinite number of common roots:

```wolfram
GroebnerBasis[{x^2+y^2+z^2-1,x y-z+2},{x,y,z}]
(* Output *)
{4-y^2+y^4-4 z+z^2+y^2 z^2,-2 x-y+y^3+x z+y z^2,2+x y-z,-1+x^2+y^2+z^2}
```

Polynomials with no common roots:

```wolfram
GroebnerBasis[{x^2+y^2+z^2-1,x y-z+2,z^2-3+x,x-y^2+1},{x,y,z}]
(* Output *)
{1}
```

Eliminate a variable:

```wolfram
GroebnerBasis[{x^2+y^2+z^2-1,x y-z+2,z^2-2x+3y},{x,y},{z}]
(* Output *)
{28+8 y+41 y^2-26 y^3+55 y^4-8 y^5+7 y^6-6 y^7+y^8,16+20 x-10 y+28 y^2-57 y^3+4 y^4-6 y^5+6 y^6-y^7}
```

A lexicographic Gröbner basis:

```wolfram
polys={x^2+y^2+z^2-1,x-z+2,z^2-x y};
```

```wolfram
GroebnerBasis[polys,{x,y,z}]
(* Output *)
{12-28 z+27 z^2-12 z^3+3 z^4,-6+4 y+11 z-6 z^2+3 z^3,2+x-z}
```

A degree reverse lexicographic Gröbner basis:

```wolfram
GroebnerBasis[polys,{x,y,z},MonomialOrder->DegreeReverseLexicographic]
(* Output *)
{2+x-z,-2 y+y z-z^2,3+y^2-4 z+2 z^2,-6+4 y+11 z-6 z^2+3 z^3}
```

### Generalizations & Extensions

Polynomial equations can be given instead of polynomials:

```wolfram
GroebnerBasis[{x^2-2y^2==1,x y==3},{x,y}]
(* Output *)
{-9+y^2+2 y^4,3 x-y-2 y^3}
```

### Options

#### CoefficientDomain

By default, Gröbner bases are computed over the field of rational numbers:

```wolfram
polys={a x^2+5x-1, 2x+3 x y+y^2};
```

```wolfram
GroebnerBasis[polys,{x,y}]
(* Output *)
{-4-12 y-19 y^2-15 y^3+a y^4,18+4 a x+27 y+45 y^2+2 a y^2-3 a y^3,2 x+3 x y+y^2}
```

This computes the strong Gröbner basis over the ring of integers:

```wolfram
GroebnerBasis[polys,{x,y},CoefficientDomain->Integers]
(* Output *)
{-4-12 y-19 y^2-15 y^3+a y^4,18+4 a x+27 y+45 y^2+2 a y^2-3 a y^3,2 x+3 x y+y^2,6+2 a x+9 y+a x y+15 y^2+a y^2-a y^3,-1+5 x+a x^2}
```

This computes the Gröbner basis over the field of rational functions ℚ(`*a*`):

```wolfram
GroebnerBasis[polys,{x,y},CoefficientDomain->RationalFunctions]
(* Output *)
{-4-12 y-19 y^2-15 y^3+a y^4,18+4 a x+27 y+(45+2 a) y^2-3 a y^3}
```

This uses approximate arithmetic:

```wolfram
GroebnerBasis[polys,{x,y},CoefficientDomain->InexactNumbers[20]]
(* Output *)
{-4.00000000000000000000000000009-12.00000000000000000000000000016 y-19.00000000000000000000000000024 y^2-15.00000000000000000000000000015 y^3+1. a y^4,4.50000000000000000000000000021+1. a x+6.75000000000000000000000000028 y+11.25000000000000000000000000049 y^2+0.50000000000000000000000000003 a y^2-0.75000000000000000000000000003 a y^3,0.66666666666666666666666666667 x+1. x y+0.33333333333333333333333333333 y^2}
```

#### Method

The [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) method setting uses `"GroebnerWalk"` for lexicographic bases over the rationals:

```wolfram
polys={x y^4+y z^4-2 x^2y-3,y^4+x y^2z+x^2-2x y+y^2+z^2,-x^3y^2+x y z^3+y^4+x y^2z-2x y};
```

```wolfram
GroebnerBasis[polys,{x,y,z}];//Timing
(* Output *)
{0.3900000000000097,Null}
```

In this case the `"Buchberger"` method is much slower than `"GroebnerWalk"`:

```wolfram
TimeConstrained[GroebnerBasis[polys,{x,y,z},Method->"Buchberger"];,60]
(* Output *)
$Aborted
```

These polynomials are close to the lexicographic Gröbner basis:

```wolfram
polys={z^19-53z^7+11z^5-21z+17,x-z^50,y-2z^47-3z^21};
```

The `"GroebnerWalk"` method computes the degree reverse lexicographic basis first:

```wolfram
GroebnerBasis[polys,{x,y,z},Method->"GroebnerWalk"]//Timing
(* Output *)
{1.4220000000000235,{17-21 z+11 z^5-53 z^7+z^19,4114+y-5082 z-39593 z^2+48909 z^3+95506 z^4-115316 z^5-38445 z^7+184657 z^9+1428 z^10-298636 z^11-748 z^14+924 z^15+3604 z^16-4452 z^17,-30634+x+75684 z-46746 z^2+2057 z^3-2541 z^4-39644 z^5+48972 z^6+143259 z^7-175636 z^8-19239 z^10+92408 z^12+714 z^13-149318 z^14-374 z^17+462 z^18}}
```

Computing the lexicographic basis directly with the `"Buchberger"` method is faster here:

```wolfram
GroebnerBasis[polys,{x,y,z},Method->"Buchberger"]//Timing
(* Output *)
{0.,{17-21 z+11 z^5-53 z^7+z^19,4114+y-5082 z-39593 z^2+48909 z^3+95506 z^4-115316 z^5-38445 z^7+184657 z^9+1428 z^10-298636 z^11-748 z^14+924 z^15+3604 z^16-4452 z^17,-30634+x+75684 z-46746 z^2+2057 z^3-2541 z^4-39644 z^5+48972 z^6+143259 z^7-175636 z^8-19239 z^10+92408 z^12+714 z^13-149318 z^14-374 z^17+462 z^18}}
```

#### Modulus

This computes the Gröbner basis over the field of integers modulo 7:

```wolfram
polys={3 x^2+y z-5x-1, 2x+3 x y+y^2,x-3y+x z-2 z^2};
```

```wolfram
GroebnerBasis[polys,{x,y,z},Modulus->7]
(* Output *)
{6+6 z+z^2+4 z^3+3 z^4+6 z^5+3 z^6+z^7,1+4 y+4 z+y z+4 z^3+z^4+z^6,1+3 y+y^2+6 z^2+3 z^3+3 z^4+3 z^5+4 z^6,1+x+y+3 z+2 z^2+2 z^3+6 z^4+z^5}
```

#### MonomialOrder

By default, [GroebnerBasis](https://reference.wolfram.com/language/ref/GroebnerBasis.html) uses the `Lexicographic` monomial order:

```wolfram
polys={ -5x^2+y z-x-1, 2x+3 x y+y^2,x-3y+x z-2 z^2};
```

```wolfram
GroebnerBasis[polys,{x,y,z}]
(* Output *)
{1540+82 z-3491 z^2+676 z^3+1807 z^4-1768 z^5+310 z^6+1296 z^7+520 z^8,16801600618344+60608534091110 y+19757744527703 z+21523317746112 z^2-16946787439081 z^3+4872263781864 z^4-16795641902050 z^5-17077664603568 z^6-5202247795160 z^7,-12397969782796+121217068182220 x+124916501217576 z+18417605948925 z^2-169792044662095 z^3+66190507027566 z^4-36996183157154 z^5-88258341985544 z^6-38225933906680 z^7}
```

This gives the Gröbner basis in the `DegreeReverseLexicographic` monomial order:

```wolfram
GroebnerBasis[polys,{x,y,z},MonomialOrder->DegreeReverseLexicographic]
(* Output *)
{x-3 y+x z-2 z^2,2 x+3 x y+y^2,1+x+5 x^2-y z,-1+27 y+5 y^2-z-29 y z+18 z^2+y z^2-20 z^3,6-156 y-20 y^2+6 z+174 y z+y^2 z-104 z^2+120 z^3,180-20 x-4185 y-559 y^2+15 y^3+162 z+4680 y z-2808 z^2+3240 z^3,4026-20 x-106386 y-17140 y^2+4086 z+114129 y z-70866 z^2+78768 z^3+1560 z^4}
```

A monomial order may be specified by giving a full rank square rational weight matrix:

```wolfram
wmat={{1,3,1},{-1,2,0},{4,-3,2}};
```

For the order to be well-founded the first nonzero entry in each column must be positive:

```wolfram
GroebnerBasis[polys,{x,y,z},MonomialOrder->wmat]
(* Output *)
{-21 x-13 x^2-15 x^3+6 z+6 x z+20 x^2 z+20 x z^2,60+39 x+287 x^2-15 x^3+6 z-14 x z+20 x^2 z+40 z^3,-x+3 y-x z+2 z^2,10-38 x+141 x^2-305 x^3+325 x^4-12 z-82 x z-50 x^2 z+20 z^2,30+913 x+189 x^2+4545 x^3+42 z+1132 x z+630 x^2 z+1950 x^3 z-200 z^2}
```

Eliminate `*z*` and return a degree reverse lexicographic basis with respect to `{*x*,*y*}`:

```wolfram
GroebnerBasis[polys,{x,y},z,MonomialOrder->EliminationOrder]
(* Output *)
{2 x+3 x y+y^2,486+1142 x+5688 x^2+5670 x^3+12150 x^4+85 y^2+723 y^3+45 y^4,2916-14944 x-14760 x^2-16200 x^3+4374 y-10388 y^2+16311 y^3+2187 y^4+1755 y^5}
```

#### ParameterVariables

Parameters are ordered lexicographically after all other variables:

```wolfram
polys={ -5x^2+y z-x-1, 2x+3 x y+y^2,x-3y+x z-2 z^2};
```

```wolfram
GroebnerBasis[polys,ParameterVariables->x]
(* Output *)
{2+8 x+113 x^2+212 x^3+938 x^4+364 x^5+1980 x^6-1000 x^7+1625 x^8,-289274-1958323 x-7458871 x^2-14140053 x^3-20245223 x^4-19904560 x^5+4016500 x^6-24451375 x^7+202126 z,-46518+1236672 x+1655875 x^2+15987939 x^3+2473650 x^4+39330525 x^5-23125875 x^6+31086250 x^7+202126 y}
```

This is an equivalent input:

```wolfram
GroebnerBasis[polys,{y,z}]
(* Output *)
{2+8 x+113 x^2+212 x^3+938 x^4+364 x^5+1980 x^6-1000 x^7+1625 x^8,-289274-1958323 x-7458871 x^2-14140053 x^3-20245223 x^4-19904560 x^5+4016500 x^6-24451375 x^7+202126 z,-46518+1236672 x+1655875 x^2+15987939 x^3+2473650 x^4+39330525 x^5-23125875 x^6+31086250 x^7+202126 y}
```

#### Sort

By default, [GroebnerBasis](https://reference.wolfram.com/language/ref/GroebnerBasis.html) is not allowed to reorder the variables:

```wolfram
polys={3 x^7+5x y z^2-10y^2z-6x z+y^3+w,-2x^2 z+3 x^3 y^2+y^4-12x z-8x z^2+3 y^2 z-11w x y^2,10 x^2w-7y z w^2-2 x z^4 w+4 x^2 y+3 x y^2-6y z^3-w+2,w^3-w x^2 y+x y z^2-2 w x z^2-3w-2 x y^2-3};
```

```wolfram
Timing[GroebnerBasis[polys,{w,x,y,z},MonomialOrder->DegreeReverseLexicographic];]
(* Output *)
{18.105999999999973,Null}
```

Reordering the variables may make computations faster; the Gröbner basis may be different:

```wolfram
Timing[GroebnerBasis[polys,{w,x,y,z},Sort->True,MonomialOrder->DegreeReverseLexicographic];]
(* Output *)
{9.914999999999996,Null}
```

#### Tolerance

Find an approximate GCD of a pair of univariate polynomials:

```wolfram
p1=x^14+3.00001x^10-7.99998x^7-25.00002x^6+3.00001x^13+9.00006x^9-3.00001x^5-2.00001x^8-6.00005x^4+16.00004x+2.00001;
p2=x^13-3.00003x^9-2.99999x^6+2.99999x^12-9.00006x^8-8.99997x^5-1.99998x^7+5.99999x^3+5.99994;
```

The polynomials are close to polynomials with integer coefficients:

```wolfram
{r1,r2}=Rationalize[{p1,p2},10^-4]
(* Output *)
{2+16 x-6 x^4-3 x^5-25 x^6-8 x^7-2 x^8+9 x^9+3 x^10+3 x^13+x^14,6+6 x^3-9 x^5-3 x^6-2 x^7-9 x^8-3 x^9+3 x^12+x^13}
```

```wolfram
PolynomialGCD[r1,r2]
(* Output *)
-2+3 x^5+x^6
```

With the default setting [Tolerance](https://reference.wolfram.com/language/ref/Tolerance.html)->0, the approximate GCD has a too low degree:

```wolfram
GroebnerBasis[{p1,p2},x,CoefficientDomain->InexactNumbers]
(* Output *)
{-0.6648618759149794+0.2210262129043201 x-0.07347709926728069 x^2+0.024430608194963576 x^3-0.00811507606219543 x^4+1. x^5}
```

With a higher setting of [Tolerance](https://reference.wolfram.com/language/ref/Tolerance.html), [GroebnerBasis](https://reference.wolfram.com/language/ref/GroebnerBasis.html) gives a better approximate GCD:

```wolfram
GroebnerBasis[{p1,p2},x,CoefficientDomain->InexactNumbers,Tolerance->10^(-2)]
(* Output *)
{-2.0001465223696937+3.0002401628321618 x^5+1. x^6}
```

### Applications

Solve a system of polynomial equations:

```wolfram
polys={x^2+y^2-1,y^3-2x y-3};
```

A Gröbner basis has the same set of roots as the input polynomials:

```wolfram
gb=GroebnerBasis[polys,{y,x}]
(* Output *)
{8+4 x-x^2-8 x^3+x^4+4 x^5+x^6,-1+2 x+2 x^2-2 x^3-x^4+3 y}
```

Solve the first polynomial of the Gröbner basis for its only variable `*x*`:

```wolfram
solx=Solve[gb[[1]]==0,x]
(* Output *)
{{x->Root},{x->Root},{x->Root},{x->Root},{x->Root},{x->Root}}
```

Solve the second polynomial of the Gröbner basis for the other variable `*y*`:

```wolfram
soly=Solve[gb[[2]]==0,y]
(* Output *)
{{y->(1)/(3) (1-2 x-2 x^2+2 x^3+x^4)}}
```

This method finds all common roots of `*polys*`:

```wolfram
RootReduce[(polys/.soly)/.solx]
(* Output *)
{{{0,0}},{{0,0}},{{0,0}},{{0,0}},{{0,0}},{{0,0}}}
```

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html) and [Solve](https://reference.wolfram.com/language/ref/Solve.html) use Gröbner bases to solve systems of equations:

```wolfram
Reduce[Thread[polys==0],{x,y}]
(* Output *)
(x==Root||x==Root||x==Root||x==Root||x==Root||x==Root)&&y==(1)/(3) (1-2 x-2 x^2+2 x^3+x^4)
```

Get a fuzzy solution to a system of overdetermined equations:

```wolfram
positiveSolve[polys_,vars_,tol_]:=Catch[Module[{gb,sols},gb=GroebnerBasis[polys,vars,CoefficientDomain->InexactNumbers,Tolerance->tol];
If[!ListQ[gb]||gb==={}||First[gb]==1,Throw[$Failed]];
sols=NSolve[gb,vars];
If[!ListQ[sols]||sols==={},Throw[$Failed]];
sols=vars/.sols;
sols=Select[sols,Apply[And,Map[Re[#]>=0&,#]]&];
Re[Apply[Plus,sols]/Length[sols]]]]
```

```wolfram
polys={-4+x^2-1.49071x y+y^2,-8+x^2-0.4x z+z^2,-4+t^2-0.894427t x+x^2, -4+y^2-1.49071y z+z^2,-8+t^2-0.666667t y+y^2,-4+t^2-0.894427t z+z^2};
```

This gives a low-precision approximate solution to this overdetermined set of polynomials:

```wolfram
positiveSolve[polys,{x,y,z,t},10^(-3)]
(* Output *)
{2.2361289413922347,3.0001393281755613,2.236170705591807,1.000022383436162}
```

### Properties & Relations

A Gröbner basis generates the same ideal as the input polynomials:

```wolfram
{p1,p2}={x^2+y^2+z^2-1,x-2y^3-3};
```

```wolfram
{g1,g2}=GroebnerBasis[{p1,p2},{x,y,z}]
(* Output *)
{8+y^2+12 y^3+4 y^6+z^2,-3+x-2 y^3}
```

Use [PolynomialReduce](https://reference.wolfram.com/language/ref/PolynomialReduce.html) to show that `p1 is in the ideal generated by `g1 and `g2:

```wolfram
{{q1,q2},r}=PolynomialReduce[p1,{g1,g2},{x,y,z}]
(* Output *)
{{1,3+x+2 y^3},0}
```

```wolfram
q1 g1+q2 g2==p1//Expand
(* Output *)
True
```

By Hilbert's Nullstellensatz, if the ideal is $\langle 1 \rangle$ then the polynomials have no common zero:

```wolfram
eqs={x+y^2==0,y-x+1==0,y^3-y==0};
```

```wolfram
GroebnerBasis[eqs,{x,y}]
(* Output *)
{1}
```

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html) or [Solve](https://reference.wolfram.com/language/ref/Solve.html) proves that there is no common solution:

```wolfram
Reduce[eqs,{x,y}]
(* Output *)
False
```

```wolfram
Solve[eqs,{x,y}]
(* Output *)
{}
```

Conversely, if the ideal is not $\langle 1 \rangle$, then there is at least one common zero:

```wolfram
eqs={x+y^2==0,y-x+1==0};
```

```wolfram
GroebnerBasis[eqs,{x,y}]
(* Output *)
{1+y+y^2,-1+x-y}
```

Use [FindInstance](https://reference.wolfram.com/language/ref/FindInstance.html) to find a solution instance:

```wolfram
FindInstance[eqs,{x,y}]
(* Output *)
{{x->(1)/(2) (1-ⅈ Sqrt[3]),y->-1+(1)/(2) (1-ⅈ Sqrt[3])}}
```

[GroebnerBasis](https://reference.wolfram.com/language/ref/GroebnerBasis.html) of univariate polynomials is equivalent to computing [PolynomialGCD](https://reference.wolfram.com/language/ref/PolynomialGCD.html):

```wolfram
GroebnerBasis[{(x-1)(x-2),(x-2)(x-3)},x]
(* Output *)
{-2+x}
```

```wolfram
PolynomialGCD[(x-1)(x-2),(x-2)(x-3)]
(* Output *)
-2+x
```

[GroebnerBasis](https://reference.wolfram.com/language/ref/GroebnerBasis.html) of linear polynomials is equivalent to a Gaussian elimination process:

```wolfram
GroebnerBasis[{x+2y+3z,4x+5y+6z,7x+8y+9z},{x,y,z}]
(* Output *)
{y+2 z,x-z}
```

```wolfram
RowReduce[{{1,2,3},{4,5,6},{7,8,9}}]
(* Output *)
{{1,0,-1},{0,1,2},{0,0,0}}
```

```wolfram
%.{x,y,z}
(* Output *)
{x-z,y+2 z,0}
```

[GroebnerBasis](https://reference.wolfram.com/language/ref/GroebnerBasis.html) is used to solve systems of polynomial equations:

```wolfram
eqs={x^2+y^2+z^2==1,x-2y^3==3};
```

```wolfram
GroebnerBasis[eqs,{x,y,z}]
(* Output *)
{8+y^2+12 y^3+4 y^6+z^2,-3+x-2 y^3}
```

Use [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) to directly solve the system:

```wolfram
Reduce[eqs,{z,y,x}]
(* Output *)
(y==Root[8+z^2+#1^2+12 #1^3+4 #1^6&,1]||y==Root[8+z^2+#1^2+12 #1^3+4 #1^6&,2]||y==Root[8+z^2+#1^2+12 #1^3+4 #1^6&,3]||y==Root[8+z^2+#1^2+12 #1^3+4 #1^6&,4]||y==Root[8+z^2+#1^2+12 #1^3+4 #1^6&,5]||y==Root[8+z^2+#1^2+12 #1^3+4 #1^6&,6])&&x==3+2 y^3
```

[Solve](https://reference.wolfram.com/language/ref/Solve.html) gives solutions in terms of replacement rules:

```wolfram
Solve[eqs,{x,y}]
(* Output *)
{{x->3+2 Root[8+z^2+#1^2+12 #1^3+4 #1^6&,1]^3,y->Root[8+z^2+#1^2+12 #1^3+4 #1^6&,1]},{x->3+2 Root[8+z^2+#1^2+12 #1^3+4 #1^6&,2]^3,y->Root[8+z^2+#1^2+12 #1^3+4 #1^6&,2]},{x->3+2 Root[8+z^2+#1^2+12 #1^3+4 #1^6&,3]^3,y->Root[8+z^2+#1^2+12 #1^3+4 #1^6&,3]},{x->3+2 Root[8+z^2+#1^2+12 #1^3+4 #1^6&,4]^3,y->Root[8+z^2+#1^2+12 #1^3+4 #1^6&,4]},{x->3+2 Root[8+z^2+#1^2+12 #1^3+4 #1^6&,5]^3,y->Root[8+z^2+#1^2+12 #1^3+4 #1^6&,5]},{x->3+2 Root[8+z^2+#1^2+12 #1^3+4 #1^6&,6]^3,y->Root[8+z^2+#1^2+12 #1^3+4 #1^6&,6]}}
```

Eliminate a variable from a system of polynomial equations:

```wolfram
polys={x^2+y^2+z^2-1,x y z-3};
```

```wolfram
GroebnerBasis[polys,{x,y},{z}]
(* Output *)
{9-x^2 y^2+x^4 y^2+x^2 y^4}
```

Eliminate a variable using [Resolve](https://reference.wolfram.com/language/ref/Resolve.html):

```wolfram
Resolve[Exists[z,And@@Thread[polys==0]]]
(* Output *)
9-x^2 y^2+x^4 y^2+x^2 y^4==0
```

Eliminate a variable using [Eliminate](https://reference.wolfram.com/language/ref/Eliminate.html):

```wolfram
Eliminate[Thread[polys==0],z]
(* Output *)
x^4 y^2+x^2 y^2 (-1+y^2)==-9
```

Eliminate a variable using [Resultant](https://reference.wolfram.com/language/ref/Resultant.html):

```wolfram
Resultant[polys[[1]],polys[[2]],z]
(* Output *)
9-x^2 y^2+x^4 y^2+x^2 y^4
```

Use [NonCommutativeGroebnerBasis](https://reference.wolfram.com/language/ref/NonCommutativeGroebnerBasis.html) to compute a Gröbner basis of noncommutative polynomials:

```wolfram
NonCommutativeGroebnerBasis[{2 x**y**x+3 x**x+1,3 y**x**x+2 y**x-x**y},{x,y}]
(* Output *)
{(1)/(3)+(2 y)/(27)+NonCommutativeMultiply,(y)/(3)+x**y,(y)/(3)+y**x,6 y+NonCommutativeMultiply}
```

## Tech Notes ▪Algebraic Operations on Polynomials ▪Polynomials Modulo Primes ▪Complex Polynomial Systems ▪Implementation notes: Algebra and Calculus

## Related Guides ▪Polynomial Systems ▪Theorem Proving ▪Polynomial Division ▪Computational Geometry ▪Polynomial Algebra

## Related Links [NKS|Online](http://www.wolframscience.com/nks/search/?q=GroebnerBasis) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1991 (2.0) | Updated in 1996 (3.0) ▪ 2007 (6.0)
