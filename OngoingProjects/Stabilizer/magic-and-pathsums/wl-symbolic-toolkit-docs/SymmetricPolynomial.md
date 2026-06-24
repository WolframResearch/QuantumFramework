# SymmetricPolynomial | [SpanFromLeft]

> [SymmetricPolynomial](https://reference.wolfram.com/language/ref/SymmetricPolynomial.html)[*k*,{*x*_1,…,*x*_*n*}]  — gives the `*k*`$^{th}$ elementary symmetric polynomial in the variables `*x*_1,…,*x*_*n*`.

## Details

A symmetric polynomial of `*n*` variables `{*x*_1,…,*x*_*n*}` is invariant under any permutation of its variables. The `*k*`$^{th}$ elementary symmetric polynomial is the sum of all square-free monomials of degree `*k*`.

The degree `*k*` must satisfy `0<=*k*<=*n*`.

The elementary symmetric polynomials form a basis for the symmetric polynomials.

Expressing a general symmetric polynomial in terms of elementary symmetric polynomials is accomplished by using [SymmetricReduction](https://reference.wolfram.com/language/ref/SymmetricReduction.html).

## Examples

### Basic Examples

The elementary symmetric polynomial of degree 3 in variables `x_1,x_2,x_3,x_4:

```wolfram
SymmetricPolynomial[3,{x_1,x_2,x_3,x_4}]
(* Output *)
x_1 x_2 x_3+x_1 x_2 x_4+x_1 x_3 x_4+x_2 x_3 x_4
```

### Scope

The zeroth elementary symmetric polynomial is defined to be 1:

```wolfram
SymmetricPolynomial[0,{x_1,x_2,x_3,x_4}]
(* Output *)
1
```

### Applications

The 2×3 matrices with entries 0 or 1:

```wolfram
tabs=Table[Partition[IntegerDigits[n,2,6],3],{n,64}];
```

Select matrices whose column sums are `1,1,1 and whose row sums are `2,1:

```wolfram
MatrixForm/@Cases[tabs,a_/;Total[a,{1}]=={1,1,1}&&Total[a,{2}]=={2,1}]
(* Output *)
{({{0, 1, 1}, {1, 0, 0}}),({{1, 0, 1}, {0, 1, 0}}),({{1, 1, 0}, {0, 0, 1}})}
```

You can also count how many such matrices there are by using [SymmetricPolynomial](https://reference.wolfram.com/language/ref/SymmetricPolynomial.html). The generating function of 2×3 matrices whose row sums are `2,1 is given by:

```wolfram
genfun=SymmetricPolynomial[2,{x_1,x_2,x_3}]SymmetricPolynomial[1,{x_1,x_2,x_3}]
(* Output *)
(x_1+x_2+x_3) (x_1 x_2+x_1 x_3+x_2 x_3)
```

The coefficient of `x_1^(1)x_2^(1)x_3^(1)` counts how many of these matrices have column sums `1,1,1:

```wolfram
Coefficient[genfun,{x_1^1x_2^1x_3^1}]
(* Output *)
{3}
```

### Properties & Relations

The `*k*`$^{th}$ elementary symmetric polynomial is the sum of all monomials constructed from `*k*`-subsets of the variables:

```wolfram
Plus@@Subsets[Times[x_1,x_2,x_3,x_4],{2}]
(* Output *)
x_1 x_2+x_1 x_3+x_2 x_3+x_1 x_4+x_2 x_4+x_3 x_4
```

```wolfram
SymmetricPolynomial[2,{x_1,x_2,x_3,x_4}]
(* Output *)
x_1 x_2+x_1 x_3+x_2 x_3+x_1 x_4+x_2 x_4+x_3 x_4
```

The generating function for the symmetric polynomials in $n$ variables is given by $\prod_{i=1}^{n}(x_{i}+1)$:

```wolfram
CoefficientList[∏_{i=1}^{4}(1+x_it),t]//TableForm
(* Output *)
{{1}, {x_1+x_2+x_3+x_4}, {x_1 x_2+x_1 x_3+x_2 x_3+x_1 x_4+x_2 x_4+x_3 x_4}, {x_1 x_2 x_3+x_1 x_2 x_4+x_1 x_3 x_4+x_2 x_3 x_4}, {x_1 x_2 x_3 x_4}}
```

Check:

```wolfram
TableForm[Table[SymmetricPolynomial[i,{x_1,x_2,x_3,x_4}],{i,0,4}]]
(* Output *)
{{1}, {x_1+x_2+x_3+x_4}, {x_1 x_2+x_1 x_3+x_2 x_3+x_1 x_4+x_2 x_4+x_3 x_4}, {x_1 x_2 x_3+x_1 x_2 x_4+x_1 x_3 x_4+x_2 x_3 x_4}, {x_1 x_2 x_3 x_4}}
```

The monic polynomial with roots $\alpha_{i}$ has coefficients that are elementary symmetric polynomials of the $\alpha_{i}$:

```wolfram
SolveAlways[x^4+a x^3+b x^2+c x+d==(x-α_1)(x-α_2)(x-α_3)(x-α_4),x][[1]]//TableForm
(* Output *)
{{a->-α_1-α_2-α_3-α_4}, {b->α_1 α_2+α_1 α_3+α_2 α_3+α_1 α_4+α_2 α_4+α_3 α_4}, {c->-α_1 α_2 α_3-α_1 α_2 α_4-α_1 α_3 α_4-α_2 α_3 α_4}, {d->α_1 α_2 α_3 α_4}}
```

```wolfram
SymmetricPolynomial[#,{-α_1,-α_2,-α_3,-α_4}]&/@Range[4]//TableForm
(* Output *)
{{-α_1-α_2-α_3-α_4}, {α_1 α_2+α_1 α_3+α_2 α_3+α_1 α_4+α_2 α_4+α_3 α_4}, {-α_1 α_2 α_3-α_1 α_2 α_4-α_1 α_3 α_4-α_2 α_3 α_4}, {α_1 α_2 α_3 α_4}}
```

The elementary symmetric polynomials `*e*_*k*=SymmetricPolynomial[*k*,{*x*_1,…,*x*_*n*}]` are related to the power sum polynomials $s_{p}=\sum_{k=1}^{n}x_{k}^{p}$ through the Newton-Girard formulas [[MathWorld](https://mathworld.wolfram.com/Newton-GirardFormulas.html)]. Generate all the Newton-Girard formulas for $n=4$:

```wolfram
n=4;
ToeplitzMatrix[Table[If[k==0,1,(-1)^ke_k],{k,0,n-1}],UnitVector[n,1]].Table[s_k,{k,n}]-Table[(-1)^(k+1)k e_k,{k,n}]==0//Thread
(* Output *)
{-e_1+s_1==0,2 e_2-e_1 s_1+s_2==0,-3 e_3+e_2 s_1-e_1 s_2+s_3==0,4 e_4-e_3 s_1+e_2 s_2-e_1 s_3+s_4==0}
```

Verify them:

```wolfram
%/.
{s_j_:>Sum[x_k^j,{k,n}],e_j_:>SymmetricPolynomial[j,{x_1,x_2,x_3,x_4}]}//Simplify
(* Output *)
{True,True,True,True}
```

The elementary symmetric polynomials can be defined in terms of the generalized Bell polynomial [BellY](https://reference.wolfram.com/language/ref/BellY.html):

```wolfram
ee[n_,vars_]:=BellY[Table[{(1)/(n!),(-1)^(k-1)(k-1)!Total[vars^k]},{k,n}]]
```

Verify for the case of five variables:

```wolfram
Table[ee[k,{x,y,z,u,v}]==SymmetricPolynomial[k,{x,y,z,u,v}]//Simplify,{k,5}]
(* Output *)
{True,True,True,True,True}
```

### Neat Examples

Find integers $a_{1},a_{2},a_{3}$ such that the roots of $x^{3}+a_{1} x^{2}+a_{2} x+a_{3}=0$ are $a_{1},a_{2},a_{3}$:

```wolfram
x^3+a_1x^2+a_2x+a_3/.{ToRules[Reduce[Table[SymmetricPolynomial[i,{a_1,a_2,a_3}]==a_i(-1)^i,{i,3}],{a_1,a_2,a_3},
Integers]]}
(* Output *)
{x^3,-2 x+x^2+x^3,-1-x+x^2+x^3}
```

Check:

```wolfram
Solve[#==0,x]&/@%
(* Output *)
{{{x->0},{x->0},{x->0}},{{x->-2},{x->0},{x->1}},{{x->-1},{x->-1},{x->1}}}
```

## Tech Notes ▪Symmetric Polynomials

## Related Guides ▪Polynomial Algebra

## History Introduced in 2007 (6.0)
