# SymmetricReduction | [SpanFromLeft]

> [SymmetricReduction](https://reference.wolfram.com/language/ref/SymmetricReduction.html)[*f*,{*x*_1,…,*x*_*n*}]  — gives a pair of polynomials $\{p,q \}$ in $x_{1},\ldots,x_{n}$ such that $f==p+q$, where $p$ is the symmetric part and $q$ is the remainder.
> [SymmetricReduction](https://reference.wolfram.com/language/ref/SymmetricReduction.html)[*f*,{*x*_1,…,*x*_*n*},{*s*_1,…,*s*_*n*}]  — gives the pair $\{p,q \}$ with the elementary symmetric polynomials in $p$ replaced by $s_{1},\ldots,s_{n}$.

## Details

If $\mathit{f}$ is a symmetric polynomial, then $p$ is the unique polynomial in elementary symmetric polynomials equal to $\mathit{f}$, and $q$ is zero.

If $\mathit{f}$ is not a symmetric polynomial, then the output $p$ is not unique, but depends on the ordering of its variables.

For a given ordering, a nonsymmetric polynomial $\mathit{f}$ can be expressed uniquely as a sum of its symmetric part $p$ and a remainder $q$ that does not contain descending monomials. A monomial $c x_{1}^{e_{1}} \ldots x_{n}^{e_{n}}$ is called descending if $e_{1}\geq \ldots \geq e_{n}$.

Changing the ordering of the variables may produce different pairs $\{p,q \}$.

[SymmetricReduction](https://reference.wolfram.com/language/ref/SymmetricReduction.html) does not check to see that $\mathit{f}$ is a polynomial, and will attempt to symmetrize the polynomial part of $\mathit{f}$.

## Examples

### Basic Examples

Write a symmetric polynomial as a sum of elementary symmetric polynomials:

```wolfram
SymmetricReduction[x^2+y^2,{x,y}]
(* Output *)
{-2 x y+(x+y)^2,0}
```

Write a nonsymmetric polynomial as a symmetric part and a remainder:

```wolfram
SymmetricReduction[x^2-y^2,{x,y}]
(* Output *)
{-2 x y+(x+y)^2,-2 y^2}
```

Name the first two elementary symmetric polynomials `s1 and `s2:

```wolfram
SymmetricReduction[x^2-y^2,{x,y},{s1,s2}]
(* Output *)
{s1^2-2 s2,-2 y^2}
```

### Scope

```wolfram
SymmetricReduction[x^10-y^10,{x,y}]
(* Output *)
{-2 x^5 y^5+25 x^4 y^4 (x+y)^2-50 x^3 y^3 (x+y)^4+35 x^2 y^2 (x+y)^6-10 x y (x+y)^8+(x+y)^10,-2 y^10}
```

```wolfram
SymmetricReduction[x^10-y^10,{x,y},{s1,s2}]
(* Output *)
{s1^10-10 s1^8 s2+35 s1^6 s2^2-50 s1^4 s2^3+25 s1^2 s2^4-2 s2^5,-2 y^10}
```

[SymmetricReduction](https://reference.wolfram.com/language/ref/SymmetricReduction.html) will reduce the polynomial part of an expression:

```wolfram
SymmetricReduction[x y+Sin[x y+x z+y z],{x,y,z},{s1,s2,s3}]
(* Output *)
{s2+Sin[x y+x z+y z],-x z-y z}
```

### Applications

Let the roots of the equation $x^{3}+a x^{2}+b x+c==0$ be $\alpha$, $\beta$, $\gamma$. The coefficients $a$, $b$, $c$ are trivially related to the symmetric polynomials of $\alpha$, $\beta$, $\gamma$:

```wolfram
poly1=x^3+a x^2+b x+c;
```

```wolfram
SolveAlways[poly1==(x-α)(x-β)(x-γ),x]
(* Output *)
{{a->-α-β-γ,b->α β+α γ+β γ,c->-α β γ}}
```

```wolfram
SymmetricPolynomial[#,{α,β,γ}]&/@Range[3]
(* Output *)
{α+β+γ,α β+α γ+β γ,α β γ}
```

A similar expression holds for the monic polynomial with roots $\alpha^{2}$, $\beta^{2}$, $\gamma^{2}$:

```wolfram
SolveAlways[x^3+a^{_}x^2+b^{_}x+c^{_}==(x-α^2)(x-β^2)(x-γ^2),x]
(* Output *)
{{a^{_}->-α^2-β^2-γ^2,b^{_}->α^2 β^2+α^2 γ^2+β^2 γ^2,c^{_}->-α^2 β^2 γ^2}}
```

Use [SymmetricReduction](https://reference.wolfram.com/language/ref/SymmetricReduction.html) to solve for `a^{_}`, `b^{_}`, `c^{_}`:

```wolfram
{a^{_},b^{_},c^{_}}=Part[SymmetricReduction[#,{α,β,γ},{-a,b,-c}]&/@{-α^2-β^2-γ^2,α^2β^2+α^2γ^2+β^2γ^2,-α^2β^2γ^2},All,1]
(* Output *)
{-a^2+2 b,b^2-2 a c,-c^2}
```

The monic polynomial with roots $\alpha^{2}$, $\beta^{2}$, $\gamma^{2}$:

```wolfram
poly2=x^3+a^{_}x^2+b^{_}x+c^{_}
(* Output *)
-c^2+(b^2-2 a c) x+(-a^2+2 b) x^2+x^3
```

Check:

```wolfram
x^2/.NSolve[poly1/.{a->3,b->5,c->7},x]//Sort
(* Output *)
{-2.8751297941627785-1.4313819935719294 ⅈ,-2.8751297941627785+1.4313819935719294 ⅈ,4.750259588325559}
```

```wolfram
x/.NSolve[poly2/.{a->3,b->5,c->7},x]//Sort
(* Output *)
{-2.875129794162779-1.4313819935719294 ⅈ,-2.875129794162779+1.4313819935719294 ⅈ,4.750259588325558}
```

Consider solving the following symmetric system of equations:

```wolfram
system={
Cos[a_1]+Cos[a_2]+Cos[a_3]-k==0,
Cos[5a_1]+Cos[5a_2]+Cos[5a_3]==0,
Cos[7a_1]+Cos[7a_2]+Cos[7a_3]==0
};
```

Use [ChebyshevT](https://reference.wolfram.com/language/ref/ChebyshevT.html) to convert to a symmetric system of polynomials:

```wolfram
eqns=system/.Cos[n_.a_i_]:>ChebyshevT[n,x_i]
(* Output *)
{-k+x_1+x_2+x_3==0,5 x_1-20 x_1^3+16 x_1^5+5 x_2-20 x_2^3+16 x_2^5+5 x_3-20 x_3^3+16 x_3^5==0,-7 x_1+56 x_1^3-112 x_1^5+64 x_1^7-7 x_2+56 x_2^3-112 x_2^5+64 x_2^7-7 x_3+56 x_3^3-112 x_3^5+64 x_3^7==0}
```

[Solve](https://reference.wolfram.com/language/ref/Solve.html) is able to solve the equations in the variables `x_1,x_2,x_3:

```wolfram
SetOptions[Solve,Cubics->False,Quartics->False];
```

```wolfram
r1=Solve[eqns,{x_1,x_2,x_3}];//Timing
(* Output *)
{0.671875,Null}
```

The leaf count of the solution is enormous:

```wolfram
LeafCount[r1]
(* Output *)
824707
```

Convert to a system of equations of symmetric polynomials $s_{1},s_{2},s_{3}$:

```wolfram
symeqns=Map[SymmetricReduction[#,{x_1,x_2,x_3},{e_1,e_2,e_3}][[1]]&,eqns,{2}]
(* Output *)
{-k+e_1==0,5 e_1-20 e_1^3+16 e_1^5+60 e_1 e_2-80 e_1^3 e_2+80 e_1 e_2^2-60 e_3+80 e_1^2 e_3-80 e_2 e_3==0,-7 e_1+56 e_1^3-112 e_1^5+64 e_1^7-168 e_1 e_2+560 e_1^3 e_2-448 e_1^5 e_2-560 e_1 e_2^2+896 e_1^3 e_2^2-448 e_1 e_2^3+168 e_3-560 e_1^2 e_3+448 e_1^4 e_3+560 e_2 e_3-1344 e_1^2 e_2 e_3+448 e_2^2 e_3+448 e_1 e_3^2==0}
```

Solve the new system of equations:

```wolfram
r2=Solve[symeqns,{e_1,e_2,e_3}];//Timing
(* Output *)
{0.03125,Null}
```

The leaf count of the symmetric solution is much smaller:

```wolfram
LeafCount[r2]
(* Output *)
2533
```

Solving for the variables `x_1,x_2,x_3 in terms of the symmetric polynomials $s_{1},s_{2},s_{3}$ is also quick:

```wolfram
Solve[Table[SymmetricPolynomial[k,{x_1,x_2,x_3}]==e_k,{k,3}],{x_1,x_2,x_3}];//Timing
(* Output *)
{0.03125,Null}
```

### Properties & Relations

The order of variables can affect the decomposition into symmetric and nonsymmetric parts:

```wolfram
SymmetricReduction[x^2-y^2,{x,y}]
(* Output *)
{-2 x y+(x+y)^2,-2 y^2}
```

```wolfram
SymmetricReduction[x^2-y^2,{y,x}]
(* Output *)
{2 x y-(x+y)^2,2 x^2}
```

Another basis for the symmetric polynomials consists of the complete symmetric polynomials. They are the sum of all monomials of a given degree, and can be defined by the generating function `[Product](https://reference.wolfram.com/language/ref/Product.html)[1-*x*_*i**t*,{*i*,*n*}]^(-1)`:

```wolfram
complete=Solve[1+h_1t+h_2t^2+h_3t^3+O[t]^4==(1)/(∏_{i=1}^{3}(1-x_it))+O[t]^4,{h_1,h_2,h_3}][[1]]
(* Output *)
{h_1->x_1+x_2+x_3,h_2->x_1^2+x_1 x_2+x_2^2+x_1 x_3+x_2 x_3+x_3^2,h_3->x_1^3+x_1^2 x_2+x_1 x_2^2+x_2^3+x_1^2 x_3+x_1 x_2 x_3+x_2^2 x_3+x_1 x_3^2+x_2 x_3^2+x_3^3}
```

A determinant formula expresses the elementary symmetric polynomials in the basis of the complete symmetric polynomials:

```wolfram
{e_1,e_2,e_3}=={h_1,Det[{{h_1, h_2}, {1, h_1}}],Det[{{h_1, h_2, h_3}, {1, h_1, h_2}, {0, 1, h_1}}]}
(* Output *)
{e_1,e_2,e_3}=={h_1,h_1^2-h_2,h_1^3-2 h_1 h_2+h_3}
```

Check:

```wolfram
%[[2]]/.complete//Expand
(* Output *)
{x_1+x_2+x_3,x_1 x_2+x_1 x_3+x_2 x_3,x_1 x_2 x_3}
```

```wolfram
SymmetricPolynomial[#,{x_1,x_2,x_3}]&/@Range[3]
(* Output *)
{x_1+x_2+x_3,x_1 x_2+x_1 x_3+x_2 x_3,x_1 x_2 x_3}
```

Any symmetric polynomial can also be expressed in terms of the complete symmetric polynomials:

```wolfram
SymmetricReduction[x y z(x+y)(x+z)(y+z),{x,y,z},{e_1,e_2,e_3}]/.{e_1->h_1,e_2->Det[{{h_1, h_2}, {1, h_1}}],e_3->Det[{{h_1, h_2, h_3}, {1, h_1, h_2}, {0, 1, h_1}}]}//Expand
(* Output *)
{h_1^4 h_2-2 h_1^2 h_2^2-h_1^3 h_3+3 h_1 h_2 h_3-h_3^2,0}
```

## Tech Notes ▪Symmetric Polynomials

## Related Guides ▪Polynomial Factoring & Decomposition ▪Polynomial Algebra

## History Introduced in 2007 (6.0)
