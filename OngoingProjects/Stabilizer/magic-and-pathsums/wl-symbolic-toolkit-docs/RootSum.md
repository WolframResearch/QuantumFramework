# RootSum | [SpanFromLeft]

> [RootSum](https://reference.wolfram.com/language/ref/RootSum.html)[*f*,*form*] — represents the sum of `*form*[*x*]` for all `*x*` that satisfy the polynomial equation `*f*[*x*]==0

## Details

`*f*` must be a [Function](https://reference.wolfram.com/language/ref/Function.html) object such as `(#^5-2#+1)&`.

`*form*` need not correspond to a polynomial function.

[Normal](https://reference.wolfram.com/language/ref/Normal.html)[*expr*] expands [RootSum](https://reference.wolfram.com/language/ref/RootSum.html) objects into explicit sums involving [Root](https://reference.wolfram.com/language/ref/Root.html) objects.

`*f*` and `*form*` can contain symbolic parameters.

[RootSum](https://reference.wolfram.com/language/ref/RootSum.html)[*f*,*form*] is automatically simplified whenever `*form*` is a rational function.

[RootSum](https://reference.wolfram.com/language/ref/RootSum.html) is often generated in computing integrals of rational functions.

## Examples

### Basic Examples

Integrating a rational function of any order:

```wolfram
Integrate[1/(x^5+11 x+1),{x,1,3}]
(* Output *)
-RootSum[1+11 #1+#1^5&,(Log[1-#1])/(11+5 #1^4)&]+RootSum[1+11 #1+#1^5&,(Log[3-#1])/(11+5 #1^4)&]
```

Evaluate numerically:

```wolfram
N[%,50]
(* Output *)
0.0512788051842869498842709401030724212861398575534094737568852992997+0`51.44057710043037 ⅈ
```

Automatic simplification of [RootSum](https://reference.wolfram.com/language/ref/RootSum.html) objects:

```wolfram
RootSum[#^5-11 #+1&,(#^2-1)/(#^3-2#+c)&]
(* Output *)
(538-88 c+396 c^2+5 c^3-5 c^4)/(97-529 c-53 c^2+88 c^3+c^5)
```

### Scope

Compute a numerical approximation of a [RootSum](https://reference.wolfram.com/language/ref/RootSum.html):

```wolfram
N[RootSum[#^5-3#-7&,Sin]]
(* Output *)
0.29218763022095295
```

Evaluate to high precision:

```wolfram
N[RootSum[#^5-3#-7&,Sin],50]
(* Output *)
0.2921876302209531616506435821830642195268316406029177210472486143553
```

Sums over roots of polynomials with inexact number coefficients:

```wolfram
RootSum[#^5-3.2#+2.1&,f]
(* Output *)
f[-1.4670004247651094]+f[-0.14580443451492514-1.3775607790522524 ⅈ]+f[-0.14580443451492514+1.3775607790522524 ⅈ]+f[0.7144017006750915]+f[1.0442075931198682]
```

Sums of numeric functions over roots of quadratics:

```wolfram
RootSum[#^2-#+a&,Sin[#]&]
(* Output *)
Sin[(1)/(2) (1-Sqrt[1-4 a])]+Sin[(1)/(2) (1+Sqrt[1-4 a])]
```

Sums of rational functions of roots:

```wolfram
RootSum[#^5-a #+b&,(#^2-1)/(#^3-2#+c)&]
(* Output *)
(-16 a+8 a^2-a^3-10 b^2+a b^2+8 a b c+8 a c^2-4 a^2 c^2-5 b c^3+5 c^4)/(-32 b+16 a b-2 a^2 b+b^3+16 a c-8 a^2 c+a^3 c-10 b^2 c+20 b c^2+3 a b c^2-8 a c^3-c^5)
```

Sums of logarithms of linear functions over roots of polynomials with rational coefficients:

```wolfram
RootSum[#^5-2 #+3/7&,Log[2#+1]&]
(* Output *)
ⅈ π+Log[(313)/(7)]
```

Sums of numeric functions over roots of polynomials with multiple factors:

```wolfram
RootSum[(#^3-a)^2(#^4-b)^3&,5Tan[#]+7&]
(* Output *)
126+5 (2 RootSum[a-#1^3&,Tan[#1]&]+3 RootSum[-b+#1^4&,Tan[#1]&])
```

Represent a [RootSum](https://reference.wolfram.com/language/ref/RootSum.html) explicitly in terms of [Root](https://reference.wolfram.com/language/ref/Root.html) objects:

```wolfram
Normal[RootSum[#^5-3#-7&,Sin]]
(* Output *)
Sin[Root]+Sin[Root]+Sin[Root]+Sin[Root]+Sin[Root]
```

Derivatives:

```wolfram
D[RootSum[#^5+11#+1&,Exp[a #]&],a]
(* Output *)
RootSum[1+11 #1+#1^5&,ℯ^(a #1) #1&]
```

```wolfram
D[RootSum[#^5+a #+1&,Exp[#]&],a]
(* Output *)
-RootSum[1+a #1+#1^5&,(ℯ^#1 #1)/(a+5 #1^4)&]
```

Integrals:

```wolfram
Integrate[RootSum[#^5+11#+1&,Exp[-a #]&],a]
(* Output *)
-RootSum[1+11 #1+#1^5&,(ℯ^(-a #1))/(#1)&]
```

```wolfram
Integrate[RootSum[#^5+11#+1&,Sin[a #]&],{a,0,1}]
(* Output *)
-11-RootSum[1+11 #1+#1^5&,(Cos[#1])/(#1)&]
```

Limits:

```wolfram
Limit[RootSum[#^5+11#+1&,Exp[a #]&],a->1]
(* Output *)
ℯ^Root+ℯ^Root+ℯ^Root+ℯ^Root+ℯ^Root
```

```wolfram
Limit[RootSum[#^5+2#^4+11#+1&,Sin[a #]/a&],a->0]
(* Output *)
-2
```

Series:

```wolfram
Series[RootSum[#^5+2#^4+11#+1&,Sin[a #]/a&],{a,0,5}]
(* Output *)
-2+(4 a^2)/(3)+(73 a^4)/(120)+O[a]^6
```

```wolfram
Series[RootSum[#^5+a #+1&,Exp[#]&],{a,0,2}]
(* Output *)
RootSum[1+#1^5&,ℯ^#1&]-(1)/(5) RootSum[1+#1^5&,(ℯ^#1)/(#1^3)&] a-(1)/(50) RootSum[1+#1^5&,(-2 ℯ^#1+ℯ^#1 #1)/(#1^2)&] a^2+O[a]^3
```

### Applications

Integrate a rational function:

```wolfram
Integrate[(x-2)/(x^5-2x+5),x]
(* Output *)
RootSum[5-2 #1+#1^5&,(-2 Log[x-#1]+Log[x-#1] #1)/(-2+5 #1^4)&]
```

Sum a rational function:

```wolfram
Sum[1/(k^3+11k+1),{k,1,n}]
(* Output *)
(1)/(5351)(-RootSum[1+11 #1+#1^3&,484 PolyGamma[0,1-#1]-9 PolyGamma[0,1-#1] #1+66 PolyGamma[0,1-#1] #1^2&]+RootSum[1+11 #1+#1^3&,484 PolyGamma[0,1+n-#1]-9 PolyGamma[0,1+n-#1] #1+66 PolyGamma[0,1+n-#1] #1^2&])
```

Matrix exponential of any order:

```wolfram
m={{1,2,3},{1,0,0},{0,1,0}};
```

```wolfram
MatrixExp[m t]
(* Output *)
{{RootSum[-3-2 #1-#1^2+#1^3&,(ℯ^(t #1) #1^2)/(-2-2 #1+3 #1^2)&],RootSum[-3-2 #1-#1^2+#1^3&,(3 ℯ^(t #1)+2 ℯ^(t #1) #1)/(-2-2 #1+3 #1^2)&],3 RootSum[-3-2 #1-#1^2+#1^3&,(ℯ^(t #1) #1)/(-2-2 #1+3 #1^2)&]},{RootSum[-3-2 #1-#1^2+#1^3&,(ℯ^(t #1) #1)/(-2-2 #1+3 #1^2)&],RootSum[-3-2 #1-#1^2+#1^3&,(-ℯ^(t #1) #1+ℯ^(t #1) #1^2)/(-2-2 #1+3 #1^2)&],3 RootSum[-3-2 #1-#1^2+#1^3&,(ℯ^(t #1))/(-2-2 #1+3 #1^2)&]},{RootSum[-3-2 #1-#1^2+#1^3&,(ℯ^(t #1))/(-2-2 #1+3 #1^2)&],RootSum[-3-2 #1-#1^2+#1^3&,(-ℯ^(t #1)+ℯ^(t #1) #1)/(-2-2 #1+3 #1^2)&],RootSum[-3-2 #1-#1^2+#1^3&,(-2 ℯ^(t #1)-ℯ^(t #1) #1+ℯ^(t #1) #1^2)/(-2-2 #1+3 #1^2)&]}}
```

### Properties & Relations

Vieta's formulas:

```wolfram
Table[RootSum[#^4+a #^3+b #^2+c #+d&,#^k&],{k,4}]
(* Output *)
{-a,a^2-2 b,-a^3+3 a b-3 c,a^4-4 a^2 b+2 b^2+4 a c-4 d}
```

```wolfram
Table[SymmetricReduction[r^k+s^k+t^k+u^k,{r,s,t,u},{-a,b,-c,d}][[1]],{k,4}]
(* Output *)
{-a,a^2-2 b,-a^3+3 a b-3 c,a^4-4 a^2 b+2 b^2+4 a c-4 d}
```

The residue theorem:

```wolfram
NIntegrate[1/(x^6-2x+4),{x,-Infinity,Infinity}]
(* Output *)
0.7154464939404108
```

```wolfram
N[RootSum[#^6-2#+4&,If[Im[#]>0,2Pi I/(6#^5-2),0]&]]
(* Output *)
0.7154464939402045+8.326672684688674×10^-17 ⅈ
```

## Tech Notes ▪Algebraic Numbers

## Related Guides ▪Polynomial Algebra

## History Introduced in 1996 (3.0)
