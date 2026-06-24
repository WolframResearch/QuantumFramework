# MultiplicativeOrder | [SpanFromLeft]

> [MultiplicativeOrder](https://reference.wolfram.com/language/ref/MultiplicativeOrder.html)[*k*,*n*] — gives the multiplicative order of `*k*` modulo `*n*`, defined as the smallest integer $m$ such that $k^{m}\equiv 1 mod n$.
> [MultiplicativeOrder](https://reference.wolfram.com/language/ref/MultiplicativeOrder.html)[*k*,*n*,{*r*_1,*r*_2,…}] — gives the generalized multiplicative order of `*k*` modulo `*n*`, defined as the smallest integer $m$ such that $k^{m}\equiv r_{i}mod n$ for some $i$.

## Details

[MultiplicativeOrder](https://reference.wolfram.com/language/ref/MultiplicativeOrder.html) is also known as modulo order or haupt[Hyphen]exponent.

Integer mathematical function, suitable for both symbolic and numerical manipulation.

Typically used in modular arithmetic and cryptography.

[MultiplicativeOrder](https://reference.wolfram.com/language/ref/MultiplicativeOrder.html)[*k*,*n*] gives the smallest positive integer `*m*` such that the remainder when dividing `*k*^(*m*)` by `*n*` is equal to 1.

[MultiplicativeOrder](https://reference.wolfram.com/language/ref/MultiplicativeOrder.html) returns unevaluated if there is no integer $m$ satisfying the necessary conditions.

$$
\text{[Graphics]}
$$

For a [FiniteFieldElement](https://reference.wolfram.com/language/ref/FiniteFieldElement.html) object `*a*`, [MultiplicativeOrder](https://reference.wolfram.com/language/ref/MultiplicativeOrder.html)[*a*] gives the multiplicative order of `*a*`, defined as the smallest positive integer `*m*` such that $a^{m}$ is the multiplicative identity of the finite field.

## Examples

### Basic Examples

The multiplicative order of 5 modulo 8:

```wolfram
MultiplicativeOrder[5,8]
(* Output *)
2
```

Plot the sequence with a fixed modulus:

```wolfram
DiscretePlot[MultiplicativeOrder[k,7],{k,1,50}]
```

*([Graphics])*

Plot the sequence, varying the modulus:

```wolfram
DiscretePlot[MultiplicativeOrder[7,n],{n,1,50}]
```

*([Graphics])*

### Scope

#### Numerical Evaluation

Compute using integers:

```wolfram
MultiplicativeOrder[5,7]
(* Output *)
6
```

```wolfram
MultiplicativeOrder[-5,7]
(* Output *)
3
```

Generalized multiplicative order:

```wolfram
MultiplicativeOrder[5,7,{3,11}]
(* Output *)
2
```

Compute using large numbers:

```wolfram
MultiplicativeOrder[10^10000,7919]
(* Output *)
3959
```

Multiplicative order of finite field elements:

```wolfram
ℱ=FiniteField[1009,7];
a=ℱ[1234]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[1009, 998, +, #, +, #, ^, 7, &, Polynomial], 2251], index -> 1234, shortIndex -> 1234, indexShortened -> True, characteristic -> 1009, shortCharacteristic -> 1009, extensionDegree -> 7, field -> FiniteField[...], fieldDisplayed -> False|>
```

```wolfram
MultiplicativeOrder[a]
(* Output *)
1064726745878753869968
```

```wolfram
a^%
(* Output *)
![image](img/image_001.png)
```

[](https://reference.wolfram.com/language/ref/.html) formatting:

```wolfram
MultiplicativeOrder[k,n]//
(* Output *)
ord_n(k)
```

#### Symbolic Manipulation

Use [Solve](https://reference.wolfram.com/language/ref/Solve.html) to find solutions of equations:

```wolfram
Solve[MultiplicativeOrder[17,n]==2&&0<n<15,n,Integers]
(* Output *)
{{n->3},{n->6},{n->9},{n->12}}
```

Use [FindInstance](https://reference.wolfram.com/language/ref/FindInstance.html) to find solutions:

```wolfram
FindInstance[MultiplicativeOrder[17,n]==MultiplicativeOrder[5,n],n,Integers]
(* Output *)
{{n->171}}
```

### Applications

#### Basic Applications

Find all primitive roots modulo 43:

```wolfram
Select[Range[43],MultiplicativeOrder[#,43]==EulerPhi[43]&]
(* Output *)
{3,5,12,18,19,20,26,28,29,30,33,34}
```

A rational number $\frac{r}{s}$ has a digit cycle of length $s-1$ if $s$ is prime and 10 is a primitive root for $s$:

```wolfram
Select[Range[2,500],Length[RealDigits[1/#,10][[1,-1]]]==#-1&]==Select[Range[2,500],MultiplicativeOrder[10,#]==#-1&]
```

```wolfram
PrimitiveRootList[11]
(* Output *)
{2,6,7,8}
```

Compute [MultiplicativeOrder](https://reference.wolfram.com/language/ref/MultiplicativeOrder.html) using [NestWhileList](https://reference.wolfram.com/language/ref/NestWhileList.html):

```wolfram
ord[n_,m_]:=NestWhileList[Mod[n #,m]&,Mod[n,m],Mod[#,m]!=1&]//Length
```

```wolfram
MultiplicativeOrder[4,13]==ord[4,13]
(* Output *)
True
```

Count number of possible multiplicative orders modulo a given prime number:

```wolfram
l1=Table[Length[DeleteDuplicates[Table[MultiplicativeOrder[i,Prime[n]],{i,1,Prime[n]-1}]]],{n,1,50}];
```

The number of divisors of $p-1$ where $p$ is prime:

```wolfram
l2=Table[Length[Divisors[Prime[n]-1]],{n,1,50}];
```

These are in fact the same list:

```wolfram
l1==l2
(* Output *)
True
```

#### Number Theory

The repetition period in Rule $90$ for odd $n$ divides `*q*[*n*]`:

```wolfram
q[n_]=2^MultiplicativeOrder[2,n,{1,-1}]-1;
```

```wolfram
Table[q[n],{n,3,50,2}]
(* Output *)
{1,3,7,7,31,63,15,15,511,63,2047,1023,511,16383,31,31,4095,262143,4095,1023,127,4095,8388607,2097151}
```

The digits of $\frac{1}{13}$ in base $4$ repeat with period $6$:

```wolfram
With[{b=4,n=13},MultiplicativeOrder[b,FixedPoint[#/GCD[#,b]&,n]]]
(* Output *)
6
```

```wolfram
With[{b=4,n=13},RealDigits[1/n,b,20]]
(* Output *)
{{1,0,3,2,3,0,1,0,3,2,3,0,1,0,3,2,3,0,1,0},-1}
```

The function `digitCycleLength` gives the digit period for any rational number $r$ in base $b$:

```wolfram
digitCycleLength[r_Rational,b_Integer?Positive]:=
MultiplicativeOrder[b,FixedPoint[Quotient[#,GCD[#,b]]&,Denominator[r]]]
```

This shows that the decimal representation of $\frac{123}{999}$ in base 10 repeats every 3 digits:

```wolfram
digitCycleLength[(123)/(999),10]
(* Output *)
3
```

```wolfram
N[(123)/(999),18]
(* Output *)
0.12312312312312312312312312312
```

Build an RSA-like toy encryption scheme:

```wolfram
{p,q}={47,59};
n=p q;
φ=EulerPhi[n];
d=NestWhile[#1+1&,Round[n/3],GCD[φ,#1]=!=1&];
e=17;
PlainText=1504;
CypherText=PowerMod[PlainText,e,n];
```

Perform a cycling attack. One of the outputs will be the plaintext:

```wolfram
Cy[1]=PowerMod[CypherText,e,n];
Cy[j_]:=PowerMod[Cy[j-1],e,n];
Table[Cy[i],{i,Divisors[MultiplicativeOrder[e,φ]]}]
(* Output *)
{470,2209,2444,1504,2209,2444}
```

### Properties & Relations

The multiplicative order of a primitive root modulo `*n*` is [EulerPhi](https://reference.wolfram.com/language/ref/EulerPhi.html)[*n*]:

```wolfram
MultiplicativeOrder[PrimitiveRoot[109],109]==EulerPhi[109]
(* Output *)
True
```

[EulerPhi](https://reference.wolfram.com/language/ref/EulerPhi.html) divides [MultiplicativeOrder](https://reference.wolfram.com/language/ref/MultiplicativeOrder.html):

```wolfram
Divisible[EulerPhi[12],MultiplicativeOrder[5,12]]
(* Output *)
True
```

The result is always positive:

```wolfram
MultiplicativeOrder[5,3]
(* Output *)
2
```

```wolfram
MultiplicativeOrder[-5,3]
(* Output *)
1
```

Find the smallest integer such that $5^{m}$≅ 2, 3, or 4 mod 7:

```wolfram
MultiplicativeOrder[5,7,{2,3,4}]
(* Output *)
2
```

Find which of the remainders satisfies $5^{2}\equiv r_{i} mod 7$:

```wolfram
PowerMod[5,2,7]
(* Output *)
4
```

Solve the discrete log problem with $5^{m}\equiv 4 mod 7$:

```wolfram
MultiplicativeOrder[5,7,{4}]
(* Output *)
2
```

### Possible Issues

For nonzero integers `*k*` and `*n*`, [MultiplicativeOrder](https://reference.wolfram.com/language/ref/MultiplicativeOrder.html)[*k*,*n*] exists if and only if `*k*` and `*n*` are coprime:

```wolfram
CoprimeQ[10,21]
(* Output *)
True
```

```wolfram
MultiplicativeOrder[10,21]
(* Output *)
6
```

```wolfram
MultiplicativeOrder[21,10]
(* Output *)
1
```

However, 10 and 22 are not coprime:

```wolfram
CoprimeQ[10,22]
(* Output *)
False
```

```wolfram
MultiplicativeOrder[10,22]
(* Output *)
MultiplicativeOrder[10,22]
```

```wolfram
MultiplicativeOrder[22,10]
(* Output *)
MultiplicativeOrder[22,10]
```

### Interactive Examples

[MultiplicativeOrder](https://reference.wolfram.com/language/ref/MultiplicativeOrder.html) of each integer below a given prime number:

```wolfram
Manipulate[
DiscretePlot[MultiplicativeOrder[k,Prime[n]],{k,1,Prime[n]-1}],
{n,10,30,1}]
```

### Neat Examples

Visualize when a number has multiplicative order modulo 12:

```wolfram
ArrayMesh[Boole[Table[NumberQ[MultiplicativeOrder[a+b^2+c^3,12]],{a,10},{b,10},{c,10}]]]
```

*([Graphics3D])*

Ulam spiral of [MultiplicativeOrder](https://reference.wolfram.com/language/ref/MultiplicativeOrder.html):

```wolfram
ulam[n_]:=Partition[Permute[Range[n^2],Accumulate[Take[Flatten[{{n^2+1}/2,Table
	[(-1)^j i,{j,n},{i,{-1,n}},{j}]}],n^2]]],n]
```

```wolfram
ArrayPlot[MapAt[MultiplicativeOrder[#,61]&,ulam[101],{All,All}],ColorFunction->"Rainbow"]
```

![image](img/image_003.png)

## Tech Notes ▪Integer and Number Theoretic Functions

## Related Guides ▪Cryptographic Number Theory ▪Number Theoretic Functions ▪Number Theory ▪Diophantine Equations ▪Mathematical Functions ▪Cryptography ▪Finite Fields

## Related Links [NKS|Online](http://www.wolframscience.com/nks/search/?q=MultiplicativeOrder) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1999 (4.0) | Updated in 2023 (13.3)
