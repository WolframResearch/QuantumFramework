# FiniteFieldIndex | [SpanFromLeft]

> [FiniteFieldIndex](https://reference.wolfram.com/language/ref/FiniteFieldIndex.html)[*u*] — gives the index of the [FiniteFieldElement](https://reference.wolfram.com/language/ref/FiniteFieldElement.html) object `*u*`.

## Details

[FiniteFieldIndex](https://reference.wolfram.com/language/ref/FiniteFieldIndex.html) has the [Listable](https://reference.wolfram.com/language/ref/Listable.html) attribute.

[FiniteFieldIndex](https://reference.wolfram.com/language/ref/FiniteFieldIndex.html)[*u*] is equivalent to `*u*["Index"]`.

If `*u*` is an element of [FiniteField](https://reference.wolfram.com/language/ref/FiniteField.html)[*p*,*f*,"Polynomial"], then $u=\sum_{i=0}^{d-1}u_{i}\alpha^{i}$, where `α` is the field generator and `*d*` is the degree of `*f*`. The index `*ind*` of `*u*` satisfies [IntegerDigits](https://reference.wolfram.com/language/ref/IntegerDigits.html)[*ind*,*p*,*d*]=={*u*_*d-1*,…,*u*_0}.

If `*u*` is a nonzero element of [FiniteField](https://reference.wolfram.com/language/ref/FiniteField.html)[*p*,*f*,"Exponential"], then $u=\alpha^{k}$, with $0 \leq k<p^{d}-2$, where `α` is the field generator and `*d*` is the degree of `*f*`. The index `*ind*` of `*u*` satisfies `*ind*==*k*+1 The index of the field zero is 0.

## Examples

### Basic Examples

Create a matrix of finite field elements:

```wolfram
(m=FromFiniteFieldIndex[{{11,12},{21,22}},FiniteField[5,2]])//MatrixForm
(* Output *)
![image](img/image_001.png)
```

Find the indices of matrix elements:

```wolfram
FiniteFieldIndex[m]//MatrixForm
(* Output *)
({{11, 12}, {21, 22}})
```

### Scope

Find the index of a finite field element:

```wolfram
FiniteFieldIndex[FiniteField[11,3][{1,2,3}]]
(* Output *)
386
```

Use a finite field in the exponential representation:

```wolfram
FiniteFieldIndex[FiniteField[11,3,"Exponential"][{1,2,3}]]
(* Output *)
80
```

Find indices of a vector and a matrix of finite field elements:

```wolfram
a=FiniteField[7,5][{0,1}];
FiniteFieldIndex[{a,a^2,a^3}]
(* Output *)
{7,49,343}
```

```wolfram
FiniteFieldIndex[{{a,2a,3a},{1/a,2/a,3/a}}]//MatrixForm
(* Output *)
({{7, 14, 21}, {12010, 7206, 2402}})
```

### Properties & Relations

For a single field element, [FiniteFieldIndex](https://reference.wolfram.com/language/ref/FiniteFieldIndex.html)[*u*] is equivalent to `*u*["Index"]`:

```wolfram
u=FiniteField[11,3][1234];
FiniteFieldIndex[u]===u["Index"]
(* Output *)
True
```

[FiniteFieldIndex](https://reference.wolfram.com/language/ref/FiniteFieldIndex.html) can be applied to lists of elements:

```wolfram
FiniteFieldIndex[{u,u^2,u^3}]
(* Output *)
{1234,1173,1323}
```

```wolfram
FiniteFieldIndex[Table[i u^j+k,{i,3},{j,3},{k,3}]]
(* Output *)
{{{1235,1236,1237},{1174,1175,1176},{1324,1325,1326}},{{1138,1139,1140},{884,885,886},{1195,1196,1197}},{{1041,1042,1043},{715,716,717},{1066,1056,1057}}}
```

Use [FromFiniteFieldIndex](https://reference.wolfram.com/language/ref/FromFiniteFieldIndex.html) to get field elements with specified indices:

```wolfram
FromFiniteFieldIndex[{123,234,345},FiniteField[19,2]]
(* Output *)
![image](img/image_003.png)
```

```wolfram
FiniteFieldIndex[%]
(* Output *)
{123,234,345}
```

Convert elements of a finite field to polynomials in a variable representing the field generator:

```wolfram
a=FiniteField[5,4][99];
FromFiniteField[{a,a^77,a^156},FiniteField[5,4],t]
(* Output *)
{4+4 t+3 t^2,3+2 t^2,3}
```

## Related Guides ▪Finite Fields

## History Introduced in 2024 (14.0)
