# ToFiniteField | [SpanFromLeft]

> [ToFiniteField](https://reference.wolfram.com/language/ref/ToFiniteField.html)[*k*,*ff*] — converts the integer `*k*` to an element of the prime subfield of the finite field `*ff*`.
> [ToFiniteField](https://reference.wolfram.com/language/ref/ToFiniteField.html)[*expr*,*ff*]  — converts the coefficients of the rational expression `*expr*` to elements of the finite field `*ff*`.
> [ToFiniteField](https://reference.wolfram.com/language/ref/ToFiniteField.html)[*expr**,**ff*,*t*]  — converts the coefficients of the rational expression `*expr*` to elements of the finite field `*ff*`, with `*t*` representing the field generator.

## Details

[ToFiniteField](https://reference.wolfram.com/language/ref/ToFiniteField.html) replaces integers `*k*` with elements `*ff*[{*k*}]` of the prime subfield of `*ff*` and replaces `*t*` with the field generator `*ff*[{0,1}]`.

[ToFiniteField](https://reference.wolfram.com/language/ref/ToFiniteField.html) goes inside [List](https://reference.wolfram.com/language/ref/List.html), [Plus](https://reference.wolfram.com/language/ref/Plus.html), [Times](https://reference.wolfram.com/language/ref/Times.html) and integer [Power](https://reference.wolfram.com/language/ref/Power.html) in `*expr*`.

## Examples

### Basic Examples

Convert an integer to an element of the prime subfield of a finite field:

```wolfram
ToFiniteField[21,FiniteField[11,2]]
(* Output *)
![image](img/image_001.png)
```

Use `*t*` to represent the field generator:

```wolfram
ToFiniteField[2+3t+4t^2,FiniteField[19,3],t]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[19, 17, +, 4, #, +, #, ^, 3, &, Polynomial], 234], index -> 1503, shortIndex -> 1503, indexShortened -> True, characteristic -> 19, shortCharacteristic -> 19, extensionDegree -> 3, field -> FiniteField[...], fieldDisplayed -> False|>
```

Convert the coefficients of a rational expression to elements in the prime subfield of a finite field:

```wolfram
ToFiniteField[x+9x/(2y+3z),FiniteField[5,2]]
(* Output *)
![image](img/image_003.png)
```

Use `*t*` to represent the field generator:

```wolfram
ToFiniteField[(2+3t)x+(4+t^2)y,FiniteField[7,3],t]
(* Output *)
![image](img/image_005.png)
```

### Scope

Convert integers and rational numbers to elements of the prime subfield of a finite field:

```wolfram
ToFiniteField[{123,456/789},FiniteField[7,2]]
(* Output *)
![image](img/image_007.png)
```

```wolfram
PolynomialMod[{123,456/789},7]
(* Output *)
{4,3}
```

Convert a polynomial in `*t*` to a a polynomial in the field generator:

```wolfram
ToFiniteField[12+34t+56t^2+78t^3,FiniteField[23,4],t]
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[23, 5, +, 19, #, +, 3, #, ^, 2, +, #, ^, 4, &, Polynomial], 1211109], index -> 115058, shortIndex -> 115058, indexShortened -> True, characteristic -> 23, shortCharacteristic -> 23, extensionDegree -> 4, field -> FiniteField[...], fieldDisplayed -> False|>
```

```wolfram
θ=FiniteField[23,4][{0,1}]
(* Output *)
![image](img/image_009.png)
```

```wolfram
12+34t+56t^2+78t^3/.t->θ
(* Output *)
<|interpretation -> FiniteFieldElement[FiniteField[23, 5, +, 19, #, +, 3, #, ^, 2, +, #, ^, 4, &, Polynomial], 1211109], index -> 115058, shortIndex -> 115058, indexShortened -> True, characteristic -> 23, shortCharacteristic -> 23, extensionDegree -> 4, field -> FiniteField[...], fieldDisplayed -> False|>
```

Convert the coefficients of a polynomial to elements in the prime subfield of a finite field:

```wolfram
ToFiniteField[123+456x+789x^2+x^3,FiniteField[29,3]]
(* Output *)
![image](img/image_011.png)
```

Convert the coefficients of a rational function, with `*t*` used to represent the field generator:

```wolfram
ToFiniteField[t((t^2+3)x+(5t+7)y)/(4t^3x^2+(3t^2+1)y^2),FiniteField[17,4],t]
(* Output *)
![image](img/image_013.png)
```

```wolfram
θ=FiniteField[17,4][{0,1}]
(* Output *)
![image](img/image_015.png)
```

```wolfram
t((t^2+3)x+(5t+7)y)/(4t^3x^2+(3t^2+1)y^2)/.t->θ
(* Output *)
![image](img/image_017.png)
```

### Properties & Relations

[FromFiniteField](https://reference.wolfram.com/language/ref/FromFiniteField.html) converts finite field elements to polynomials in the field generator:

```wolfram
ToFiniteField[(2+3t)x+(4+t^2)y,FiniteField[7,3],t]
(* Output *)
![image](img/image_019.png)
```

```wolfram
FromFiniteField[%,FiniteField[7,3],t]
(* Output *)
(2+3 t) x+(4+t^2) y
```

[FromFiniteFieldIndex](https://reference.wolfram.com/language/ref/FromFiniteFieldIndex.html) gives finite field elements with specified indices:

```wolfram
FromFiniteFieldIndex[{123,456,789},FiniteField[5,7]]
(* Output *)
![image](img/image_021.png)
```

[ToFiniteField](https://reference.wolfram.com/language/ref/ToFiniteField.html) converts integers to elements of the prime subfield:

```wolfram
ToFiniteField[{123,456,789},FiniteField[5,7]]
(* Output *)
![image](img/image_023.png)
```

## Related Guides ▪Finite Fields

## History Introduced in 2024 (14.0)
