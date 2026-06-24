# FromFiniteField | [SpanFromLeft]

> [FromFiniteField](https://reference.wolfram.com/language/ref/FromFiniteField.html)[*a*,*ff*] — converts the element `*a*` of the prime subfield of the finite field `*ff*` to an integer.
> [FromFiniteField](https://reference.wolfram.com/language/ref/FromFiniteField.html)[*expr*,*ff*,*t*] — converts the elements of the finite field `*ff*` in the coefficients of the rational expression `*expr*` to polynomials in `*t*`, where `*t*` represents the field generator.

## Details

[FromFiniteField](https://reference.wolfram.com/language/ref/FromFiniteField.html) replaces elements $u=\sum_{i=0}^{d-1}u_{i}\alpha^{i}$ of `*ff*`, where `α=*ff*[{0,1}]` is the field generator, with polynomials $\sum_{i=0}^{d-1}u_{i}t^{i}$.

[FromFiniteField](https://reference.wolfram.com/language/ref/FromFiniteField.html) goes inside [List](https://reference.wolfram.com/language/ref/List.html), [Plus](https://reference.wolfram.com/language/ref/Plus.html), [Times](https://reference.wolfram.com/language/ref/Times.html) and integer [Power](https://reference.wolfram.com/language/ref/Power.html) in `*expr*`.

## Examples

### Basic Examples

Convert an element of the prime subfield of a finite field to an integer:

```wolfram
ff=FiniteField[5,2];
a=ff[3]
(* Output *)
![image](img/image_001.png)
```

```wolfram
FromFiniteField[a,ff]
(* Output *)
3
```

Convert an element of a finite field to a polynomial in a variable representing the field generator:

```wolfram
ff=FiniteField[7,3];
a=ff[123]
(* Output *)
![image](img/image_003.png)
```

```wolfram
FromFiniteField[a,ff,t]
(* Output *)
4+3 t+2 t^2
```

Convert prime field coefficients in a rational expression to integers:

```wolfram
ff=FiniteField[19];
expr=ff[7]x+ff[12]y/(ff[17]+ff[15]x)
(* Output *)
![image](img/image_005.png)
```

```wolfram
FromFiniteField[expr,ff]
(* Output *)
7 x+(12 y)/(17+15 x)
```

Convert finite field coefficients in a rational expression to polynomials in the field generator:

```wolfram
ff=FiniteField[11,3];
expr=ff[123]x+ff[234]/(ff[345]x+ff[456]y)
(* Output *)
![image](img/image_007.png)
```

```wolfram
FromFiniteField[expr,ff,t]
(* Output *)
(2+t^2) x+(3+10 t+t^2)/((4+9 t+2 t^2) x+(5+8 t+3 t^2) y)
```

### Scope

Convert an element of the prime subfield of a finite field to an integer:

```wolfram
ff=FiniteField[19,4];
a=ff[16]
(* Output *)
![image](img/image_009.png)
```

```wolfram
FromFiniteField[a,ff]
(* Output *)
16
```

`*b*` is not an element of the prime subfield:

```wolfram
b=ff[123]
(* Output *)
![image](img/image_011.png)
```

```wolfram
FromFiniteField[b,ff]
(* Output *)
FromFiniteField
(* Output *)
FromFiniteField[<|interpretation -> FiniteFieldElement[FiniteField[19, 2, +, 11, #, +, 2, #, ^, 2, +, #, ^, 4, &, Polynomial], 96], index -> 123, shortIndex -> 123, indexShortened -> True, characteristic -> 19, shortCharacteristic -> 19, extensionDegree -> 4, field -> FiniteField[...], fieldDisplayed -> False|>,FiniteField[...]]
```

Convert an element of a finite field to a polynomial in a variable representing the field generator:

```wolfram
ff=FiniteField[5,5];
a=ff[1234]
(* Output *)
![image](img/image_013.png)
```

```wolfram
FromFiniteField[a,ff,t]
(* Output *)
4+t+4 t^2+4 t^3+t^4
```

```wolfram
%/.t->ff[{0,1}]
(* Output *)
![image](img/image_015.png)
```

Convert finite field coefficients in a rational expression to polynomials in the field generator:

```wolfram
ff=FiniteField[29,3];
expr=ff[123](ff[234]x^2+ff[345]y^2)/(ff[456]+ff[678]x^3+ff[789]y^3)
(* Output *)
![image](img/image_017.png)
```

```wolfram
FromFiniteField[expr,ff,t]
(* Output *)
((7+4 t) ((2+8 t) x^2+(26+11 t) y^2))/(21+15 t+(11+23 t) x^3+(6+27 t) y^3)
```

```wolfram
%/.t->ff[{0,1}]
(* Output *)
![image](img/image_019.png)
```

### Properties & Relations

[ToFiniteField](https://reference.wolfram.com/language/ref/ToFiniteField.html) converts coefficients to finite field elements, with `*t*` representing the field generator:

```wolfram
ToFiniteField[(2+3t)x+(4+t^2)y,FiniteField[7,3],t]
(* Output *)
![image](img/image_021.png)
```

```wolfram
FromFiniteField[%,FiniteField[7,3],t]
(* Output *)
(2+3 t) x+(4+t^2) y
```

[FiniteFieldIndex](https://reference.wolfram.com/language/ref/FiniteFieldIndex.html) gives indices of field elements:

```wolfram
a=FiniteField[5,4][99];
FiniteFieldIndex[{a,a^77,a^156}]
(* Output *)
{99,53,3}
```

[FromFiniteField](https://reference.wolfram.com/language/ref/FromFiniteField.html) gives integers only for elements of the prime subfield:

```wolfram
FromFiniteField[{a,a^77,a^156},FiniteField[5,4],t]
(* Output *)
{4+4 t+3 t^2,3+2 t^2,3}
```

## Related Guides ▪Finite Fields

## History Introduced in 2024 (14.0)
