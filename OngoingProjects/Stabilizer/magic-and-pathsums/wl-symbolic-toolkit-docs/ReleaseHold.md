# ReleaseHold | [SpanFromLeft]

> [ReleaseHold](https://reference.wolfram.com/language/ref/ReleaseHold.html)[*expr*] — removes [Hold](https://reference.wolfram.com/language/ref/Hold.html), [HoldForm](https://reference.wolfram.com/language/ref/HoldForm.html), [HoldPattern](https://reference.wolfram.com/language/ref/HoldPattern.html), [HoldComplete](https://reference.wolfram.com/language/ref/HoldComplete.html) and [HoldCompleteForm](https://reference.wolfram.com/language/ref/HoldCompleteForm.html) in `*expr*`.

## Details

[ReleaseHold](https://reference.wolfram.com/language/ref/ReleaseHold.html) removes only one layer of [Hold](https://reference.wolfram.com/language/ref/Hold.html) etc.; it does not remove inner occurrences in nested [Hold](https://reference.wolfram.com/language/ref/Hold.html) etc. functions.

## Examples

### Basic Examples

```wolfram
Hold[1+1]
(* Output *)
Hold[1+1]
```

```wolfram
ReleaseHold[%]
(* Output *)
2
```

### Scope

[ReleaseHold](https://reference.wolfram.com/language/ref/ReleaseHold.html) removes all standard unevaluated containers:

```wolfram
ReleaseHold/@{Hold[1+2],HoldForm[2+3],HoldComplete[3+4],HoldCompleteForm[4+5],HoldPattern[_*_]}
(* Output *)
{3,5,7,9,_^2}
```

### Properties & Relations

[ReleaseHold](https://reference.wolfram.com/language/ref/ReleaseHold.html) removes only the outermost layer of [Hold](https://reference.wolfram.com/language/ref/Hold.html):

```wolfram
ReleaseHold[f[Hold[1+2]]]
(* Output *)
f[3]
```

```wolfram
ReleaseHold[f[Hold[1+g[Hold[2+3]]]]]
(* Output *)
f[1+g[Hold[2+3]]]
```

## Tech Notes ▪Non[Hyphen]Standard Evaluation

## Related Guides ▪Evaluation Control

## History Introduced in 1991 (2.0) | Updated in 1996 (3.0) ▪ 2025 (14.2)
