# Unevaluated | [SpanFromLeft]

> [Unevaluated](https://reference.wolfram.com/language/ref/Unevaluated.html)[*expr*] — represents the unevaluated form of `*expr*` when it appears as the argument to a function.

## Details

`*f*[Unevaluated[*expr*]]` effectively works by temporarily setting attributes so that `*f*` holds its argument unevaluated, then evaluating `*f*[*expr*]`.

## Examples

### Basic Examples

Feed an unevaluated expression to [Length](https://reference.wolfram.com/language/ref/Length.html):

```wolfram
Length[Unevaluated[5+6+7+8]]
(* Output *)
4
```

### Applications

Use [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html) and [Unevaluated](https://reference.wolfram.com/language/ref/Unevaluated.html) to suppress evaluation of symbols wherever they would occur:

```wolfram
SetAttributes[symbolLength,HoldAll];
symbolLength[s_Symbol]:=StringLength[SymbolName[Unevaluated[s]]]
```

Find the length of a symbol's name even if it has a value:

```wolfram
xyzzy=42;
```

```wolfram
symbolLength[xyzzy]
(* Output *)
5
```

### Properties & Relations

[Unevaluated](https://reference.wolfram.com/language/ref/Unevaluated.html) works only where it appears; it is not propagated:

```wolfram
f[x_]:=g[x]
```

```wolfram
f[Unevaluated[1+1]]
(* Output *)
g[2]
```

[Unevaluated](https://reference.wolfram.com/language/ref/Unevaluated.html) stops [Evaluate](https://reference.wolfram.com/language/ref/Evaluate.html):

```wolfram
Hold[Evaluate[Unevaluated[1+2]]]
(* Output *)
Hold[Unevaluated[1+2]]
```

[Unevaluated](https://reference.wolfram.com/language/ref/Unevaluated.html) inside a held function remains:

```wolfram
SetAttributes[f,HoldAll]
```

```wolfram
f[Unevaluated[1+2]]
(* Output *)
f[Unevaluated[1+2]]
```

## Tech Notes ▪Evaluation ▪Non[Hyphen]Standard Evaluation

## Related Guides ▪Evaluation Control

## History Introduced in 1991 (2.0)
