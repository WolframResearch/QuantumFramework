# HoldPattern | [SpanFromLeft]

> [HoldPattern](https://reference.wolfram.com/language/ref/HoldPattern.html)[*expr*] — is equivalent to `*expr*` for pattern matching, but maintains `*expr*` in an unevaluated form.

## Details

[HoldPattern](https://reference.wolfram.com/language/ref/HoldPattern.html) has attribute [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html).

The left[Hyphen]hand sides of rules are usually evaluated, as are parts of the left[Hyphen]hand sides of assignments. You can use [HoldPattern](https://reference.wolfram.com/language/ref/HoldPattern.html) to stop any part from being evaluated.

## Examples

### Basic Examples

Set up a pattern whose left-hand side is kept unevaluated:

```wolfram
HoldPattern[_+_]->0
(* Output *)
HoldPattern[_+_]->0
```

Use the pattern:

```wolfram
a+b/.%
(* Output *)
0
```

Make a definition without the argument of `f` being evaluated:

```wolfram
f[HoldPattern[_+_]]:=0
```

```wolfram
f[a+b]
(* Output *)
0
```

[Log](https://reference.wolfram.com/language/ref/Log.html)[*a*,*b*] autoevaluates to [Log](https://reference.wolfram.com/language/ref/Log.html)[*b*]/[Log](https://reference.wolfram.com/language/ref/Log.html)[*a*], so there is a match:

```wolfram
MatchQ[Log[a,b],HoldPattern[Log[_]/Log[_]]]
(* Output *)
True
```

[Cases](https://reference.wolfram.com/language/ref/Cases.html)[*e*,*patt*->*rhs*] finds elements that match `*patt*`; use [HoldPattern](https://reference.wolfram.com/language/ref/HoldPattern.html) to find rules:

```wolfram
Cases[{a->b,c->d},HoldPattern[a->_]]
(* Output *)
{a->b}
```

## Tech Notes ▪Patterns and Transformation Rules ▪Evaluation ▪Evaluation in Patterns, Rules and Definitions

## Related Guides ▪Evaluation Control ▪Patterns

## History Introduced in 1996 (3.0)
