# SequenceHold | [SpanFromLeft]

> [SequenceHold](https://reference.wolfram.com/language/ref/SequenceHold.html) — is an attribute that specifies that [Sequence](https://reference.wolfram.com/language/ref/Sequence.html) objects appearing in the arguments of a function should not automatically be flattened out.

## Details

The attribute [SequenceHold](https://reference.wolfram.com/language/ref/SequenceHold.html) prevents [Sequence](https://reference.wolfram.com/language/ref/Sequence.html) objects from being flattened out.

## Examples

### Basic Examples

Add [SequenceHold](https://reference.wolfram.com/language/ref/SequenceHold.html) to the list of attributes of the symbol `*h*`:

```wolfram
SetAttributes[h,SequenceHold];
```

Then [Sequence](https://reference.wolfram.com/language/ref/Sequence.html) objects are not automatically spliced into expressions with head `*h*`:

```wolfram
h[a,Sequence[b,c]]
(* Output *)
h[a,Sequence[b,c]]
```

### Properties & Relations

The attribute [HoldAll](https://reference.wolfram.com/language/ref/HoldAll.html) does not prevent splicing of sequences:

```wolfram
Attributes[Hold]
(* Output *)
{HoldAll,Protected}
```

```wolfram
Hold[a,Sequence[b,c]]
(* Output *)
Hold[a,b,c]
```

[SequenceHold](https://reference.wolfram.com/language/ref/SequenceHold.html) is also necessary:

```wolfram
SetAttributes[h,{HoldAll,SequenceHold}];
```

```wolfram
h[a,Sequence[b,c]]
(* Output *)
h[a,Sequence[b,c]]
```

The attribute [HoldAllComplete](https://reference.wolfram.com/language/ref/HoldAllComplete.html) implies [SequenceHold](https://reference.wolfram.com/language/ref/SequenceHold.html):

```wolfram
Attributes[HoldComplete]
(* Output *)
{HoldAllComplete,Protected}
```

```wolfram
HoldComplete[a,Sequence[b,c]]
(* Output *)
HoldComplete[a,Sequence[b,c]]
```

Assignment operators are [SequenceHold](https://reference.wolfram.com/language/ref/SequenceHold.html), so that sequences can be returned as results:

```wolfram
splice[x_]:=Sequence[x,x,x]
```

```wolfram
{a,splice[b],c}
(* Output *)
{a,b,b,b,c}
```

Rules have the same property:

```wolfram
{f[1],g[2],h[3]}/.g[x_]:>Sequence[x,x]
(* Output *)
{f[1],2,2,h[3]}
```

## Tech Notes ▪Attributes ▪Non[Hyphen]Standard Evaluation

## Related Guides ▪Evaluation Control ▪Attributes

[Graphics] | Related Workflows
▪Handle Code Symbolically

## History Introduced in 1996 (3.0)
