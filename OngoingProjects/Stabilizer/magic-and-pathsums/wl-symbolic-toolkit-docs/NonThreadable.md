[EXPERIMENTAL]

# NonThreadable

> [NonThreadable](https://reference.wolfram.com/language/ref/NonThreadable.html)  — is an attribute that can be assigned to a symbol `*f*` to indicate that `*f*` and `f[*arg*_1,*arg*_2,…]` should not combine with other list arguments in arithmetic and many other functions that work with lists.

## Details

The [NonThreadable](https://reference.wolfram.com/language/ref/NonThreadable.html) attribute is typically used for symbols that represent non-scalar quantities.

[VectorSymbol](https://reference.wolfram.com/language/ref/VectorSymbol.html), [MatrixSymbol](https://reference.wolfram.com/language/ref/MatrixSymbol.html) and [ArraySymbol](https://reference.wolfram.com/language/ref/ArraySymbol.html) have the [NonThreadable](https://reference.wolfram.com/language/ref/NonThreadable.html) attribute.

Most standard built[Hyphen]in functions that produce array results have the [NonThreadable](https://reference.wolfram.com/language/ref/NonThreadable.html) attribute.

Expressions whose heads are non-threadable objects are non-threadable.

Expressions obtained by applying [Listable](https://reference.wolfram.com/language/ref/Listable.html) functions to non-threadable objects are non-threadable.

## Examples

### Basic Examples

Give the [NonThreadable](https://reference.wolfram.com/language/ref/NonThreadable.html) attribute to symbol `*x*`:

```wolfram
SetAttributes[x,NonThreadable]
```

Arithmetic operations do not combine `*x*` with list elements:

```wolfram
x+y+{{1,2},{3,4}}
(* Output *)
x+{{1+y,2+y},{3+y,4+y}}
```

Expressions with non-threadable heads are non-threadable:

```wolfram
x[1]y[2]{1,2,3}
(* Output *)
{y[2],2 y[2],3 y[2]} x[1]
```

Expressions obtained by applying [Listable](https://reference.wolfram.com/language/ref/Listable.html) functions to non-threadable objects are non-threadable:

```wolfram
Cos[x[1,2]]+Cos[y[1,2]]+{1,2,3}
(* Output *)
Cos[x[1,2]]+{1+Cos[y[1,2]],2+Cos[y[1,2]],3+Cos[y[1,2]]}
```

### Scope

[ArraySymbol](https://reference.wolfram.com/language/ref/ArraySymbol.html) has the [NonThreadable](https://reference.wolfram.com/language/ref/NonThreadable.html) attribute:

```wolfram
Attributes[ArraySymbol]
(* Output *)
{NHoldAll,NonThreadable,Protected,ReadProtected}
```

Arithmetic operations do not combine symbolic arrays with list elements:

```wolfram
ArraySymbol[a,{p,q,r}]+{1,2,3}
(* Output *)
a+{1,2,3}
```

[SymbolicOnesArray](https://reference.wolfram.com/language/ref/SymbolicOnesArray.html) has the [NonThreadable](https://reference.wolfram.com/language/ref/NonThreadable.html) attribute:

```wolfram
Attributes[SymbolicOnesArray]
(* Output *)
{NHoldAll,NonThreadable,Protected,ReadProtected}
```

A symbolic $2 \times 2$ matrix of ones is not added to matrix elements:

```wolfram
SymbolicOnesArray[{2,2}]+{{1,2},{3,4}}
(* Output *)
{{1,2},{3,4}}+2,2
```

When the symbolic matrix is converted to an explicit matrix, the matrices are added correctly:

```wolfram
Normal[%]
(* Output *)
{{2,3},{4,5}}
```

[Inverse](https://reference.wolfram.com/language/ref/Inverse.html) has the [NonThreadable](https://reference.wolfram.com/language/ref/NonThreadable.html) attribute:

```wolfram
Attributes[Inverse]
(* Output *)
{NonThreadable,Protected,ReadProtected}
```

Arithmetic operations do not combine expressions with head [Inverse](https://reference.wolfram.com/language/ref/Inverse.html) with list elements:

```wolfram
Inverse[a]+{{1,2},{3,4}}
(* Output *)
a+{{1,2},{3,4}}
```

[Dot](https://reference.wolfram.com/language/ref/Dot.html) does not have the [NonThreadable](https://reference.wolfram.com/language/ref/NonThreadable.html) attribute:

```wolfram
Attributes[Dot]
(* Output *)
{Flat,OneIdentity,Protected,ReadProtected}
```

Expressions with head [Dot](https://reference.wolfram.com/language/ref/Dot.html) are threadable, unless they are known to be nonscalars:

```wolfram
x.y+{1,2}
(* Output *)
{1+x.y,2+x.y}
```

```wolfram
VectorSymbol[u,n].VectorSymbol[v,n]+{1,2}
(* Output *)
{1+u.v,2+u.v}
```

```wolfram
MatrixSymbol[a,{n,n}].VectorSymbol[v,n]+{1,2}
(* Output *)
a.v+{1,2}
```

## Related Guides ▪Symbolic Vectors, Matrices and Arrays

## History Introduced in 2024 (14.1)
