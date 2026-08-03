# Join | [SpanFromLeft]

> [Join](https://reference.wolfram.com/language/ref/Join.html)[*list*_1,*list*_2,…] — concatenates lists or other expressions that share the same head.
> [Join](https://reference.wolfram.com/language/ref/Join.html)[*list*_1,*list*_2,…,*n*] — joins the objects at level `*n*` in each of the `*list*_*i*`.

## Details

The `*list*_*i*` do not need to have head [List](https://reference.wolfram.com/language/ref/List.html), but must all have the same head.

[Join](https://reference.wolfram.com/language/ref/Join.html) works on [Association](https://reference.wolfram.com/language/ref/Association.html) objects, keeping the last value associated with any given key.

[Join](https://reference.wolfram.com/language/ref/Join.html) works on [SparseArray](https://reference.wolfram.com/language/ref/SparseArray.html) objects by effectively concatenating the corresponding ordinary lists.

[Join](https://reference.wolfram.com/language/ref/Join.html)[*list*_1,*list*_2,…,*n*] handles ragged arrays by effectively concatenating all successive elements at level `*n*` in each of the `*list*_*i*`.

## Examples

### Basic Examples

```wolfram
Join[{a,b,c},{x,y},{u,v,w}]
(* Output *)
{a,b,c,x,y,u,v,w}
```

Infix syntax:

```wolfram
{a,b,c}~Join~{x,y}~Join~{u,v,w}
(* Output *)
{a,b,c,x,y,u,v,w}
```

Join two associations:

```wolfram
Join[<|a->b|>,<|c->d,a->e|>]
(* Output *)
<|a->e,c->d|>
```

### Scope

Join two matrices to make longer columns:

```wolfram
Join[{{a,b},{c,d}}, {{1,2},{3,4}}]
(* Output *)
{{a,b},{c,d},{1,2},{3,4}}
```

Join columns of two matrices to make longer rows:

```wolfram
Join[{{a,b},{c,d}}, {{1,2},{3,4}}, 2]
(* Output *)
{{a,b,1,2},{c,d,3,4}}
```

With ragged arrays, successive elements are effectively concatenated:

```wolfram
Join[{{1},{5, 6}},{{2, 3},{7}},{{4},{8}},2]
(* Output *)
{{1,2,3,4},{5,6,7,8}}
```

The second row comes from the concatenation of nothing with `{3,4}`:

```wolfram
Join[{{x}}, {{1,2},{3,4}},2]
(* Output *)
{{x,1,2},{3,4}}
```

Join depth 3 arrays at different levels:

```wolfram
aa = Array[a_#&, {2,2,2}];
bb = Array[b_#&, {2,2,2}];
```

```wolfram
Join[aa, bb]
(* Output *)
{{{a_1,a_1},{a_1,a_1}},{{a_2,a_2},{a_2,a_2}},{{b_1,b_1},{b_1,b_1}},{{b_2,b_2},{b_2,b_2}}}
```

```wolfram
Join[aa, bb, 2]
(* Output *)
{{{a_1,a_1},{a_1,a_1},{b_1,b_1},{b_1,b_1}},{{a_2,a_2},{a_2,a_2},{b_2,b_2},{b_2,b_2}}}
```

```wolfram
Join[aa, bb, 3]
(* Output *)
{{{a_1,a_1,b_1,b_1},{a_1,a_1,b_1,b_1}},{{a_2,a_2,b_2,b_2},{a_2,a_2,b_2,b_2}}}
```

### Generalizations & Extensions

Join expressions with any head:

```wolfram
Join[f[a,b,c],f[x,y],f[u,v,w]]
(* Output *)
f[a,b,c,x,y,u,v,w]
```

[Join](https://reference.wolfram.com/language/ref/Join.html) works with [SparseArray](https://reference.wolfram.com/language/ref/SparseArray.html) objects:

```wolfram
SparseArray[Range[5]]
(* Output *)
SparseArray[...]
```

```wolfram
Join[%,%,%]
(* Output *)
SparseArray[...]
```

### Applications

Augment a matrix by adding a row:

```wolfram
Join[IdentityMatrix[3], {{1,2,3}}]//MatrixForm
(* Output *)
({{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 2, 3}})
```

Augment by a column:

```wolfram
Join[IdentityMatrix[3], Transpose[{{1,2,3}}],2]//MatrixForm
(* Output *)
({{1, 0, 0, 1}, {0, 1, 0, 2}, {0, 0, 1, 3}})
```

Make a block matrix:

```wolfram
Join[Join[({{a_1, a_1, a_1}, {a_2, a_2, a_2}, {a_3, a_3, a_3}}), ({{b_1, b_1, b_1}, {b_2, b_2, b_2}, {b_3, b_3, b_3}}),2], Join[({{c_1, c_1, c_1}, {c_2, c_2, c_2}, {c_3, c_3, c_3}}),({{d_1, d_1, d_1}, {d_2, d_2, d_2}, {d_3, d_3, d_3}}),2]]//MatrixForm
(* Output *)
({{a_1, a_1, a_1, b_1, b_1, b_1}, {a_2, a_2, a_2, b_2, b_2, b_2}, {a_3, a_3, a_3, b_3, b_3, b_3}, {c_1, c_1, c_1, d_1, d_1, d_1}, {c_2, c_2, c_2, d_2, d_2, d_2}, {c_3, c_3, c_3, d_3, d_3, d_3}})
```

This can also be done with [ArrayFlatten](https://reference.wolfram.com/language/ref/ArrayFlatten.html):

```wolfram
ArrayFlatten[({{({{a_1, a_1, a_1}, {a_2, a_2, a_2}, {a_3, a_3, a_3}}), ({{b_1, b_1, b_1}, {b_2, b_2, b_2}, {b_3, b_3, b_3}})}, {({{c_1, c_1, c_1}, {c_2, c_2, c_2}, {c_3, c_3, c_3}}), ({{d_1, d_1, d_1}, {d_2, d_2, d_2}, {d_3, d_3, d_3}})}})]//MatrixForm
(* Output *)
({{a_1, a_1, a_1, b_1, b_1, b_1}, {a_2, a_2, a_2, b_2, b_2, b_2}, {a_3, a_3, a_3, b_3, b_3, b_3}, {c_1, c_1, c_1, d_1, d_1, d_1}, {c_2, c_2, c_2, d_2, d_2, d_2}, {c_3, c_3, c_3, d_3, d_3, d_3}})
```

### Properties & Relations

[Join](https://reference.wolfram.com/language/ref/Join.html)[*list*_1,*list*_2,…] is equivalent to [Flatten](https://reference.wolfram.com/language/ref/Flatten.html)[{*list*_1,*list*_2,…},1]:

```wolfram
Join[{1,2},{{a,b},{c,d}}, {3,4,5}]
(* Output *)
{1,2,{a,b},{c,d},3,4,5}
```

```wolfram
Flatten[{{1,2},{{a,b},{c,d}}, {3,4,5}},1]
(* Output *)
{1,2,{a,b},{c,d},3,4,5}
```

### Neat Examples

Successively double a list by joining to itself:

```wolfram
NestList[Join[#,#]&,{x},4]
(* Output *)
{{x},{x,x},{x,x,x,x},{x,x,x,x,x,x,x,x},{x,x,x,x,x,x,x,x,x,x,x,x,x,x,x,x}}
```

Build up the Thue-Morse sequence [[more info](http://www.wolframscience.com/nksonline/page-889c-text)]:

```wolfram
NestList[Join[#,1-#]&,{1},5]
(* Output *)
{{1},{1,0},{1,0,0,1},{1,0,0,1,0,1,1,0},{1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1},{1,0,0,1,0,1,1,0,0,1,1,0,1,0,0,1,0,1,1,0,1,0,0,1,1,0,0,1,0,1,1,0}}
```

## Tech Notes ▪Rearranging Nested Lists ▪Structural Operations ▪Combining Lists

## Related Guides ▪Structural Operations on Expressions ▪Rearranging & Restructuring Lists ▪Constructing Matrices ▪Tabular Transformation ▪Parts of Matrices ▪List Manipulation ▪Handling Arrays of Data

## Related Links [An Elementary Introduction to the Wolfram Language](https://www.wolfram.com/language/elementary-introduction/03-first-look-at-lists.html)[: First Look at Lists](https://www.wolfram.com/language/elementary-introduction/03-first-look-at-lists.html) [An Elementary Introduction to the Wolfram Language](https://www.wolfram.com/language/elementary-introduction/27-applying-functions-repeatedly.html)[: Applying Functions Repeatedly](https://www.wolfram.com/language/elementary-introduction/27-applying-functions-repeatedly.html) [An Elementary Introduction to the Wolfram Language](https://www.wolfram.com/language/elementary-introduction/38-assigning-names-to-things.html)[: Assigning Names to Things](https://www.wolfram.com/language/elementary-introduction/38-assigning-names-to-things.html) [NKS|Online](http://www.wolframscience.com/nks/search/?q=Join) ([A New Kind of Science](http://www.wolframscience.com/nks/))

## History Introduced in 1988 (1.0) | Updated in 2003 (5.0) ▪ 2007 (6.0) ▪ 2014 (10.0)
