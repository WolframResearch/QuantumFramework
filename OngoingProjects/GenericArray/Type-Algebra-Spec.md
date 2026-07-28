# Type algebra for mixed array containers - specification

Status: specification, not yet implemented. Written against `Wolfram/Arrays` as of the
container/lazy-registry work (kernel files `Classification.wl`, `Shape.wl`,
`Accessors.wl`, `Structural.wl`, `Vector.wl`, `ArrayObject.wl`).

The requirement: when containers of different kinds meet in one operation, the result
must still be a container, its tier and element domain follow from the operands by a
stated rule, and the operands are coerced to a common representation rather than one
silently winning. `ArrayObject` wraps such a result like any other container.

## 1. What exists today

- `arrayTier` (`ArrayObject.wl:109`) computes `"Explicit" | "Lazy" | "Symbolic"`, but is
  **private to `ArrayObject.wl`** and reachable only as the `"Tier"` property of a
  handle. Nothing dispatches on it.
- `arrayElementType` (`ArrayObject.wl:117-121`) covers `NumericArray` and
  `TabularColumn` and answers `Missing["NotApplicable"]` for every other container. It
  knows nothing of symbolic domains.
- There are **no binary or n-ary operations** over containers. The only site where
  containers meet is `ArrayContract[{a, b, ...}, c]`, whose list argument is exactly the
  mixing point, and the tensor-network contraction that calls into it.
- Consequently there is no rule for what a mixed operation returns, and no coercion.

## 2. Axis 1 - tier

A total order, most specific to most general:

```
Explicit  <  Lazy  <  Symbolic
```

The join of a set of operands is the maximum. This gives exactly the required cases:

| operands | result tier |
|---|---|
| explicit, explicit | Explicit |
| lazy, explicit | Lazy |
| symbolic, lazy | Symbolic |
| symbolic, explicit | Symbolic |

The order is not arbitrary: binding a lazy container's parameters yields an explicit
one, and specializing a symbolic container yields a lazy or explicit one, so generality
increases to the right. A result may be *more* specific than the join only when an
operation genuinely collapses generality (a symbolic array multiplied by a zero array,
say); such collapses are opportunistic simplifications, never the contract.

## 3. Axis 2 - element domain

Numeric container types and symbolic domains are two spellings of one lattice, and must
unify:

```
Integers  <  Rationals  <  Algebraics  <  Reals  <  Complexes
```

Numeric types map into it, retaining precision as a secondary component:

| container type | domain | precision |
|---|---|---|
| `"Integer8" ... "Integer64"`, `"UnsignedInteger*"` | Integers | the type's bit width |
| `"Real32"`, `"Real64"` | Reals | 32, 64 |
| `"ComplexReal32"`, `"ComplexReal64"` | Complexes | 32, 64 |
| exact `List`/`SparseArray` entries | inferred from the values | exact |
| `QuantityArray` | domain of its magnitudes | as magnitudes, unit is separate metadata |
| symbolic (`ArraySymbol` etc.) | its declared `$Assumptions` domain, default Complexes | symbolic |

Join rules:

- Within the same tier and both numeric: widen the domain, then take the wider
  precision. `Integer64` with `Real32` gives `Real64`, not `Real32`: widening the domain
  must not silently narrow precision below what an operand already carried.
- Exact with machine: the result is machine only where the operation is machine; an
  exact operand does not force exactness on a machine one, but `ArrayPack`'s fidelity
  rule still applies to any packing decision (exact values that do not survive at
  machine precision are not coerced).
- Across tiers: the symbolic domain absorbs. `Real64` joined with a symbolic array over
  `Reals` is `Reals` (symbolic), not `Real64` - the symbolic operand cannot be narrowed
  to a machine type without materializing it, which the tier join already forbids.
- Unknown or inapplicable domain on one side (`Missing["NotApplicable"]`) does not
  poison the join: it contributes nothing, and the result takes the other operand's
  domain. Only if every operand is unknown is the result unknown.

## 4. Coercion

Two entry points, both operating on raw containers (never on `ArrayObject`):

- `ArrayCoerce[a, spec]` converts one container toward a target, where `spec` names a
  tier, an element domain or type, or both. Coercion may only move *along* the
  lattices: an explicit array coerces to lazy or symbolic freely, while coercing a lazy
  container down to explicit is materialization and must be requested explicitly, since
  it can be arbitrarily expensive and, for an unbound-parameter container, impossible.
- `ArrayUnify[{a, b, ...}]` computes the joined tier and domain and returns the operands
  coerced to it, so a caller can then apply a single-representation operation. This is
  the function a mixing operation calls first.

Failure is a clean refusal: operands whose domains have no join in the lattice, or a
requested narrowing that would lose values, fail with a message naming both operands'
types, never a silent cast.

## 5. Public surface

Promote the two private helpers and add the algebra:

```
ArrayTier[a]                  "Explicit" | "Lazy" | "Symbolic"
ArrayElementDomain[a]         the unified domain (Integers ... Complexes), or Missing
ArrayElementType[a]           the concrete container type where one exists (Real64 ...)
ArrayUnify[{a, b, ...}]       joined tier and domain, plus the coerced operands
ArrayCoerce[a, spec]          coerce one container toward a tier and/or domain
```

`ArrayObject` gains `"Domain"` beside its existing `"Tier"` and `"ElementType"`, and its
summary box shows the domain where the element type is inapplicable, so a symbolic array
over `Reals` and a `Real64` `NumericArray` describe themselves in comparable terms.

## 6. Where the algebra must be applied

- `ArrayContract[{a, b, ...}, c]` - the existing mixing point.
- Any binary operation added later (addition, product, `Dot`).
- The tensor-network contraction, whose leaf tensors may be of different containers:
  with the algebra in place, a network mixing a `SparseArray` leaf, a lazy leaf and a
  symbolic leaf contracts to a symbolic container, and one mixing a `SparseArray` with
  net-backed leaves contracts to a net.

## 7. Tests to pin

- The full tier join table (3 x 3), including reflexive cases.
- Domain join across every adjacent pair in the chain, plus the precision rule
  (`Integer64` with `Real32` gives `Real64`).
- Symbolic absorption: a symbolic array over `Reals` joined with each explicit numeric
  type gives symbolic over `Reals`.
- `Missing` domain contributes nothing rather than poisoning.
- Coercion refuses narrowing that loses values, with a message.
- A mixed `ArrayContract` returns a container satisfying `ArrayContainerQ` whose tier is
  the join, and `ArrayObject` wraps that result and reports the joined tier and domain.
