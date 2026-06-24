# Eliminate | [SpanFromLeft]

> [Eliminate](https://reference.wolfram.com/language/ref/Eliminate.html)[*eqns*,*vars*] — eliminates variables between a set of simultaneous equations.

## Details and Options

Equations are given in the form `*lhs*==*rhs*`.

Simultaneous equations can be combined either in a list or with `&&`.

A single variable or a list of variables can be specified.

Variables can be any expressions.

[Eliminate](https://reference.wolfram.com/language/ref/Eliminate.html) works primarily with linear and polynomial equations.

## Examples

### Basic Examples

Eliminate the variable $y$ between two equations:

```wolfram
Eliminate[{x==2+y,y==z},y]
(* Output *)
2+z==x
```

Eliminate multiple variables in a system of equations:

```wolfram
Eliminate[{f==x^5+y^5,a==x+y,b==x y},{x,y}]
(* Output *)
f==a^5-5 a^3 b+5 a b^2
```

### Scope

A system of linear equations:

```wolfram
Eliminate[2x+3y+4z==1&&9x+8y+7z==2,z]
(* Output *)
1-11 y==22 x
```

A system of polynomial equations:

```wolfram
Eliminate[x^2+y^2+z^2==1&&x-y+z==2&&x^3-y^2==z+1,z]
(* Output *)
-18 x+4 x^2-28 x^3+8 x^4+4 x^5+4 x^6==-27&&y==-12-2 x-5 x^2+8 x^3+4 x^4+2 x^5
```

Eliminate two variables:

```wolfram
Eliminate[x^2+y^2+z^2==1&&x-y+z==2&&x^3-y^2==z+1,{y,z}]
(* Output *)
-18 x+4 x^2-28 x^3+8 x^4+4 x^5+4 x^6==-27
```

A system of equations involving radicals:

```wolfram
Eliminate[Sqrt[x]+Sqrt[y-1]==1&&2x^(1/3)+3 y^2==2,y]
(* Output *)
-2864000 x+4956544 x^2-3109824 x^3+4114800 x^4-1219104 x^5+819936 x^6-310608 x^7+121500 x^8-23328 x^9+14580 x^10+729 x^12==-1000000
```

A system of transcendental equations:

```wolfram
Eliminate[Sin[x+y]==1&&Cos[x-y]==2,y]
(* Output *)
Eliminate
(* Output *)
2 x-ArcCos[2]==(π)/(2)||2 x+ArcCos[2]==-(π)/(2)
```

A system of modular equations:

```wolfram
Eliminate[x^2+y^2+z^2==1&&x-y+z==2&&x^3-y^2==z+1&&Modulus==5,z]
(* Output *)
Modulus==5&&-2 x+x^2-2 x^3+2 x^4+x^5+x^6==2&&y==-2-2 x-2 x^3-x^4+2 x^5
```

Find a modulus for which a system of equations has a solution and eliminate a variable:

```wolfram
Eliminate[y==5x+1&&2x+3y==5&&x y==7,y,Mode->Modular]
(* Output *)
Modulus==1969&&x==-579
```

### Options

#### InverseFunctions

By default, [Eliminate](https://reference.wolfram.com/language/ref/Eliminate.html) uses inverse functions but prints warning messages:

```wolfram
Eliminate[f[x-y]==0 && g[x+y]==0,y]
(* Output *)
InverseFunction
(* Output *)
InverseFunction
(* Output *)
g^((-1))[0]==2 x-f^((-1))[0]
```

With [InverseFunctions](https://reference.wolfram.com/language/ref/InverseFunctions.html)->[True](https://reference.wolfram.com/language/ref/True.html), [Eliminate](https://reference.wolfram.com/language/ref/Eliminate.html) does not print inverse function messages:

```wolfram
Eliminate[f[x-y]==0 && g[x+y]==0,y,InverseFunctions->True]
(* Output *)
g^((-1))[0]==2 x-f^((-1))[0]
```

```wolfram
Eliminate[Sin[x^2-y]==1&&Cos[x-y^2]==2,y,InverseFunctions->True]
(* Output *)
x-x^4-ArcCos[2]==(π^2)/(4)-π x^2||x-x^4+ArcCos[2]==(π^2)/(4)+π x^2
```

With [InverseFunctions](https://reference.wolfram.com/language/ref/InverseFunctions.html)->[False](https://reference.wolfram.com/language/ref/False.html), [Eliminate](https://reference.wolfram.com/language/ref/Eliminate.html) does not use inverse functions:

```wolfram
Eliminate[f[x-y]==0 && g[x+y]==0,y,InverseFunctions->False]
(* Output *)
Eliminate
(* Output *)
f[x-y]==0&&g[x+y]==0
```

```wolfram
Eliminate[Sin[x^2-y]==1&&Cos[x-y^2]==2,y,InverseFunctions->False]
(* Output *)
Eliminate
(* Output *)
Eliminate[Sin[x^2-y]==1&&Cos[x-y^2]==2,y,InverseFunctions->False]
```

Eliminating variables from algebraic equations does not require using inverse functions:

```wolfram
Eliminate[x^2+y^2+z^2==1&&x-y+z==2&&x^3-y^2==z+1,z,InverseFunctions->False]
(* Output *)
-18 x+4 x^2-28 x^3+8 x^4+4 x^5+4 x^6==-27&&y==-12-2 x-5 x^2+8 x^3+4 x^4+2 x^5
```

#### Mode

Find a modulus for which a system of equations has a solution and eliminate a variable:

```wolfram
Eliminate[x^2+y^2==1&&2x+3y==5&&x^3-x y==7,y,Mode->Modular]
(* Output *)
Modulus==191721&&x==47377
```

#### WorkingPrecision

By default, [Eliminate](https://reference.wolfram.com/language/ref/Eliminate.html) computes with exact coefficients:

```wolfram
Eliminate[x^2+y^2+z^2==E&&x-y+z==2&&x^3-y^2==z+Pi,z]
(* Output *)
(-16+2 ℯ-4 π) x+(12-4 ℯ-4 π) x^2+(-16-4 ℯ-8 π) x^3+8 x^4+4 x^5+4 x^6==-8-2 ℯ-ℯ^2-8 π-4 ℯ π-4 π^2&&(-2+ℯ+2 π) y==-4-2 ℯ-6 π+2 x-2 ℯ x-2 π x-2 x^2-ℯ x^2-2 π x^2+8 x^3+4 x^4+2 x^5&&(-2+2 x) y==-ℯ-2 π-2 x+2 x^2+2 x^3&&y+y^2==-2-π+x+x^3
```

This performs the elimination using 20-digit approximate number coefficients:

```wolfram
Eliminate[x^2+y^2+z^2==E&&x-y+z==2&&x^3-y^2==z+Pi,z,WorkingPrecision->20]
(* Output *)
-23.12980695744108248312999859164 x-11.43949792819535389529172341976 x^2-52.0058685425545268491422969541 x^3+8. x^4+4. x^5+4. x^6==-119.59571547961878934284431696075&&7.00146713563863171228557423852 y==-28.28611957845684990149643524422-9.71974896409767694764586170988 x-11.00146713563863171228557423853 x^2+8. x^3+4. x^4+2. x^5&&(-2.+2. x) y==-9.00146713563863171228557423852-2. x+2. x^2+2. x^3&&y+y^2==-5.14159265358979323846264338358+x+x^3
```

### Applications

Rewrite $f$ in terms of $a$ and $b$:

```wolfram
Eliminate[{f==x^5+y^5,a==x+y,b==x y},{x,y}]
(* Output *)
f==a^5-5 a^3 b+5 a b^2
```

Find a condition for two polynomials to have a common root:

```wolfram
f=a x^3+(a^2-1)x^2+(2a-3)x+4;
g=(1-a^2)x^4-2a x+7a^2-1;
Eliminate[f==0&&g==0,x]
(* Output *)
-102 a-2270 a^2+4602 a^3+2549 a^4-9336 a^5+1749 a^6+7464 a^7-4122 a^8-3102 a^9+2372 a^10+546 a^11-623 a^12+49 a^14==0
```

This solves the same problem using [Resolve](https://reference.wolfram.com/language/ref/Resolve.html):

```wolfram
Resolve[Exists[x,f==0&&g==0]]
(* Output *)
a==0||-102-2270 a+4602 a^2+2549 a^3-9336 a^4+1749 a^5+7464 a^6-4122 a^7-3102 a^8+2372 a^9+546 a^10-623 a^11+49 a^13==0
```

The condition is equivalent to [Resultant](https://reference.wolfram.com/language/ref/Resultant.html) of the polynomials being zero:

```wolfram
Resultant[f,g,x]
(* Output *)
102 a+2270 a^2-4602 a^3-2549 a^4+9336 a^5-1749 a^6-7464 a^7+4122 a^8+3102 a^9-2372 a^10-546 a^11+623 a^12-49 a^14
```

### Properties & Relations

Equations returned by [Eliminate](https://reference.wolfram.com/language/ref/Eliminate.html) do not contain the elimination variables:

```wolfram
Eliminate[x^2+y^2+z^2==1&&2x-3y+5z==7,z]
(* Output *)
-24-42 y-34 y^2==29 x^2+x (-28-12 y)
```

Equations returned by [Eliminate](https://reference.wolfram.com/language/ref/Eliminate.html) are implied by the input equations:

```wolfram
inp=x^2+y^2+z^2==1&&x-y+z==2&&x^3-y^2==z+1;
out=Eliminate[x^2+y^2+z^2==1&&x-y+z==2&&x^3-y^2==z+1,z]
(* Output *)
-18 x+4 x^2-28 x^3+8 x^4+4 x^5+4 x^6==-27&&y==-12-2 x-5 x^2+8 x^3+4 x^4+2 x^5
```

Use [Resolve](https://reference.wolfram.com/language/ref/Resolve.html) to check this property:

```wolfram
Resolve[ForAll[{x,y,z},Implies[inp,out]]]
(* Output *)
True
```

Use [Resolve](https://reference.wolfram.com/language/ref/Resolve.html) to eliminate an existential quantifier:

```wolfram
eqns=x z==y && x^2+x^2z^2==1;
Resolve[Exists[z,eqns]]
(* Output *)
-1+x^2+y^2==0&&x≠0
```

[Eliminate](https://reference.wolfram.com/language/ref/Eliminate.html) gives the same set of equations, but does not give inequations:

```wolfram
Eliminate[eqns,z]
(* Output *)
y^2==1-x^2
```

Eliminate a variable using [GroebnerBasis](https://reference.wolfram.com/language/ref/GroebnerBasis.html):

```wolfram
GroebnerBasis[eqns,{x,y},{z}]
(* Output *)
{-1+x^2+y^2}
```

Eliminate a variable from a pair of polynomials using [Resultant](https://reference.wolfram.com/language/ref/Resultant.html):

```wolfram
Resultant[eqns[[1,1]]-eqns[[1,2]],eqns[[2,1]]-eqns[[2,2]],z]//Factor
(* Output *)
x^2 (-1+x^2+y^2)
```

Use [Reduce](https://reference.wolfram.com/language/ref/Reduce.html) to eliminate an existential quantifier and solve the resulting system:

```wolfram
Reduce[Exists[z,eqns],{x,y}]
(* Output *)
(y==-Sqrt[1-x^2]||y==Sqrt[1-x^2])&&x≠0
```

### Possible Issues

When the input contains only equations, [Eliminate](https://reference.wolfram.com/language/ref/Eliminate.html) returns only equations:

```wolfram
eqns=x z==y && x^2+x^2z^2==1
(* Output *)
x z==y&&x^2+x^2 z^2==1
```

The zero set of the result is the Zariski closure of the projection of the zero set of `*eqns*`:

```wolfram
Eliminate[eqns,z]
(* Output *)
y^2==1-x^2
```

Use [Resolve](https://reference.wolfram.com/language/ref/Resolve.html) to get a result with zero set equal to the projection of the zero set of `*eqns*`:

```wolfram
Resolve[Exists[z,eqns]]
(* Output *)
-1+x^2+y^2==0&&x≠0
```

When the input contains inequations, [Eliminate](https://reference.wolfram.com/language/ref/Eliminate.html) returns equations and inequations:

```wolfram
eqin=eqns&&t≠0
(* Output *)
x z==y&&x^2+x^2 z^2==1&&t≠0
```

The zero set of the result is equal to the projection of the zero set of `*eqin*`:

```wolfram
Eliminate[eqin,z]
(* Output *)
y^2==1-x^2&&t≠0&&x≠0
```

## Tech Notes ▪Solving Equations ▪Eliminating Variables ▪Solving Logical Combinations of Equations

## Related Guides ▪Polynomial Systems ▪Polynomial Algebra

## History Introduced in 1988 (1.0)
