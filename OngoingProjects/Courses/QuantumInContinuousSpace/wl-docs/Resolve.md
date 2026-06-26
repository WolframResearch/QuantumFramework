# Resolve | [SpanFromLeft]

> [Resolve](https://reference.wolfram.com/language/ref/Resolve.html)[*expr*] — attempts to resolve `*expr*` into a form that eliminates [ForAll](https://reference.wolfram.com/language/ref/ForAll.html) and [Exists](https://reference.wolfram.com/language/ref/Exists.html) quantifiers.
> [Resolve](https://reference.wolfram.com/language/ref/Resolve.html)[*expr*,*dom*] — works over the domain `*dom*`. Common choices of `*dom*` are [Complexes](https://reference.wolfram.com/language/ref/Complexes.html), [Reals](https://reference.wolfram.com/language/ref/Reals.html), and [Booleans](https://reference.wolfram.com/language/ref/Booleans.html).

## Details and Options

[Resolve](https://reference.wolfram.com/language/ref/Resolve.html) is in effect automatically applied by [Reduce](https://reference.wolfram.com/language/ref/Reduce.html).

`*expr*` can contain equations, inequalities, domain specifications, and quantifiers, in the same form as in [Reduce](https://reference.wolfram.com/language/ref/Reduce.html).

The statement `*expr*` can be any logical combination of:

*lhs*==*rhs* | equations
*lhs*!=*rhs* | inequations
`*lhs*>*rhs*` or `*lhs*>=*rhs*` | inequalities
*expr*∈*dom* | domain specifications
{*x*,*y*,…}∈*reg* | region specification
[ForAll](https://reference.wolfram.com/language/ref/ForAll.html)[*x*,*cond*,*expr*] | universal quantifiers
[Exists](https://reference.wolfram.com/language/ref/Exists.html)[*x*,*cond*,*expr*] | existential quantifiers

The result of [Resolve](https://reference.wolfram.com/language/ref/Resolve.html)[*expr*] always describes exactly the same mathematical set as `*expr*`, but without quantifiers.

[Resolve](https://reference.wolfram.com/language/ref/Resolve.html)[*expr*] assumes by default that quantities appearing algebraically in inequalities are real, while all other quantities are complex.

When a quantifier such as [ForAll](https://reference.wolfram.com/language/ref/ForAll.html)[*x*,…] is eliminated, the result will contain no mention of the localized variable `*x*`.

[Resolve](https://reference.wolfram.com/language/ref/Resolve.html)[*expr*] can in principle always eliminate quantifiers if `*expr*` contains only polynomial equations and inequalities over the reals or complexes.

[Resolve](https://reference.wolfram.com/language/ref/Resolve.html)[*expr*] can in principle always eliminate quantifiers for any Boolean expression `*expr*`.

## Examples

### Basic Examples

Prove that the unit disk is nonempty:

```wolfram
Resolve[Exists[{x,y},x^2+y^2<1]]
(* Output *)
True
```

Find the conditions for a quadratic form over the reals to be positive:

```wolfram
Resolve[ForAll[x,ax^2+bx+c>0],Reals]
(* Output *)
(a==0&&b==0&&c>0)||(a>=0&&b==0&&c>0&&-b^2+4 a c>0)||(a>0&&-b^2+4 a c>0)
```

Find conditions for a quadratic to have at least two distinct complex roots:

```wolfram
Resolve[Exists[{x,y},a x^2+b x+c==0 && a y^2+b y+c==0 && x≠y]]
(* Output *)
(a==0&&b==0&&c==0)||(a≠0&&-b^2+4 a c≠0)
```

Find the projection of a geometric region:

```wolfram
Resolve[Exists[z,{x,y,z}∈Cylinder[{{0,0,0},{1,1,1}},2]],Reals]
(* Output *)
x^2+x y+y^2<=2||-3 x+x^2-3 y+x y+y^2<=-1||(x+y<=2&&x+y>=0&&x^2-2 x y+y^2<=8)
```

```wolfram
RegionPlot[%,{x,-2,3},{y,-2,3}]
```

*([Graphics])*

### Scope

#### Complex Domain

Decide the existence of solutions of a univariate polynomial equation:

```wolfram
Resolve[Exists[x, x^7+3 x-11==0]]
(* Output *)
True
```

Decide the existence of solutions of a multivariate polynomial system:

```wolfram
Resolve[Exists[{x,y},x^2+y^2==1&&x^3-y^2==2&&x y≠3]]
(* Output *)
True
```

Decide the truth value of fully quantified polynomial formulas:

```wolfram
Resolve[ForAll[{a,b,c},Exists[x,a x^2+b x+c==0]]]
(* Output *)
False
```

```wolfram
Resolve[ForAll[{b,c},Exists[x,x^2+b x+c==0]]]
(* Output *)
True
```

Find conditions under which a polynomial equation has solutions:

```wolfram
Resolve[Exists[x, (a^2-b^3+1)x^7+(a^2-a b+1) x^2-3==0]]
(* Output *)
1+a^2-a b≠0||1+a^2-b^3≠0
```

Find conditions under which a polynomial system has solutions:

```wolfram
Resolve[Exists[{x,y},a x^2+b y^2==1&&b x-a y==2&&x≠y]]
(* Output *)
(a≠0&&b==0)||(a≠0&&a+b==0&&-4 a^2-4 a b-a b^2+b^3≠0)||(a≠0&&a^2-a b+b^2==0&&-4 a^2-4 a b-a b^2+b^3≠0)||(b≠0&&-a+b≠0&&a+b≠0&&a^2+b^2≠0&&a^2-a b+b^2≠0)||(b≠0&&a+b≠0&&a^2-a b+b^2≠0&&8 a^2-4 a b+a^2 b+4 b^2-2 a b^2+b^3≠0)
```

Find conditions under which a quantified polynomial formula is true:

```wolfram
Resolve[ForAll[c,Exists[x,a x^2+b x+c==0]]]
(* Output *)
a≠0||b≠0
```

#### Real Domain

Decide the existence of solutions of a univariate polynomial equation:

```wolfram
Resolve[Exists[x, x^7+3 x-11==0],Reals]
(* Output *)
True
```

Decide the existence of solutions of a univariate polynomial inequality:

```wolfram
Resolve[Exists[x, x^4+3 x+11<=0],Reals]
(* Output *)
False
```

Decide the existence of solutions of a multivariate polynomial system:

```wolfram
Resolve[Exists[{x,y},x^2+y^2<=11&&x^3-y^2==2&&x y>3],Reals]
(* Output *)
True
```

Decide the truth value of fully quantified polynomial formulas:

```wolfram
Resolve[ForAll[{a,b,c},Exists[x,a x^2+b x+c>=0]],Reals]
(* Output *)
False
```

```wolfram
Resolve[ForAll[{b,c},Exists[x,x^2+b x+c>=0]],Reals]
(* Output *)
True
```

Decide the existence of solutions of an exp-log equation:

```wolfram
Resolve[Exists[x,Exp[x]-2x Log[x]+3x^x==4],Reals]
(* Output *)
True
```

Decide the existence of solutions of an exp-log inequality:

```wolfram
Resolve[Exists[x,x^Sqrt[Pi]-2Log[1/(x^2+1)]+E^x<1],Reals]
(* Output *)
False
```

Decide the existence of solutions of an elementary function equation in a bounded interval:

```wolfram
Resolve[Exists[x,-1<x<1,Sin[Sin[Sin[x]]]-x==1/7],Reals]
(* Output *)
True
```

Decide the existence of solutions of a holomorphic function equation in a bounded interval:

```wolfram
Resolve[Exists[x,1<x<2,Gamma[x]-BesselJ[2,x]==x+1],Reals]
(* Output *)
Resolve
(* Output *)
False
```

Decide the existence of solutions of a periodic elementary function equation:

```wolfram
Resolve[Exists[x,Exp[Sin[x]]-Sin[3 Cos[x]]==Tan[x]],Reals]
(* Output *)
True
```

Fully quantified formulas exp-log in the first variable and polynomial in the other variables:

```wolfram
Resolve[Exists[x,ForAll[y,E^x y^3+Log[x]y==1&&x y+E^x/x>=2]],Reals]
(* Output *)
False
```

```wolfram
Resolve[Exists[{x,y},y^2 E^E^x-y E^x-x==1],Reals]
(* Output *)
True
```

Fully quantified formulas elementary and bounded in the first variable:

```wolfram
Resolve[ForAll[x,-1<x<0,Exists[y,Sin[x-Cos[x]] y^3-x==1]],Reals]
(* Output *)
True
```

```wolfram
Resolve[Exists[{x,y},Sin[x-Exp[x]] y^3-x==2&&x^2+y^2<=1],Reals]
(* Output *)
False
```

Fully quantified formulas holomorphic and bounded in the first variable:

```wolfram
Resolve[ForAll[x,-1<x<1,Exists[y,y^3-BesselJ[2,x+2] y-y-3 x==-2]],Reals]
(* Output *)
True
```

```wolfram
Resolve[Exists[{x,y},1<x<2&&Gamma[x] y^2-x==-5],Reals]
(* Output *)
False
```

Find conditions under which a linear system has solutions:

```wolfram
Resolve[Exists[{x,y},a x+b y>=1&&(a^2-b)x-(a-b^2)y==2],Reals]
(* Output *)
(a^2-b<0&&a^2+a^2 b-b^2-a b^2<0)||(a^2-b<0&&-a^2-a^2 b+b^2+a b^2<0)||(-a^2+b<0&&a^2+a^2 b-b^2-a b^2<0)||(-a^2+b<0&&-a^2-a^2 b+b^2+a b^2<0)||(-a<0&&a^2-b==0&&a-b^2≠0)||(a<0&&a^2-b==0&&a-b^2≠0)||(a^2-b<0&&2 a-a^2+b<=0&&a^2+a^2 b-b^2-a b^2==0)||(-a^2+b<0&&-2 a+a^2-b<=0&&a^2+a^2 b-b^2-a b^2==0)||(a==0&&a^2-b==0&&-b<0&&a+2 b-b^2==0)||(a==0&&a^2-b==0&&b<0&&a+2 b-b^2==0)||(a==0&&a^2-b==0&&a-b^2<0&&-a-2 b+b^2<=0)||(a==0&&a^2-b==0&&-a+b^2<0&&a+2 b-b^2<=0)
```

Find conditions under which a quadratic system has solutions:

```wolfram
Resolve[Exists[x,a x^2+b x+c>=0&&c x^2+a^3-b^3==2],Reals]
(* Output *)
(a==0&&b≠0&&-2 b^2+a^3 b^2-b^5+c^3==0)||(a==0&&b<=0&&a^3-b^3==2&&c==0)||(a≠0&&-b^2+4 a c<=0&&4 a^2-4 a^5+a^8+4 a^2 b^3-2 a^5 b^3+a^2 b^6-2 b^2 c+a^3 b^2 c-b^5 c+4 a c^2-2 a^4 c^2+2 a b^3 c^2+c^4==0)||(a>0&&a^3-b^3==2&&c==0)||(c≠0&&-2 c+a^3 c-b^3 c<=0&&-2 a c+a^4 c-a b^3 c-c^3<=0)||(c≠0&&-2 c+a^3 c-b^3 c<=0&&4 a^2 c^2-4 a^5 c^2+a^8 c^2+4 a^2 b^3 c^2-2 a^5 b^3 c^2+a^2 b^6 c^2-2 b^2 c^3+a^3 b^2 c^3-b^5 c^3+4 a c^4-2 a^4 c^4+2 a b^3 c^4+c^6<=0)
```

Find conditions under which a polynomial system has solutions:

```wolfram
Resolve[Exists[{x,y},a x^4+b y^4==1&&x-a y^3>=2],Reals]
(* Output *)
(a<=0&&b>0)||0<a<=(1)/(16)||(a>(1)/(16)&&b<=Root[-531441 a^20+(-314928 a^15+25509168 a^16) #1^3+(-69984 a^10-57946752 a^11-408146688 a^12) #1^6+(-6912 a^5+5038848 a^6-216670464 a^7+2176782336 a^8) #1^9+(-256+4096 a) #1^12&,2])
```

Find conditions under which a formula linear in quantified variables is true:

```wolfram
Resolve[ForAll[c,Exists[{x,y},a x+b y+c==0&&c x+a y>=b]],Reals]
(* Output *)
(a==0&&b<0)||(a==0&&b≠0&&b^2<=0)||(a≠0&&b==0)||(a≠0&&a^4+a b^3==0)||(a<0&&a^4+a b^3>0)||(a>0&&a^4+a b^3<0)
```

Find conditions under which a formula quadratic in quantified variables is true:

```wolfram
Resolve[ForAll[c,c^2<a^2-b^3,Exists[x,a x^2+b x+c==0]],Reals]
(* Output *)
(a==0&&b≠0)||(a==0&&a^2-b^3<=0)||(a≠0&&b^2>=0&&a b^2==0&&16 a^4-16 a^2 b^3-b^4<=0)||(a≠0&&b^2>=0&&16 a^4-16 a^2 b^3-b^4<0)||(a<0&&b^2>=0&&a b^2<0&&16 a^4-16 a^2 b^3-b^4<=0)||(a>0&&b^2>=0&&a b^2>0&&16 a^4-16 a^2 b^3-b^4<=0)||(a b^2==0&&a^2-b^3<=0&&16 a^4-16 a^2 b^3-b^4<=0)||(b^2>=0&&a^2-b^3<=0&&16 a^4-16 a^2 b^3-b^4<0)||(a>=0&&a b^2>=0&&a^2-b^3<=0&&16 a^4-16 a^2 b^3-b^4<=0)||(a^2-b^3<=0&&16 a^4-16 a^2 b^3-b^4<0)||(a<=0&&a b^2<=0&&a^2-b^3<=0&&16 a^4-16 a^2 b^3-b^4<=0)
```

Find conditions under which a quantified polynomial formula is true:

```wolfram
Resolve[ForAll[y,Exists[x,(a^2-1)x y-(2 a-1) y^4+a^3>=2]],Reals]
(* Output *)
a>=Root
```

#### Integer Domain

Decide the existence of solutions of a linear system of equations:

```wolfram
Resolve[Exists[{x,y,z}, 2 x+3y-5z==1&&3x-4y+7z==3],Integers]
(* Output *)
True
```

Decide the existence of solutions of a linear system of equations and inequalities:

```wolfram
Resolve[Exists[{x,y,z},2 x+3y==4&&3x-4y<=5&&x-2y>-21],Integers]
(* Output *)
True
```

Decide the existence of solutions of a univariate polynomial equation:

```wolfram
Resolve[Exists[x,x^1000-2x^777+1==0],Integers]
(* Output *)
True
```

Decide the existence of solutions of a univariate polynomial inequality:

```wolfram
Resolve[Exists[x, x^5-2x+1<0],Integers]
(* Output *)
True
```

Decide the existence of solutions of Frobenius equations:

```wolfram
Resolve[Exists[{x,y,z,t},1221x+3434y+5566z+8778t==90909&&x>=0&&y>=0&&z>=0&&t>=0],Integers]
(* Output *)
False
```

```wolfram
Resolve[Exists[{x,y,z,t},1221x+3434y+5566z+8778t==909090&&x>=0&&y>=0&&z>=0&&t>=0],Integers]
(* Output *)
True
```

Decide the existence of solutions of binary quadratic equations:

```wolfram
Resolve[Exists[{x,y},x^2+x y+y^2==107],Integers]
(* Output *)
False
```

```wolfram
Resolve[Exists[{x,y},x^2-3y^2==22&&x>0&&y>0],Integers]
(* Output *)
True
```

Decide the existence of solutions of a Thue equation:

```wolfram
Resolve[Exists[{x,y},x^3-2x^2 y+y^3==2],Integers]
(* Output *)
True
```

Decide the existence of solutions of a sum of squares equation:

```wolfram
Resolve[Exists[{x,y,z,t},x^2+4y^2+9z^2+16t^2==354],Integers]
(* Output *)
True
```

Decide the existence of solutions of a bounded system of equations and inequalities:

```wolfram
Resolve[Exists[{x,y,z},x^4+y^4+z^4<=500&&x+y^2+z^3==32],Integers]
(* Output *)
True
```

Decide the existence of solutions of a system of congruences:

```wolfram
Resolve[Exists[{x,y},Mod[x^2+y^2,2]==1&&Mod[x-2y,3]==2],Integers]
(* Output *)
True
```

#### Boolean Domain

Decide the satisfiability of a Boolean formula:

```wolfram
Resolve[Exists[{a,b,c,d},Xor[a,b,c,d]&&(a||b)&&(c||d)],Booleans]
(* Output *)
True
```

Find conditions under which a quantified Boolean formula is true:

```wolfram
Resolve[Exists[a,ForAll[b,Nand[a,b,c,d]&&(a||b)&&(c||d)]],Booleans]
(* Output *)
(c&&!d)||(!c&&d)
```

#### Finite Field Domains

Decide the existence of solutions of univariate equations:

```wolfram
ℱ=FiniteField[53,4];
Resolve[Exists[x,x^5+ℱ[123]x==ℱ[234]]]
(* Output *)
True
```

```wolfram
Resolve[Exists[x,x^9+2 x+3==0],ℱ]
(* Output *)
True
```

Verify that a univariate equation is satisfied by all field elements:

```wolfram
ℱ=FiniteField[3,2];
Resolve[ForAll[x,x^8==x],ℱ]
(* Output *)
True
```

Decide the existence of solutions of systems of linear equations:

```wolfram
ℱ=FiniteField[71,2];
Resolve[Exists[{x,y},ℱ[123]x+ℱ[234]y==ℱ[345]&&ℱ[321]x+ℱ[432]y==ℱ[543]]]
(* Output *)
True
```

```wolfram
Resolve[Exists[{x,y,z},ℱ[1234]x+ℱ[2345]y+ℱ[3456]z==ℱ[4567]&&ℱ[1]x+ℱ[2]y+ℱ[3]z==ℱ[4]&&ℱ[2471]x+ℱ[4696]y+ℱ[1809]z==ℱ[123]]]
(* Output *)
False
```

Decide the existence of solutions of systems of polynomial equations:

```wolfram
ℱ=FiniteField[7,5];
Resolve[Exists[{x,y,z},x^2+y^2+z^2==21],ℱ]
(* Output *)
True
```

```wolfram
Resolve[Exists[{x,y,z},ℱ[321]x^3+ℱ[432]y^3+ℱ[543]z^3==ℱ[654]&&x^2==ℱ[333]y z+ℱ[111]]]
(* Output *)
True
```

Eliminate quantifiers:

```wolfram
ℱ=FiniteField[2,5];
Resolve[Exists[z,ℱ[1]x+ℱ[3]y+ℱ[5]z==ℱ[7]&&ℱ[21]x+ℱ[23]y+ℱ[25]z==ℱ[27]]]
(* Output *)
![image](img/image_001.png)
```

```wolfram
Resolve[Exists[{y,z},ℱ[1]x^2+ℱ[2]y^3+ℱ[3]z^4==ℱ[4]&&ℱ[5]x^4+ℱ[6]y^3+ℱ[7]z^2==ℱ[8]&&x y z!=ℱ[0]]]
(* Output *)
![image](img/image_003.png)
```

#### Mixed Domains

Decide the existence of solutions of an equation involving a real and a complex variable:

```wolfram
Resolve[Exists[{x,y},Element[x,Reals],x^2+y^2==-1]]
(* Output *)
True
```

Decide the existence of solutions of an inequality involving [Abs](https://reference.wolfram.com/language/ref/Abs.html)[*x*]:

```wolfram
Resolve[Exists[x,Abs[x]^2-3Abs[x]<1]]
(* Output *)
True
```

Find under what conditions a fourth power of a complex number is real:

```wolfram
Resolve[Exists[x,Element[x,Reals],y^4==x]]
(* Output *)
Im[y]^3 Re[y]-Im[y] Re[y]^3==0
```

#### Geometric Regions

Test $\mathcal{R}_{1}\subseteq \mathcal{R}_{2}$:

```wolfram
ℛ_1=Rectangle[];
ℛ_2=Disk[{0,0},2];
```

```wolfram
Resolve[Subscript[∀, {x,y},{x,y}∈ℛ_1]{x,y}∈ℛ_2,Reals]
(* Output *)
True
```

Get conditions for $\mathcal{R}_{1}\subseteq \mathcal{R}_{2}$:

```wolfram
ℛ_1=Rectangle[];
ℛ_2=Disk[{0,0},r];
```

```wolfram
Resolve[Subscript[∀, {x,y},{x,y}∈ℛ_1]{x,y}∈ℛ_2,Reals]
(* Output *)
r>0&&r^2>=2
```

Project a cone to the $x$-$y$ plane:

```wolfram
ℛ=Cone[{{3,2,1},{1,2,3}},1];
```

```wolfram
Resolve[Exists[z,{x,y,z}∈ℛ],Reals]
(* Output *)
-12 x+2 x^2-4 y+y^2<=-21||-4 x+2 x^2-4 y+y^2<=-6||(x<=(11)/(4)&&x>=1&&-4 x+2 x^2+28 y-7 y^2>=26)
```

Plot it:

```wolfram
RegionPlot[%,{x,0,4},{y,0,4}]
```

*([Graphics])*

An implicitly defined region:

```wolfram
ℛ=ImplicitRegion[a+2 b-3 c>=1&&a b c==7,{a,b,c}];
```

```wolfram
Resolve[∃_z,{x,y,z}∈ℛx^2+y z==1,Reals]
(* Output *)
(x≠0&&7-x+x^3==0&&y≠0&&x y<0&&-21-x y+x^2 y+2 x y^2<=0)||(x≠0&&7-x+x^3==0&&y≠0&&-x y<0&&21+x y-x^2 y-2 x y^2<=0)
```

A parametrically defined region:

```wolfram
ℛ=ParametricRegion[{s^2,t^2,s t},{s,t}];
```

```wolfram
Resolve[∃_{y,z}(x^2+y^2+z^2==1&&{x,y,z}∈ℛ),Reals]
(* Output *)
x==0||x==1||(x>0&&x^2<2&&-x^2+x^3<=0&&-4 x^2+3 x^4<=0)
```

Derived regions:

```wolfram
ℛ=RegionIntersection[Ball[{0,0,0},2],Ball[{1,1,1},2]];
```

```wolfram
Resolve[∃_z{x,y,z}∈ℛ,Reals]
(* Output *)
(x+y>=(3)/(2)&&x^2+y^2<=4)||(x+y<=(1)/(2)&&-2 x+x^2-2 y+y^2<=2)||(x^2+y^2<=4&&-3 x+2 x^2-3 y+2 x y+2 y^2<=(7)/(4))||(-2 x+x^2-2 y+y^2<=2&&-3 x+2 x^2-3 y+2 x y+2 y^2<=(7)/(4))
```

Plot it:

```wolfram
RegionPlot[%,{x,-2,3},{y,-2,3}]
```

*([Graphics])*

Regions dependent on parameters:

```wolfram
ℛ_1=InfiniteLine[{{2,0},{0,t}}];
ℛ_2=Circle[];
```

```wolfram
Resolve[Subscript[∃, {x,y},{x,y}∈ℛ_1]{x,y}∈ℛ_2,Reals]
(* Output *)
t==0||-4 t^2+3 t^4<=0
```

The conditions on $t$ indicate when the line intersects the circle:

```wolfram
Graphics[{Table[{If[%,Green,Red],ℛ_1},{t,-2,2,0.5}],{Blue,ℛ_2}}]
```

*([Graphics])*

A condition for $\mathcal{R}_{1}\subseteq \mathcal{R}_{2}$:

```wolfram
ℛ_1=Disk[{a,b}];
ℛ_2=Triangle[{{0,0},{0,5},{7,0}}];
```

```wolfram
Resolve[Subscript[∀, {x,y},{x,y}∈ℛ_1]{x,y}∈ℛ_2,Reals]
(* Output *)
(a==1&&b==1&&5 a+7 b<=35&&-350 a+25 a^2-490 b+70 a b+49 b^2>=-1151)||(a==1&&b>=1&&5 a+7 b<=35&&-350 a+25 a^2-490 b+70 a b+49 b^2>=-1151)||(a>=1&&b==1&&5 a+7 b<=35&&-350 a+25 a^2-490 b+70 a b+49 b^2>=-1151)||(a>=1&&b>=1&&5 a+7 b<=35&&-350 a+25 a^2-490 b+70 a b+49 b^2>=-1151)
```

The condition tells when $\mathcal{R}_{1}\subseteq \mathcal{R}_{2}$:

```wolfram
Graphics[{{Blue,ℛ_2},Table[If[%,{Opacity[0.8],Green,ℛ_1},{Opacity[0.2],Red,ℛ_1}],{a,0,6},{b,0,4}]}]
```

*([Graphics])*

Vector variables:

```wolfram
ℛ_1=Sphere[];
ℛ_2=InfinitePlane[{{1,0,0},{0,1,0},{0,0,1}}];
```

```wolfram
Resolve[Subscript[∃, x,x∈ℛ_1](y∈ℛ_2&&x.y==1),Reals]
(* Output *)
(y≠0&&y^2>=1&&y==0&&y+y+y==1)||(y≠0&&y≠0&&y+y+y==1&&y^2+y^2+y^2==0)||(y≠0&&y^2+y^2>=1&&y==0&&y+y+y==1)||(y≠0&&y≠0&&y+y+y==1&&y^2+y^2==0)||(y≠0&&y+y+y==1&&y^2+y^2>=1&&y^2+y^2==0)||(y≠0&&y+y+y==1&&-y^4+y^6-2 y^2 y^2+3 y^4 y^2-y^4+3 y^2 y^4+y^6+y^4 y^2+2 y^2 y^2 y^2+y^4 y^2>=0)
```

### Options

#### Backsubstitution

Here the solutions for `*x*` are expressed in terms of `*y*`:

```wolfram
Resolve[Exists[z,x y+z^4-z==1&&x^2+y^2+z^3==10&&a(x^2-x y+y)^2<=0&&a>0],Reals]
(* Output *)
a>0&&((y==Root&&x==(y)/(2)+(1)/(2) Sqrt[-4 y+y^2])||(y==Root&&x==(y)/(2)+(1)/(2) Sqrt[-4 y+y^2]))
```

With `Backsubstitution->[True](https://reference.wolfram.com/language/ref/True.html)`, [Resolve](https://reference.wolfram.com/language/ref/Resolve.html) gives explicit numeric values for `*x*`:

```wolfram
Resolve[Exists[z,x y+z^4-z==1&&x^2+y^2+z^3==10&&a(x^2-x y+y)^2<=0&&a>0],Reals,Backsubstitution->True]
(* Output *)
(a>0&&y==Root&&x==(1)/(2) (Root+Sqrt[(-4+Root) Root]))||(a>0&&y==Root&&x==(1)/(2) (Root+Sqrt[(-4+Root) Root]))
```

#### Cubics

By default, [Resolve](https://reference.wolfram.com/language/ref/Resolve.html) does not use general formulas for solving cubics in radicals:

```wolfram
Resolve[Exists[x,x>a&&x^3+x+1==0],Reals]
(* Output *)
a<Root
```

With [Cubics](https://reference.wolfram.com/language/ref/Cubics.html)->[True](https://reference.wolfram.com/language/ref/True.html), [Resolve](https://reference.wolfram.com/language/ref/Resolve.html) expresses roots of cubics in terms of radicals:

```wolfram
Resolve[Exists[x,x>a&&x^3+x+1==0],Reals,Cubics->True]
(* Output *)
a<-((2)/(3 (-9+Sqrt[93])))^(1/3)+(((1)/(2) (-9+Sqrt[93]))^(1/3))/(3^(2/3))
```

#### Method

This locally sets system options in ["InequalitySolvingOptions"](https://reference.wolfram.com/language/tutorial/RealPolynomialSystems.html#343227045) and ["ReduceOptions"](https://reference.wolfram.com/language/tutorial/RealPolynomialSystems.html#227757976) groups:

```wolfram
Resolve[∃_x(a x^2+b x+c) (c x^2+b x+a)<=0,Reals,Method->{"QuadraticQE"->False}]//AbsoluteTiming
(* Output *)
{0.0463146,(c<0&&a>=(b^2)/(4 c))||c==0||(c>0&&a<=(b^2)/(4 c))}
```

By default [Resolve](https://reference.wolfram.com/language/ref/Resolve.html) uses the quadratic case of virtual substitution algorithm here:

```wolfram
Resolve[∃_x(a x^2+b x+c) (c x^2+b x+a)<=0,Reals]//AbsoluteTiming
(* Output *)
{0.0081424,(a==0&&b≠0)||(a≠0&&-b^2+4 a c<=0)||(b≠0&&c==0)||(c≠0&&-b^2+4 a c<=0)||(a==0&&b==0&&a^2+b^2+c^2==0)||(a==0&&c==0&&a^2+b^2+c^2==0)||(a==0&&b≠0&&a b+b c>0)||(b==0&&c==0&&a^2+b^2+c^2==0)||(b==0&&-b^2+4 a c<=0&&a^2+b^2+c^2==0)||(b≠0&&c==0&&a b+b c>0)||(b≠0&&-b^2+4 a c<=0&&a b+b c>0)||a c<0}
```

#### Quartics

By default, [Resolve](https://reference.wolfram.com/language/ref/Resolve.html) does not use general formulas for solving quartics in radicals:

```wolfram
Resolve[Exists[x,x>a&&x^4+x-9==0],Reals]
(* Output *)
a<Root
```

With [Quartics](https://reference.wolfram.com/language/ref/Quartics.html)->[True](https://reference.wolfram.com/language/ref/True.html), [Resolve](https://reference.wolfram.com/language/ref/Resolve.html) expresses roots of quartics in terms of radicals:

```wolfram
Resolve[Exists[x,x>a&&x^4+x-9==0],Reals,Quartics->True]
(* Output *)
a<-(1)/(2) Sqrt[-12 ((2)/(1+Sqrt[6913]))^(1/3)+((1)/(2) (1+Sqrt[6913]))^(1/3)]+(1)/(2) Sqrt(12 ((2)/(1+Sqrt[6913]))^(1/3)-((1)/(2) (1+Sqrt[6913]))^(1/3)+(2)/(Sqrt[-12 ((2)/(1+Sqrt[6913]))^(1/3)+((1)/(2) (1+Sqrt[6913]))^(1/3)]))
```

#### WorkingPrecision

This computation takes a long time due to high degrees of algebraic numbers involved:

```wolfram
Resolve[Exists[{x,y, z}, x^8 + 2y^8 + 3z^8 <=1 && x^3- 9y^3 + 7z^3-x y z>=10],Reals]//Timing
(* Output *)
{98.20299999999999,False}
```

With [WorkingPrecision](https://reference.wolfram.com/language/ref/WorkingPrecision.html)->100 you get an answer faster, but it may be incorrect:

```wolfram
Resolve[Exists[{x,y, z}, x^8 + 2y^8 + 3z^8 <=1 && x^3- 9y^3 + 7z^3>=10],Reals,WorkingPrecision->100]//Timing
(* Output *)
{26.453000000000017,False}
```

### Applications

#### Polynomials

Find conditions for a quintic to have all roots equal:

```wolfram
f[x_]:= x^5+a x^4 + b x^3 + c x^2 + d x + e
```

```wolfram
Resolve[Exists[x,f[x]==0,ForAll[y,f[y]==0, x==y]]]
(* Output *)
2 a^2-5 b==0&&a b-5 c==0&&b^2-20 d==0&&a c-10 d==0&&b c-100 e==0&&a d-25 e==0&&c^2-20 a e==0&&b d-10 a e==0&&c d-5 b e==0&&2 d^2-5 c e==0
```

Find conditions for a quadratic to be always positive:

```wolfram
Resolve[ForAll[x,a x^2+b x+c>0],Reals]
(* Output *)
(a==0&&b==0&&c>0)||(a>=0&&b==0&&c>0&&-a b^2+4 a^2 c>0)||(a>0&&-a b^2+4 a^2 c>0)
```

#### Theorem Proving

Prove the inequality between the arithmetic mean and the geometric mean:

```wolfram
Resolve[ForAll[{a,b},a>0&&b>0,(a+b)/2>=Sqrt[a b]],Reals]
(* Output *)
True
```

Prove a special case of Hölder's inequality:

```wolfram
Resolve[ForAll[{a,b,c,d},a>0&&b>0&&c>0&&d>0,a c+b d<=Sqrt[a^2+b^2]Sqrt[c^2+d^2]],Reals]
(* Output *)
True
```

Prove a special case of Minkowski's inequality:

```wolfram
Resolve[ForAll[{a,b,c,d},Sqrt[(a+c)^2+(b+d)^2]<=Sqrt[a^2+b^2]+Sqrt[c^2+d^2]],Reals]
(* Output *)
True
```

#### Geometry

The region $\mathcal{R}$ is a subset of $\mathcal{S}$, $\mathcal{R}\subseteq \mathcal{S}$ if $\forall_{\{x,y,\ldots \}}\{x,y,\ldots \}\in \mathcal{R}\Rightarrow \{x,y,\ldots \}\in \mathcal{S}$ is true. Show that [Disk](https://reference.wolfram.com/language/ref/Disk.html)[{0,0},{2,1}] is a subset of [Rectangle](https://reference.wolfram.com/language/ref/Rectangle.html)[{-2,-1},{2,1}]:

```wolfram
ℛ=Disk[{0,0},{2,1}];
𝒮=Rectangle[{-2,-1},{2,1}];
```

```wolfram
Resolve[∀_{x,y}{x,y}∈ℛ⟹{x,y}∈𝒮,Reals]
(* Output *)
True
```

Plot it:

```wolfram
Graphics[{{StandardBlue,EdgeForm[Thick],𝒮},{StandardYellow,EdgeForm[Thick],ℛ}}]
```

*([Graphics])*

Show that [Cylinder](https://reference.wolfram.com/language/ref/Cylinder.html)[]⊆[Ball](https://reference.wolfram.com/language/ref/Ball.html)[{0,0,0},2]:

```wolfram
ℛ=Cylinder[];
𝒮=Ball[{0,0,0},2];
```

```wolfram
Resolve[∀_{x,y,z}{x,y,z}∈ℛ⟹{x,y,z}∈𝒮,Reals]
(* Output *)
True
```

Plot it:

```wolfram
Graphics3D[{{Opacity[0.3],𝒮},{LightBlue,EdgeForm[Gray],ℛ}}]
```

*([Graphics3D])*

A region $\mathcal{R}$ is disjoint from $\mathcal{S}$ if $\mathcal{R}\cap \mathcal{S}=\emptyset$. Show that [Circle](https://reference.wolfram.com/language/ref/Circle.html)[{0,0},2] and [Disk](https://reference.wolfram.com/language/ref/Disk.html)[{0,0},1] are disjoint:

```wolfram
ℛ=Circle[{0,0},2];
𝒮=Disk[{0,0},1];
```

There are no points in the intersection, so they are disjoint:

```wolfram
With[{p={x,y}},Resolve[∃_p(p∈ℛ&&p∈𝒮),Reals]]
(* Output *)
False
```

Plot it:

```wolfram
Graphics[{{StandardBlue,EdgeForm[Thick],𝒮},{Thick,ℛ}}]
```

*([Graphics])*

A region $\mathcal{R}$ intersects $\mathcal{S}$ if $\mathcal{R}\cap \mathcal{S}\neq \emptyset$. Show that [Circle](https://reference.wolfram.com/language/ref/Circle.html)[{0,0},1] intersects [Disk](https://reference.wolfram.com/language/ref/Disk.html)[{1/2,0},1]:

```wolfram
ℛ=Circle[{0,0},1];
𝒮=Disk[{1/2,0},1];
```

There are points in the intersection:

```wolfram
With[{p={x,y}},Resolve[∃_p(p∈ℛ&&p∈𝒮),Reals]]
(* Output *)
True
```

Plot it:

```wolfram
Graphics[{{StandardBlue,EdgeForm[Thick],𝒮},{Thick,ℛ}}]
```

*([Graphics])*

### Properties & Relations

For fully quantified systems of equations and inequalities, [Resolve](https://reference.wolfram.com/language/ref/Resolve.html) is equivalent to [Reduce](https://reference.wolfram.com/language/ref/Reduce.html):

```wolfram
Resolve[Exists[{x,y},x^2+y^2<=3&&x^3-y^2==1]]
(* Output *)
True
```

```wolfram
Reduce[Exists[{x,y},x^2+y^2<=3&&x^3-y^2==1]]
(* Output *)
True
```

A solution instance can be found with [FindInstance](https://reference.wolfram.com/language/ref/FindInstance.html):

```wolfram
FindInstance[x^2+y^2<=3&&x^3-y^2==1,{x,y}]
(* Output *)
{{x->(5)/(4),y->(Sqrt[61])/(8)}}
```

For systems with free variables, [Resolve](https://reference.wolfram.com/language/ref/Resolve.html) may return an unsolved system:

```wolfram
Resolve[Exists[x,x^2-2a x+b<0&&2x>a+b^2-1],Reals]
(* Output *)
-2 a+3 a^2-4 b+2 b^2+2 a b^2-b^4>1||(a-b^2>-1&&-2 a+3 a^2-4 b+2 b^2+2 a b^2-b^4==1)||(a^2-b>0&&a-b^2>-1)
```

[Reduce](https://reference.wolfram.com/language/ref/Reduce.html) eliminates quantifiers and solves the resulting system:

```wolfram
Reduce[Exists[x,x^2-2a x+b<0&&2x>a+b^2-1],{a,b},Reals]
(* Output *)
(a<=Root&&Root[1+2 a-3 a^2+4 #1+(-2-2 a) #1^2+#1^4&,1]<b<Root[1+2 a-3 a^2+4 #1+(-2-2 a) #1^2+#1^4&,2])||(Root<a<=Root&&Root[1+2 a-3 a^2+4 #1+(-2-2 a) #1^2+#1^4&,1]<b<a^2)||(a>Root&&Root[1+2 a-3 a^2+4 #1+(-2-2 a) #1^2+#1^4&,1]<b<Root[1+2 a-3 a^2+4 #1+(-2-2 a) #1^2+#1^4&,2])
```

[Eliminate](https://reference.wolfram.com/language/ref/Eliminate.html) can be used to eliminate variables from systems of complex polynomial equations:

```wolfram
eqns=x z==y && x^2+x^2z^2==1;
Eliminate[eqns,z]
(* Output *)
y^2==1-x^2
```

[Resolve](https://reference.wolfram.com/language/ref/Resolve.html) gives the same equations, but may also give inequations:

```wolfram
Resolve[Exists[z,eqns]]
(* Output *)
-1+x^2+y^2==0&&x≠0
```

Find a formula description of a [TransformedRegion](https://reference.wolfram.com/language/ref/TransformedRegion.html):

```wolfram
f=Evaluate[{Indexed[#,1] Indexed[#,2],Indexed[#,1]+Indexed[#,2]}]&;
ℛ_1=Disk[{0,0}];
ℛ_2=TransformedRegion[ℛ_1,f]
(* Output *)
TransformedRegion[Disk[{0,0}],{#1 #1,#1+#1}&]
```

```wolfram
c_1=Resolve[Exists[{u,v},{u,v}∈ℛ_1,And@@Thread[{x,y}==f[{u,v}]]],Reals]
(* Output *)
(x==0&&y^2<=1)||(x≠0&&4 x-y^2<=0&&2 x+4 x^2-y^2-4 x y^2+y^4<=0)||(x≠0&&2 x-y^2==-1&&4 x-y^2<=0)||(y<0&&4 x-y^2<=0&&2 x+4 x^2-y^2-4 x y^2+y^4<=0&&x^2+4 x^3+4 x^4-2 x^2 y^2-4 x^3 y^2+x^2 y^4>=0)||(y<0&&2 x-y^2>=-1&&4 x-y^2<=0&&x+2 x^2-x y^2==0)||(y>0&&2 x-y^2>=-1&&4 x-y^2<=0&&x+2 x^2-x y^2==0)||(y>0&&4 x-y^2<=0&&2 x+4 x^2-y^2-4 x y^2+y^4<=0&&x^2+4 x^3+4 x^4-2 x^2 y^2-4 x^3 y^2+x^2 y^4>=0)
```

Compute a formula description of $\mathcal{R}_{2}$ using [RegionMember](https://reference.wolfram.com/language/ref/RegionMember.html):

```wolfram
c_2=RegionMember[ℛ_2,{x,y}]
(* Output *)
(x|y)∈Reals&&((x==0&&-1<=y<=1)||(x≠0&&2 x-y^2==-1&&4 x-y^2<=0)||(x≠0&&4 x-y^2<=0&&2 x+4 x^2-y^2-4 x y^2+y^4<=0)||(y>0&&2 x-y^2>=-1&&4 x-y^2<=0&&x+2 x^2-x y^2==0)||(y>0&&4 x-y^2<=0&&2 x+4 x^2-y^2-4 x y^2+y^4<=0&&x^2+4 x^3+4 x^4-2 x^2 y^2-4 x^3 y^2+x^2 y^4>=0)||(y<0&&2 x-y^2>=-1&&4 x-y^2<=0&&x+2 x^2-x y^2==0)||(y<0&&4 x-y^2<=0&&2 x+4 x^2-y^2-4 x y^2+y^4<=0&&x^2+4 x^3+4 x^4-2 x^2 y^2-4 x^3 y^2+x^2 y^4>=0))
```

Check that the formulas are equivalent:

```wolfram
Resolve[∀_{x,y}(c_1⟺c_2),Reals]
(* Output *)
True
```

[Resolve](https://reference.wolfram.com/language/ref/Resolve.html) shows that the polynomial $f$ is non-negative:

```wolfram
f=x^2+2y^2+3 z^2-x y-x z-y z;
```

```wolfram
Resolve[ForAll[{x,y,z},f>=0],Reals]
(* Output *)
True
```

Use [PolynomialSumOfSquaresList](https://reference.wolfram.com/language/ref/PolynomialSumOfSquaresList.html) to represent $f$ as a sum of squares:

```wolfram
PolynomialSumOfSquaresList[f, {x,y,z}]
(* Output *)
{-(x)/(2 Sqrt[3])-(y)/(2 Sqrt[3])+Sqrt[3] z,-(7 x)/(2 Sqrt[69])+(1)/(2) Sqrt[(23)/(3)] y,Sqrt[(17)/(23)] x}
```

```wolfram
f-%.%//Expand
(* Output *)
0
```

The Motzkin polynomial is non-negative, but is not a sum of squares:

```wolfram
m=x^4 y^2+x^2 y^4-3 x^2 y^2+1;
```

```wolfram
Resolve[ForAll[{x,y},m>=0],Reals]
(* Output *)
True
```

```wolfram
PolynomialSumOfSquaresList[m, {x,y}]
(* Output *)
PolynomialSumOfSquaresList[1-3 x^2 y^2+x^4 y^2+x^2 y^4,{x,y}]
```

### Possible Issues

Because `*x*` appears in an inequality, it is assumed to be real:

```wolfram
Resolve[Exists[x,x^2<-1]]
(* Output *)
False
```

This allows complex values of `*x*` for which both sides of the inequality are real:

```wolfram
Resolve[Exists[x,x^2<-1],Complexes]
(* Output *)
True
```

## Tech Notes ▪Quantifiers ▪Complex Polynomial Systems ▪Real Polynomial Systems ▪Implementation notes: Algebra and Calculus

## Related Guides ▪Polynomial Systems ▪Theorem Proving ▪Polynomial Algebra ▪Boolean Computation ▪Logic & Boolean Algebra ▪Solvers over Regions ▪Operations on Sets ▪Finite Fields

## History Introduced in 2003 (5.0) | Updated in 2014 (10.0) ▪ 2024 (14.0)
