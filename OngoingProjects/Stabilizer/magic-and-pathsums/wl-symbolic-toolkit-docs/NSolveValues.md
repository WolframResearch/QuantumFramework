# NSolveValues | [SpanFromLeft]

> [NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html)[*expr*,*vars*]  — attempts to find numerical approximations to the values of `*vars*` determined by the solutions of the system `*expr*`.
> [NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html)[*expr*,*vars*,[Reals](https://reference.wolfram.com/language/ref/Reals.html)] — finds solutions over the domain of real numbers.

## Details and Options

The system `*expr*` can be any logical combination of:

*lhs*==*rhs* | equations
*lhs*!=*rhs* | inequations
`*lhs*>*rhs*` or `*lhs*>=*rhs*` | inequalities
*expr*∈*dom* | domain specifications
{*x*,*y*,…}∈*reg* | region specification
[ForAll](https://reference.wolfram.com/language/ref/ForAll.html)[*x*,*cond*,*expr*] | universal quantifiers
[Exists](https://reference.wolfram.com/language/ref/Exists.html)[*x*,*cond*,*expr*] | existential quantifiers

[NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html)[{*expr*_1,*expr*_2,…},*vars*] is equivalent to [NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html)[*expr*_1&&*expr*_2&&…,*vars*].

A single variable or a list of variables can be specified.

If a single variable is specified, the result is a list of values of the variable for which `*expr*` is [True](https://reference.wolfram.com/language/ref/True.html).

If a list of variables is specified, the result is a list of lists of values for the variables for which `*expr*` is [True](https://reference.wolfram.com/language/ref/True.html).

When a single variable is specified and a particular root of an equation has multiplicity greater than one, [NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html) gives several copies of the corresponding solution.

[NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html)[*expr*,*vars*] assumes by default that quantities appearing algebraically in inequalities are real, while all other quantities are complex.

In [NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html)[*expr*,*vars*,[Reals](https://reference.wolfram.com/language/ref/Reals.html)] all variables, parameters, constants and function values are restricted to be real.

[NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html)[*expr*&&*vars*∈[Reals](https://reference.wolfram.com/language/ref/Reals.html),*vars*,[Complexes](https://reference.wolfram.com/language/ref/Complexes.html)] solves for real values of variables, but function values are allowed to be complex.

[NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html)[…,*x*∈*reg*,[Reals](https://reference.wolfram.com/language/ref/Reals.html)] constrains `*x*` to be in the region `*reg*`. The different coordinates for `*x*` can be referred to using [Indexed](https://reference.wolfram.com/language/ref/Indexed.html)[*x*,*i*].

[NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html) deals primarily with linear and polynomial equations.

The following options can be given:

| [MaxRoots](https://reference.wolfram.com/language/ref/MaxRoots.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | maximum number of roots returned |
| --- | --- | --- |
| [Method](https://reference.wolfram.com/language/ref/Method.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | what method should be used |
| [RandomSeeding](https://reference.wolfram.com/language/ref/RandomSeeding.html) | 1234 | the seeding of pseudorandom generators |
| [VerifySolutions](https://reference.wolfram.com/language/ref/VerifySolutions.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | whether to verify solutions |
| [WorkingPrecision](https://reference.wolfram.com/language/ref/WorkingPrecision.html) | [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) | precision to be used in computations |

Possible [Method](https://reference.wolfram.com/language/ref/Method.html) settings include `"EndomorphismMatrix"`,  `"Homotopy"`, `"Monodromy"`, and `"Symbolic"`.

## Examples

### Basic Examples

Approximate solutions to a polynomial equation:

```wolfram
NSolveValues[x^5-2x+3==0,x]
(* Output *)
{-1.4236058485523315,-0.24672925691056408-1.3208163474502472 ⅈ,-0.24672925691056408+1.3208163474502472 ⅈ,0.9585321811867297-0.4984277790318459 ⅈ,0.9585321811867297+0.4984277790318459 ⅈ}
```

Approximate real solutions to a polynomial equation:

```wolfram
NSolveValues[x^5-2x+3==0,x,Reals]
(* Output *)
{-1.4236058485523317}
```

Find three approximate solutions to a polynomial equation of a high degree:

```wolfram
NSolveValues[x^1000000-2x+3==0,x,MaxRoots->3]
(* Output *)
{0.9998709946689086+0.016062240400403986 ⅈ,0.9998708937279452+0.016068522762298948 ⅈ,0.999870792747509+0.01607480512355958 ⅈ}
```

Approximate solutions to a system of polynomial equations:

```wolfram
NSolveValues[{x^2+y^3==1,2x+3y==4},{x,y}]
(* Output *)
{{7.9364095819855995,-3.9576063879903995},{0.7192952090072029-0.25567875198601236 ⅈ,0.8538031939951981+0.17045250132400824 ⅈ},{0.7192952090072029+0.25567875198601236 ⅈ,0.8538031939951981-0.17045250132400824 ⅈ}}
```

Approximate real solutions to a system of polynomial equations:

```wolfram
NSolveValues[{x^2+y^3==1,2x+3y==4},{x,y},Reals]
(* Output *)
{{7.9364095819855995,-3.9576063879903995}}
```

Solve equations in a geometric region:

```wolfram
NSolveValues[{x,y}∈InfiniteLine[{{0,0},{2,1}}]&&{x,y}∈Circle[],{x,y}]
(* Output *)
{{-0.8944271909999159,-0.4472135954999579},{0.894427190999916,0.447213595499958}}
```

```wolfram
Graphics[{{Blue,InfiniteLine[{{0,0},{2,1}}],Circle[]},{Red,Point[%]}}]
```

*([Graphics])*

### Scope

#### Complex Equations in One Variable

Univariate polynomial equations:

```wolfram
NSolveValues[x^5-2x+17==0,x]
(* Output *)
{-1.8325135554153524,-0.48796127686225305-1.7207007442321953 ⅈ,-0.48796127686225305+1.7207007442321953 ⅈ,1.4042180545699292-0.963419154203649 ⅈ,1.4042180545699292+0.963419154203649 ⅈ}
```

Polynomial equations with inexact coefficients:

```wolfram
NSolveValues[x^3-1.234x+5.678==0,x]
(* Output *)
{-2.013460015723761,1.0067300078618806-1.3440669351593868 ⅈ,1.0067300078618806+1.3440669351593868 ⅈ}
```

Polynomial equations with multiple roots:

```wolfram
NSolveValues[(x^2-1)(x^4-1)==0,x]
(* Output *)
{-1.,-1.,0.-1. ⅈ,0.+1. ⅈ,1.,1.}
```

Find five roots of a polynomial of a high degree:

```wolfram
NSolveValues[x^1234567+9x^2+7x-1==0,x,MaxRoots->5]
(* Output *)
{0.1233080254051604,0.7427959949631621-0.6695210601471077 ⅈ,0.7427994024079472-0.6695172797598118 ⅈ,0.7428028098334925-0.6695134993551742 ⅈ,0.7428062172397976-0.6695097189331949 ⅈ}
```

Algebraic equations:

```wolfram
NSolveValues[Sqrt[x]+3x^(1/3)==5,x]
(* Output *)
{1.8086326695171464}
```

Transcendental equations:

```wolfram
NSolveValues[Sin[x]==1/3,x]
(* Output *)
{0.3398369094541219,2.8017557441356713,-3.4814295630439154,6.623022216633708,-5.943348397725464,9.084941051315258,-9.764614870223502,12.906207523813295,-12.226533704905052,15.368126358494845,-16.047800177403087,19.18939283099288,-18.509719012084638,21.65131166567443,-22.330985484582676,25.47257813817247}
```

Find all solutions:

```wolfram
NSolveValues[Sin[x]==1/3,x,MaxRoots->Infinity]
(* Output *)
{1. (0.3398369094541219+6.283185307179586 1),1. (2.8017557441356713+6.283185307179586 1)}
```

Specify the number of solutions returned:

```wolfram
NSolveValues[5==x 2^(x^2),x,MaxRoots->3]
(* Output *)
{1.3675873372820502,2.1052847120151+1.9011827939980588 ⅈ,2.9400810435860656+2.892459320555923 ⅈ}
```

Univariate elementary function equations over bounded regions:

```wolfram
NSolveValues[Sin[E^x]-Cos[2 x]==1&&-1<=Re[x]<=1&&-1<=Im[x]<=1,x]
(* Output *)
{0.9152269083602887+0.20167483871637634 ⅈ,0.9152269083602887-0.20167483871637634 ⅈ}
```

Univariate holomorphic function equations over bounded regions:

```wolfram
NSolveValues[Gamma[x]-Log[x]==I/2&&Abs[x-2]<3/2,x]
(* Output *)
{1.49328747884441-0.9848501122597985 ⅈ,2.671009023704611+0.8195706035753547 ⅈ}
```

Equation with a purely imaginary period over a vertical stripe in the complex plane:

```wolfram
NSolveValues[Cos[Exp[x]]==3 Exp[-x]+1&&0<=Re[x]<=1,x]
(* Output *)
{0.5465407839974525+1.1111486424967023 ⅈ,0.7185533775573053+3.141592653589793 ⅈ,0.5465407839974525+5.172036664682884 ⅈ,0.5465407839974525-1.1111486424967023 ⅈ}
```

Find all solutions:

```wolfram
NSolveValues[Cos[Exp[x]]==3 Exp[-x]+1&&0<=Re[x]<=1,x,MaxRoots->Infinity]
(* Output *)
{(0.5465407839974525+5.172036664682884 ⅈ)+(0.+6.283185307179586 ⅈ) 1,(0.5465407839974525+1.1111486424967023 ⅈ)+(0.+6.283185307179586 ⅈ) 1,(0.7185533775573053+3.141592653589793 ⅈ)+(0.+6.283185307179586 ⅈ) 1}
```

Unrestricted transcendental function equations:

```wolfram
NSolveValues[BesselJ[2,x]x-2^x==3,x]
(* Output *)
{-6.890303031977378-0.9074308637760139 ⅈ,-1.1046221329230759-2.321653496321717 ⅈ,2.600791627339534-2.4846514034482636 ⅈ,2.600791627339534+2.4846514034482636 ⅈ,-1.1046221329230759+2.321653496321717 ⅈ,-6.890303031977378+0.9074308637760139 ⅈ,-13.248469225030814+0.2392314083387334 ⅈ,-19.001428003356665,-13.24846922503092-0.23923140833781917 ⅈ,6.990042059422345-4.791810246920608 ⅈ,-20.122201476016354,-25.119401219070067,-26.603052974984244,-31.306673474721176,-33.00254507825711,-37.524098509537616}
```

#### Systems of Complex Equations in Several Variables

Systems of linear equations:

```wolfram
NSolveValues[{x+2y+3z==4,3x+4y+5z==6,7x+9y+8z==10},{x,y,z}]
(* Output *)
{{-1.0000000000000009,1.0000000000000007,0.9999999999999998}}
```

Linear equations with inexact coefficients:

```wolfram
NSolveValues[1.23 x+2.34 y+3.45 z==4.56 &&5.67x+6.78y+7.89z==8.9&&9.01x+0.12y+1.23z==2.34,{x,y,z}]
(* Output *)
{{-0.0049999999999998665,-1.0600024348672987,2.0424799123447763}}
```

Underdetermined systems of linear equations:

```wolfram
NSolveValues[x+2y+3z==4&&3x+4y+5z==6&&6x+7y+8z==9,{x,y,z}]
(* Output *)
{{-2.+1. z,3.-2. z,z}}
```

Linear equations with no solutions:

```wolfram
NSolveValues[x+2y+3z==4&&3x+4y+5z==6&&6x+7y+8z==0,{x,y,z}]
(* Output *)
{}
```

Systems of polynomial equations:

```wolfram
NSolveValues[x^2+y^2==1&&x+2y==3,{x,y}]
(* Output *)
{{0.5999999999999996-0.8000000000000003 ⅈ,1.2000000000000002+0.40000000000000013 ⅈ},{0.5999999999999996+0.8000000000000003 ⅈ,1.2000000000000002-0.40000000000000013 ⅈ}}
```

Find five out of a trillion roots of a polynomial system:

```wolfram
NSolveValues[x^10000==y^2+2y+1&&y^10000==z^2+2z+1&& z^10000==x^2+2 x+1,{x,y,z},MaxRoots->5]
(* Output *)
{{-0.9986965182709251-0.0005377715054718102 ⅈ,-0.9986737839010562-0.0006393697568540476 ⅈ,-0.9986880333934715+0.00007815562728396285 ⅈ},{-0.9987330143392446-0.0011090523316215804 ⅈ,-0.9986808891052007-0.0011826143719502681 ⅈ,-0.9987233601689871-0.0004838946994072922 ⅈ},{-0.9987720466249359-0.00169588465563482 ⅈ,-0.998715555793081-0.0017401411636269817 ⅈ,-0.998766434037532-0.0010665017406316093 ⅈ},{-0.9988100033667913-0.002294546288706712 ⅈ,-0.9987598403740191-0.0023205557171392115 ⅈ,-0.9988076990039676-0.001664497333761613 ⅈ},{-0.9988474335970712-0.0029040879240262926 ⅈ,-0.9987544212332762-0.002946179253627027 ⅈ,-0.998845611637445-0.0016444595649509639 ⅈ}}
```

Underdetermined systems of polynomial equations:

```wolfram
NSolveValues[x^2+y^2+z^2==1&&x y z==2,{x,y,z}]
(* Output *)
NSolveValues
(* Output *)
{{0.3257881805145985-1.043807816388064 ⅈ,0.3226483801300118+1.0204471405784412 ⅈ,1.7089901266964773+0.006328083599261027 ⅈ},{0.3257881805145985+1.043807816388064 ⅈ,0.3226483801300118-1.0204471405784412 ⅈ,1.7089901266964773-0.006328083599261027 ⅈ},{0.439870637792798-0.9660453744355169 ⅈ,-1.6309961781361495+0.03582094291924612 ⅈ,-0.4554102597234386-1.0613700628567102 ⅈ},{0.439870637792798+0.9660453744355169 ⅈ,-1.6309961781361495-0.03582094291924612 ⅈ,-0.4554102597234386+1.0613700628567102 ⅈ},{-1.647999877905746+0.04243749789045728 ⅈ,0.41639313921064797-0.9579706895720473 ⅈ,-0.4354162981085605-1.0767383216507473 ⅈ},{-1.647999877905746-0.04243749789045728 ⅈ,0.41639313921064797+0.9579706895720473 ⅈ,-0.4354162981085605+1.0767383216507473 ⅈ}}
```

Algebraic equations:

```wolfram
NSolveValues[Sqrt[x y]==Sqrt[x+y]&&Sqrt[x]-y^(1/3)==1,{x,y}]
(* Output *)
{{4.370091649877744,1.2967278352908524},{4.370091649877743,1.2967278352908522}}
```

Transcendental equations:

```wolfram
NSolveValues[Sin[x+y]==x y+1&&Cos[x-y]==AiryAi[x y]+2,{x,y}]
(* Output *)
{{0.5979497570909443-1.2705733856994925 ⅈ,0.530877134382463+0.16012535979852152 ⅈ},{5.490657859723929-1.757840928114804 ⅈ,-0.7906866811232162-0.24491347288881815 ⅈ},{11.757817159423988-1.2314777055486594 ⅈ,-0.22233802508658002-0.10031850110086371 ⅈ},{15.578720088747158+2.233921600120761 ⅈ,-0.05227456682338649-0.22388368971659106 ⅈ},{13.037717826107372+1.4274948874781008 ⅈ,0.03389135111634696+0.15226193453068027 ⅈ},{17.015038669364753-1.93816303872407 ⅈ,-0.23793968348802577+0.06290055911159166 ⅈ},{23.62582244782386-2.504257662126576 ⅈ,-0.29255015584005417+0.025787702299229692 ⅈ},{29.99366186875182-2.8113751096293127 ⅈ,-0.30570569875364867+0.01398437381762144 ⅈ},{36.313388651461-3.0219339296123495 ⅈ,-0.3071324225793418+0.008487243071173593 ⅈ},{42.615112210550016-3.1826101903246196 ⅈ,-0.30426503191307275+0.005496177459698996 ⅈ},{48.90831510310423-3.3126319039213845 ⅈ,-0.29969722806444354+0.0037092933589936874 ⅈ},{55.196947175698-3.421844852833278 ⅈ,-0.2944966076677351+0.002570839779905514 ⅈ},{61.482928139131644-3.5159883873861335 ⅈ,-0.2891434524926428+0.0018102551783918247 ⅈ},{67.76728892196952-3.598708950634197 ⅈ,-0.28386287336217364+0.0012832373814148207 ⅈ},{74.05062354114726-3.672470468157755 ⅈ,-0.27875996872825565+0.0009073813612056942 ⅈ},{80.33329266804478-3.739016921688897 ⅈ,-0.27388038773796297+0.0006330387880293382 ⅈ}}
```

Specify the number of solutions returned:

```wolfram
NSolveValues[3^x^2==7&&x^2-y^3==4,{x,y},MaxRoots->5]
(* Output *)
{{1.3308808170386341,-1.3062336217467332},{1.3308808170386341,0.6531168108733666-1.1312314997100241 ⅈ},{1.3308808170386341,0.6531168108733666+1.1312314997100241 ⅈ},{1.9695741298473832+1.4518879102061062 ⅈ,-1.6866143945936043+0.7125611497138453 ⅈ},{1.9695741298473832+1.4518879102061062 ⅈ,0.22621113989476538-1.8169314869634954 ⅈ}}
```

#### Real Equations in One Variable

Polynomial equations:

```wolfram
NSolveValues[x^5-2x+1==0,x,Reals]
(* Output *)
{-1.2906488013467097,0.5187900636758842,1.}
```

Polynomial equations with multiple roots:

```wolfram
NSolveValues[(x^2-1)(x^4-1)==0,x,Reals]
(* Output *)
{-1.,-1.,1.,1.}
```

Algebraic equations:

```wolfram
NSolveValues[Sqrt[x]+3x^(1/3)==5,x,Reals]
(* Output *)
{1.8086326695171466}
```

Piecewise equations:

```wolfram
NSolveValues[RealAbs[x]^2-x+UnitStep[x]==9,x,Reals]
(* Output *)
{3.3722813232690143,-2.5413812651491097}
```

Transcendental equations, solvable using inverse functions:

```wolfram
NSolveValues[E^x-x==7,x,Reals]
(* Output *)
{-6.999087285366495,2.22154230138681}
```

```wolfram
NSolveValues[ (27^(2x-1))^(1/x)==Sqrt[9^(2x-1)], x , Reals]
(* Output *)
{0.5,3.}
```

Transcendental equations, solvable using special function zeros:

```wolfram
NSolveValues[AiryBi[1-x^2]==0&&2<x<3,x,Reals]
(* Output *)
{2.066662358208605,2.4146920800926184,2.677657955809564,2.8942636506316703}
```

Transcendental inequalities, solvable using special function zeros:

```wolfram
NSolveValues[900<AiryAiZero[2t+1]^2<1000,t,Reals]
(* Output *)
{17.5,18.}
```

Exp-log equations:

```wolfram
NSolveValues[E^(2E^x)-Log[x^2+1]-20x==11,x,Reals]
(* Output *)
{-0.3516271290675965,0.38318251567919076}
```

High-degree sparse polynomial equations:

```wolfram
NSolveValues[x^1000000-2x^777777+3x^12345+9x^67-10==0,x,Reals]
(* Output *)
{-1.0000023076691011,0.9999150468764029,1.0000018393911132,1.0000020798897298}
```

Algebraic equations involving high-degree radicals:

```wolfram
NSolveValues[2x^(123451/67890)-x^2+4Sqrt[x]-4x-9/8==0,x,Reals]
(* Output *)
{0.3002931677856876,0.6146644976191242,4.670464740861935,20.75164149004275}
```

Equations involving irrational real powers:

```wolfram
NSolveValues[x^Pi-x^x^Sqrt[2]-Sqrt[3]x+2^(1/3)==0,x,Reals]
(* Output *)
{0.2650798144381758,1.403117910902082,2.0492718931323264}
```

Tame elementary function equations:

```wolfram
NSolveValues[10Sin[Tan[E^-x^2]]-x==3,x,Reals]
(* Output *)
{-2.99875666305434,-1.3422490524391641,0.973157151243307}
```

Elementary function equations in bounded intervals:

```wolfram
NSolveValues[2Sin[Exp[x]]-Cos[Pi x]==3/2&&-1<x<1,x,Reals]
(* Output *)
{-0.6764666779144783,0.342857185985482}
```

Holomorphic function equations in bounded intervals:

```wolfram
NSolveValues[Cos[x]-BesselJ[5,x]==1/2&&0<=x<=10,x]
(* Output *)
{1.0468364058841102,5.712182959861007,6.814812041680406}
```

#### Systems of Real Equations and Inequalities in Several Variables

Linear systems:

```wolfram
NSolveValues[2 x+3y-5z==1&&3x-4y+7z==3&&x+y-z==8,{x,y,z},Reals]
(* Output *)
{{0.09090909090909073,19.36363636363637,11.454545454545459}}
```

Polynomial systems:

```wolfram
NSolveValues[x y==z^2-x&&x y z==2&&x^2+y^2+z^2==5,{x,y,z},Reals]
(* Output *)
{{1.0460001409568287,1.247144736559103,1.5331385166354763},{1.1530813961467938,1.1115709908406748,1.560388806093419}}
```

Quantified polynomial systems:

```wolfram
NSolveValues[Exists[x, x^2+a x+b==0&&2 x+a==0&&a^2+x b^2==1],{a,b},Reals]
(* Output *)
{{3.0572858304191923,2.3367491622204883},{1.016844173210528,0.2584930181480509},{-0.985377547246573,0.2427422276544179}}
```

```wolfram
NSolveValues[ForAll[x,Exists[y,a x^2+b y^2-3y==1&&y<0&&a y-y==b+1]],{a,b},Reals]
(* Output *)
{{0,-0.6388969194713526},{1.,-1.}}
```

Algebraic systems:

```wolfram
NSolveValues[Sqrt[x+2y]-3x+4y==5&&x+y^(1/3)==2,{x,y},Reals]
(* Output *)
{{0.8749864865272123,1.4238794345833141}}
```

Piecewise systems:

```wolfram
NSolveValues[Max[x^2,1-y^3]==Min[y^2,x]+7y&&x+UnitStep[4-y]==5,{x,y},Reals]
(* Output *)
{{4.,1.8150729063673248}}
```

Transcendental systems, solvable using inverse functions:

```wolfram
NSolveValues[Sin[x+y]==1/2&&E^x-y==1&&x^2<=10,{x,y},Reals]
(* Output *)
{{-2.730385576667298,-0.9348058525207941},{0.24542910765918458,0.2781696679391143},{0.9727596256144202,1.6452342523770742},{1.7939038704158041,5.012880212362081},{2.059458309698381,6.8417208754727},{2.454093108207093,10.63587628175038},{2.608307710444274,12.576056781906393},{2.8627928133479914,16.510361883789066},{2.9702747164747634,18.49727508305549},{3.156973465012508,22.499366539304138}}
```

```wolfram
NSolveValues[ 3^x-2^2y==77 && Sqrt[3^x]-2^y==7, {x, y}, Reals]
(* Output *)
{{4.,1.}}
```

Systems exp-log in the first variable and polynomial in the other variables:

```wolfram
NSolveValues[E^x y^3+Log[x]y==1&&x y+E^x/x==10,{x,y},Reals]
(* Output *)
{{0.11037344375004365,-1.0638862099946733},{0.11109415744014833,-0.5314045787609156},{0.11415282270288646,1.5812311488285882},{3.4373769561435754,0.27662750335046726}}
```

Quantified system:

```wolfram
NSolveValues[Exists[a,a x^2+Sinh[x^2+1]a^2==1&&x^2-a^2==1],x,Reals]
(* Output *)
{-1.1602987294098774,1.1602987294098774,-1.066777618050499,1.066777618050499}
```

Systems elementary and bounded in the first variable and polynomial in the other variables:

```wolfram
NSolveValues[Sin[x-Cos[x]] y^3-x==1&&x^2+y^2==1,{x,y},Reals]
(* Output *)
{{-1.,0},{-0.8383977531025116,-0.5450589028652408},{-0.12223534929360626,-0.9925011432653718}}
```

Quantified system:

```wolfram
NSolveValues[Exists[y,y^3-Cos[x] y+2 x^2 Sin[x^2-1]==0&&x^2+y^2==2],x,Reals]
(* Output *)
{-1.0722719806684207,1.0722719806684207}
```

Systems holomorphic and bounded in the first variable and polynomial in the other variables:

```wolfram
NSolveValues[y^3-BesselJ[2,x+2] y-y-3 x==-2&&x^2<=1&&x y-x^2+y^2==3,{x,y},Reals]
(* Output *)
{{-0.14480600520889547,-1.6671978061820676}}
```

Quantified system:

```wolfram
NSolveValues[Exists[y,y^4-Gamma[x+2] y-y-3 ArcSin[x/3]==1&&x^2+y^3==1&&0<x<2],x,Reals]
(* Output *)
{1.1119416736295196}
```

#### Systems with Mixed Variable Domains

Mixed real and complex variables:

```wolfram
NSolveValues[x^2+y^2==1&&x^4+ y^4==2&&Element[x,Reals],{x,y}]
(* Output *)
{{1.1687708944803692,0.-0.6050003337060545 ⅈ},{1.1687708944803692,0.+0.6050003337060545 ⅈ},{-1.16877089448037,0.-0.6050003337060559 ⅈ},{-1.16877089448037,0.+0.6050003337060559 ⅈ}}
```

#### Geometric Regions

Solve over special regions in 2D:

```wolfram
ℛ_1=Circle[];
ℛ_2=Line[{{-2,1},{1,-2}}];
```

```wolfram
NSolveValues[{x,y}∈ℛ_1&&{x,y}∈ℛ_2,{x,y}]
(* Output *)
{{-1.,0.},{0.,-1.0000000000000002}}
```

Plot it:

```wolfram
Graphics[{{Blue,ℛ_1,ℛ_2},{Red,Point[%]}}]
```

*([Graphics])*

Solve over special regions in 3D:

```wolfram
ℛ_1=Sphere[];
ℛ_2=InfinitePlane[{{0,0,0},{0,1,0},{1,0,1}}];
```

```wolfram
NSolveValues[2 x y==z^2&&{x,y,z}∈ℛ_1&&{x,y,z}∈ℛ_2,{x,y,z},Reals]
(* Output *)
{{0.,-0.9999999999999997,0.},{0.,1.0000000000000002,0.},{0.6666666666666665,0.333333333333334,0.6666666666666665},{-0.6666666666666665,-0.33333333333333276,-0.6666666666666665}}
```

Plot it:

```wolfram
Show[{ContourPlot3D[2 x y==z^2,{x,-1.2,1.2},{y,-1.2,1.2},{z,-1.2,1.2},Mesh->None,ContourStyle->Opacity[0.5]],Graphics3D[{{Opacity[0.5],Green,ℛ_1},{Opacity[0.5],Yellow,ℛ_2},{PointSize[Large],Red,Point[%]}}]}]
```

*([Graphics3D])*

An implicitly defined region:

```wolfram
ℛ=ImplicitRegion[a+2 b-3 c>=1&&a b c==7,{a,b,c}];
```

```wolfram
NSolveValues[y^2+x z==1&&x+y==z^2&&{x,y,z}∈ℛ,{x,y,z},Reals]
(* Output *)
{{3.191502458676469,-2.086745339882667,-1.0510742689238448}}
```

A parametrically defined region:

```wolfram
ℛ=ParametricRegion[{s+t,s-t,s t},{s,t}];
```

```wolfram
NSolveValues[x y==z&&x+2 y+3 z==1&&{x,y,z}∈ℛ,{x,y,z},Reals]
(* Output *)
{{-0.38196601125010654,1.618033988749893,-0.6180339887498931},{-0.20601132958329937,0.8726779962499646,-0.1797815543055432},{-2.6180339887498754,-0.6180339887498878,1.6180339887498836},{0.5393446629166322,0.12732200375003497,0.06867044319443263}}
```

Derived regions:

```wolfram
ℛ_1=Disk[{0,0},2];
ℛ_2=Circle[{1,1},2];
ℛ_3=RegionIntersection[ℛ_1,ℛ_2];
```

```wolfram
NSolveValues[x^2==x y+1&&{x,y}∈ℛ_3,{x,y},Reals]
(* Output *)
{{0.6278006311513903,-0.9650617369000644},{-0.8677206193059607,0.28472404750378727}}
```

Plot it:

```wolfram
Show[{ContourPlot[x^2==x y+1,{x,-2,3},{y,-2,3}],Graphics[{{Opacity[0.5],Yellow,ℛ_1},{Green,ℛ_2},{Red,Point[%]}}]}]
```

*([Graphics])*

Regions dependent on parameters:

```wolfram
ℛ_1=Circle[];
ℛ_2=InfiniteLine[{{-2,0},{0,t}}];
ℛ_3=Circle[{1,1},{2t,t}];
ℛ_4=RegionIntersection[ℛ_1,ℛ_2,ℛ_3];
```

```wolfram
NSolveValues[{x,y}∈ℛ_4,{t,x,y},Reals]
(* Output *)
{{0.4207205988612039,0.8070444370749472,0.5904907082980915},{0.9911578098013685,-0.8062242239561307,0.5916100917887853}}
```

Plot it:

```wolfram
Graphics[{{Blue,ℛ_1,ℛ_2,ℛ_3},{Red,Point[{x,y}]}}/.(Thread[{t,x,y}->#]&/@%)]
```

*([Graphics])*

Find values of parameters $a$, $b$ and $r$ for which the circles contain the given points:

```wolfram
ℛ_1=Circle[{a,b},r];
ℛ_2=Circle[{a+1,b},r];
ℛ_3=RegionIntersection[ℛ_1,ℛ_2];
```

```wolfram
NSolveValues[({0,1}|{0,-1})∈ℛ_3,{a,b,r},Reals]
(* Output *)
{{-0.5,0.,1.118033988749895}}
```

Plot it:

```wolfram
Show[{Graphics[{{Blue,ℛ_1},{Green,ℛ_2},{Red,Point[{{0,1},{0,-1}}]}}/.Thread[{a,b,r}->%[[1]]]]}]
```

*([Graphics])*

Use $x \in \mathcal{R}$ to specify that $x$ is a vector in $\mathbb{R}^{2}$:

```wolfram
ℛ=RegionIntersection[Circle[],Line[{{-2,-1},{1,2}}]];
```

```wolfram
NSolveValues[x∈ℛ,x]
(* Output *)
{{-1.,0.},{0.,1.}}
```

In this case, $x$ is a vector in $\mathbb{R}^{3}$:

```wolfram
ℛ=Sphere[];
```

```wolfram
NSolveValues[x.{1,2,3}==0&&x.{-3,-2,-1}==0&&x∈ℛ,x]
(* Output *)
{{-0.408248290463863,0.816496580927726,-0.408248290463863},{0.408248290463863,-0.816496580927726,0.408248290463863}}
```

### Generalizations & Extensions

Working precision can be given as the last argument:

```wolfram
NSolveValues[x^5-x+2,x,Reals,30]
(* Output *)
{-1.2671683045421243172528914279777347048339936021970811826775}
```

### Options

#### MaxRoots

Find $3$ out of $12345$ roots of a polynomial:

```wolfram
NSolveValues[x^12345+x+1==0,x,MaxRoots->3]
(* Output *)
{-0.999399332443321,-0.9994161447953653-0.00045506980321965355 ⅈ,-0.9994467756614861-0.0009335099314335696 ⅈ}
```

Find $3$ out of $1000000000$ roots of a polynomial system:

```wolfram
NSolveValues[x^1000==y+1&&y^1000==z+1&&z^1000==x+1,{x,y,z},MaxRoots->3]
(* Output *)
{{-0.9946935031356067-0.006417854589319048 ⅈ,-0.995078441016325+0.0008388304820450095 ⅈ,-0.995208728843315-0.005377442735937451 ⅈ},{-0.9949406114995196-0.011653666589919393 ⅈ,-0.9955893269596023-0.005061292258790306 ⅈ,-0.9955789655530017-0.011355235526799212 ⅈ},{-0.9954241601716333-0.017549317351262986 ⅈ,-0.995925125175366-0.011183874167711702 ⅈ,-0.9958450775469632-0.017462760235084338 ⅈ}}
```

Find $5$ out of infinitely many solutions of a transcendental system:

```wolfram
NSolveValues[Sin[x+y]==2&&Cos[x-y]==3,{x,y},MaxRoots->5]
(* Output *)
{{0.7853981633974483-1.5398525354819514 ⅈ,0.7853981633974483+0.22289463855713468 ⅈ},{3.9269908169872414-1.5398525354819514 ⅈ,-2.356194490192345+0.22289463855713468 ⅈ},{7.0685834705770345-1.5398525354819514 ⅈ,-5.497787143782138+0.22289463855713468 ⅈ},{10.210176124166829-1.5398525354819514 ⅈ,-8.639379797371932+0.22289463855713468 ⅈ},{13.351768777756622-1.5398525354819514 ⅈ,-11.780972450961725+0.22289463855713468 ⅈ}}
```

With the default [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) setting [NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html) may not give all solutions:

```wolfram
NSolveValues[Sin[x]==1/2,x]
(* Output *)
{0.5235987755982989,2.6179938779914944,-3.6651914291880923,6.806784082777885,-5.759586531581288,8.901179185171081,-9.94837673636768,13.089969389957473,-12.042771838760874,15.184364492350667,-16.231562043547264,19.373154697137057,-18.32595714594046,21.467549799530254,-22.51474735072685,25.656340004316643}
```

With [MaxRoots](https://reference.wolfram.com/language/ref/MaxRoots.html)->[Infinity](https://reference.wolfram.com/language/ref/Infinity.html),  [NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html) attempts to find all solutions:

```wolfram
NSolveValues[Sin[x]==1/2,x,MaxRoots->Infinity]
(* Output *)
{1. (0.5235987755982989+6.283185307179586 1),1. (2.617993877991494+6.283185307179586 1)}
```

#### Method

Solve a square polynomial system using the automatically chosen method:

```wolfram
NSolveValues[x^2+2y^2==3&&x^3-4x y==5,{x,y}]
(* Output *)
{{-0.3692815807235537-3.2866081970081487 ⅈ,-2.624155084978011+0.23125231377479902 ⅈ},{-0.3692815807235537+3.2866081970081487 ⅈ,-2.624155084978011-0.23125231377479902 ⅈ},{1.7316064067187487,0.027742135918021837},{-0.9965216226358309-0.5717177306118897 ⅈ,1.1102840170189954-0.2565690723571353 ⅈ},{-0.9965216226358309+0.5717177306118897 ⅈ,1.1102840170189954+0.2565690723571353 ⅈ},{1.,-0.9999999999999994}}
```

Use the `"EndomorphismMatrix"` method:

```wolfram
NSolveValues[x^2+2y^2==3&&x^3-4x y==5,{x,y},Method->"EndomorphismMatrix"]
(* Output *)
{{-0.3692815807235537-3.2866081970081487 ⅈ,-2.624155084978011+0.23125231377479902 ⅈ},{-0.3692815807235537+3.2866081970081487 ⅈ,-2.624155084978011-0.23125231377479902 ⅈ},{1.7316064067187487,0.027742135918021837},{-0.9965216226358309-0.5717177306118897 ⅈ,1.1102840170189954-0.2565690723571353 ⅈ},{-0.9965216226358309+0.5717177306118897 ⅈ,1.1102840170189954+0.2565690723571353 ⅈ},{1.,-0.9999999999999994}}
```

Use the `"Homotopy"` method:

```wolfram
NSolveValues[x^2+2y^2==3&&x^3-4x y==5,{x,y},Method->"Homotopy"]
(* Output *)
{{-0.3692815807235462+3.2866081970081433 ⅈ,-2.6241550849780046-0.231252313774798 ⅈ},{-0.36928158072354544-3.286608197008143 ⅈ,-2.624155084978004+0.23125231377479816 ⅈ},{-0.996521622635828-0.571717730611889 ⅈ,1.1102840170189938-0.2565690723571348 ⅈ},{-0.9965216226358279+0.5717177306118892 ⅈ,1.1102840170189938+0.25656907235713483 ⅈ},{1.7316064067187473,0.027742135918021427},{1.,-1.}}
```

Use the `"Monodromy"` method:

```wolfram
NSolveValues[x^2+2y^2==3&&x^3-4x y==5,{x,y},Method->"Monodromy"]
(* Output *)
{{1.7316064067187473,0.027742135918021344},{-0.996521622635828+0.5717177306118894 ⅈ,1.1102840170189938+0.2565690723571349 ⅈ},{-0.996521622635828-0.5717177306118894 ⅈ,1.1102840170189938-0.2565690723571349 ⅈ},{1.,-1.},{-0.36928158072354567+3.2866081970081438 ⅈ,-2.6241550849780046-0.23125231377479788 ⅈ},{-0.36928158072354567-3.2866081970081438 ⅈ,-2.6241550849780046+0.23125231377479788 ⅈ}}
```

Use the `"Symbolic"` method:

```wolfram
NSolveValues[x^2+2y^2==3&&x^3-4x y==5,{x,y},Method->"Symbolic"]
(* Output *)
{{1.,-1.},{1.7316064067187473,0.02774213591802142},{-0.996521622635828-0.5717177306118894 ⅈ,1.1102840170189938-0.256569072357135 ⅈ},{-0.996521622635828+0.5717177306118894 ⅈ,1.1102840170189938+0.256569072357135 ⅈ},{-0.36928158072354567-3.2866081970081438 ⅈ,-2.6241550849780064+0.23125231377479588 ⅈ},{-0.36928158072354567+3.2866081970081438 ⅈ,-2.6241550849780064-0.23125231377479588 ⅈ}}
```

This system has $40$ roots, which is strictly less than the bound of $64$ provided by the Bernstein-Khovanskii-Kushnirenko theorem:

```wolfram
syst=t1^2+ t2^2+ t3^2-12t1-68==0&&t4^2+ t5^2+ t6^2-12t5-68==0&&t7^2+ t8^2+ t9^2-24 t8-12t9+100==0&&t1 t4+t2 t5+t3 t6-6 t1-6t5-52==0&&t1 t7+t2 t8+t3 t9-6 t1-12t8-6t9+64==0&&t4 t7+t5 t8+t6 t9-6t5-12t8-6t9+32==0&&2t2+2t3-t4-t5-2 t6-t7-t9+18==0&&t1+t2+2t3+2t4+2t6- 2t7+t8- t9- 38==0&&t1+t3- 2t4+t5- t6+2t7- 2t8+8==0;
vars={t1,t2,t3,t4,t5,t6,t7,t8,t9};
```

The `"Homotopy"` method, used by default, returns multiple copies of some of the roots:

```wolfram
(sols=NSolveValues[syst,vars])//Length//AbsoluteTiming
(* Output *)
{4.7082557,64}
```

Remove the multiple copies:

```wolfram
Union[sols,SameTest->(Norm[#1-#2]<10^-6&)]//Length
(* Output *)
40
```

The `"Monodromy"` method runs faster here and does not produce multiple copies of roots:

```wolfram
NSolveValues[syst,vars,Method->"Monodromy"]//Length//AbsoluteTiming
(* Output *)
{2.6839248,40}
```

The `"Monodromy"` method returns finitely many solutions of transcendental systems:

```wolfram
NSolveValues[Sin[x+y]+Cos[x-y]==1&&Sin[x-y]+Cos[x+y]==2,{x,y},Method->"Monodromy"]
(* Output *)
{{0.7853981633974483-0.48121182505960347 ⅈ,-0.3217505543966422},{3.9269908169872414-0.48121182505960347 ⅈ,-3.4633432079864352},{7.0685834705770345-0.48121182505960347 ⅈ,-6.604935861576228},{10.210176124166829-0.48121182505960347 ⅈ,-9.746528515166021},{13.351768777756622-0.48121182505960347 ⅈ,-12.888121168755815},{16.493361431346415-0.48121182505960347 ⅈ,-16.02971382234561},{19.634954084936208-0.48121182505960347 ⅈ,-19.171306475935403},{22.776546738526-0.48121182505960347 ⅈ,-22.312899129525196},{25.918139392115794-0.48121182505960347 ⅈ,-25.45449178311499},{29.059732045705587-0.48121182505960347 ⅈ,-28.596084436704782},{32.201324699295384-0.48121182505960347 ⅈ,-31.737677090294575},{35.34291735288517-0.48121182505960347 ⅈ,-34.879269743884365},{38.48451000647497-0.48121182505960347 ⅈ,-38.02086239747416},{41.62610266006476-0.48121182505960347 ⅈ,-41.16245505106395},{44.767695313654556-0.48121182505960347 ⅈ,-44.30404770465375},{47.909287967244346-0.48121182505960347 ⅈ,-47.445640358243544}}
```

Use the `"Symbolic"` method to obtain all solutions:

```wolfram
NSolveValues[Sin[x+y]+Cos[x-y]==1&&Sin[x-y]+Cos[x+y]==2,{x,y},Method->"Symbolic"]
(* Output *)
{{1. ((0.7853981633974482+0.4812118250596033 ⅈ)+6.283185307179586 1),-1. (0.3217505543966422+6.283185307179586 2)},{1. ((-2.356194490192345+0.4812118250596033 ⅈ)+6.283185307179586 1),-1. (-2.819842099193151+6.283185307179586 2)},{1. ((0.7853981633974482-0.4812118250596033 ⅈ)+6.283185307179586 1),-1. (0.3217505543966422+6.283185307179586 2)},{1. ((-2.356194490192345-0.4812118250596033 ⅈ)+6.283185307179586 1),-1. (-2.819842099193151+6.283185307179586 2)}}
```

The [Method](https://reference.wolfram.com/language/ref/Method.html) option may also be used to locally set system options from the `"NSolveOptions"` group:

```wolfram
"NSolveOptions"/.SystemOptions[]
(* Output *)
{"ComplexEquationMethod"->Automatic,"EndomorphismMatrixMaxVolume"->4096,"HomotopyMaxVolume"->131072,"MonomialOrder"->Automatic,"ReorderVariables"->True,"SelectCriterion"->(True&),"Tolerance"->0,"UseSlicingHyperplanes"->True}
```

By default [NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html) introduces slicing hyperplanes for underdetermined complex systems:

```wolfram
NSolveValues[x^2+y^2==1,{x,y}]
(* Output *)
NSolveValues
(* Output *)
{{-0.9977144816511953,0.06757080067223352},{-0.33598731323105446,-0.9418665114270585}}
```

With [Method](https://reference.wolfram.com/language/ref/Method.html)->{"UseSlicingHyperplanes"->[False](https://reference.wolfram.com/language/ref/False.html)}, [NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html) gives parametric solutions:

```wolfram
NSolveValues[x^2+y^2==1,{x,y},Method->{"UseSlicingHyperplanes"->False}]
(* Output *)
Solve
(* Output *)
{{x,-1. Sqrt[1.-1. x^2]},{x,Sqrt[1.-1. x^2]}}
```

#### VerifySolutions

[NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html) verifies solutions obtained using non-equivalent transformations:

```wolfram
funs={Sqrt[x+Sqrt[x]]-2,y^9-y-2 x^(1/7)-3};
eqns=And@@Thread[funs==0];
sol1=NSolveValues[eqns,{x,y}];
```

With [VerifySolutions](https://reference.wolfram.com/language/ref/VerifySolutions.html)->[False](https://reference.wolfram.com/language/ref/False.html), [NSolve](https://reference.wolfram.com/language/ref/NSolve.html) does not verify the solutions:

```wolfram
sol2=NSolveValues[eqns,{x,y},VerifySolutions->False];
```

Some of the solutions returned with [VerifySolutions](https://reference.wolfram.com/language/ref/VerifySolutions.html)->[False](https://reference.wolfram.com/language/ref/False.html) are not correct:

```wolfram
Length/@{sol1, sol2}
(* Output *)
{9,18}
```

This uses a fast numeric test in an attempt to select correct solutions:

```wolfram
sol3=Select[sol2,Chop[funs/.Thread[{x,y}->#],10^-8]=={0,0}&];
```

In this case, the simple numeric verification gives the correct solution set:

```wolfram
sol3===sol1
(* Output *)
True
```

#### WorkingPrecision

By default, [NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html) finds solutions of exact equations using machine-precision computations:

```wolfram
NSolveValues[x^5+2 x+7==0,x]
(* Output *)
{-1.3399779064938688,-0.5649399330245253-1.3374320953062897 ⅈ,-0.5649399330245253+1.3374320953062897 ⅈ,1.2349288862714596-0.9763463870286491 ⅈ,1.2349288862714596+0.9763463870286491 ⅈ}
```

This computes the solutions using 50-digit precision:

```wolfram
NSolveValues[x^5+2 x+7==0,x,WorkingPrecision->50]
(* Output *)
{-1.33997790649386879667997942271842329495403945969056267772238947526151957779219,-0.56493993302452527131791933695282521424843129696002478343220291034475585405606-1.33743209530628955835549902563516070893076035471422720389066995905102539786482 ⅈ,-0.56493993302452527131791933695282521424843129696002478343220291034475585405606+1.33743209530628955835549902563516070893076035471422720389066995905102539786482 ⅈ,1.23492888627145966965790904831203686172545102680530612229347730257107226556601-0.97634638702864905304299216265590087317416581856953055816606714302388594290732 ⅈ,1.23492888627145966965790904831203686172545102680530612229347730257107226556601+0.97634638702864905304299216265590087317416581856953055816606714302388594290732 ⅈ}
```

### Applications

#### Geometry

Find intersection points of a circle and a parabola:

```wolfram
pts=NSolveValues[x^2+y^2==1&&y-2x^2+3/2==0,{x,y}]
(* Output *)
{{-0.5877852522924726,-0.8090169943749473},{0.9510565162951534,0.30901699437494756},{0.5877852522924757,-0.8090169943749489},{-0.951056516295154,0.30901699437494784}}
```

```wolfram
Show[{ContourPlot[{x^2+y^2==1,y-2x^2+3/2==0},{x,-1.5,1.5},{y,-1.5,1.5}],Graphics[{Red,PointSize[Medium],Point[pts]}]}]
```

*([Graphics])*

Find the intersection of [InfiniteLine](https://reference.wolfram.com/language/ref/InfiniteLine.html)[{0,0},{1,1}] and [InfiniteLine](https://reference.wolfram.com/language/ref/InfiniteLine.html)[{{0,1},{1,0}}]:

```wolfram
line1=InfiniteLine[{0,0},{1,1}];
line2=InfiniteLine[{{0,1},{1,0}}];
```

```wolfram
sol=NSolveValues[{x,y}∈line1&&{x,y}∈line2,{x,y}]
(* Output *)
{{0.5,0.5}}
```

Plot it:

```wolfram
Graphics[{{Lighter[Blue,0.5],line1,line2},{Red,Point[sol]}},Axes->True,ImageSize->Small]
```

*([Graphics])*

Find the intersections of [InfiniteLine](https://reference.wolfram.com/language/ref/InfiniteLine.html)[{0,0},{1,1}] and [Circle](https://reference.wolfram.com/language/ref/Circle.html)[{0,0},1]:

```wolfram
line=InfiniteLine[{0,0},{1,1}];
circle=Circle[{0,0},1];
```

```wolfram
sol=NSolveValues[{x,y}∈line&&{x,y}∈circle,{x,y}]
(* Output *)
{{-0.7071067811865476,-0.7071067811865476},{0.7071067811865476,0.7071067811865476}}
```

Plot it:

```wolfram
Graphics[{{Lighter[Blue,0.5],line,circle},{Red,Point[sol]}},ImageSize->Small]
```

*([Graphics])*

Find all pairwise intersections between five random lines:

```wolfram
lines=InfiniteLine/@RandomReal[1,{5,2,2}];
```

Use [BooleanCountingFunction](https://reference.wolfram.com/language/ref/BooleanCountingFunction.html) to express that exactly two conditions are true:

```wolfram
sol=NSolveValues[{x,y}∈BooleanRegion[BooleanCountingFunction[{2},5],lines],{x,y}]
(* Output *)
{{0.7110237406474238,0.876221057974256},{0.47235056870279113,0.8051090504973661},{-1.0221680390646353,0.3598213897854922},{0.8426183725145454,0.9154293123792554},{0.7164587080519104,0.7120450894476968},{0.7264752767307693,0.40947108269156157},{0.7150456939628657,0.7547285021244847},{1.4557977919876302,0.4301789300752427},{0.6893621237183204,0.7223754105739162},{0.4343784830097863,0.4011775010468124}}
```

Plot it:

```wolfram
Graphics[{{Lighter[Blue,0.5],lines},{Red,PointSize[Medium],Point[sol]}}]
```

*([Graphics])*

Find the pairwise intersections of the circles [Circle](https://reference.wolfram.com/language/ref/Circle.html)[{1/3 [Cos](https://reference.wolfram.com/language/ref/Cos.html)[*k* 2π/5],1/3 [Sin](https://reference.wolfram.com/language/ref/Sin.html)[*k* 2π/5]}] for `*k*=0,…,4:

```wolfram
circles=Table[Circle[{1/3 Cos[k 2π/5],1/3 Sin[k 2π/5]}],{k,0,4}];
```

```wolfram
sol=NSolveValues[{x,y}∈BooleanRegion[BooleanCountingFunction[{2},5],circles],{x,y}]
(* Output *)
{{0.324908164063325,0.9999645076315318},{-0.26124716218830746,-0.8040360902007072},{0.324908164063325,-0.9999645076315318},{-0.26124716218830746,0.8040360902007072},{1.0115062882979888,0.7349023357933402},{-0.5751672901730063,-0.417883497028289},{1.0115062882979886,-0.7349023357933401},{-0.5751672901730065,0.4178834970282891},{-1.2502905320048137,0.},{0.7109458690881828,0.},{-0.38636102229558167,1.1890969577253128},{0.21969435562891512,-0.6761497015294375},{-0.8506206167401118,-0.61801205325984},{0.6839539500734455,0.4969216319256135},{-0.8506206167401121,0.6180120532598402},{0.6839539500734455,-0.4969216319256135},{-0.38636102229558167,-1.1890969577253128},{0.2196943556289151,0.6761497015294374},{1.051424905353574,0.},{-0.8454135757702758,0.}}
```

Plot it:

```wolfram
Graphics[{{Lighter[Blue,0.5],circles},{Red,PointSize[Medium],Point[sol]}}]
```

*([Graphics])*

Find the intersection of [InfiniteLine](https://reference.wolfram.com/language/ref/InfiniteLine.html)[{{-1,1,1},{1,1,1}}] and [InfinitePlane](https://reference.wolfram.com/language/ref/InfinitePlane.html)[{{2,0,0},{0,2,0},{0,0,2}}]:

```wolfram
line=InfiniteLine[{{-1,1,1},{1,1,1}}];
plane=InfinitePlane[{{2,0,0},{0,2,0},{0,0,2}}];
```

```wolfram
sol=NSolveValues[{{x,y,z}∈line,{x,y,z}∈plane},{x,y,z}]
(* Output *)
{{0.,1.,1.}}
```

Plot it:

```wolfram
Graphics3D[{{Opacity[0.5],line,plane},{Red,PointSize[Large],Point[sol]}},Axes->True,PlotRange->{{-2,2},{-2,2},{-2,2}}]
```

*([Graphics3D])*

Find the intersections of [InfiniteLine](https://reference.wolfram.com/language/ref/InfiniteLine.html)[{{-1,1,1},{1,1,1}}] and [Sphere](https://reference.wolfram.com/language/ref/Sphere.html)[{0,0,0},3]:

```wolfram
line=InfiniteLine[{{-1,1,1},{1,1,1}}];
sphere=Sphere[{0,0,0},3];
```

```wolfram
sol=NSolveValues[{{x,y,z}∈line,{x,y,z}∈sphere},{x,y,z}]
(* Output *)
{{-2.6457513110645907,1.,1.},{2.6457513110645907,1.,1.}}
```

Plot it:

```wolfram
Graphics3D[{{Opacity[0.5],line,sphere},{Red,PointSize[Large],Point[sol]}},PlotRangePadding->0.5]
```

*([Graphics3D])*

Find the intersections of [InfiniteLine](https://reference.wolfram.com/language/ref/InfiniteLine.html)[{{-1,1/3,1/2},{1,1/3,1/2}}] and the boundary of [Tetrahedron](https://reference.wolfram.com/language/ref/Tetrahedron.html)[{{0,0,0},{1,0,0},{0,1,0},{0,0,1}}]:

```wolfram
line=InfiniteLine[{{-1,1/3,1/2},{1,1/3,1/2}}];
tet=RegionBoundary@Tetrahedron[{{0,0,0},{1,0,0},{0,1,0},{0,0,1}}];
```

```wolfram
sol=NSolveValues[{x,y,z}∈RegionIntersection[line,tet],{x,y,z}]
(* Output *)
{{0.,0.3333333333333333,0.5},{0.16666666666666666,0.3333333333333333,0.5}}
```

Plot it:

```wolfram
Graphics3D[{{Opacity[0.5],line,Tetrahedron[{{0,0,0},{1,0,0},{0,1,0},{0,0,1}}]},{Red,PointSize[Medium],Point[sol]}},PlotRangePadding->0.2]
```

*([Graphics3D])*

Find the intersection for three random planes:

```wolfram
planes=InfinitePlane/@RandomReal[1,{3,3,3}];
```

```wolfram
sol=NSolveValues[{x,y,z}∈RegionIntersection@@planes,{x,y,z}]
(* Output *)
{{0.4098644591774831,0.6616584491523385,0.1039625578597415}}
```

Plot it:

```wolfram
Graphics3D[{{Opacity[0.5],planes},{Red,PointSize[Medium],Point[sol]}},PlotRangePadding->0.1]
```

*([Graphics3D])*

Find the intersections of the spheres [Sphere](https://reference.wolfram.com/language/ref/Sphere.html)[{1/3 [Cos](https://reference.wolfram.com/language/ref/Cos.html)[*k* 2π/3],1/3 [Sin](https://reference.wolfram.com/language/ref/Sin.html)[*k* 2π/3],0}] for `*k*=0,1,2:

```wolfram
spheres=Table[Sphere[{1/3 Cos[k 2π/3],1/3 Sin[k 2π/3],0}],{k,0,2}];
```

```wolfram
sol=NSolveValues[{x,y,z}∈RegionIntersection@@spheres,{x,y,z}]
(* Output *)
{{0.,0.,-0.9428090415820632},{0.,0.,0.9428090415820635}}
```

Plot it:

```wolfram
Graphics3D[{{Opacity[0.5],spheres},{Red,PointSize[Large],Point[sol]}}]
```

*([Graphics3D])*

Find all the intersections of exactly three planes among 10 random planes:

```wolfram
planes=InfinitePlane/@RandomReal[{-1,1},{10,3,3}];
```

Use [BooleanCountingFunction](https://reference.wolfram.com/language/ref/BooleanCountingFunction.html) to find the condition of exactly three things being true:

```wolfram
sol=NSolveValues[{x,y,z}∈BooleanRegion[BooleanCountingFunction[{3},10],planes],{x,y,z}];//Quiet
```

Plot it:

```wolfram
Graphics3D[{{Opacity[0.3],planes},{Red,PointSize[Medium],Point[sol]}},PlotRange->{{-1,1},{-1,1},{-1,1}}]
```

*([Graphics3D])*

#### Chemistry

A polynomial system for the equilibria of a certain chemical reaction network:

```wolfram
reactionpolys={50+x_1-x_4-(4 x_2 x_4)/(125)-(x_2 x_5)/(50)-(99 x_6)/(100)+x_7,-(9)/(200) x_3 x_5+15 (50-x_4-x_6)+x_8,-(4)/(125) x_2 x_4-(1)/(50) x_4 (350-x_2-x_3+x_4+x_5)+16 (50-x_4-x_6)+(101 x_6)/(100),x_1-(x_2 x_5)/(50)-(9 x_3 x_5)/(200)-(11 x_5 (350-x_2-x_3+x_4+x_5))/(10000)+x_7+(43)/(500) (100-x_1-x_5-x_7-x_8)+x_8,(1)/(50) x_4 (350-x_2-x_3+x_4+x_5)-(101 x_6)/(100),(9 x_3 x_5)/(200)-(273 x_8)/(250),(x_2 x_5)/(100)-x_7+(23 x_8)/(250),-(3 x_1)/(2)+(x_2 x_5)/(100)};
vars=Variables[reactionpolys];
```

Since the variables represent chemical species quantities, we are interested in real-valued solutions with all components nonnegative::

```wolfram
nonnegSolutions[polys_,vars_,meth_:Automatic]:=Module[
{solns,realsolns},
solns=NSolveValues[polys,vars,Method->meth];
realsolns=Select[solns,FreeQ[#,Complex]&];
Select[realsolns,Apply[And,Map[#>=0&,#]]&]
]
```

Find the solutions:

```wolfram
nonnegsolns=nonnegSolutions[reactionpolys,vars]
(* Output *)
{{0.25622064316116705,6.97674545003859,367.56996902631346,36.67719296367814,5.508742829934045,12.811032158062812,8.060954138627576,83.4415562378894},{0.2562206431612769,6.976745450041114,367.5699690262705,36.67719296366311,5.508742829934475,12.81103215807785,8.06095413862747,83.44155623788646},{0.701280751097668,14.6721264880418,234.97359857167103,14.510172279042063,7.169520570203967,35.06403755487517,7.438773617888077,69.42230968740843},{0.8628563267160883,9.496206390697303,37.10128286324027,6.729376566694116,13.629489891269026,43.14281633578979,3.211390952815757,20.83811372545244}}
```

The apparent multiplicity is an artifact of the [Automatic](https://reference.wolfram.com/language/ref/Automatic.html) method, which tends to be fast but can sometimes overstate multiplicity:

```wolfram
nonnegSolnsE=nonnegSolutions[reactionpolys,vars,"EndomorphismMatrix"]
(* Output *)
{{0.2562206431611096,6.976745450037269,367.5699690263356,36.67719296368546,5.508742829933561,12.811032158055479,8.060954138627627,83.44155623789091},{0.7012807510973713,14.67212648803852,234.9735985717492,14.510172279048598,7.1695205702039955,35.06403755486857,7.438773617888614,69.4223096874191},{0.8628563267158091,9.496206390692775,37.10128286321059,6.729376566693518,13.629489891271184,43.14281633579046,3.2113909528140976,20.83811372543896}}
```

#### Mechanics

Direct kinematics of a [Gough-Stewart](https://link.springer.com/article/10.1007/BF02142677) parallel 6-degree of freedom platform:

```wolfram
eqns={t1^2+ t2^2+ t3^2-12t1-68,t4^2+ t5^2+ t6^2-12t5-68,t7^2+ t8^2+ t9^2-24 t8-12t9+100,t1 t4+t2 t5+t3 t6-6 t1-6t5-52,t1 t7+t2 t8+t3 t9-6 t1-12t8-6t9+64,t4 t7+t5 t8+t6 t9-6t5-12t8-6t9+32,2t2+2t3-t4-t5-2 t6-t7-t9+18,t1+t2+2t3+2t4+2t6- 2t7+t8- t9- 38,t1+t3- 2t4+t5- t6+2t7- 2t8+8};
vars={t1,t2,t3,t4,t5,t6,t7,t8,t9};
Length[solns=NSolveValues[eqns,vars]]
(* Output *)
64
```

This solution set claims multiple solutions; remove the multiple copies:

```wolfram
Union[solns,SameTest->(Norm[#1-#2]<10^-6&)]//Length
(* Output *)
40
```

A different method does not produce multiple copies of solutions:

```wolfram
Length[NSolveValues[eqns,vars,Method->"Monodromy"]]
(* Output *)
40
```

Set up an overdetermined system of six equations in four variables that arises from a camera pose estimation procedure:

```wolfram
coords={{1,2,1.49071,4},{1,3,.400000,8},{1,4,.894427,4},{2,3,1.49071,4},{2,4,.666667,8},{3,4,.894427,4}};
vars=Array[x,4];
polys=MapThread[x[#1]^2+x[#2]^2-#3*x[#1]*x[#2]-#4&,Transpose[coords]]
(* Output *)
{-4+x[1]^2-1.49071 x[1] x[2]+x[2]^2,-8+x[1]^2-0.4 x[1] x[3]+x[3]^2,-4+x[1]^2-0.894427 x[1] x[4]+x[4]^2,-4+x[2]^2-1.49071 x[2] x[3]+x[3]^2,-8+x[2]^2-0.666667 x[2] x[4]+x[4]^2,-4+x[3]^2-0.894427 x[3] x[4]+x[4]^2}
```

Use the first four polynomials so the subsystem will be exactly determined:

```wolfram
realsols=NSolveValues[polys[[1;;4]],vars,Reals]
(* Output *)
{{2.2250807284985488,2.999950067090309,2.246974836013612,0.7970654888514712},{2.2250807284985488,2.999950067090309,2.246974836013612,1.1931067918973002},{-2.225080728487858,-2.999950067090227,-2.246974836024156,-1.193106791988619},{-2.225080728487858,-2.999950067090227,-2.246974836024156,-0.7970654887505906}}
```

For each solution the first four residuals from plugging in a solution are small whereas the last two are non-negligible:

```wolfram
residuals=(polys/.Thread[vars->#])&/@realsols
(* Output *)
{{-8.881784197001252×10^-16,0.,-1.1102230246251565×10^-16,-2.6645352591003757×10^-15,0.040908556997401146,0.08230298621629151},{-8.881784197001252×10^-16,0.,0.,-2.6645352591003757×10^-15,0.03702916188754379,0.07454743575458389},{1.509903313490213×10^-14,3.375077994860476×10^-14,-2.220446049250313×10^-16,1.3322676295501878×10^-14,0.03702916192238859,0.07454743582509438},{1.509903313490213×10^-14,3.375077994860476×10^-14,-1.1102230246251565×10^-16,1.3322676295501878×10^-14,0.040908557037893645,0.08230298629808719}}
```

Now polish each solution by using it as a starting point for minimizing the sum of squares of residuals:

```wolfram
newsols=Map[FindMinimum[polys.polys,Transpose[{vars,#}]]&,realsols]
(* Output *)
{{2.099282803988261×10^-10,{x[1]->2.2360674544354247,x[2]->2.99999781255834,x[3]->2.2360674544354247,x[4]->0.9977225194226561}},{2.0995801785780716×10^-10,{x[1]->2.236067454875393,x[2]->2.9999978123321838,x[3]->2.236067454875393,x[4]->1.002276484194443}},{2.0995801785780716×10^-10,{x[1]->-2.236067454875393,x[2]->-2.9999978123321838,x[3]->-2.236067454875393,x[4]->-1.00227648419443}},{2.099282803988261×10^-10,{x[1]->-2.2360674544354247,x[2]->-2.99999781255834,x[3]->-2.2360674544354247,x[4]->-0.9977225194226546}}}
```

A system of equations for static equilibria of a weight suspended by two massless cables in a plane:

```wolfram
polys={-1+(3 ct1)/(4)+(7 ct2)/(5)-(9 ct3)/(10),(1)/(8)+(3 st1)/(4)+(7 st2)/(5)-(9 st3)/(10),(3 ct1)/(4)+(3 ct2)/(5)+(st2)/(5)-x,-(ct2)/(5)+(3 st1)/(4)+(3 st2)/(5)-y,ct1 f1+ct3 f3,-1+f1 st1+f3 st3,ct1 t1-x,st1 t1-y0,1+ct3 t3-x,-(1)/(8)+st3 t3-y0,-1+ct1^2+st1^2,-1+ct2^2+st2^2,-1+ct3^2+st3^2};
vars=Variables[polys];
Length[solns=NSolveValues[polys==0,vars]]
(* Output *)
20
```

Each real-valued solution corresponds to a distinct equilibrium position, and there are six such:

```wolfram
realsols=Union[Cases[solns,{_Real..}],SameTest->(Norm[#1-#2]<=10^-9&)]
(* Output *)
{{-0.29281420382404505,0.4387188013100926,-0.6726714789265602,-0.7823097040366003,0.3405397735858377,-0.9561693584501078,0.8986244006129788,0.7399412689117659,-0.7627550341609386,1.1545821642367786,0.22334550804061748,-0.26569613873181197,0.7293229916682546},{-0.20790094624983305,0.9999850088505232,0.2711703363370643,0.5827153692957288,0.44675637574279997,0.9781498845005422,0.005475588938079906,0.9625313754319094,-2.14121398406532,-2.0460924822535915,0.44516041341055507,0.53690076496815,-2.0944282112044386},{-0.19189808534297012,0.982802649001525,0.2577779384387856,-0.5879973464645284,-0.43772390165221625,-0.9814148586819451,-0.1846590185059623,-0.9662041887997844,-2.1304340841224114,-2.2933451244427996,0.40882622169249483,-1.0434170849153412,2.0908396656001953},{0.06744505250647932,0.14900491386727416,-0.8231214791177298,1.0513111821690044,0.08614249497903796,0.9977229900590627,-0.988836455458336,-0.5678653278859702,-0.8567055883882784,1.285084377247979,-0.05778055339144323,0.12518938649584055,-0.8547548612470611},{0.7142123856253132,0.3748611065982416,0.06718315384058107,-0.08844285466479979,0.9402205554616933,-0.699929045132004,0.9270809839274786,0.9977406596105186,1.3245249858488282,-0.8038897692229479,0.9459921499634256,-0.04367041481216427,-0.9270735085986515},{0.9265404773993471,0.47450969622505457,0.39913214751620707,0.3992557520268481,-0.9268274113960465,0.3761950873424454,-0.8802502758809033,-0.9168934119182559,0.867270389377733,-0.49216501505597526,0.8035611206083623,-0.3409057892667189,0.32626285988147286}}
```

#### Economics

A reduced 8-dimensional system arising in [economics](http://homepages.math.uic.edu/~jan/Demo/redeco8.html); we use a non-default method to get better accuracy:

```wolfram
eqns={x1+x1*x2+x2*x3+x3*x4+x4*x5+x5*x6+x6*x7-1*u8,
x2+x1*x3+x2*x4+x3*x5+x4*x6+x5*x7-2*u8,
x3+x1*x4+x2*x5+x3*x6+x4*x7-3*u8,
x4+x1*x5+x2*x6+x3*x7-4*u8,
x5+x1*x6+x2*x7-5*u8,
x6+x1*x7-6*u8,
x7-7*u8,
x1+x2+x3+x4+x5+x6+x7+1};
vars={x1,x2,x3,x4,x5,x6,x7,u8};
Length[solns=NSolveValues[eqns,vars,Method->"Monodromy"]]
(* Output *)
64
```

Find the real-valued solutions:

```wolfram
Cases[solns,{_Real..}]
(* Output *)
{{0.2518736446040636,1.0912208114082294,-1.1142421117940848,-0.6371519417366686,0.2911452761960025,-0.3328783138382896,-0.5499673648392522,-0.07856676640560746},{-0.7422463999814358,0.6020750900222276,0.7366348698613789,0.06997887540047072,-0.5553830480454723,-0.6836284496967303,-0.4274309375604388,-0.06106156250863412},{-0.3475156125202294,1.59955565082714,-0.7141951080636515,1.3220970404197319,0.7678488168410543,-1.9822793228378468,-1.6455114646661986,-0.23507306638088551},{0.6052692125387936,-0.5293864596512925,1.1585268189920672,2.026014965669394,-1.9841555720804958,-0.4579792560558256,-1.8182897094126411,-0.25975567277323447},{1.,1.,1.,1.,1.,1.,-7.,-1.},{-0.14285714285714285,-0.14285714285714285,-0.14285714285714285,-0.14285714285714285,-0.14285714285714285,-0.14285714285714285,-0.14285714285714285,-0.02040816326530612},{1.2046584696630864,-0.4666323105787761,-0.8034565962066553,0.4340262121531496,-0.972072018441767,0.21119002813555293,-0.6077137847245906,-0.0868162549606558},{1.599389257124293,1.299351542532976,-0.16371972464107304,-1.7234009172703337,-1.4085903408362763,1.736529424420716,-2.339559241330302,-0.33422274876147173}}
```

#### Difference Equations

A nonlinear system of polynomial difference equations:

```wolfram
des={a[n+1]==b[n]^2*a[n]-3*c[n]-1,b[n+1]==a[n]*c[n]^2-2*b[n-1]+1,c[n+1]==a[n-1]-b[n-2]};
```

[RSolve](https://reference.wolfram.com/language/ref/RSolve.html) cannot find an analytic form for the solution:

```wolfram
RSolve[des,{a,b,c},n]
(* Output *)
RSolve[{a[1+n]==-1+a[n] b[n]^2-3 c[n],b[1+n]==1-2 b[-1+n]+a[n] c[n]^2,c[1+n]==a[-1+n]-b[-2+n]},{a,b,c},n]
```

Set up a polynomial system to find the asymptotic values:

```wolfram
rep={a[_]:>a,b[_]:>b,c[_]:>c};
asymp=des/.rep
(* Output *)
{a==-1+a b^2-3 c,b==1-2 b+a c^2,c==a-b}
```

Solve it and give asymptotic possible real asymptotic values (they depend on initial conditions):

```wolfram
desolns=NSolveValues[asymp,{a,b,c}];
Select[desolns,MatchQ[#,{_Real..}]&]
(* Output *)
{{-2.827425691765565,-1.4513424802950081,-1.376083211470557},{2.793407474408426,1.6185565498521337,1.1748509245562924},{0.,0.33333333333333304,-0.3333333333333326}}
```

### Properties & Relations

Solutions approximately satisfy the equations:

```wolfram
polys={x^2+y^2-1,2x+3y-4};
NSolveValues[Thread[polys==0],{x,y}]
(* Output *)
{{0.615384615384615-0.39970403251589404 ⅈ,0.9230769230769234+0.26646935501059604 ⅈ},{0.615384615384615+0.39970403251589404 ⅈ,0.9230769230769234-0.26646935501059604 ⅈ}}
```

```wolfram
polys/.(Thread[{x,y}->#]&/@%)
(* Output *)
{{6.661338147750939×10^-16+4.996003610813204×10^-16 ⅈ,0.+0. ⅈ},{6.661338147750939×10^-16-4.996003610813204×10^-16 ⅈ,0.+0. ⅈ}}
```

For univariate equations, [NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html) repeats solutions according to their multiplicity:

```wolfram
NSolveValues[(x-1)^2(x-2)^3==0,x]
(* Output *)
{1.,1.,2.,2.,2.}
```

Find solutions over specified domains:

```wolfram
NSolveValues[(x^4-1)(x^4-4)==0,x]
(* Output *)
{-1.4142135623730951,-1.,0.-1. ⅈ,0.+1. ⅈ,0.-1.4142135623730951 ⅈ,0.+1.4142135623730951 ⅈ,1.,1.4142135623730951}
```

```wolfram
NSolveValues[(x^4-1)(x^4-4)==0,x,Reals]
(* Output *)
{-1.4142135623730951,-1.,1.,1.4142135623730951}
```

[NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html) gives values of the solutions:

```wolfram
NSolveValues[x^4==4,x]
(* Output *)
{-1.4142135623730951,0.-1.4142135623730951 ⅈ,0.+1.4142135623730951 ⅈ,1.4142135623730951}
```

[NSolve](https://reference.wolfram.com/language/ref/NSolve.html) represents solutions in terms of replacement rules:

```wolfram
NSolve[x^4==4,x]
(* Output *)
{{x->-1.4142135623730951},{x->0.-1.4142135623730951 ⅈ},{x->0.+1.4142135623730951 ⅈ},{x->1.4142135623730951}}
```

[NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html) is a global equation solver:

```wolfram
NSolveValues[{x^2+y^3==1,x^4+y^4==2},{x,y}]
(* Output *)
{{0.4286528668755243-1.425668037428391 ⅈ,-0.5523289324629057-1.3494243351694064 ⅈ},{0.4286528668755243+1.425668037428391 ⅈ,-0.5523289324629057+1.3494243351694064 ⅈ},{-1.1536260276627142,-0.6916372321061667},{1.2009746404322088-0.05002165858925961 ⅈ,0.32414717979453855-0.6980877340703181 ⅈ},{1.2009746404322088+0.05002165858925961 ⅈ,0.32414717979453855+0.6980877340703181 ⅈ},{-1.2009746404322121+0.05002165858926033 ⅈ,0.32414717979453883-0.6980877340703202 ⅈ},{-1.2009746404322121-0.05002165858926033 ⅈ,0.32414717979453883+0.6980877340703202 ⅈ},{0.-0.7162099605832936 ⅈ,1.148000737442902},{0.+0.7162099605832936 ⅈ,1.148000737442902},{-0.4286528668755351-1.4256680374283892 ⅈ,-0.5523289324628892+1.3494243351693982 ⅈ},{-0.4286528668755351+1.4256680374283892 ⅈ,-0.5523289324628892-1.3494243351693982 ⅈ},{1.1536260276627133,-0.6916372321061641}}
```

[FindRoot](https://reference.wolfram.com/language/ref/FindRoot.html) is a local equation solver:

```wolfram
FindRoot[{x^2+y^3==1,x^4+y^4==2},{{x,1},{y,-1}}]
(* Output *)
{x->1.1536260276627142,y->-0.6916372321061672}
```

[NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html) gives approximate results:

```wolfram
NSolveValues[{x^2+y^2==1,x+y^2==2},{x,y}]
(* Output *)
{{0.5000000000000003-0.86602540378444 ⅈ,1.2712298784187062+0.3406250193166061 ⅈ},{0.5000000000000003+0.86602540378444 ⅈ,1.2712298784187062-0.3406250193166061 ⅈ},{0.49999999999999983-0.8660254037844385 ⅈ,-1.2712298784187057-0.3406250193166066 ⅈ},{0.49999999999999983+0.8660254037844385 ⅈ,-1.2712298784187057+0.3406250193166066 ⅈ}}
```

Use [SolveValues](https://reference.wolfram.com/language/ref/SolveValues.html) to get exact solutions:

```wolfram
SolveValues[{x^2+y^2==1,x+y^2==2},{x,y}]
(* Output *)
{{(1)/(2)+(ⅈ Sqrt[3])/(2),-Sqrt[(3)/(2)-(ⅈ Sqrt[3])/(2)]},{(1)/(2)+(ⅈ Sqrt[3])/(2),Sqrt[(3)/(2)-(ⅈ Sqrt[3])/(2)]},{(1)/(2)-(ⅈ Sqrt[3])/(2),-Sqrt[(3)/(2)+(ⅈ Sqrt[3])/(2)]},{(1)/(2)-(ⅈ Sqrt[3])/(2),Sqrt[(3)/(2)+(ⅈ Sqrt[3])/(2)]}}
```

Use [FindInstance](https://reference.wolfram.com/language/ref/FindInstance.html) to get exact solution instances:

```wolfram
FindInstance[{x^2+y^2==1,x+y^2==2},{x,y}]
(* Output *)
{{x->(1)/(2) (1-ⅈ Sqrt[3]),y->-Sqrt[2+(1)/(2) (-1+ⅈ Sqrt[3])]}}
```

Use [NDSolveValue](https://reference.wolfram.com/language/ref/NDSolveValue.html) to solve differential equations numerically:

```wolfram
NDSolveValue[{y''[x]==-y[x],y[0]==0,y'[0]==1},y[x],{x,0,7}]
(* Output *)
InterpolatingFunction[...][x]
```

```wolfram
Plot[%,{x,0,7}]
```

*([Graphics])*

### Possible Issues

Solutions obtained with machine-precision numeric computations may not be accurate:

```wolfram
poly=Expand[Product[x-i,{i,30}]];
```

```wolfram
NSolveValues[poly==0,x]
(* Output *)
{1.0000000000000029,1.9999999999939861,3.000000001935441,3.999999871226411,5.00000537723431,5.999777466514384,6.9930312038646525,8.140798407562976,8.385018146901867,9.129441151541732,9.492712058236105,10.602895782736502,11.105371028608976,12.108098380383126,12.610381108137158,14.392260909180427,14.81771412476811,16.811709162215493,17.394067761945177,19.791753137723486,20.459384326489396,23.185503980922576,23.660474438981687,26.53025161910415,27.02089570041308,28.97052366327305,29.571955737119755,31.337461687278243,31.79221781358922,31.971239859831275}
```

With higher [WorkingPrecision](https://reference.wolfram.com/language/ref/WorkingPrecision.html), more accurate results are produced:

```wolfram
NSolveValues[poly==0,x,WorkingPrecision->20]
(* Output *)
{1.,2.,3.,4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.,21.,22.,23.,24.,25.,26.,27.,28.,29.,30.}
```

Approximate solutions may not satisfy the equations due to numeric errors:

```wolfram
eqns={x^2+y^2==1,x^3-2x y==2};
sols=NSolveValues[eqns,{x,y}]
(* Output *)
{{-0.31643436536437397-2.2779040662116814 ⅈ,-2.484529090574299+0.29011820802870963 ⅈ},{-0.31643436536437397+2.2779040662116814 ⅈ,-2.484529090574299-0.29011820802870963 ⅈ},{-0.6710195936385324-0.5384250269557445 ⅈ,0.9867593241523941-0.3661417064419465 ⅈ},{-0.6710195936385324+0.5384250269557445 ⅈ,0.9867593241523941+0.3661417064419465 ⅈ},{0.9874539590029026-0.21613625162125424 ⅈ,-0.5022302335780996-0.42495370266128146 ⅈ},{0.9874539590029026+0.21613625162125424 ⅈ,-0.5022302335780996+0.42495370266128146 ⅈ}}
```

```wolfram
eqns/.(Thread[{x,y}->#]&/@sols)
(* Output *)
{{False,False},{False,False},{False,False},{False,False},{False,False},{False,False}}
```

The equations are satisfied up to a certain tolerance:

```wolfram
Function[{x,y},Norm/@{x^2+y^2-1,x^3-2x y-2}]@@@sols
(* Output *)
{{1.1520670436369372×10^-14,4.572393917918694×10^-14},{1.1520670436369372×10^-14,4.572393917918694×10^-14},{3.9983445677392×10^-15,5.816763469748014×10^-15},{3.9983445677392×10^-15,5.816763469748014×10^-15},{2.7761126175751524×10^-15,8.042797661766892×10^-15},{2.7761126175751524×10^-15,8.042797661766892×10^-15}}
```

Using higher [WorkingPrecision](https://reference.wolfram.com/language/ref/WorkingPrecision.html) will give solutions with smaller tolerance:

```wolfram
sols=NSolveValues[eqns,{x,y},WorkingPrecision->25]
(* Output *)
{{-0.31643436536436919281338185525367994749-2.27790406621167937399455326792316736456 ⅈ,-2.48452909057429637214674356387776638743+0.29011820802870734342693208222921302971 ⅈ},{-0.31643436536436919281338185525367994749+2.27790406621167937399455326792316736456 ⅈ,-2.48452909057429637214674356387776638743-0.29011820802870734342693208222921302971 ⅈ},{-0.67101959363853400283420150812620330829-0.53842502695574370694151331995655004928 ⅈ,0.98675932415239468841956351583305747345-0.3661417064419467372185090734054666597 ⅈ},{-0.67101959363853400283420150812620330829+0.53842502695574370694151331995655004928 ⅈ,0.98675932415239468841956351583305747345+0.3661417064419467372185090734054666597 ⅈ},{0.98745395900290319564758336337988325687-0.21613625162125236231744611152141729119 ⅈ,-0.50223023357809831627281995195529108363-0.42495370266128180206759636311303434353 ⅈ},{0.98745395900290319564758336337988325687+0.21613625162125236231744611152141729119 ⅈ,-0.50223023357809831627281995195529108363+0.42495370266128180206759636311303434353 ⅈ}}
```

```wolfram
Function[{x,y},Norm/@{x^2+y^2-1,x^3-2x y-2}]@@@sols
(* Output *)
{{0`22.877214315977312,0`22.49963005646669},{0`22.877214315977312,0`22.49963005646669},{0`23.66058977613498,0`23.455889898053915},{0`23.66058977613498,0`23.455889898053915},{0`23.753664040791534,0`23.472135880195257},{0`23.753664040791534,0`23.472135880195257}}
```

If the solutions set is infinite, [NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html) gives its intersection with random hyperplanes:

```wolfram
NSolveValues[x^2+y^2==1,{x,y}]
(* Output *)
NSolveValues
(* Output *)
{{-0.9977144816511953,0.06757080067223352},{-0.33598731323105446,-0.9418665114270585}}
```

```wolfram
NSolveValues[x^2+y^2+z^2==1,{x,y,z}]
(* Output *)
NSolveValues
(* Output *)
NSolveValues
(* Output *)
{{-11.230368419360744+6.317399394968035 ⅈ,0.7553697879857248-12.367421869931137 ⅈ,-10.982516189508647-7.31059240962063 ⅈ},{-11.230368419360744-6.317399394968035 ⅈ,0.7553697879857248+12.367421869931137 ⅈ,-10.982516189508647+7.31059240962063 ⅈ}}
```

Use [ContourPlot](https://reference.wolfram.com/language/ref/ContourPlot.html) and [ContourPlot3D](https://reference.wolfram.com/language/ref/ContourPlot3D.html) to view the real part of solutions:

```wolfram
ContourPlot[x^2+y^2==1,{x,-1,1},{y,-1,1}]
```

*([Graphics])*

```wolfram
ContourPlot3D[x^2+y^2+z^2==1,{x,-1,1},{y,-1,1},{z,-1,1}]
(* Output *)
![image](img/image_001.png)
```

[NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html) by default can claim multiple solutions when the actual count should be smaller:

```wolfram
polys={9*y2^2-5.656854249492381*y2+z2,x3^2+y3^2+z3^2-1,x4^2+y4^2+z4^2-1,y5^2+z5^2-0.888888888888889,x3-2.828427124746190*y2*x3+y2*y3+z2*z3-1/3,x3*x4+y3*y4+z3*z4-1/3,1/3*x4+y4*y5+z4*z5-1/3,8/3-2.828427124746190*y2+x3+x4,y2+y3+y4+y5+0.8888888888888889,z2+z3+z4+z5};
vars={x3,x4,y2,y3,y4,y5,z2,z3,z4,z5};
Length[solns=NSolveValues[polys,vars]]
(* Output *)
122
```

```wolfram
Length[Union[solns,SameTest->(Norm[#1-#2]<10^-9&)]]
(* Output *)
40
```

Compare to a nondefault method:

```wolfram
Length[NSolveValues[polys,vars,Method->"Monodromy"]]
(* Output *)
40
```

Validate against another nondefault method:

```wolfram
Length[NSolveValues[polys,vars,Method->"EndomorphismMatrix"]]
(* Output *)
40
```

Use a particular method to solve a polynomial system with large coefficients:

```wolfram
polys={-1+a+b+c+d,-2772 a^2-4620 a b-1980 b^2-3960 a c-3465 b c-1540 c^2-3465 a d-3080 b d-2772 c d-1260 d^2+13860 r,926640 a^3+2432430 a^2 b+2116400 a b^2+617760 b^3+2230800 a^2 c+3848130 a b c+1684800 b^2 c+1731600 a c^2+1526525 b c^2+462000 c^3+2108106 a^2 d+3603600 a b d+1573000 b^2 d+3211065 a c d+2836680 b c d+1287000 c^2 d+1474200 a d^2+1310400 b d^2+1192464 c d^2+368550 d^3+16216200 a r+14414400 b r+13513500 c r+12972960 d r,213840 a^2+374220 a b+162800 b^2+343200 a c+296010 b c+133200 c^2+324324 a d+277200 b d+247005 c d+113400 d^2+332640 a lambda+277200 b lambda+237600 c lambda+207900 d lambda+831600 mu+1247400 r,486486 a^2+846560 a b+370656 b^2+769626 a c+673920 b c+305305 c^2+720720 a d+629200 b d+567336 c d+262080 d^2+720720 a lambda+617760 b lambda+540540 c lambda+480480 d lambda+2162160 mu+2882880 r,2230800 a^2+3848130 a b+1684800 b^2+3463200 a c+3053050 b c+1386000 c^2+3211065 a d+2836680 b d+2574000 c d+1192464 d^2+3088800 a lambda+2702700 b lambda+2402400 c lambda+2162160 d lambda+10810800 mu+13513500 r,2108106 a^2+3603600 a b+1573000 b^2+3211065 a c+2836680 b c+1287000 c^2+2948400 a d+2620800 b d+2384928 c d+1105650 d^2+2702700 a lambda+2402400 b lambda+2162160 c lambda+1965600 d lambda+10810800 mu+12972960 r};
vars={r,a,b,c,d,lambda,mu};
Length[solmach=NSolveValues[polys,{r,a,b,c,d,lambda,mu},Method->"EndomorphismMatrix"]]
(* Output *)
21
```

Use the same method but with higher precision than default:

```wolfram
Length[solbig=NSolveValues[polys,{r,a,b,c,d,lambda,mu},Method->"EndomorphismMatrix",WorkingPrecision->100]]
(* Output *)
21
```

At machine precision some residuals are not so small:

```wolfram
Map[Max,Abs[(polys/.Thread[vars->#]&)/@solmach[[1;;4]]]]
(* Output *)
{0.18948332555679537,0.18948332555679537,0.03128431343679951,0.03128431343679951}
```

The solution using higher precision, not surprisingly, gives smaller residuals:

```wolfram
Map[Max,Abs[(polys/.Thread[vars->#]&)/@solbig[[1;;4]]]]
(* Output *)
{0`73.97122773936566,0`73.97122773936566,0`74.36587395368046,0`74.36587395368046}
```

Despite the disparities in residuals, the two solutions agree to all digits of [MachinePrecision](https://reference.wolfram.com/language/ref/MachinePrecision.html):

```wolfram
solmach[[1;;4]]-solbig[[1;;4]]
(* Output *)
{{0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ},{0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ},{0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ},{0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ,0.+0. ⅈ}}
```

With domain [Reals](https://reference.wolfram.com/language/ref/Reals.html) specified, [NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html) may not find solutions for which the function is
not real valued in any neighborhood of the solution:

```wolfram
f=(x^2-7 x+11)^(x^2-13 x+42)-1;
```

```wolfram
NSolveValues[f==0,x,Reals]
(* Output *)
{2.,5.,6.,7.}
```

```wolfram
Plot[{Abs[f],Im[f]},{x,1.5,7.5},PlotLegends->"Expressions"]
(* Output *)
![image](img/image_003.png)
```

This gives all real solutions:

```wolfram
NSolveValues[f==0&&Element[x,Reals],x]
(* Output *)
{2.,3.,4.,5.,6.,7.}
```

[NSolveValues](https://reference.wolfram.com/language/ref/NSolveValues.html) may not give all solutions:

```wolfram
NSolveValues[Cos[x]==x,x]
(* Output *)
{0.7390851332151607,-2.4868856989085604+1.8093613412957033 ⅈ,-9.109987453936563+2.950170861699437 ⅈ,-15.487957788766868+3.4566149353000477 ⅈ,-21.818976539497065+3.7903125588244637 ⅈ,-28.131612972193157+4.039950881461229 ⅈ,-34.434967951887934+4.239540266390788 ⅈ,-40.73292708439513+4.4058542622751045 ⅈ,-47.02744973977256+4.548424206260375 ⅈ,-53.31963845956909+4.673192083315507 ⅈ,-59.610164228277455+4.7841145036634956 ⅈ,-65.8994602437461+4.883959600084061 ⅈ,-72.18781942101943+4.97474014452765 ⅈ,-78.4754473341386+5.057965670789603 ⅈ,-84.76249275138699+5.134797433655855 ⅈ,-91.04906613536808+5.206147951828253 ⅈ}
```

Get $100$ solutions:

```wolfram
NSolveValues[Cos[x]==x,x,MaxRoots->100]//Length
(* Output *)
100
```

### Neat Examples

Solve the equation $sin(z+sin(z+sin(z)))=cos(z+cos(z+cos(z)))$:

```wolfram
roots=NSolveValues[Sin[z+Sin[z+Sin[z]]]==Cos[z+Cos[z+Cos[z]]]&&-3<Re[z]<3&&-3<Im[z]<3,z,MaxRoots->Infinity];
```

```wolfram
ListPlot[ReIm[roots],PlotLabel->[Sin[z+Sin[z+Sin[z]]]==Cos[z+Cos[z+Cos[z]]]]]
```

*([Graphics])*

## Tech Notes ▪Numerical Mathematics: Basic Operations ▪Numerical Equation Solving ▪Equations in One Variable ▪Numerical Mathematics ▪Numerical Solution of Polynomial Equations ▪Implementation notes: Numerical and Related Functions

## Related Guides ▪Equation Solving ▪Polynomial Algebra ▪Polynomial Equations ▪Precollege Education ▪Symbolic Vectors, Matrices and Arrays ▪Solvers over Regions

## History Introduced in 2021 (12.3) | Updated in 2024 (14.0)
