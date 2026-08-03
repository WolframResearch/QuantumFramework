# FourierParameters | [SpanFromLeft]

> [FourierParameters](https://reference.wolfram.com/language/ref/FourierParameters.html)  — is an option to [Fourier](https://reference.wolfram.com/language/ref/Fourier.html) and related functions that specifies the conventions to use in computing Fourier transforms.

## Details

A typical setting is [FourierParameters](https://reference.wolfram.com/language/ref/FourierParameters.html)->{*a*,*b*}.

Some common choices for `{*a*,*b*}` are `{0,1}` (default), `{-1,1}` (data analysis), `{1,-1}` (signal processing).

## Examples

### Basic Examples

Use a nondefault definition of the discrete Fourier transform:

```wolfram
Fourier[{1,2,3,4,5,6},FourierParameters->{1,1}]
(* Output *)
{21.+0. ⅈ,-3.0000000000000018-5.196152422706633 ⅈ,-3.-1.7320508075688783 ⅈ,-3.+0. ⅈ,-3.+1.7320508075688783 ⅈ,-3.0000000000000018+5.196152422706633 ⅈ}
```

Use the same definition to get the inverse:

```wolfram
InverseFourier[%,FourierParameters->{1,1}]
(* Output *)
{0.9999999999999993,1.9999999999999991,3.0000000000000004,4.,5.,6.000000000000001}
```

A nondefault definition used for the continuous Fourier transform:

```wolfram
FourierTransform[Exp[-t^2],t,ω,FourierParameters->{1,1}]
(* Output *)
ℯ^(-(ω^2)/(4)) Sqrt[π]
```

```wolfram
InverseFourierTransform[%,ω,t,FourierParameters->{1,1}]
(* Output *)
ℯ^(-t^2)
```

### Scope

A typical pure mathematics or systems-engineering definition of Fourier transform:

```wolfram
FourierTransform[Exp[-Abs[t]],t,ω,FourierParameters->{1,-1}]
(* Output *)
(2)/(1+ω^2)
```

Use the same definition for the inverse transform:

```wolfram
InverseFourierTransform[%,ω,t,FourierParameters->{1,-1}]
(* Output *)
ℯ^(-Abs[t])
```

A typical data-analysis definition of discrete Fourier transform:

```wolfram
Fourier[Range[6],FourierParameters->{-1,1}]
(* Output *)
{3.5+0. ⅈ,-0.5000000000000002-0.8660254037844388 ⅈ,-0.5-0.28867513459481303 ⅈ,-0.5+0. ⅈ,-0.5+0.28867513459481303 ⅈ,-0.5000000000000002+0.8660254037844388 ⅈ}
```

Use the same definition to get the correct inverse:

```wolfram
InverseFourier[%,FourierParameters->{-1,1}]
(* Output *)
{0.9999999999999996,1.9999999999999991,3.,4.000000000000001,5.,6.000000000000001}
```

A common signal-processing definition of Fourier transform:

```wolfram
FourierTransform[Exp[-Abs[t]],t,ω,FourierParameters->{1,-1}]
(* Output *)
(2)/(1+ω^2)
```

Discrete Fourier transform:

```wolfram
Fourier[Range[6],FourierParameters->{1,-1}]
(* Output *)
{21.+0. ⅈ,-3.+5.196152422706632 ⅈ,-3.+1.7320508075688772 ⅈ,-3.+0. ⅈ,-3.-1.7320508075688772 ⅈ,-3.-5.196152422706632 ⅈ}
```

Discrete-time Fourier transform:

```wolfram
FourierSequenceTransform[(1/2)^n UnitStep[n],n,ω,FourierParameters->{1,1}]
(* Output *)
(2 ℯ^(ⅈ ω))/(-1+2 ℯ^(ⅈ ω))
```

### Possible Issues

The same [FourierParameters](https://reference.wolfram.com/language/ref/FourierParameters.html) values need to be used for both forward and inverse transforms:

```wolfram
Fourier[{1,2,3,4,5,6},FourierParameters->{-1,1}]
(* Output *)
{3.5+0. ⅈ,-0.5000000000000002-0.8660254037844388 ⅈ,-0.5-0.28867513459481303 ⅈ,-0.5+0. ⅈ,-0.5+0.28867513459481303 ⅈ,-0.5000000000000002+0.8660254037844388 ⅈ}
```

Here the inverse uses a different choice of [FourierParameters](https://reference.wolfram.com/language/ref/FourierParameters.html):

```wolfram
InverseFourier[%]
(* Output *)
{0.4082482904638629,0.8164965809277258,1.2247448713915892,1.6329931618554527,2.041241452319315,2.4494897427831788}
```

The second parameter needs to be relatively prime to the data length to guarantee invertibility:

```wolfram
Fourier[{1,0,0,0,0,0},FourierParameters->{-1,3}]
(* Output *)
{0.16666666666666666,0.16666666666666663,0.16666666666666663,0.16666666666666666,0.16666666666666666,0.16666666666666663}
```

```wolfram
InverseFourier[%,FourierParameters->{-1,3}]//Chop
(* Output *)
InverseFourier
(* Output *)
{0.9999999999999998,0,0.9999999999999998,0,0.9999999999999998,0}
```

## Related Guides ▪Fourier Analysis

## History Introduced in 1999 (4.0)
