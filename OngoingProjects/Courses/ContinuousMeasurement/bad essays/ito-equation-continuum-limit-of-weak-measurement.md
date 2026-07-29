---
Template: Default
Title: The Itô Equation: The Continuum Limit of Weak Measurement
Author: Mads Bahrami
---

# The Itô Equation: The Continuum Limit of Weak Measurement

We derive the continuous equations for the measured position record and for the quantum state conditioned on that record. I want each nontrivial algebraic or probabilistic claim to be reproducible, so each one is followed by a focused Wolfram Language computation.

The detector is ideal and efficient: all measurement backaction is associated with its observed output, so a conditioned pure state remains pure. It is also memoryless: once the state at the beginning of a time bin is known from the outputs already observed, the next output is governed by the same measurement rule.

Let

$$
t_k=k\Delta t,
\qquad
|\psi_k\rangle=|\psi(t_k)\rangle .
$$

During the bin $[t_k,t_{k+1})$, the detector measures one real number $\bar x_k$. This is the only newly observed random quantity in the derivation. The sequence $\{\bar x_k\}$ is the finite-resolution position record.

The code uses descriptive names for the symbols in the equations:

```wl
ClearAll[
    strength,
    binWidth,
    position,
    reading,
    meanPosition,
    positionVariance,
    responsePDF,
    conditionalReadingDistribution,
    conditionalRecordMoments,
    bornMomentRules,
    recordMoments,
    leadingAccumulatedVariance,
    recordIncrementMoments,
    centeredValue,
    conditionalCenteredDistribution,
    standardizedConditionalDistribution,
    bimodalAssumptions,
    separation,
    stateWidth,
    bimodalPositionDistribution,
    bimodalReadingDensity,
    bimodalReadingDistribution,
    bimodalStandardizedDistribution,
    scale,
    power,
    totalTime,
    maximumPositionVariance,
    standardWiener,
    w,
    t,
    centeredOperator,
    centeredMean,
    centeredSecondMoment,
    hamiltonian,
    hbar,
    measurementExponent,
    measurementFactor,
    normChange,
    normExpectation,
    unitaryFactor,
    itoRules
];

assumptions = And[
    strength > 0,
    binWidth > 0,
    positionVariance >= 0,
    maximumPositionVariance >= 0,
    totalTime > 0,
    hbar > 0,
    Element[
        {
            position,
            reading,
            meanPosition,
            centeredValue,
            scale,
            positionVariance,
            maximumPositionVariance,
            totalTime,
            hbar
        },
        Reals
    ]
]
```

Thus `strength`, `binWidth`, `reading`, and `meanPosition` represent $\lambda$, $\Delta t$, $\bar x_k$, and $\langle\hat x\rangle_k$, respectively.

## One weak position reading

The measurement operator associated with the observed value $\bar x_k$ is

$$
\hat M(\bar x_k)
=
\left(\frac{2\lambda\Delta t}{\pi}\right)^{1/4}
\exp\!\left[-\lambda\Delta t(\hat x-\bar x_k)^2\right],
\qquad
\lambda>0 .
$$

The limit is taken with $\lambda$ fixed, so $\lambda\Delta t\to0$. Each measurement becomes weak, but the measurement strength per unit time remains finite.

To check the normalization, work in the position basis, where $\hat x|x\rangle=x|x\rangle$. The diagonal element of $\hat M^\dagger(\bar x_k)\hat M(\bar x_k)$ is the Gaussian response

$$
g_{\Delta t}(\bar x_k\mid x)
=
\sqrt{\frac{2\lambda\Delta t}{\pi}}
\exp\!\left[-2\lambda\Delta t(\bar x_k-x)^2\right].
$$

Its normalization, center, and variance follow from one integral:

```wl
responsePDF[reading_, position_] := Sqrt[2 strength binWidth / Pi] Exp[
    -2 strength binWidth (reading - position)^2
];

FullSimplify[
    Integrate[
        {1, reading, (reading - position)^2} responsePDF[reading, position],
        {reading, -Infinity, Infinity}
    ],
    assumptions
]
```

The result is

$$
\left\{1,\ x,\ \frac{1}{4\lambda\Delta t}\right\}.
$$

The first entry proves

$$
\int_{-\infty}^{\infty}
\hat M^\dagger(\bar x_k)\hat M(\bar x_k)\,d\bar x_k
=\hat I .
$$

The same response can be represented directly as a Wolfram Language distribution:

```wl
conditionalReadingDistribution[position_] := NormalDistribution[
    position,
    1 / (2 Sqrt[strength binWidth])
];

FullSimplify[
    {
        PDF[conditionalReadingDistribution[position], reading] == responsePDF[reading, position],
        Mean[conditionalReadingDistribution[position]],
        StandardDeviation[conditionalReadingDistribution[position]],
        Variance[conditionalReadingDistribution[position]]
    },
    assumptions
]
```

The returned list confirms the same density, mean $x$, standard deviation $1/(2\sqrt{\lambda\Delta t})$, and variance $1/(4\lambda\Delta t)$. In both the mathematics and Wolfram Language, $\mathcal N(m,s)$ or `NormalDistribution[m,s]` means a normal Gaussian distribution with mean $m$ and standard deviation $s$, not variance $s$.

The probability density for the actual detector output is

$$
p(\bar x_k\mid\psi_k)
=
\bigl\|\hat M(\bar x_k)|\psi_k\rangle\bigr\|^2
=
\int_{-\infty}^{\infty}
|\psi_k(x)|^2
g_{\Delta t}(\bar x_k\mid x)\,dx .
$$

Here $x$ is only the integration variable labeling position eigenvalues. It is not a second measured value and it is not a hidden particle position. The formula says that the detector reading is a position measurement broadened by the Gaussian response $g_{\Delta t}$.

Once $\bar x_0,\ldots,\bar x_{k-1}$ have been observed, they determine $|\psi_k\rangle$. The density above then gives the probabilities for the next, not yet observed value $\bar x_k$.

The inner Gaussian averages needed for the first two moments are:

```wl
conditionalRecordMoments = FullSimplify[
    Integrate[
        {
            reading,
            (reading - meanPosition)^2
        } responsePDF[reading, position],
        {reading, -Infinity, Infinity}
    ],
    assumptions
]
```

The result is

$$
\left\{
x,\ 
(x-\langle\hat x\rangle_k)^2
+\frac{1}{4\lambda\Delta t}
\right\}.
$$

The outer Born average needs only normalization, $\langle x\rangle_k=\langle\hat x\rangle_k$, and $\langle x^2\rangle_k=\operatorname{Var}_{\psi_k}(\hat x)+\langle\hat x\rangle_k^2$. Because the two inner moments are at most quadratic in $x$, these three identities complete the average without assuming a shape for $|\psi_k(x)|^2$:

```wl
bornMomentRules = {
    position^2 -> positionVariance + meanPosition^2,
    position -> meanPosition
};

recordMoments = FullSimplify[
    Expand[conditionalRecordMoments] /. bornMomentRules,
    assumptions
]
```

The constant terms are unchanged because the Born density integrates to one. The returned pair gives

$$
\mathbb E[\bar x_k\mid\psi_k]
=
\int x|\psi_k(x)|^2dx
=
\langle\hat x\rangle_k ,
$$

and

$$
\operatorname{Var}(\bar x_k\mid\psi_k)
=
\left\langle
(\hat x-\langle\hat x\rangle_k)^2
\right\rangle_k
+\frac{1}{4\lambda\Delta t}.
$$

No probability distribution for the state position has been assumed beyond the Born density $|\psi_k(x)|^2$.

Given the observed $\bar x_k$, the normalized state update is

$$
|\psi_{k+1}\rangle
=
\frac{\hat U\,\hat M(\bar x_k)|\psi_k\rangle}
{\bigl\|\hat U\,\hat M(\bar x_k)|\psi_k\rangle\bigr\|},
\qquad
\hat U
=
\exp\!\left(-\frac{i}{\hbar}\hat H\Delta t\right),
\qquad
\hat H=\hat H^\dagger .
$$

## The record variable that has a continuous limit

The variance of $\bar x_k$ diverges like $1/\Delta t$, so the individual readings do not approach finite values as $\Delta t\to0$. We therefore seek a scaled increment whose variance is proportional to $\Delta t$. Such increments vanish individually but accumulate a finite variance over a fixed duration.

The power of $\Delta t$ is fixed by this requirement. Keeping only the leading detector variance, scaling each reading by $\Delta t^\alpha$ gives accumulated variance

$$
\frac{T}{\Delta t}
\times
\frac{\Delta t^{2\alpha}}{4\lambda\Delta t}
=
\frac{T}{4\lambda}\Delta t^{2\alpha-2}.
$$

The three relevant cases are computed directly:

```wl
leadingAccumulatedVariance[power_] := (totalTime / binWidth) binWidth^(2 power) / (4 strength binWidth);

FullSimplify[
    {
        Limit[
            leadingAccumulatedVariance[1 / 2],
            binWidth -> 0,
            Direction -> "FromAbove"
        ],
        leadingAccumulatedVariance[1],
        Limit[
            leadingAccumulatedVariance[3 / 2],
            binWidth -> 0,
            Direction -> "FromAbove"
        ]
    },
    assumptions
]
```

The result is $\{\infty,T/(4\lambda),0\}$. A smaller power makes the accumulated variance diverge, a larger power removes the noise, and $\alpha=1$ is the unique finite nonzero choice.

The remaining positive scale factor is chosen so that the leading variance over time $T$ is exactly $T$:

```wl
FullSimplify[
    Solve[
        scale^2 / (4 strength) == 1 && scale > 0,
        scale,
        Reals
    ],
    assumptions
]
```

The result is $2\sqrt{\lambda}$. This motivates defining the accumulated record $Y$ by

$$
\Delta Y_k
\equiv
Y_{t_{k+1}}-Y_{t_k}
=
2\sqrt{\lambda}\,\bar x_k\Delta t,
\qquad
Y_0=0 .
$$

Only $\bar x_k$ is measured. $Y$ is calculated from the measured readings, and the calculation is exactly invertible:

$$
\bar x_k
=
\frac{\Delta Y_k}{2\sqrt{\lambda}\,\Delta t}.
$$

The mean and variance of one record increment follow from the moments of $\bar x_k$:

```wl
recordIncrementMoments = FullSimplify[
    {
        2 Sqrt[strength] binWidth recordMoments[[1]],
        (2 Sqrt[strength] binWidth)^2 recordMoments[[2]]
    },
    assumptions
]
```

Thus

$$
\mathbb E[\Delta Y_k\mid\psi_k]
=
2\sqrt{\lambda}\,
\langle\hat x\rangle_k\Delta t ,
$$

and

$$
\operatorname{Var}(\Delta Y_k\mid\psi_k)
=
\Delta t
+4\lambda\Delta t^2
\left\langle
(\hat x-\langle\hat x\rangle_k)^2
\right\rangle_k .
$$

## The centered record fluctuation

The deterministic part of $\Delta Y_k$ is its conditional mean. Subtracting it isolates the unpredictable part of the same measured record:

$$
\begin{aligned}
\Delta W_k
&\equiv
\Delta Y_k
-2\sqrt{\lambda}\,
\langle\hat x\rangle_k\Delta t\\
&=
2\sqrt{\lambda}\,\Delta t
\bigl(\bar x_k-\langle\hat x\rangle_k\bigr).
\end{aligned}
$$

This is a definition, not a new measurement and not an independent noise added afterward. After $\bar x_k$ is observed, $\Delta W_k$ is calculated from it. Before the reading occurs, the definition is a one-to-one change of random variable:

$$
\bar x_k
=
\langle\hat x\rangle_k
+\frac{\Delta W_k}{2\sqrt{\lambda}\,\Delta t}.
$$

At fixed position eigenvalue $x$, `TransformedDistribution` performs precisely this change of variable:

```wl
conditionalCenteredDistribution[position_] := TransformedDistribution[
    2 Sqrt[strength] binWidth (reading - meanPosition),
    Distributed[
        reading,
        conditionalReadingDistribution[position]
    ],
    Assumptions -> assumptions
];

FullSimplify[
    conditionalCenteredDistribution[position],
    assumptions
]
```

It returns

$$
\mathcal N\!\left(
2\sqrt{\lambda}\,\Delta t
(x-\langle\hat x\rangle_k),
\sqrt{\Delta t}
\right).
$$

Equivalently, the inverse map has Jacobian

$$
\left|\frac{d\bar x_k}{d(\Delta W_k)}\right|
=
\frac{1}{2\sqrt{\lambda}\,\Delta t}.
$$

The transformed conditional density is computed by:

```wl
FullSimplify[
    PDF[
        conditionalCenteredDistribution[position],
        centeredValue
    ],
    assumptions
]
```

Writing the returned expression with $w$ for a possible value of $\Delta W_k$ gives

$$
\frac{1}{\sqrt{2\pi\Delta t}}
\exp\!\left[
-\frac{
\bigl(
w-2\sqrt{\lambda}\,\Delta t
(x-\langle\hat x\rangle_k)
\bigr)^2
}{2\Delta t}
\right].
$$

Averaging this conditional density with the Born density gives the exact finite-bin distribution

$$
p_{\Delta W_k}(w\mid\psi_k)
=
\frac{1}{\sqrt{2\pi\Delta t}}
\int_{-\infty}^{\infty}
|\psi_k(x)|^2
\exp\!\left[
-\frac{
\bigl(
w-2\sqrt{\lambda}\,\Delta t
(x-\langle\hat x\rangle_k)
\bigr)^2
}{2\Delta t}
\right]dx .
$$

To see its limiting form, divide $\Delta W_k$ by $\sqrt{\Delta t}$. The conditional distribution at fixed $x$ becomes:

```wl
standardizedConditionalDistribution[position_] := TransformedDistribution[
    centeredValue / Sqrt[binWidth],
    Distributed[
        centeredValue,
        conditionalCenteredDistribution[position]
    ],
    Assumptions -> assumptions
];

FullSimplify[
    {
        standardizedConditionalDistribution[position],
        Limit[
            4 strength binWidth positionVariance,
            binWidth -> 0,
            Direction -> "FromAbove"
        ]
    },
    assumptions
]
```

The first result is

$$
\mathcal N\!\left(
2\sqrt{\lambda\Delta t}\,
(x-\langle\hat x\rangle_k),
1
\right).
$$

The mean-square size of its state-dependent center is

$$
4\lambda\Delta t
\left\langle
(\hat x-\langle\hat x\rangle_k)^2
\right\rangle_k ,
$$

and the second result confirms that it tends to zero when the position variance remains finite. Therefore

$$
\frac{\Delta W_k}{\sqrt{\Delta t}}
\sim
\mathcal N(0,1)
\qquad
\text{to leading order},
$$

where $\sim$ means "is distributed as."

The exact first two finite-bin moments use the previously computed record moments:

```wl
FullSimplify[
    {
        recordIncrementMoments[[1]] - 2 Sqrt[strength] binWidth meanPosition,
        recordIncrementMoments[[2]]
    },
    assumptions
]
```

Thus

$$
\mathbb E[\Delta W_k\mid\psi_k]=0,
\qquad
\operatorname{Var}(\Delta W_k\mid\psi_k)
=
\Delta t
+4\lambda\Delta t^2
\left\langle
(\hat x-\langle\hat x\rangle_k)^2
\right\rangle_k .
$$

The derivation did not assume a Gaussian quantum state. To test that point, take a symmetric bimodal Born density, the equal mixture of $\mathcal N(-a,s)$ and $\mathcal N(a,s)$ with $a,s>0$. It has mean zero and variance $a^2+s^2$, but it is not Gaussian.

The next cell integrates this Born density against the detector response. It then verifies the resulting measured-output density, transforms the measured reading to $\Delta W_k/\sqrt{\Delta t}$, and computes its variance and excess kurtosis:

```wl
bimodalAssumptions = And[
    assumptions,
    separation > 0,
    stateWidth > 0,
    Element[{separation, stateWidth}, Reals]
];

bimodalPositionDistribution = MixtureDistribution[
    {1 / 2, 1 / 2},
    {
        NormalDistribution[-separation, stateWidth],
        NormalDistribution[separation, stateWidth]
    }
];

bimodalReadingDensity = FullSimplify[
    Integrate[
        PDF[bimodalPositionDistribution, position] responsePDF[reading, position],
        {position, -Infinity, Infinity}
    ],
    bimodalAssumptions
];

bimodalReadingDistribution = MixtureDistribution[
    {1 / 2, 1 / 2},
    {
        NormalDistribution[
            -separation,
            Sqrt[stateWidth^2 + 1 / (4 strength binWidth)]
        ],
        NormalDistribution[
            separation,
            Sqrt[stateWidth^2 + 1 / (4 strength binWidth)]
        ]
    }
];

bimodalStandardizedDistribution = TransformedDistribution[
    2 Sqrt[strength binWidth] reading,
    Distributed[reading, bimodalReadingDistribution],
    Assumptions -> bimodalAssumptions
];

FullSimplify[
    {
        PDF[bimodalReadingDistribution, reading] == bimodalReadingDensity,
        Variance[bimodalStandardizedDistribution],
        Kurtosis[bimodalStandardizedDistribution] - 3,
        Limit[
            Kurtosis[bimodalStandardizedDistribution],
            binWidth -> 0,
            Direction -> "FromAbove"
        ]
    },
    bimodalAssumptions
]
```

The first result is `True`. The other results are

$$
\operatorname{Var}\!\left(
\frac{\Delta W_k}{\sqrt{\Delta t}}
\right)
=
1+4\lambda\Delta t(a^2+s^2),
$$

$$
\operatorname{Kurt}\!\left(
\frac{\Delta W_k}{\sqrt{\Delta t}}
\right)-3
=
-\frac{
32\lambda^2\Delta t^2a^4
}{
\left[1+4\lambda\Delta t(a^2+s^2)\right]^2
},
\qquad
\lim_{\Delta t\to0}\operatorname{Kurt}=3 .
$$

The nonzero excess kurtosis shows that the exact finite-bin record need not be Gaussian. Its vanishing limit agrees with the general Wiener limit.

A numerical check should start from simulated detector readings, not from the target normal distribution. Here the first returned row contains the sampled mean, variance, and kurtosis of $\Delta W_k/\sqrt{\Delta t}$ for $\lambda=1$, $\Delta t=1/10$, $a=2$, and $s=1/2$. The second row contains their exact values:

```wl
BlockRandom[
    SeedRandom[20260727];
    With[{
        parameters = {
            strength -> 1,
            binWidth -> 1 / 10,
            separation -> 2,
            stateWidth -> 1 / 2
        }
    },
        With[{
            measuredReadings = RandomVariate[
                bimodalReadingDistribution /. parameters,
                50000
            ]
        },
            With[{
                standardizedRecord = 2 Sqrt[(strength binWidth) /. parameters] measuredReadings
            },
                N @ {
                    {
                        Mean[standardizedRecord],
                        Variance[standardizedRecord],
                        Kurtosis[standardizedRecord]
                    },
                    {0, 27 / 10, 1675 / 729}
                }
            ]
        ]
    ]
]
```

The sampled row is approximately $\{0.00078,2.7229,2.2967\}$, close to the exact row $\{0,2.7,2.2977\}$. This checks the path from measured readings to the centered record; it is not used in the derivation.

## From finite-bin fluctuations to a Wiener process

Define the accumulated centered record at finite resolution by

$$
W_{t_n}^{(\Delta t)}
\equiv
\sum_{k=0}^{n-1}\Delta W_k
=
Y_{t_n}
-2\sqrt{\lambda}
\sum_{k=0}^{n-1}
\langle\hat x\rangle_k\Delta t .
$$

Between grid points, take $W_t^{(\Delta t)}$ to remain constant until the next detector reading. This makes $W^{(\Delta t)}$ a defined finite-resolution process, not only a list of increments.

At bin $k$, the previously observed readings determine $|\psi_k\rangle$, and the next increment has conditional mean zero. The increments need not be independent at finite $\Delta t$ because each new reading changes the state used in the next bin.

Over a fixed interval $T=n\Delta t$, the leading conditional variances sum to $T$. If the position variance is bounded above by `maximumPositionVariance` on that interval, the correction is bounded by

$$
4\lambda T\Delta t\,
\text{maximumPositionVariance}.
$$

Its limit is:

```wl
FullSimplify[
    Limit[
        4 strength totalTime binWidth
            maximumPositionVariance,
        binWidth -> 0,
        Direction -> "FromAbove"
    ],
    assumptions
]
```

The result is zero. To control rare large increments as well, assume that the conditional position moments of some order $2+\epsilon$ remain uniformly bounded on every finite time interval, for some $\epsilon>0$. The martingale functional central-limit theorem then gives

$$
W^{(\Delta t)}\Longrightarrow W
\qquad
(\Delta t\to0),
$$

where the arrow denotes convergence in distribution of the accumulated process and $W_t$ is a standard Wiener process. This theorem is the one probabilistic limit used here; Wolfram Language checks the finite-bin inputs and the properties of the limiting process, not the theorem itself.

For comparison, the built-in standard Wiener process gives the distribution of an increment of duration $\Delta t$ directly:

```wl
standardWiener = WienerProcess[];

FullSimplify[
    TransformedDistribution[
        w[t + binWidth] - w[t],
        Distributed[w, standardWiener],
        Assumptions -> t >= 0 && binWidth > 0
    ],
    assumptions && t >= 0
]
```

Inside this cell, `w` is only the local process symbol required by `Distributed`. The result is `NormalDistribution[0,Sqrt[binWidth]]`, exactly the limiting distribution $\mathcal N(0,\sqrt{\Delta t})$.

We now make the single notational transition from finite bins to differentials:

$$
\Delta t\longrightarrow dt,
\qquad
\Delta Y_k\longrightarrow dY_t,
\qquad
\Delta W_k\longrightarrow dW_t .
$$

Because $W_t$ is a Wiener process, its Itô increments satisfy

$$
(dW_t)^2=dt,
\qquad
dt\,dW_t=0,
\qquad
dt^2=0 .
$$

These rules specify the Itô convention. No Stratonovich equation is being used. The exact finite-bin definition of $\Delta W_k$ now becomes the record SDE

$$
\boxed{
dY_t
=
2\sqrt{\lambda}\,
\langle\hat x\rangle_t\,dt
+dW_t
}.
$$

The measured finite-bin values are still recovered from the accumulated record by

$$
\bar x_k
=
\frac{Y_{t_{k+1}}-Y_{t_k}}
{2\sqrt{\lambda}\,\Delta t}.
$$

Thus the experiment measures the sequence $\{\bar x_k\}$, $Y$ is the same data in the scaling that has a continuous limit, and $dW_t$ is the centered random part of that limiting record.

## The conditioned state equation

For the operator expansions below, assume that the states considered lie in a common invariant analytic domain for $\hat H$ and $\hat x$. This is the regularity needed to interpret the uppercase-$O$ remainders when these operators are unbounded.

Return to the exact state update and substitute

$$
\bar x_k
=
\langle\hat x\rangle_k
+\frac{\Delta W_k}
{2\sqrt{\lambda}\,\Delta t}.
$$

Define

$$
\hat A_k
\equiv
\hat x-\langle\hat x\rangle_k .
$$

The measurement exponent becomes

$$
\begin{aligned}
-\lambda\Delta t(\hat x-\bar x_k)^2
&=
-\lambda\Delta t
\left(
\hat A_k
-\frac{\Delta W_k}
{2\sqrt{\lambda}\,\Delta t}
\right)^2\\
&=
\sqrt{\lambda}\,\hat A_k\Delta W_k
-\lambda\hat A_k^2\Delta t
-\frac{(\Delta W_k)^2}{4\Delta t}.
\end{aligned}
$$

The identity is purely algebraic. Using `centeredOperator` for $\hat A_k$ and `centeredValue` for $\Delta W_k$, Wolfram Language verifies it exactly:

```wl
FullSimplify[
    -strength binWidth (centeredOperator - centeredValue / (2 Sqrt[strength] binWidth))^2 == Sqrt[strength] centeredOperator centeredValue - strength centeredOperator^2 binWidth - centeredValue^2 / (4 binWidth),
    assumptions
]
```

The result is `True`. The last term is a scalar, so it multiplies the whole unnormalized state by the positive number

$$
\exp\!\left[-\frac{(\Delta W_k)^2}{4\Delta t}\right].
$$

The same factor appears in the norm and cancels exactly. The scalar ratio is:

```wl
FullSimplify[
    Exp[-centeredValue^2 / (4 binWidth)] / Sqrt[Exp[-centeredValue^2 / (2 binWidth)]],
    assumptions
]
```

The result is $1$. Therefore the state-changing part of the measurement operator is

$$
\exp\!\left(
\sqrt{\lambda}\,\hat A_k\Delta W_k
-\lambda\hat A_k^2\Delta t
\right).
$$

Since $\Delta W_k$ is of order $\sqrt{\Delta t}$, its square contributes at order $\Delta t$. Expand the exponential through second order and then apply the Itô-order rules:

```wl
measurementExponent = Sqrt[strength] centeredOperator centeredValue - strength centeredOperator^2 binWidth;

itoRules = {
    centeredValue^2 -> binWidth,
    centeredValue binWidth -> 0,
    binWidth^2 -> 0
};

measurementFactor = Expand[
    Total[{1, measurementExponent, measurementExponent^2 / 2}]
] /. itoRules // Expand
```

The returned factor is

$$
1
+\sqrt{\lambda}\,\hat A_k\Delta W_k
-\frac{\lambda}{2}\hat A_k^2\Delta t .
$$

The coefficient $-1/2$ comes from combining the explicit term $-\lambda\hat A_k^2\Delta t$ with the quadratic contribution $+\lambda\hat A_k^2(\Delta W_k)^2/2$ and then using $(\Delta W_k)^2=\Delta t$ at Itô order.

This factor preserves the norm to the same order. Squaring it, applying the Itô rules, and then replacing powers of `centeredOperator` by their state expectations gives:

```wl
normChange = Expand[measurementFactor^2 - 1] /. itoRules // Expand;

normExpectation = normChange /. {
    centeredOperator^2 -> centeredSecondMoment,
    centeredOperator -> centeredMean
};

{
    normExpectation,
    normExpectation /. centeredMean -> 0
}
```

The first expression is $2\sqrt{\lambda}\,\langle\hat A_k\rangle_k\Delta W_k$, and the second is zero because

$$
\langle\hat A_k\rangle_k
=
\langle\hat x-\langle\hat x\rangle_k\rangle_k
=0 .
$$

The unitary factor is

$$
\hat U
=
\hat I
-\frac{i}{\hbar}\hat H\Delta t
+O(\Delta t^2).
$$

Products of its order-$\Delta t$ term with the order-$\sqrt{\Delta t}$ measurement term are $O(\Delta t^{3/2})$ and do not survive. The complete retained factor is computed by:

```wl
unitaryFactor = 1 - I hamiltonian binWidth / hbar;

Expand[
    unitaryFactor measurementFactor
] /. itoRules // Expand
```

The result is

$$
\hat I
-\frac{i}{\hbar}\hat H\Delta t
-\frac{\lambda}{2}\hat A_k^2\Delta t
+\sqrt{\lambda}\,\hat A_k\Delta W_k .
$$

Replacing the finite increments by the differentials already defined gives the conditioned state SDE:

$$
\boxed{
d|\psi_t\rangle
=
\left[
-\frac{i}{\hbar}\hat H\,dt
-\frac{\lambda}{2}
\bigl(\hat x-\langle\hat x\rangle_t\bigr)^2dt
+\sqrt{\lambda}
\bigl(\hat x-\langle\hat x\rangle_t\bigr)dW_t
\right]|\psi_t\rangle
}.
$$

The same $dW_t$ appears in both boxed equations because both are derived from the same observed reading $\bar x_k$. The record equation describes what the detector data become in the continuous limit. The state equation describes how that same data condition the quantum state.
