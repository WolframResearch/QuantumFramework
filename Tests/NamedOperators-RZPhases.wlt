(* R_Z(theta) = diag(e^(-i theta/2), e^(i theta/2)). Its phases are passed through
   FullSimplify only where that can change their form: an exact angle, or a symbolic one
   holding a number or a numeric constant such as Pi. A plain symbolic angle or an inexact
   number skips it, so a new symbolic R_Z costs no more than a general 2x2 operator, and
   every angle keeps the matrix it had.

   References are independent of the route under test: the phases passed through
   FullSimplify, which is the stored form for every angle, and MatrixExp of
   -i theta/2 sigma_z.

   Regimes covered:
     general symbolic     theta, sums, multiples, theta[1], f[t], theta^2, Sqrt[2] theta;
                          substituting a value afterwards, the derivative in theta
     exactly solvable     0, 1/3, Sqrt[2], Pi/2, Pi/7, ArcCos[1/3], the default Pi/2
     mixed                theta + Pi/2, theta + 2 Pi, i theta, Degree theta, 0.5 theta,
                          0.3 + theta
     numerical reference  0.3, 1.5 + 0.5 i, a 30-digit 1/3, against MatrixExp
     scaling              a distinct symbolic R_Z costs less than one general 2x2
                          operator

   Qudit R_Z is built by another route; RotationWires-RZ-qudit-* in
   NamedOperators-RotationWires.wlt checks it. *)

rzMatrix[args___] := Normal[QuantumOperator["RZ"[args], {2}]["Matrix"]]

simplifiedForm[angle_] := DiagonalMatrix[FullSimplify[{Exp[- I angle / 2], Exp[I angle / 2]}]]

zeroQ[m_] := AllTrue[Flatten[FullSimplify[m]], # === 0 &]

numZeroQ[m_] := Max[Abs[N[Flatten[m]]]] < 10^-12

$angles = {
    \[Theta], -\[Theta], 2 \[Theta], \[Theta] / 3, \[Theta] + \[Phi], \[Theta] - \[Phi], a \[Theta] + b,
    \[Theta][1], \[Theta][1] + \[Theta][2], f[t], \[Theta]^2, Sqrt[2] \[Theta],
    0, 1 / 3, Sqrt[2], Pi / 2, Pi / 7, ArcCos[1 / 3],
    \[Theta] + Pi / 2, \[Theta] + 2 Pi, I \[Theta], Degree \[Theta], 0.5 \[Theta], 0.3 + \[Theta],
    0.3, 1.5 + 0.5 I, N[1 / 3, 30]
};

Scan[
    Function[angle,
        VerificationTest[
            {
                rzMatrix[angle] === simplifiedForm[angle],
                Normal[QuantumCircuitOperator[{"RZ"[angle] -> 3}]["Elements"][[1]]["Matrix"]] === simplifiedForm[angle],
                If[ InexactNumberQ[angle] || ! FreeQ[angle, _Real],
                    numZeroQ[rzMatrix[angle] - MatrixExp[- I angle / 2 PauliMatrix[3]] /. {\[Theta] -> 0.7, \[Phi] -> -0.2}],
                    zeroQ[rzMatrix[angle] - MatrixExp[- I angle / 2 PauliMatrix[3]]]
                ]
            },
            {True, True, True},
            TestID -> "RZPhases-form-and-matrix-" <> ToString[angle, InputForm]
        ]
    ],
    $angles
]

VerificationTest[
    {rzMatrix[] === simplifiedForm[Pi / 2], Precision[rzMatrix[N[1 / 3, 30]]] >= 30},
    {True, True},
    TestID -> "RZPhases-default-angle-and-precision"
]

VerificationTest[
    {
        numZeroQ[(rzMatrix[\[Theta]] /. \[Theta] -> 0.7) - rzMatrix[0.7]],
        (rzMatrix[\[Theta]] /. \[Theta] -> Pi / 2) === rzMatrix[Pi / 2] || zeroQ[(rzMatrix[\[Theta]] /. \[Theta] -> Pi / 2) - rzMatrix[Pi / 2]],
        zeroQ[D[rzMatrix[\[Theta]], \[Theta]] - (- I / 2) PauliMatrix[3] . rzMatrix[\[Theta]]],
        zeroQ[rzMatrix[\[Theta]] . rzMatrix[\[Phi]] - rzMatrix[\[Theta] + \[Phi]]]
    },
    {True, True, True, True},
    TestID -> "RZPhases-substitution-derivative-group-law"
]

(* A distinct symbolic R_Z is one diagonal 2x2 matrix placed on one wire, so it costs less
   than a general 2x2 operator built through the generic constructor. With its phases
   passed through FullSimplify it cost about twice as much. *)
VerificationTest[
    With[{k = 300},
        With[
            {
                general = Min @ Table[ClearSystemCache[]; First @ AbsoluteTiming[Table[QuantumOperator[{{Cos[Unique[]], 0}, {0, 1}}], k]], 3],
                rotation = Min @ Table[ClearSystemCache[]; First @ AbsoluteTiming[QuantumCircuitOperator[Table["RZ"[Unique[]] -> Mod[j, 10] + 1, {j, k}]]], 3]
            },
            rotation / general < 1
        ]
    ],
    True,
    TestID -> "RZPhases-distinct-rotation-cost"
]
