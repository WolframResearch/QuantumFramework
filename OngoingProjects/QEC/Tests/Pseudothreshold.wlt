(* ::Package:: *)

(* ============================================================================
   Tests/QEC/Pseudothreshold.wlt

   The pseudothreshold (Got26 sec. 14.7.6) as an exact root, and as an object
   that will not be read as a bare threshold.

   The load-bearing checks: the repetition code under bit flips crosses at
   exactly 1/2 (its rate is the binomial tail 3p^2 - 2p^3); below a crossing the
   encoded rate is lower and above it higher; and a code that does not help at
   low noise -- a phase flip on the bit-flip repetition code, or the bare-ancilla
   circuit whose rate starts at 32p/15 -- reports that, instead of handing back
   the crossing where encoding starts to help as if it were where it stops.
   ============================================================================ *)

Needs["Wolfram`QuantumFramework`"];
Needs["Wolfram`QuantumFramework`PackageScope`"];

Get[FileNameJoin[{
    ParentDirectory[DirectoryName[FindFile["Wolfram`QuantumFramework`"], 2]],
    "OngoingProjects", "QEC", "QECCore", "QECCore.wl"
}]];

QECClearCache[];

steane = QECCode["SteaneCode"];
rep3 = QECCode["Repetition", 3];

steanePs = QECPseudothreshold[steane, QECNoiseModel["Depolarizing", qecP]];


VerificationTest[
    QECPseudothreshold[rep3, QECNoiseModel["BitFlip", qecP]]["Value"],
    1/2,
    TestID -> "QEC-Pseudothreshold-repetition-code-crosses-at-one-half"
]

VerificationTest[
    {Head[steanePs["Value"]], 0.08 < N[steanePs["Value"]] < 0.082},
    {Root, True},
    TestID -> "QEC-Pseudothreshold-Steane-is-an-exact-root"
]

VerificationTest[
    With[{v = N[steanePs["Value"]], rate = steanePs["LogicalErrorRate"]},
        {(rate /. qecP -> v / 2) < v / 2, (rate /. qecP -> 2 v) > 2 v, Abs[(rate /. qecP -> v) - v] < 10^-12}
    ],
    {True, True, True},
    TestID -> "QEC-Pseudothreshold-encoding-helps-below-and-hurts-above"
]

VerificationTest[
    Simplify[steanePs["LogicalErrorRate"] - QECLogicalErrorRate[steane, QECNoiseModel["Depolarizing", qecP]]],
    0,
    TestID -> "QEC-Pseudothreshold-carries-the-rate-it-came-from"
]

VerificationTest[
    QECPseudothreshold[rep3, QECNoiseModel["PhaseFlip", qecP]]["Value"],
    Missing["EncodingDoesNotHelp"],
    TestID -> "QEC-Pseudothreshold-no-value-when-encoding-does-not-help"
]

(* The bare-ancilla circuit: the rate starts at 32p/15, so encoding hurts at low
   noise.  The polynomial does cross p, further up -- that is reported as a
   crossing, not as the pseudothreshold. *)
VerificationTest[
    With[{ps = QECPseudothreshold[QECCode["BitFlipCode"], QECNoiseModel["Circuit", qecP], "Rounds" -> 1]},
        {ps["Value"], ps["LowNoiseRatio"], Length[ps["Crossings"]], ps["Rounds"], ps["Extraction"]}
    ],
    {Missing["EncodingDoesNotHelp"], 32/15, 1, 1, "BareAncilla"},
    TestID -> "QEC-Pseudothreshold-the-hook-error-has-no-pseudothreshold"
]

VerificationTest[
    {Length[steanePs["Caveats"]] >= 2, StringContainsQ[First[steanePs["Caveats"]], "not a threshold"]},
    {True, True},
    TestID -> "QEC-Pseudothreshold-says-it-is-not-a-threshold"
]

VerificationTest[
    {steanePs["Level"], steanePs["Rounds"], steanePs["Parameter"], Head[steanePs["Code"]], Head[steanePs["Noise"]]},
    {"CodeCapacity", None, qecP, QECCode, QECNoiseModel},
    TestID -> "QEC-Pseudothreshold-carries-its-context"
]

VerificationTest[
    QECPseudothreshold[rep3, QuantumChannel["BitFlip"[qecQ]]]["Value"],
    1/2,
    TestID -> "QEC-Pseudothreshold-accepts-a-QuantumChannel"
]

VerificationTest[
    QECPseudothreshold[steane, QECNoiseModel["Depolarizing", 1/10]],
    $Failed,
    {QECPseudothreshold::params},
    TestID -> "QEC-Pseudothreshold-needs-a-symbolic-rate"
]
