(* Baseline times for Tests/QuantumEntanglementTiming.wlt, recorded 2026-10-18 on 96f48907 (paclet 2.1.3, WL 15.0.1, this Mac),
   with no other kernel running: minimum of 7 runs, ClearSystemCache[] before each, QF property cache off.
   "BaselineAnswer" is "correct", "wrong" (fast but wrong) or "unfinished" (no answer within 300 s).
   "Limit" (seconds) applies only to wrong or unfinished cases. G10-G13: 1 s each, approved 2026-10-21 after their fixes
   (measured 0.02-0.1 s); G14 by the same rule once fixed the same day. G15: PROPOSED, pending the exact-input fix. *)

<|
    "G1" -> <|"BaselineAnswer" -> "correct", "Seconds" -> 43.606073, "Answer" -> "67 -> 0"|>,
    "G2" -> <|"BaselineAnswer" -> "correct", "Seconds" -> 0.3264, "Answer" -> "1.3485152445312587"|>,
    "G3" -> <|"BaselineAnswer" -> "correct", "Seconds" -> 0.911521, "Answer" -> "1.3688374590941177"|>,
    "G4" -> <|"BaselineAnswer" -> "correct", "Seconds" -> 0.325909, "Answer" -> "0.4455240984845671"|>,
    "G5" -> <|"BaselineAnswer" -> "correct", "Seconds" -> 0.301719, "Answer" -> "{Quantity[-1/2*(((1 + Abs[Sqrt[2/15] - Sqrt[3/10]/2 - 1/(2*Sqrt[30])]^2)^(-1) - Sqrt[221]/(15*(1 + Abs[Sqrt[2/15] - Sqrt[3/10]/2 - 1/(2*Sqrt[30])]^2)))*Log[((1 + Abs[Sqrt[2/15] - Sqrt[3/10]/2 - 1/(2*S"|>,
    "G6" -> <|"BaselineAnswer" -> "correct", "Seconds" -> 0.020921, "Answer" -> "{True, True}"|>,
    "G7" -> <|"BaselineAnswer" -> "correct", "Seconds" -> 0.251259, "Answer" -> "{True, True, True, 0.04671624156386511, 0.08638311263645726}"|>,
    "G8" -> <|"BaselineAnswer" -> "correct", "Seconds" -> 0.046834, "Answer" -> "0.004262206555527337"|>,
    "G9a" -> <|"BaselineAnswer" -> "correct", "Seconds" -> 0.436661, "Answer" -> "False"|>,
    "G9b" -> <|"BaselineAnswer" -> "correct", "Seconds" -> 2.868305, "Answer" -> "True"|>,
    "G9c" -> <|"BaselineAnswer" -> "correct", "Seconds" -> 1.014617, "Answer" -> "True"|>,
    "G10" -> <|"BaselineAnswer" -> "wrong", "Seconds" -> 0.069619, "Answer" -> "{1, (2*Sqrt[2])/3}", "Limit" -> 1.|>,
    "G11" -> <|"BaselineAnswer" -> "unfinished", "Seconds" -> Infinity, "Answer" -> "$TimedOut", "Limit" -> 1.|>,
    "G12" -> <|"BaselineAnswer" -> "wrong", "Seconds" -> 0.042674, "Answer" -> "Sqrt[2]*Sqrt[1 - (a^2*Conjugate[a]^2 + b^2*Conjugate[b]^2 + 2*a*c*Conjugate[a]*Conjugate[c] + c^2*Conjugate[c]^2)/(a*Conjugate[a] + b*Conjugate[b] + c*Conjugate[c])^2]", "Limit" -> 1.|>,
    "G13" -> <|"BaselineAnswer" -> "unfinished", "Seconds" -> Infinity, "Answer" -> "$TimedOut", "Limit" -> 1.|>,
    "G14" -> <|"BaselineAnswer" -> "unfinished", "Seconds" -> Infinity, "Answer" -> "$TimedOut", "Limit" -> 1.|>,
    "G15" -> <|"BaselineAnswer" -> "unfinished", "Seconds" -> Infinity, "Answer" -> "$TimedOut", "Limit" -> 10.|>
|>
