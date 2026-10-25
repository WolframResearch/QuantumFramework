BeginPackage["DFTHCC30`"];
FourAddress::usage = "FourAddress[rows] reconstructs four rows of signed ten-bit words using three compact differences. Invalid input returns Failure.";
SymbolicChecks::usage = "SymbolicChecks[] proves the Boolean, sign-gauge and phase-repair identities used by the C30 XOR loader.";
Begin["`Private`"];

signedWordQ[x_] := IntegerQ[x] && 0 <= x < 1024 && BitGet[x, 9] == BitGet[x, 8];
compact[x_] := BitAnd[x, 511];
expand[x_] := BitOr[x, 512 BitGet[x, 8]];
FourAddress[rows_?MatrixQ] /; Length[rows] == 4 && Length[First[rows]] > 0 && AllTrue[Flatten[rows], signedWordQ] :=
  With[{base = rows[[1]], d1 = compact[BitXor[rows[[1]], rows[[2]]]],
      d2 = compact[BitXor[rows[[1]], rows[[3]]]],
      d3 = compact[BitXor[rows[[1]], rows[[2]], rows[[3]], rows[[4]]]]},
    {base, BitXor[base, expand[d1]], BitXor[base, expand[d2]],
      BitXor[base, expand[BitXor[d1, d2, d3]]]}];
FourAddress[___] := Failure["InvalidSignedTable", <||>];

SymbolicChecks[] := Module[{f0, f1, f2, f3, a, b, m, f, x, y, z, t, theta,
    interpolation, cases, signed, vectors, xor},
  interpolation = f0 + a (f0 + f1) + b (f0 + f2) + a b (f0 + f1 + f2 + f3);
  cases = MapThread[PolynomialMod[(interpolation /. Thread[{a,b} -> #1]) - #2, 2] &,
    {{{0,0},{1,0},{0,1},{1,1}}, {f0,f1,f2,f3}}];
  signed = Mod[Range[-256,255],1024];
  vectors = {Cos[theta], Sin[theta]};
  xor = x + y - 2 x y;
  <|"four_address_arbitrary_coefficients" -> (cases === {0,0,0,0}),
    "hmr_phase_repaired" -> TrueQ[FullSimplify[(-1)^(2 m f) == 1, Element[{m,f},Integers]]],
    "anf_phase_linearity" -> (PolynomialMod[m (f0 + a f1 + b f2 + a b f3) -
      (m f0 + m a f1 + m b f2 + m a b f3), 2] === 0),
    "dirty_target_fanout" -> (PolynomialMod[(t + z) + (z + x y) - (t + x y), 2] === 0),
    "xor_involution" -> (PolynomialMod[(xor /. y -> xor) - y, {x^2-x,y^2-y,2}] === 0),
    "row_sign_projector" -> (Outer[Times,-{x,y,z},-{x,y,z}] === Outer[Times,{x,y,z},{x,y,z}]),
    "givens_unit_norm" -> TrueQ[FullSimplify[vectors.vectors == 1, Element[theta,Reals]]],
    "gauge_endpoint_sine" -> (Sin[-Pi/2] === -1 && Sin[Pi/2] === 1),
    "signed_all_512_words" -> (BitGet[signed,9] === BitGet[signed,8]),
    "compact_expand_identity" -> (expand[compact[signed]] === signed),
    "zero_table_limit" -> (FourAddress[ConstantArray[0,{4,3}]] === ConstantArray[0,{4,3}]),
    "constant_table_limit" -> (FourAddress[ConstantArray[{255,768,1023},4]] === ConstantArray[{255,768,1023},4]),
    "nontrivial_exact_table" -> With[{rows={{0,255,768},{1023,0,255},{255,768,1023},{768,1023,0}}}, FourAddress[rows] === rows],
    "reject_unsigned_gap" -> FailureQ[FourAddress[ConstantArray[256,{4,3}]]],
    "reject_empty_table" -> FailureQ[FourAddress[{}]],
    "reject_float_table" -> FailureQ[FourAddress[ConstantArray[0.,{4,3}]]],
    "reject_negative_word" -> FailureQ[FourAddress[ConstantArray[-1,{4,3}]]]|>
];
End[];
EndPackage[];
