BeginTestSection["IBMJob"]

(* Decoding of IBM Quantum SamplerV2 /results payloads into measurement counts.

   IBM returns each job's results as a RuntimeEncoder-serialized PrimitiveResult: a nested
   __type__/__value__ envelope whose measured bits live in a BitArray, packed as
   base64(zlib(.npy uint8)). The decode (ibmSamples / ibmCounts) is pure Wolfram Language,
   so this suite needs no qiskit/Python session and runs anywhere.

   Fixtures b3 and b12 are REAL payloads captured from ibm_fez (a 3-qubit GHZ+Fourier job
   and a 12-qubit deterministic-X job); their expected histograms were cross-checked against
   the IBM console and an independent numpy decode. ba/bb are a synthetic two-register
   payload with known values, exercising the multi-register packing that QuantumCircuitOperator
   itself never emits (it always builds a single classical register) but external blobs may. *)

(* ---- fixture builders: assemble the exact /results envelope shape ---- *)
makeBA[b64_, nb_] := <|"__type__" -> "BitArray", "__value__" -> <|
  "array" -> <|"__type__" -> "ndarray", "__value__" -> b64|>, "num_bits" -> nb|>|>;
makeResults[names_, fields_] := <|"__type__" -> "PrimitiveResult", "__value__" -> <|
  "pub_results" -> {<|"__type__" -> "SamplerPubResult", "__value__" -> <|
    "data" -> <|"__type__" -> "DataBin", "__value__" -> <|
      "field_names" -> names, "fields" -> fields|>|>,
    "metadata" -> <|"circuit_metadata" -> <||>|>|>|>},
  "metadata" -> <||>, "version" -> 2|>|>;
makeJob[results_, mq_, status_ : "Completed"] := IBMJob[<|"ID" -> "fixture", "Status" -> status,
  "Primitive" -> "sampler", "MeasuredQubits" -> mq, "Raw" -> <|"Results" -> results|>|>];

(* ============================================================================
   Single classical register, 1 byte/shot (3 qubits). Real ibm_fez payload.
   ============================================================================ *)
b3 = "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";
job3 = makeJob[makeResults[{"c"}, <|"c" -> makeBA[b3, 3]|>], {0, 1, 2}];

VerificationTest[
  KeySort @ Counts[job3["Samples"]],
  <|0 -> 846, 1 -> 114, 2 -> 549, 3 -> 537, 4 -> 806, 5 -> 234, 6 -> 246, 7 -> 764|>,
  TestID -> "IBMJob-3q-raw-sample-histogram"
];

VerificationTest[ Total[job3["Counts"]], 4096, TestID -> "IBMJob-3q-counts-total" ];
(* full outcome map: each integer v -> bit-list Reverse[IntegerDigits[v, 2, 3]] (qubit 0 = LSB) *)
VerificationTest[
  KeySort @ job3["Counts"],
  KeySort @ <|{0, 0, 0} -> 846, {1, 0, 0} -> 114, {0, 1, 0} -> 549, {1, 1, 0} -> 537,
    {0, 0, 1} -> 806, {1, 0, 1} -> 234, {0, 1, 1} -> 246, {1, 1, 1} -> 764|>,
  TestID -> "IBMJob-3q-counts-bitlists"
];
VerificationTest[ {job3["NumBits"], job3["Shots"]}, {3, 4096}, TestID -> "IBMJob-3q-numbits-shots" ];
VerificationTest[ MatchQ[job3["Measurement"], _QuantumMeasurement], True, TestID -> "IBMJob-3q-measurement-head" ];

(* ============================================================================
   Single classical register, 2 bytes/shot (12 qubits). Real ibm_fez payload.
   X on qubits 1,9,12 -> the mode is the prepared state, and the high byte (qubits 9 and 12)
   is non-zero in essentially every shot, so this pins the multi-byte big-endian combine.
   ============================================================================ *)
b12 = "eJydVbFOwzAQvZMVKSbq0D9Ip4DUgVZMbCxsIBYGJlTRVAyIohRYAImMzPwwTmI79t25IHpqYl+e3909X5zvy+uLqxuEV3ir1vXurqlOy+r9ZVHNy2qzbZ6b1ePttlnXnf989bCrjX93v3qqzfxwcbw8mZfLo3n5Uf7vd6BRI6L2Fo6NKeab6RgB6NdZnOrvrTGPmqBOmUEj86YzQjLia7vYQsbZuI4yIuYSixa5uYH5T30uZ1JORayCjvQyljsenx/hGFTWgkpqQKOmMRPaSLzWhF0CokQUhfIPcwj8BcPSumQuaf+tKm26LqobAvFrBHlll4fti1bOz9uUslMWkV+xyro4mXsWdwIGuvEaeYQsoZmrSCVzsmu69wW4fiMbCEr/Yq7iscPAXlWo7B6mRG/w7ITRePqwXYrPjURfijqLyFmMHjXkV9HakCeIInUq0XLPqRpz9W8tPeNQPAlRFcTHVXbjgTET+oarB+J+ctWTnd96bSHW2fNoGrnHfIVnnLQbngm4njxTVplyDH/4MmUocwRWWMQk4rL2Gcy0UyXCgMw/fh1yh2+H+Q9G7kie";
job12 = makeJob[makeResults[{"c"}, <|"c" -> makeBA[b12, 12]|>], Range[0, 11]];

VerificationTest[
  KeySort @ Counts[job12["Samples"]],
  <|1 -> 3, 256 -> 4, 257 -> 135, 259 -> 1, 265 -> 3, 385 -> 1, 393 -> 1, 769 -> 1,
    1281 -> 2, 2049 -> 6, 2304 -> 24, 2305 -> 789, 2307 -> 10, 2309 -> 3, 2313 -> 8,
    2321 -> 2, 2337 -> 4, 2369 -> 1, 2432 -> 1, 2433 -> 13, 2817 -> 6, 3329 -> 6|>,
  TestID -> "IBMJob-12q-raw-sample-histogram"
];

VerificationTest[ {Total[job12["Counts"]], job12["NumBits"]}, {1024, 12}, TestID -> "IBMJob-12q-total-numbits" ];
VerificationTest[ Max[job12["Samples"]] >= 256, True, TestID -> "IBMJob-12q-high-byte-exercised" ];
(* end-to-end: the most frequent outcome is the prepared state (X on qubits 1, 9, 12) *)
VerificationTest[
  First @ Keys @ ReverseSort @ job12["Counts"],
  {1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1},
  TestID -> "IBMJob-12q-mode-is-prepared-state"
];

(* ============================================================================
   Two classical registers (10-bit + 5-bit), synthetic known values. field 0 ("a") occupies
   the low bits, field 1 ("b") the high bits: outcome integer = a + b*2^10.
   ============================================================================ *)
ba = "eJyb7BfqGxDJyFDGUK2eklqcXKRupaBeU2qorqOgnpZfVFKUmBefX5SSChJ3S8wpTgWKF2ckFqQC+RoWOgpGmjoKtQpkAy4GBgZG5v9MDAzsjAzMLxiYAUGkHbg=";
bb = "eJyb7BfqGxDJyFDGUK2eklqcXKRupaBeU2qorqOgnpZfVFKUmBefX5SSChJ3S8wpTgWKF2ckFqQC+RoWOgqGmjoKtQpkAy4GeUYBVg4meQBWqxwa";
jobS = makeJob[makeResults[{"a", "b"}, <|"a" -> makeBA[ba, 10], "b" -> makeBA[bb, 5]|>], Range[0, 14]];

VerificationTest[
  jobS["Samples"],
  {0, 31745, 2047, 16896, 5127, 8448, 3048, 31747},
  TestID -> "IBMJob-2register-packing"
];
VerificationTest[ jobS["NumBits"], 15, TestID -> "IBMJob-2register-total-width" ];

(* ============================================================================
   Status gating and graceful failure.
   ============================================================================ *)
VerificationTest[
  makeJob[makeResults[{"c"}, <|"c" -> makeBA[b3, 3]|>], {0, 1, 2}, "Running"]["Counts"],
  Missing["JobNotComplete", "Running"],
  TestID -> "IBMJob-not-complete-gate"
];
VerificationTest[
  makeJob[<||>, {0, 1, 2}]["Counts"],
  Missing["NoSamples"],
  TestID -> "IBMJob-malformed-results-nosamples"
];

EndTestSection[]
