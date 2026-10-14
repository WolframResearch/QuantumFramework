# QuantumFramework through the local Wolfram MCP server: 15-minute demo

Dry run 2026-10-06 on QF 2.1.3 (loaded from the QF folder), Wolfram/Arrays 1.4.1, Wolfram 15.0.1,
Wolfram/AgentTools 2.2.0. Every call below went through the installed server over stdio, the same
JSON-RPC traffic Claude Code sends. Times are per call.

## 0. Setup (before the audience arrives)

```wl
Needs["Wolfram`AgentTools`"];
InstallMCPServer["ClaudeCode", "WolframLanguage"]
```

This writes a `"Wolfram"` entry under `mcpServers` in `~/.claude.json` (user scope, every project).
Start a new Claude Code session afterwards: a session only sees MCP servers that existed when it started.
The tools then appear as `mcp__Wolfram__WolframLanguageEvaluator`, `...TestReport`, `...WriteNotebook`, and so on.
Undo with `UninstallMCPServer`; the pre-install config is at `~/.claude.json.bak-2026-10-06`.

Server start: about 2 s. The evaluator returns a session id ("Pass session=... in future calls");
definitions persist only when the model passes it back. All sessions share one kernel, so a loaded
package stays loaded across sessions while user definitions stay isolated.

## 1. Load QF (about 3.4 s): its own call, nothing else in it

```wl
PacletDirectoryLoad["/Users/mohammadb/Documents/GitHub/QuantumFramework/QuantumFramework"];
Needs["Wolfram`QuantumFramework`"];
{PacletObject["Wolfram/QuantumFramework"]["Version"], PacletObject["Wolfram/Arrays"]["Version"]}
```
Out: `{"2.1.3", "1.4.1"}`

QF is not installed as a paclet on this machine, so a bare `Needs` fails. Do not put any QF symbol in
this call: the evaluator parses the whole call first, so `QuantumCircuitOperator` is created in the
session context before `Needs` runs, shadows the real symbol (`::shdw`), and every later call in that
session returns unevaluated.

## 2. Bell state (0.4 s)

```wl
bell = QuantumCircuitOperator[{"H", "CNOT"}][];
{Values @ bell["Amplitudes"], QuantumEntangledQ[bell], QuantumEntanglementMonotone[bell, "Concurrence"]}
```
Out: `{{1/Sqrt[2], 0, 0, 1/Sqrt[2]}, True, 1}`

Avoid `bell["Formula"]` here: it comes back as `RawBoxes[...]`, which reads badly in a chat client.

## 3. Symbolic angle (each under 1 s)

```wl
psi = QuantumCircuitOperator[{"RY"[θ] -> 1, "CNOT" -> {1, 2}}][];
Simplify[psi["ProbabilitiesList"], θ ∈ Reals]
```
Out: `{Cos[θ/2]^2, 0, 0, Sin[θ/2]^2}`

```wl
Simplify[Normal @ QuantumPartialTrace[psi, {2}]["DensityMatrix"], θ ∈ Reals]
```
Out: `{{Cos[θ/2]^2, 0}, {0, Sin[θ/2]^2}}`

```wl
FullSimplify[QuantumEntanglementMonotone[psi, "EntanglementEntropy"], 0 < θ < Pi]
```
Out: `Quantity[(Log[2 Csc[θ]] + Cos[θ] Log[Tan[θ/2]])/Log[2], "Bits"]`

The reduced state is diagonal with the Born weights, so the entanglement entropy is the binary entropy
of cos²(θ/2): zero at θ = 0, one bit at θ = π/2. `"Concurrence"` gives `Sin[θ]` correctly but prints
three `PolynomialExpressionQ::ivar` messages in the MCP client, so it is left out of the live run.

## 4. Amplitude damping, γ and t symbolic (0.2 s)

```wl
rho = QuantumEvolve[QuantumOperator["Z"], {QuantumOperator[{{0, 1}, {0, 0}}]} -> {γ}, QuantumState["1"], t];
Simplify[Normal @ rho["DensityMatrix"], γ > 0]
```
Out: `{{1 - E^(-t γ), 0}, {0, E^(-t γ)}}`

The excited population decays as e^(-γt) and the Hamiltonian Z, diagonal in the same basis, leaves the
populations alone. The MCP client shows four `Part::pkspec1` messages from inside `DSolve` before the
answer; they appear with numeric γ too and do not change the result. `wolframscript` does not print them.

Warn the audience before running this step. The messages come out without their text, followed by a
general warning, and on screen they look like a failure:

```
Part::pkspec1: -- Message text not found -- (#1)
Part::pkspec1: -- Message text not found -- (Internal`ProcessEquations`ProcessEquationsDump`i)
Part::pkspec1: -- Message text not found -- (#1)
Part::pkspec1: -- Message text not found -- (Internal`ProcessEquations`ProcessEquationsDump`i)
General::messages: Messages were generated which may indicate errors.
```

The correct density matrix follows right after (live run 2026-10-06, 0.18 s in the kernel).

`AnnihilationOperator[2]` gives the same result only after
``Needs["Wolfram`QuantumFramework`SecondQuantization`"]``; step 1 does not load that subpackage. Without
it, `AnnihilationOperator` is an undefined symbol and `QuantumEvolve` recurses instead of failing
(QF-B86 in `QF-Defect-Ledger.md`), so the evaluator burns its whole time constraint and returns
`EvaluationTimeExceeded`. The plain σ₋ matrix avoids the extra load.

## 5. Isolation (0.0 s)

A call without the session id: `{MemberQ[$Packages, "Wolfram`QuantumFramework`"], ValueQ[bell]}`
gives `{True, False}`. The package is shared, the definitions are not.

## 6. Notebook round trip

`WriteNotebook` (2.8 s) turns markdown with a ```` ```wl ```` block into a `.nb`; `ReadNotebook` (0.2 s)
reads it back as the same markdown.

## 7. Tests (5.6 s)

`TestReport` on a `.wlt` that loads QF with the two lines of step 1 and checks
`FullSimplify[QuantumEntanglementMonotone[psi, "Concurrence"], 0 < θ < Pi] === Sin[θ]` and
`QuantumEntangledQ[QuantumState["00"]] === False`: 2/2 passed, returned as a markdown table.
