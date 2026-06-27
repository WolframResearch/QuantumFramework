# Handoff: fix the Python / Colab edition of "Qubit Explained" (QIS course)

## Context
Single-qubit quantum-information tutorial ("Qubit Explained"), authored as literate markdown and built to several outputs. Working directory (all paths absolute):

`/Users/mohammadb/Documents/GitHub/QuantumFramework/OngoingProjects/Courses/IntroductionToQuantumInformationAndScience/`

Files and roles:
- `introduction-to-QIS.md` — original English source.
- `introduction-to-QIS-revised.md` — **revised English = source of truth. CURRENT** (all fixes applied).
- `introduction-to-QIS-persian.md` — Persian translation. CURRENT.
- `introduction-to-QIS-persian.nb` / `introduction-to-QIS-persian-rtl.nb` — Persian Mathematica notebooks, built via MarkdownToNotebook ("md2nb") at `/Users/mohammadb/Documents/GitHub/MarkdownToNotebook`. `-rtl.nb` is the right-to-left-display version. CURRENT.
- `introduction-to-QIS-python.md` — **Python/wolframclient edition** (WL embedded in `wl()` / `wlshow()`), the source for the Colab notebook. **STALE.**
- `introduction-to-QIS-colab.ipynb` — Google Colab notebook generated from `python.md` (640 cells, markdown+code). **STALE.**
- `introduction-to-QIS-persian-preview.html` — RTL HTML preview (pandoc + MathJax).

## Status (verified by grep counts)
`revised.md`, `persian.md`, `persian.nb`/`-rtl.nb` are current. `python.md` and `colab.ipynb` predate the fixes (built Jun 26 11:57; fixes landed 15:29), so they still contain every issue below.

| issue | revised.md | python.md | colab.ipynb |
|---|---|---|---|
| broken `\begin{matrix}` (manual fixed-size parens) | 0 | 6 | 6 |
| verbose `ResourceObject[<\|"Name"...` | 0 | 27 | 27 |
| md2nb `<code>[Name]()[args]` refs in prose | 0 | 4 | 4 |
| literal PlotLabel `Subscript[\[Sigma], x]` | 0 | 6 | (in code) |

The screenshots of broken `Tr[ρ.σ_x σ_x]` and `det[□]` are the **Colab** file, not `revised.md` (which renders these correctly).

---

## TASK: bring `introduction-to-QIS-python.md` current, then regenerate the Colab

### 1. Text patches to `python.md` (Python regex; verify the replacement count each time)

(a) **Matrix → pmatrix** — fixes the Itô SDE, Cartesian SDE, and Gaussian 2x2 matrices that render with tiny non-stretching parentheses. Cause: they use `\begin{matrix}` (no delimiters) wrapped in hand-typed `(\, ... \, )`, which do not grow to the matrix height. Fix uses auto-sized `pmatrix`. Regex with DOTALL:
  - find: `\(\\,\s*\\begin\{matrix\}(.*?)\\end\{matrix\}\\,\s*\)`
  - replace: `\begin{pmatrix}\1\end{pmatrix}`
  - expect **6**; afterward 0 `\begin{matrix}` remain.

(b) **Collapse verbose ResourceFunction** to the named short form (`ResourceFunction["BlochSpherePlot"]` etc.). Regex with DOTALL:
  - find: `ResourceFunction\[\s*ResourceObject\[<\|"Name" -> "(LindbladSolve|ZassenhausTerms|BlochSpherePlot)".*?\]\]`
  - replace: `ResourceFunction["\1"]`
  - expect **27**; afterward 0 `ResourceObject[<|"Name"` remain.

(c) **PlotLabel literal Subscript → box markup** so the sigmas typeset:
  - `Subscript[\[Sigma], x]` -> `\!\(\*SubscriptBox[\(\[Sigma]\), \(x\)]\)` (and `y`, `z`). Expect **6**.

(d) **Stratonovich split** (if present): replace `[\sigma_{z},\rho ]$∘d$W_{t}$` -> `[\sigma_{z},\rho ]\circ dW_{t}$`.

### 2. Translate md2nb prose conventions to Jupyter MathJax (the main Colab-specific bug)
The prose uses md2nb's inline code-reference convention `<code>[Name]()[args]</code>`, which Jupyter/Colab cannot interpret, so it shows as the garbled `Tr[ρ.σ_x σ_x]` (the `$...$` inside the code span double-renders). Convert each to plain `$...$` math. Exact map (identical to what was applied to the Persian edition):
- `<code>[Tr]()[ρ.$\sigma _{x}$]</code>` -> `$\mathrm{Tr}[\rho .\sigma_{x}]$`
- `<code>[Tr]()[ρ.$\sigma _{y}$]</code>` -> `$\mathrm{Tr}[\rho .\sigma_{y}]$`
- `<code>[Tr]()[ρ.$\sigma _{z}$]</code>` -> `$\mathrm{Tr}[\rho .\sigma_{z}]$`
- `<code>[Tr]()[$A^{\dagger }$B]</code>` -> `$\mathrm{Tr}[A^{\dagger }B]$`
- `<code>[Cross]()[v1,v2]</code>` -> `$\mathrm{Cross}[v1,v2]$`
- `<code>[TensorProduct]()[[LeviCivitaTensor]()[3],v1,v2]</code>` -> `$\mathrm{TensorProduct}[LeviCivitaTensor[3],v1,v2]$`
- `<code>[Erf]()[x]</code>` -> `$\mathrm{Erf}[x]$`
- `<code>[Cos]()[$\omega _{0}$t+ϕ(t)]</code>` -> `$\mathrm{Cos}[\omega_{0}t+\phi (t)]$`
- `<code>[Cov]()[…]</code>` -> `$\mathrm{Cov}[\ldots ]$`
- `<code>[ItoProcess]()[{a,b,$r_{t}$},{r,$r_{0}$},t]</code>` -> `$\mathrm{ItoProcess}[\{a,b,r_{t}\},\{r,r_{0}\},t]$`

Then grep for any remaining `<code>[`...`]()` and any `\!\(\*...\)` box markup sitting in **prose** (markdown) cells and rewrite as `$...$` MathJax. (`\!\(\*...\)` inside actual WL **code** cells is fine; only fix prose.)

`det[□]`: the Gram-matrix paragraph uses `$det[\mathcal{G}]$`; MathJax renders `\mathcal{G}` fine, so confirm the `$...$` delimiters are intact in the markdown cell. If the box still fails in Colab, replace `\mathcal{G}` with the literal 𝒢 in that prose math, or check the cell wasn't split mid-`$...$`.

### 3. Regenerate `introduction-to-QIS-colab.ipynb` from the fixed `python.md`
- First determine the build step: there is **no converter script in this folder** — confirm how the `.ipynb` was produced from `python.md` (a `md -> ipynb` transform, or a script elsewhere; check git history / ask the user).
- If regeneration executes the `wolframclient` cells, it needs a Wolfram Engine and, for the cloud ResourceFunctions (`BlochSpherePlot`, `LindbladSolve`, `ZassenhausTerms`), a Wolfram Cloud connection.
- After regen, re-verify in the `.ipynb`: 0 `\begin{matrix}`, 0 verbose `ResourceObject`, 0 `<code>[`...`]()` in markdown cells.

---

## Optional / deferred (also applies to the `.nb` path)
**PlotLabel code-prettiness**: two `ListLinePlot` cells build their label as a string containing raw `\!\(\*SubscriptBox...\)`. The label renders correctly but the Input-cell code looks ugly. Optional: rewrite as WL expressions, e.g. `Row[{OverVector["r"], "=", "{", AngleBracket[Subscript["\[Sigma]","x"]], ...}]`. Render-test before committing — a wrong expression breaks the plot.

## Reference: Persian RTL notebook (already done, for context)
Mathematica 15 has no bidirectional-text support, so `introduction-to-QIS-persian-rtl.nb` is produced by: pull each Persian text run, swap inline math/links for Private-Use-Area placeholders, run `arabic_reshaper` (joined glyph forms) + `python-bidi` `get_display` (visual order) per width-wrapped line, reinsert the math, and pin each cell's `PageWidth` to its measured widest line (measured at `PageWidth->Infinity` with left alignment; capped below the window). The `.md` stays the editable source; `-rtl.nb` is display-only (not searchable/editable in the notebook). Converter: `/tmp/convert_rtl2.wls` + `/tmp/bidi_helper.py` (may need recreating in a new session).
