# Template: WFR Documentation Page Structure

This is the canonical structure for a Wolfram Function Repository documentation page.
Every section must be present in the output. Leave the body empty if there is nothing to say.

## Formatting rules

- **Voice**: terse, Wolfram-documentation style. One-liner descriptions. No filler.
- **Code blocks**: `` ```wl `` fenced blocks. Show In/Out pairs where useful.
- **Tables**: all alignment headers must be explicit: `|:---|` or `|:---:|` or `|---:|`.
- **`---`** (horizontal rule) only use at the beginning of a new example.

---

## Page structure

What follows is the page structure. Inline comments (in parentheses) are guidance — do not include them in the output.

# FunctionName

One-line description explaining the function's basic purpose.

## Definition

(The function definition. All definitions, including dependencies, will be included in the generated resource function.)

```wl
FunctionName[args_] := implementation
```

## Documentation

### Usage

(Document input usage cases. Each case is a calling pattern followed by a brief explanation. Every defined case should be demonstrated in the Examples section.)

``FunctionName[arg]``
	explanation of what this calling pattern does.

``FunctionName[arg1, arg2]``
	explanation of the second calling pattern.

### Details & Options

(Detailed explanation of how the function is used and configured: acceptable input types, result formats, option specifications, background information. May include multiple cells, bullet lists, tables, hyperlinks. Include all relevant background or contextual information.)

* Detail about usage, input types, or behavior.
* Information about options and their defaults.

## Examples

(Demonstrate the function's usage, starting with the most basic case. Each example has a preceding text caption. Use page breaks between examples within a group.)

### Basic Examples

(Most basic function usage — 2–3 examples giving immediate understanding.)

Caption ending with colon:

```wl
In[1]:= FunctionName[x, y]

Out[1]= x y
```

### Scope

(Input and display conventions, standard computational attributes like threading over lists, supported data types.)

Caption:

```wl
In[2]:= FunctionName[x, y, z]

Out[2]= x y z
```

### Options

(Available options and parameters. One subsection per option if there are many.)

### Applications

(Standard industry or academic applications. Each tells a mini-story where appropriate:)
1. State the problem
2. Show the data
3. Apply the function
4. Visualize or interpret the result

Caption:

```wl
Code
```

### Properties and Relations

(How the function relates to other functions. Each entry is one fact demonstrated with code.)

### Possible Issues

(Limitations or unexpected behavior a user might experience. Be generous — include all reasonable gotchas.)

### Neat Examples

(Particularly interesting, unconventional, or otherwise unique usage. Only when genuinely impressive. Leave empty otherwise.)

## Source & Additional Information

### Contributed By

(Name of the person, people, or organization that should be publicly credited.)

Author Name

### Keywords

(Comprehensive list of relevant terms: functional areas, algorithm names, related concepts. Used for search indexing.)

* keyword 1
* keyword 2

### Categories

(Choose categories that best represent the function. Options:)

[ ] Cloud & Deployment
[ ] Data Manipulation & Analysis
[ ] External Interfaces & Connections
[ ] Geographic Data & Computation
[ ] Graphs & Networks
[ ] Images
[ ] Knowledge Representation & Natural Language
[ ] Notebook Documents & Presentation
[ ] Mathematical Computation
[ ] Just For Fun
[ ] Machine Learning
[ ] Programming Utilities
[ ] Scientific and Medical Data & Computation
[ ] Sound & Video
[ ] Symbolic & Numeric Computation
[ ] Time-Related Computation
[ ] Visualization & Graphics

### Related Symbols

(Up to twenty documented, system-level Wolfram Language symbols related to the function.)

```wl
SymbolName (documented Wolfram Language symbol)
```

### Related Resource Objects

(Names of published resource objects from any Wolfram repository that are related.)

* Resource Name

### Source/Reference Citation

(Bibliographic-style citation for the original source: published paper, algorithm, code repository.)

### Links

(Additional URLs or hyperlinks for external information.)

* Link to related material

### Tests

(Verification tests for confirming the function works properly. Use Input/Output pairs or VerificationTest expressions.)

```wl
In[1]:= FunctionName[x, y]

Out[1]= x y
```

### Compatibility

#### Wolfram Language Version

(e.g. 14.3+)

#### Operating System

[\[Checkmark]] Windows
[\[Checkmark]] Mac
[\[Checkmark]] Linux

#### Environments

[\[Checkmark]] Session
[\[Checkmark]] WebEvaluation
[\[Checkmark]] BatchJob
[\[Checkmark]] Script
[\[Checkmark]] WebAPI
[\[Checkmark]] Subkernel
[\[Checkmark]] Scheduled

#### Cloud Support

[\[Checkmark]] Supported in cloud

#### Required Features

[ ] Notebooks
[ ] Parallel Kernels
[ ] Cloud Access

## Author Notes

(Background, possible improvements, additional information, implementation details beyond the scope of the documentation. Appears near the bottom of the published resource.)

## Submission Notes

(Additional information for the reviewer. Not included in the published resource.)
