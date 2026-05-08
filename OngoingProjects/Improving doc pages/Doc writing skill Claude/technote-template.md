# Template: Tech Note (Tutorial) Page Structure

This is the canonical structure for a Wolfram paclet Tech Note page.
A tech note is a narrative tutorial that explains concepts, workflows, or algorithms with prose, math, and code examples. It teaches — unlike a reference page (which catalogs) or a guide page (which indexes).

## Formatting rules

- **Voice**: expository but concise. Explain concepts clearly, avoid filler. Wolfram documentation style — technical prose, not conversational.
- **Code blocks**: `` ```wl `` fenced blocks with In/Out pairs.
- **Math**: use Wolfram Language notation in code blocks for complex formulas. Use prose with subscript notation for inline math.
- **Tables**: used for definition tables (key terms, notation). All alignment headers must be explicit.
- **`---`** (horizontal rule) not typically used in tech notes — sections are separated by `##` headings.
- **References**: numbered citations `[n]` with a References section at the end. Use italic for book/journal titles.

### Paclet link format

Links in documentation must include the full paclet qualifier for non-system symbols. The format depends on whether you're linking to system documentation or paclet documentation:

- **System symbols**: `paclet:ref/Symbol`, `paclet:guide/GuideName`, `paclet:tutorial/TutorialName`
- **Paclet symbols**: `paclet:PublisherID/PacletName/ref/Symbol`
- **Paclet guides**: `paclet:PublisherID/PacletName/guide/GuideName`
- **Paclet tutorials**: `paclet:PublisherID/PacletName/tutorial/TechNoteName`

The skill resolves `PublisherID` and `PacletName` during the gather-context step. Use `paclet:ref/...` (short form) only for built-in Wolfram Language symbols.

## Content patterns

A tech note follows a narrative arc:

1. **Introduction** — what the reader will learn, why it matters
2. **Definition / setup** — key terms, notation, mathematical framework
3. **Core sections** — progressive explanation with interleaved prose, math, and code
4. **Advanced topics** — deeper material for readers who want more
5. **References** — cited works

Each section should be self-contained enough to be useful on its own, but sections should build on each other. Code examples demonstrate concepts introduced in the prose — they are not the primary vehicle for explanation (unlike reference pages where examples are central).

---

## Page structure

What follows is the page structure. Inline comments (in parentheses) are guidance — do not include them in the output.

# Tech Note Title

(Introductory paragraph. What this tech note covers, what the reader will learn. Set context and motivation. Can be 1–3 paragraphs.)

(Optional: definition table for key terms / notation used throughout.)

|              |                                          |
| ------------ | ---------------------------------------- |
| Term1        | definition or description                |
| Term2        | definition or description                |
| Term3        | definition or description                |

*Caption for the definition table.*

## First Concept Section

(Prose explaining the concept. Can include multiple paragraphs, bullet lists, and inline math notation.)

This shows how to use the concept:

```wl
In[1]:= CodeExample[x, y]

Out[1]= result
```

(More prose building on the example result.)

## Second Concept Section

(Continue building the narrative. Each section introduces a new concept or builds on the previous one.)

### Subsection

(Use `###` for subsections within a concept when needed.)

```wl
In[2]:= AnotherExample[args]

Out[2]= result
```

(Interpret the result, connect it to the concept.)

## Advanced Topics

(Optional. Deeper material — algorithmic details, mathematical proofs, performance characteristics, edge cases.)

### Algorithm Details

(Detailed explanation with math notation in code blocks for complex formulas.)

```wl
formula or algorithmic expression
```

(Explain the formula components.)

## References

(Numbered list of cited works. Use standard bibliographic format.)

[1] Author, A. B. *Title of Book*. Publisher, Year.

[2] Author, C. D. "Title of Paper." *Journal Name* Volume, no. Issue (Year): Pages.

## Related Guides

(Guide pages that index the symbols discussed in this tech note. Use paclet-qualified links for paclet guides, short form for system guides.)

* [PacletGuide](paclet:PublisherID/PacletName/guide/PacletGuide)
* [SystemGuide](paclet:guide/SystemGuide)

## Related Tech Notes

(Other tech notes that complement this one. Use paclet-qualified links for paclet tech notes.)

* [OtherTechNote](paclet:PublisherID/PacletName/tutorial/OtherTechNote)
