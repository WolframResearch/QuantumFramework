#!/usr/bin/env python3
"""
outcell_lint.py - deterministic pre-check for Out-cell Associations in a
"learning by computing" essay (.md with ```wl fences).

It enforces the one rule that is mechanical, not aesthetic:

    The Out cell holds only the value. A string-keyed Association is allowed
    ONLY when every key is a bare canonical name/symbol (or a bare parameter
    value) so the keys read as a legend. It is a violation when the cell
    staples two facts together (dict-of-dicts, or one value is a transform of
    a sibling value), which turns the Out cell into a labeled report.

Design point: the check keys on STRUCTURE, not on a keyword blocklist. The
descriptive-token filter ("grid-substituted", "purity at eta = 1") is a trap:
it passes role-word keys like "errors"/"orders" that clear the low lexical bar
while still stapling sub-tables. So the HARD rules here are structural:

  A. STAPLE-TRANSFORM: one top-level value expression is a subterm of another
     top-level value expression (a quantity shown next to a transform of it).
  B. NESTED-DASHBOARD: a top-level value is itself an Association / list of
     Rules / a symbol assigned (or a helper returning) an Association.

Both make the Out cell a dashboard; each table belongs in its own cell under
its own bridge. A SOFT lexical advisory (D) flags keys that are not obviously
bare canonical tokens, for a human/verifier to confirm - never a hard gate.

Usage:  python3 outcell_lint.py FILE.md [FILE2.md ...]
Exit:   nonzero if any HARD violation is found.
No Wolfram kernel needed; pure static analysis.
"""

import re, sys

# ---------- bracket/string/comment-aware top-level splitter ----------

OPENERS = {'[': ']', '(': ')', '{': '}'}
CLOSERS = {']', ')', '}'}

def _iter_tokens(s):
    """Yield (index, char) skipping over strings and (* comments *)."""
    i, n = 0, len(s)
    while i < n:
        c = s[i]
        if c == '"':                      # skip string literal
            j = i + 1
            while j < n:
                if s[j] == '\\':
                    j += 2; continue
                if s[j] == '"':
                    break
                j += 1
            yield ('str', i, j)
            i = j + 1; continue
        if c == '(' and i + 1 < n and s[i+1] == '*':   # skip (* comment *)
            j = s.find('*)', i + 2)
            j = n if j == -1 else j + 2
            yield ('cmt', i, j - 1)
            i = j; continue
        yield ('chr', i, i)
        i += 1

def depth_profile(s):
    """Return list of (index, char, depth_before_char) for structural chars,
    counting [](){} and <|...|> nesting, ignoring strings/comments."""
    prof = []
    depth = 0
    toks = list(_iter_tokens(s))
    k = 0
    while k < len(toks):
        kind, a, b = toks[k]
        if kind != 'chr':
            k += 1; continue
        c = s[a]
        # detect two-char association delimiters
        if c == '<' and a + 1 < len(s) and s[a+1] == '|':
            prof.append((a, '<|', depth)); depth += 1
            k += 1  # consume the following '|' token too
            if k < len(toks): k += 1
            continue
        if c == '|' and a + 1 < len(s) and s[a+1] == '>':
            depth -= 1; prof.append((a, '|>', depth))
            k += 1
            if k < len(toks): k += 1
            continue
        if c in OPENERS:
            prof.append((a, c, depth)); depth += 1
        elif c in CLOSERS:
            depth -= 1; prof.append((a, c, depth))
        else:
            prof.append((a, c, depth))
        k += 1
    return prof

def split_top(s, seps):
    """Split s on any separator in `seps` (1- or 2-char) at bracket depth 0."""
    prof = depth_profile(s)
    # map index -> depth for quick lookup of structural chars
    parts, last = [], 0
    i = 0
    # Build a depth lookup by walking chars with the same logic
    # Simpler: reconstruct split points from prof for single-char seps and
    # handle '->' specially.
    idx_depth = {}
    for (a, c, d) in prof:
        idx_depth[a] = d
    while i < len(s):
        matched = None
        for sep in seps:
            if s.startswith(sep, i):
                # depth at this position: use idx_depth of first char if structural,
                # else compute by scanning prof for nearest structural char <= i
                d = _depth_at(prof, i)
                if d == 0:
                    matched = sep; break
        if matched:
            parts.append(s[last:i]); i += len(matched); last = i
        else:
            i += 1
    parts.append(s[last:])
    return [p.strip() for p in parts]

def _depth_at(prof, i):
    """Depth of the code at character index i (0 = top level)."""
    d = 0
    for (a, c, dep) in prof:
        if a < i:
            # depth AFTER processing char a
            if c in ('<|',) or c in OPENERS:
                d = dep + 1
            elif c in ('|>',) or c in CLOSERS:
                d = dep
            else:
                d = dep
        else:
            break
    return d

# ---------- WL cell / statement extraction ----------

def wl_cells(md):
    """Return list of (start_line, code) for each ```wl fenced block."""
    cells, lines = [], md.splitlines()
    i = 0
    while i < len(lines):
        if lines[i].strip().startswith('```wl'):
            start = i + 1
            body = []
            i += 1
            while i < len(lines) and not lines[i].strip().startswith('```'):
                body.append(lines[i]); i += 1
            cells.append((start + 1, '\n'.join(body)))
        i += 1
    return cells

def top_statements(code):
    """Split a cell into top-level statements on ';' at depth 0."""
    stmts = split_top(code, [';'])
    # split_top keeps a trailing empty piece if code ends with ';'
    return stmts

def displayed_out(code):
    """The expression the cell displays: the last top-level statement that is
    not terminated by ';' and is not a SetDelayed definition. Returns the
    expression text or None."""
    # Does the code end (ignoring trailing ws/comments) with ';'?
    stripped = re.sub(r'\(\*.*?\*\)', '', code, flags=re.S).rstrip()
    stmts = [s for s in top_statements(code)]
    # Remove trailing empty
    while stmts and stmts[-1].strip() == '':
        stmts.pop()
    if not stmts:
        return None
    if stripped.endswith(';'):
        return None                      # output suppressed
    last = stmts[-1].strip()
    # SetDelayed definition displays Null
    if re.match(r'^[\w$]+\s*(\[.*\])?\s*:=', last, flags=re.S):
        return None
    return last

# ---------- association analysis ----------

def assoc_pairs(expr):
    """If expr is an <|...|> literal, return list of (key_text, value_text).
    Else None."""
    e = expr.strip()
    if not (e.startswith('<|') and e.endswith('|>')):
        return None
    inner = e[2:-2]
    pairs = []
    for chunk in split_top(inner, [',']):
        if not chunk.strip():
            continue
        kv = split_top(chunk, ['->', ':>'])
        if len(kv) != 2:
            return ('MALFORMED', chunk)
        pairs.append((kv[0].strip(), kv[1].strip()))
    return pairs

def assocthread_firstarg(expr):
    """If expr is AssociationThread[a, b], return (a_text, b_text) else None."""
    e = expr.strip()
    m = re.match(r'^AssociationThread\s*\[', e)
    if not m:
        return None
    inner = e[e.index('[')+1:-1]
    args = split_top(inner, [','])
    if len(args) < 2:
        return None
    return (args[0].strip(), args[1].strip())

# symbols assigned association-typed RHS, or helpers returning associations
def assoc_typed_symbols(all_code):
    syms = set()
    # name = <|...  OR  name = AssociationThread[...
    for m in re.finditer(r'([\w$]+)\s*=\s*(<\||AssociationThread\[)', all_code):
        syms.add(m.group(1))
    # helper[...] := ... AssociationThread[...]  or ends returning <|...|>
    for m in re.finditer(r'([\w$]+)\s*\[[^\]]*\]\s*:=', all_code):
        name = m.group(1)
        # crude: does its definition body contain AssociationThread or <| as the
        # returned head? mark if AssociationThread/<| appears in the RHS.
        start = m.end()
        body = all_code[start:start+800]
        if 'AssociationThread[' in body or '<|' in body:
            syms.add(name)
    return syms

def is_numeric_list(text):
    t = text.strip()
    if not (t.startswith('{') and t.endswith('}')):
        return False
    items = split_top(t[1:-1], [','])
    for it in items:
        it = it.strip()
        if not re.match(r'^[-+]?(\d+\.?\d*|\.\d+)$', it):
            return False
    return True

def value_is_assoc(value_text, assoc_syms):
    v = value_text.strip()
    if '<|' in v or 'AssociationThread[' in v:
        return True
    head = re.match(r'^([\w$]+)', v)
    if head and head.group(1) in assoc_syms:
        return True
    # value is `helper[sym]` where helper returns assoc
    call = re.match(r'^([\w$]+)\s*\[', v)
    if call and call.group(1) in assoc_syms:
        return True
    return False

CANON_STAT = {'mean', 'variance', 'stddev', 'standarddeviation', 'median',
              'min', 'max', 'norm', 'trace', 'purity', 'correlation',
              'e+', 'e-', 'real', 'imaginary'}

def key_embeds_condition(key_text):
    """HARD-decidable subset of 'key describes the value's condition/method':
    a string key that embeds a relational operator or a condition clause binding
    a parameter (e.g. "purity at eta = 1", "error when q<0"). This is never a
    bare canonical name or a bare parameter value, so it is a certain violation.
    Genuinely-semantic 'is this multiword name canonical?' stays SOFT."""
    k = key_text.strip()
    if not (k.startswith('"') and k.endswith('"')):
        return False
    s = k[1:-1]
    # strip typeset boxes / named chars so a ket's internals don't trip it
    s = re.sub(r'\\!\\\(\\\*.*?\\\)', ' ', s)
    s = re.sub(r'\\\[[A-Z][A-Za-z]+\]', ' ', s)
    if re.search(r'(==|<=|>=|=|≤|≥)', s):
        return True
    if re.search(r'\b(at|when|given|where|for)\b\s+\S', s, flags=re.I) and ' ' in s.strip():
        # condition word followed by something, in a multi-token key
        return True
    return False

def key_is_bare_canonical(key_text):
    """Heuristic soft check. A bare canonical key: a quoted single short token
    (symbol name, ket, stat name) with no descriptive whitespace/verb, or a
    non-string literal (number, symbol)."""
    k = key_text.strip()
    if not (k.startswith('"') and k.endswith('"')):
        return True     # non-string key (number/symbol) treated as legend value
    s = k[1:-1]
    low = s.lower().strip()
    if low in CANON_STAT:
        return True
    # typeset box (\!\(\*...\)) or a \[NamedCharacter] is a bare symbol/ket
    if r'\!\(\*' in s or re.search(r'\\\[[A-Z][A-Za-z]+\]', s):
        return True
    # single short token, no spaces, letters/greek/sub only -> likely a symbol
    if ' ' not in s and len(s) <= 6:
        return True
    return False

# ---------- main lint ----------

def lint(path):
    md = open(path, encoding='utf-8').read()
    cells = wl_cells(md)
    all_code = '\n'.join(c for _, c in cells)
    assoc_syms = assoc_typed_symbols(all_code)
    hard, soft = [], []
    for (line0, code) in cells:
        out = displayed_out(code)
        if out is None:
            continue
        # locate approx line of the Out expression for reporting
        outline = line0 + code[:code.rfind(out)].count('\n') if out in code else line0

        pairs = assoc_pairs(out)
        if pairs and pairs[0] == 'MALFORMED':
            soft.append((outline, "could not parse Association keys", out[:70]))
            continue
        if pairs is not None:
            keys = [k for k, _ in pairs]
            vals = [v for _, v in pairs]
            # Rule A: staple-transform (one value is a subterm of another)
            for i, vi in enumerate(vals):
                for j, vj in enumerate(vals):
                    if i != j and _is_subterm(vi, vj):
                        hard.append((outline, "STAPLE-TRANSFORM",
                            f'value of {keys[j]} contains value of {keys[i]} '
                            f'({vi!r} ⊂ {vj!r}); a quantity shown beside a '
                            f'transform of itself. Split into two cells.'))
                        break
            # Rule B: nested dashboard (a value is itself an association/table)
            for k, v in pairs:
                if value_is_assoc(v, assoc_syms):
                    hard.append((outline, "NESTED-DASHBOARD",
                        f'value of key {k} is an Association/table ({v!r}); '
                        f'the Out cell is a dict of dicts. Show each table in '
                        f'its own cell under its own bridge.'))
            # Rule C: key embeds a condition/relation (hard-decidable)
            for k in keys:
                if key_embeds_condition(k):
                    hard.append((outline, "CONDITION-IN-KEY",
                        f'key {k} embeds a condition/relation; it describes the '
                        f"value's condition, not its name. Move the condition to "
                        f'the Text bridge; the key must be a bare name/symbol.'))
            # Rule D: soft key lexicon
            for k in keys:
                if key_embeds_condition(k):
                    continue
                if not key_is_bare_canonical(k):
                    soft.append((outline, "KEY-LEXICON",
                        f'key {k} is not obviously a bare canonical name/symbol; '
                        f'confirm it reads as a legend, not an annotation.'))
            continue

        at = assocthread_firstarg(out)
        if at is not None:
            keyarg, valarg = at
            if is_numeric_list(keyarg):
                continue                 # parametric legend: keys are param values
            # string keys in AssociationThread: soft-check each
            if keyarg.strip().startswith('{'):
                for k in split_top(keyarg.strip()[1:-1], [',']):
                    if key_embeds_condition(k):
                        hard.append((outline, "CONDITION-IN-KEY",
                            f'AssociationThread key {k.strip()} embeds a '
                            f'condition/relation; move it to the Text bridge.'))
                        continue
                    if not key_is_bare_canonical(k):
                        soft.append((outline, "KEY-LEXICON",
                            f'AssociationThread key {k.strip()} is not obviously '
                            f'a bare canonical name/symbol; confirm it is a legend.'))
            # nested-table values inside AssociationThread
            if 'AssociationThread[' in valarg or '<|' in valarg:
                # legitimate when it's a flat legend of scalars; flag only if the
                # value clearly builds an association per key
                pass
    return hard, soft

def _is_subterm(small, big):
    """True if expression `small` occurs as a bracketed subterm of `big` and
    they are not identical. Token-boundary aware to avoid substring accidents."""
    small, big = small.strip(), big.strip()
    if small == big or not small:
        return False
    # require small to appear delimited by non-word chars (or bracket) in big
    pat = r'(?<![\w$])' + re.escape(small) + r'(?![\w$])'
    return re.search(pat, big) is not None

def main():
    files = sys.argv[1:]
    if not files:
        print("usage: outcell_lint.py FILE.md ...", file=sys.stderr)
        sys.exit(2)
    total_hard = 0
    for f in files:
        hard, soft = lint(f)
        print(f"\n=== {f} ===")
        if not hard and not soft:
            print("  clean: no Out-cell Association violations.")
        for (ln, kind, msg) in hard:
            print(f"  [HARD] line ~{ln}: {kind}\n         {msg}")
        for (ln, kind, msg) in soft:
            print(f"  [warn] line ~{ln}: {kind}\n         {msg}")
        print(f"  -- {len(hard)} hard, {len(soft)} soft")
        total_hard += len(hard)
    print(f"\nTOTAL HARD VIOLATIONS: {total_hard}")
    sys.exit(1 if total_hard else 0)

if __name__ == '__main__':
    main()
