#!/usr/bin/env python3
"""Extract every WL fence verbatim and build a loud, cell-attributed audit runner."""

from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path


def escape_wl_string(text: str) -> str:
    return text.rstrip().replace("\\", "\\\\").replace('"', '\\"')


VISUAL_HEADS = (
    "Plot[",
    "ListPlot[",
    "ListLinePlot[",
    "ListLogPlot[",
    "ListLogLogPlot[",
    "BarChart[",
    "Histogram[",
    "DensityPlot[",
    "ListDensityPlot[",
    "Animate[",
    "Graphics[",
    "Graphics3D[",
    "GraphicsRow[",
    "Legended[",
    "Show[",
    "blochPlot[",
    "wignerPlot[",
)


def visual_evaluation_cell(body: str) -> bool:
    """Identify rendering calls, but retain helper definitions whose RHS contains graphics."""
    if ":=" in body:
        return False
    return any(head in body for head in VISUAL_HEADS)


def main() -> int:
    if len(sys.argv) not in {3, 4}:
        print("usage: extract-and-run.py CATALOG OUTPUT [ASSERTIONS]", file=sys.stderr)
        return 2

    source_path, output_path = map(Path, sys.argv[1:3])
    assertions_path = Path(sys.argv[3]).resolve() if len(sys.argv) == 4 else None
    source = source_path.read_text(encoding="utf-8")
    headings = [(match.start(), match.group(2)) for match in re.finditer(r"^(#{2,4})\s+(.+)$", source, re.M)]
    cells = [(match.start(), match.group(1)) for match in re.finditer(r"^```wl\s*$\n(.*?)^```\s*$", source, re.M | re.S)]

    def heading_at(position: int) -> str:
        label = "preamble"
        for start, title in headings:
            if start >= position:
                break
            label = title
        return label.replace('"', "'")

    digest = hashlib.sha256(source.encode("utf-8")).hexdigest()
    lines = [
        "(* Generated audit runner: exact WL fence bytes, in document order. *)",
        "$HistoryLength = 0; issueCount = 0;",
        f'WriteString[$Output, "SOURCE SHA256: {digest}\\n"];',
        f'WriteString[$Output, "CELL COUNT: {len(cells)}\\n"];',
    ]
    for index, (position, body) in enumerate(cells, start=1):
        label = heading_at(position)
        escaped = escape_wl_string(body)
        if visual_evaluation_cell(body):
            lines.extend(
                [
                    f'cellLabel = "cell {index}: {label}";',
                    'WriteString[$Output, "PARSE GRAPHICS: ", cellLabel, "\\n"];',
                    "cellMessages = {};",
                    "Block[{$MessageList = {}},",
                    f'  cellValue = ToExpression["{escaped}", InputForm, HoldComplete];',
                    "  cellMessages = $MessageList];",
                    "If[cellMessages =!= {} || cellValue === $Failed,",
                    "  issueCount++; WriteString[$Output, \"GRAPHICS PARSE FAILURE: \", cellLabel, \" => \", ToString[cellMessages, InputForm], \"\\n\"]];",
                    "Clear[cellValue, cellMessages, cellLabel];",
                ]
            )
            continue
        lines.extend(
            [
                f'cellLabel = "cell {index}: {label}";',
                'WriteString[$Output, "START: ", cellLabel, "\\n"];',
                "cellMessages = {}; cellTiming = AbsoluteTiming[",
                "  Block[{$MessageList = {}},",
                f'    cellValue = CheckAbort[TimeConstrained[ToExpression["{escaped}", InputForm], 120, $TimedOut], $Aborted];',
                "    cellMessages = $MessageList]][[1]];",
                "If[cellMessages =!= {},",
                "  issueCount++; WriteString[$Output, \"MESSAGES: \", cellLabel, \" => \", ToString[cellMessages, InputForm], \"\\n\"]];",
                "If[MemberQ[{$Failed, $Aborted, $TimedOut}, cellValue],",
                "  issueCount++; WriteString[$Output, \"FAILURE: \", cellLabel, \" => \", ToString[cellValue, InputForm], \"\\n\"]];",
                "If[cellTiming > 30, WriteString[$Output, \"SLOW: \", cellLabel, \" => \", ToString[cellTiming, InputForm], \" seconds\\n\"]];",
                'WriteString[$Output, "DONE: ", cellLabel, "\\n"];',
                "Clear[cellValue, cellMessages, cellTiming, cellLabel];",
            ]
        )
    if assertions_path is not None:
        lines.append(f'Get["{assertions_path}"];')
    lines.extend(
        [
            'WriteString[$Output, "ALL CELLS COMPLETED\\n"];',
            'WriteString[$Output, "OPEN ISSUES: ", issueCount, "\\n"];',
            "If[issueCount > 0, Exit[1], Exit[0]];",
        ]
    )
    output_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"{len(cells)} cells -> {output_path}; sha256={digest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
