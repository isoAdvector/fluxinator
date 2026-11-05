#!/usr/bin/env python3
"""Apply the OpenFOAM Journal writing style to LaTeX sections of a paper.

    python3 papers/journalStyle.py paper2/paper2/*.tex

- References: Fig./Figs., Tab./Tabs., Eqn./Eqns. with \\ref inside a sentence, the
  full word at the start of a sentence (WritingStyle.pdf of the journal).
- Table captions above the table.
- OpenFOAM names in \\texttt (the list WORDS).

Running it again changes nothing, so it can be rerun after merging new text.
"""
import re
import sys

WORDS = ("interFoam interIsoFoam interIsoGoam ddtCorr interfaceUpwind limitedLinearV "
         "limitedLinear laxWendroff plicRDF isoAlpha setAlphaField nLayers surfCellTol "
         "constantHeight").split()
# A bare equation number directly after one of these needs no parentheses.
LINKING = ("of in by from with to and or is if so that exactly see on into using through than "
           "as between for reproduces satisfies").split()
# Their arguments are not running text.
SKIPCMD = ("cite ref eqref label texttt datasource datanote datafn figdir nolinkurl url "
           "input includegraphics begin end bibliography bibliographystyle").split()
MATHENV = "equation align gather multline tabular".split()


def sentence_start(before):
    """True if the text after `before` starts a sentence or a paragraph."""
    stripped = before.rstrip()
    if not stripped or "\n\n" in before[len(stripped):]:
        return True
    if re.search(r"\\(caption|footnote|item)\s*(\[[^\]]*\])?\{?$", stripped):
        return True
    if re.search(r"\\parhead\{[^{}]*\}$", stripped):
        return True
    if re.search(r"(e\.g|i\.e|cf|et al|Fig|Figs|Tab|Tabs|Eqn|Eqns|Sec|Sect)\.$", stripped):
        return False
    if stripped[-1] in ".?!":
        return True
    # A footnote after a full stop
    return bool(re.search(r"[.?!]\\(footnote|datafn)\{[^\n]*\}$", stripped))


def _word(pattern, start, inside, text):
    out, pos = [], 0
    for m in re.finditer(pattern, text):
        out.append(text[pos:m.start()])
        done = "".join(out)
        out = [done]
        name = start if sentence_start(done) else inside
        out.append(name + "~\\ref{" + m.group(1) + "}")
        pos = m.end()
    out.append(text[pos:])
    return "".join(out)


def references(text):
    for words, start, inside in ((r"Figures|Figs\.", "Figures", "Figs."), (r"Figure|Fig\.", "Figure", "Fig."),
                                 (r"Tables|Tabs\.", "Tables", "Tabs."), (r"Table|Tab\.", "Table", "Tab."),
                                 (r"Equations|Eqs\.|Eqns\.", "Equations", "Eqns."),
                                 (r"Equation|Eq\.|Eqn\.", "Equation", "Eqn.")):
        text = _word(r"(?<![A-Za-z])(?:" + words + r")(?:~|\s+)\\(?:eq)?ref\{([^}]*)\}", start, inside, text)
    # Lists of bare equation numbers
    text = re.sub(r"(Eqns\.~\\ref\{[^}]*\}(?:,| and)\s*)\\eqref", r"\1\\ref", text)
    text = re.sub(r"(Equations~\\ref\{[^}]*\}(?:,| and)\s*)\\eqref", r"\1\\ref", text)
    text = re.sub(r"\\eqref\{([^}]*)\}((?:,| and|--)\s*)\\eqref\{([^}]*)\}",
                  r"Eqns.~\\ref{\1}\2\\ref{\3}", text)
    # Bare equation numbers
    out, pos = [], 0
    for m in re.finditer(r"\\eqref\{([^}]*)\}", text):
        out.append(text[pos:m.start()])
        done = "".join(out)
        out = [done]
        words = re.findall(r"[A-Za-z]+", done[-30:])
        ref = "~\\ref{" + m.group(1) + "}"
        if sentence_start(done):
            out.append("Equation" + ref)
        elif done.rstrip().endswith("(") or (words and words[-1].lower() in LINKING and done.rstrip()[-1].isalpha()):
            out.append("Eqn." + ref)
        else:
            out.append("(Eqn." + ref + ")")
        pos = m.end()
    out.append(text[pos:])
    return "".join(out)


def _brace_end(text, i):
    """Index just after the group that opens at text[i] == '{'."""
    depth = 0
    while i < len(text):
        if text[i] == "\\":
            i += 2
            continue
        depth += (text[i] == "{") - (text[i] == "}")
        i += 1
        if depth == 0:
            return i
    raise ValueError("unbalanced braces")


def table_captions(text):
    def move(m):
        body = m.group(0)
        tab, cap = body.find("\\begin{tabular"), body.find("\\caption")
        if tab < 0 or cap < 0 or cap < tab:
            return body
        end = _brace_end(body, body.index("{", cap))
        label = re.match(r"\s*\\label\{[^}]*\}", body[end:])
        if label:
            end += label.end()
        block = body[cap:end]
        body = body[:cap].rstrip() + body[end:]
        at = re.search(r"\\begin\{table\}(\[[^\]]*\])?\s*(\\centering\s*)?", body).end()
        return body[:at] + block + "\n" + body[at:]
    return re.sub(r"\\begin\{table\}.*?\\end\{table\}", move, text, flags=re.S)


def typewriter(text):
    skip = re.compile(r"(?<!\\)%[^\n]*|(?<!\\)\$\$?.*?(?<!\\)\$\$?|\\\[.*?\\\]"
                      r"|\\begin\{(" + "|".join(MATHENV) + r")\*?\}.*?\\end\{\1\*?\}"
                      r"|\\(?:" + "|".join(SKIPCMD) + r")\b\*?(?:\[[^\]]*\])?(?=\{)", re.S)
    word = re.compile(r"(?<![A-Za-z\\_])(" + "|".join(sorted(WORDS, key=len, reverse=True)) + r")(?![A-Za-z_])")
    out, pos = [], 0
    while True:
        m = skip.search(text, pos)
        stop = m.start() if m else len(text)
        out.append(word.sub(r"\\texttt{\1}", text[pos:stop]))
        if not m:
            return "".join(out)
        end = m.end()
        if text[end:end + 1] == "{" and m.group(0).startswith("\\") and not m.group(0).startswith(("\\[", "\\begin{")):
            end = _brace_end(text, end)
        out.append(text[m.start():end])
        pos = end


for path in sys.argv[1:]:
    with open(path) as f:
        old = f.read()
    new = typewriter(table_captions(references(old)))
    if new != old:
        with open(path, "w") as f:
            f.write(new)
        print("changed", path)
