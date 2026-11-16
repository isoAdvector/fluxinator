#!/bin/bash
# Copyright (C) 2026 Johan Roenby
# SPDX-License-Identifier: GPL-3.0-or-later
# Makes the anonymised code snapshot for double-blind review from the committed state: this
# repository without its git history, and in thirdParty/ copies of the four separately installed
# tools (streamFunctionWave, surfacePosition without its tests, dictator and cfdTools), with
# names and repository addresses removed.
# Usage: papers/makeReviewSnapshot.sh <output directory> [<directory with the four tools>]
set -eu
root=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
out=$(readlink -f "$1")
deps=$(readlink -f "${2:-$root/..}")
name=consistentInterfacialFlowNumerics_code
work=$out/$name
rm -rf "$work"; mkdir -p "$work/thirdParty"

git -C "$root" archive HEAD | tar -x -C "$work" --exclude=paper --exclude=paper2 --exclude=.gitmodules --exclude=.github \
    --exclude=papers/makeReviewSnapshot.sh
for d in streamFunctionWave surfacePosition dictator cfdTools; do
    mkdir "$work/thirdParty/$d"
    git -C "$deps/$d" archive HEAD | tar -x -C "$work/thirdParty/$d" --exclude=test --exclude=.github
done

cat > "$work/README_REVIEW.md" <<'EOD'
# Code snapshot for review

Solver, test cases and scripts of the paper series. Names and repository addresses are
removed for double-blind review (ANONYMISED). `README.md` describes the contents.
`thirdParty/` holds copies of the four tools it refers to, which are installed separately:
streamFunctionWave (initial and exact wave of the wave cases), surfacePosition (wave
gauge), dictator (sets up the cases) and cfdTools (parameter scans).

```shell
source <OpenFOAM-v2606>/etc/bashrc     # or v2306 and later
source thirdParty/dictator/dictator
./Allwmake -j
thirdParty/streamFunctionWave/Allwmake -j
thirdParty/surfacePosition/Allwmake
export CFDTOOLS=$PWD/thirdParty/cfdTools/cfdTools
run/quickTest/Allrun                   # about a minute
```

Each case in `run/` has a README and an `Allrun` (`run/README.md`); each figure and
table of the papers has a directory in `papers/` with the script that makes it.
EOD

python3 - "$work" <<'EOD'
import os
import re
import sys

WORDS = [(r"github\.com[:/]isoAdvector", "github.com/ANONYMISED"), (r"github\.com[:/]InterFlowers", "github.com/ANONYMISED"),
         (r"InterFlowers", "ANONYMISED"), (r"johan@stromning\.com", "ANONYMISED"), (r"Johan\s+Roenby", "ANONYMISED"),
         (r"Henning\s+Scheufler", "ANONYMISED"), (r"Konstantinos\s+Missios", "ANONYMISED"), (r"Sithik\s+Aliyar", "ANONYMISED"),
         (r"Kasper\s+M(ø|oe?)ller", "ANONYMISED"), (r"/home/roenby", "/home/anonymised"), (r"(?i)roenby", "ANONYMISED"),
         (r"(?i)stromning", "ANONYMISED"), (r"Scheufler", "ANONYMISED"), (r"Missios", "ANONYMISED"), (r"Aliyar", "ANONYMISED"),
         (r"\bDLR\b", "ANONYMISED"), (r"\bJohan\b", "ANONYMISED"), (r"Roskilde", "ANONYMISED"),
         (r"M(ø|oe)ller", "ANONYMISED"), (r"\bDHI\b", "ANONYMISED")]
for top, _, files in os.walk(sys.argv[1]):
    for name in files:
        path = os.path.join(top, name)
        with open(path, "rb") as f:
            data = f.read()
        if b"\0" in data:
            continue
        text = new = data.decode("utf-8", "surrogateescape")
        for pattern, replacement in WORDS:
            new = re.sub(pattern, replacement, new)
        if new != text:
            with open(path, "wb") as f:
                f.write(new.encode("utf-8", "surrogateescape"))
EOD

echo "Identifying words left (none expected):"
grep -r -l -a -i -E 'roenby|stromning|scheufler|missios|aliyar|m(ø|oe)ller|\bdhi\b|interflowers|github.com.isoadvector|roskilde' "$work" || true
tar -C "$out" -czf "$out/$name.tar.gz" "$name"
du -sh "$out/$name.tar.gz"
