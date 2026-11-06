#!/bin/bash
# Makes the anonymised code snapshot for double-blind review from the committed state: this
# repository without the papers and git history, surfacePosition without its large tests,
# and copies of streamFunctionWave and cfdTools, with names and repository addresses removed.
# Usage: papers/makeReviewSnapshot.sh <output directory> [<directory with streamFunctionWave and cfdTools>]
set -eu
root=$(git -C "$(dirname "$0")" rev-parse --show-toplevel)
out=$(readlink -f "$1")
deps=$(readlink -f "${2:-$root/..}")
name=consistentInterfacialFlowNumerics_code
work=$out/$name
rm -rf "$work"; mkdir -p "$work/thirdParty"

git -C "$root" archive HEAD | tar -x -C "$work" --exclude=paper --exclude=paper2 --exclude=.gitmodules --exclude=.github \
    --exclude=papers/makeReviewSnapshot.sh
git -C "$root/applications/utilities/surfacePosition" archive HEAD | \
    tar -x -C "$work/applications/utilities/surfacePosition" --exclude=test
for d in streamFunctionWave cfdTools; do
    mkdir "$work/thirdParty/$d"
    git -C "$deps/$d" archive HEAD | tar -x -C "$work/thirdParty/$d"
done

cat > "$work/README_REVIEW.md" <<'EOD'
# Code snapshot for review

Solver, test cases and scripts of the two papers. Names and repository addresses are
removed for double-blind review (ANONYMISED). `README.md` describes the contents.
`thirdParty/` holds copies of the two tools it refers to: streamFunctionWave, needed
for the wave cases, and cfdTools, needed for the parameter scans.

```shell
source <OpenFOAM-v2606>/etc/bashrc
./Allwmake -j
thirdParty/streamFunctionWave/Allwmake
export CFDTOOLS=$PWD/thirdParty/cfdTools/cfdTools
run/quickTest/Allrun
```
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
         (r"\bDLR\b", "ANONYMISED"), (r"\bJohan\b", "ANONYMISED"), (r"Roskilde", "ANONYMISED")]
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
grep -r -l -a -i -E 'roenby|stromning|scheufler|missios|aliyar|interflowers|github.com.isoadvector' "$work" || true
tar -C "$out" -czf "$out/$name.tar.gz" "$name"
du -sh "$out/$name.tar.gz"
