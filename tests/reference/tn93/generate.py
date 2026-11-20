#!/usr/bin/env python3
"""Writes tests/reference/tn93_distances.tsv: HIV-TRACE's tn93 on the FASTA files in this
directory, every pair, in every ambiguity mode. src/reference_tests.rs reads it.

    python3 tests/reference/tn93/generate.py            # builds tn93 in a temporary directory
    TN93=/path/to/tn93 python3 tests/reference/tn93/generate.py

tn93 (github.com/veg/tn93) is the reference for what these modes mean. It is built at commit
076eb0f (v1.1.0) with CMake, after changing `%g` to `%.17g` in src/TN93.cpp so it prints full
precision (this leaves the computation alone); building needs git, cmake and a C++ compiler. A
binary given in $TN93 must be built the same way. It runs as

    tn93 -q -t 1000 -l 1 -f csv -a MODE FILE
    tn93 -q -t 1000 -l 1 -f csv -a resolve -g FRACTION FILE

for every mode, and for -g (what HIV-TRACE runs) at 0.015 (HIV-TRACE's default), 0.05, and the
edge cases 0 (every pair averaged) and 1 (never).

The threshold 1000 is tn93's value for a saturated model (TN93_MAX_DIST), so those pairs are
printed too; the library gives Infinity for them. tn93 leaves out pairs that compare at most one
site (the smallest overlap -l accepts); they are NA here.

The files:
- distances.fasta, ambiguity.fasta, hivtrace-g.fasta: MicrobeTrace's hand-checkable golden data
  (MicrobeTrace test/golden/README.md describes every sequence).
- n_codes.fasta: N and ? facing bases, partial codes, N and gaps.
- random.fasta: written by this script from a fixed seed: clusters of related sequences with
  ragged ends, gap runs (some longer than 17, which take tn93's gap-jump path), N, ?, partial
  ambiguity codes at varying rates, and lowercase.

tn93's reader drops characters outside ACGTURYSWKMBDHVN?- (upper-cased), which shifts the rest of
the sequence; an aligner's '.' would not survive. The files use only those characters, and no U:
tn93 reads U as a separate code (skipped in skip mode), while the library reads it as T.
"""

import os
import random
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "tn93_distances.tsv")
COMMIT = "076eb0f"
FILES = ["distances.fasta", "ambiguity.fasta", "hivtrace-g.fasta", "n_codes.fasta", "random.fasta"]
# (mode, -g fraction or None)
CONFIGS = [("skip", None), ("resolve", None), ("average", None), ("gapmm", None)] + [
    ("resolve", f) for f in (0.015, 0.05, 0.0, 1.0)
]


def write_random(path, seed=20261120):
    rnd = random.Random(seed)
    length = 237  # crosses 64-site words, not a multiple of 64
    root = [rnd.choice("ACGT") for _ in range(length)]

    def mutate(seq, rate):
        return [rnd.choice("ACGT".replace(c, "")) if rnd.random() < rate else c for c in seq]

    out = []
    for c in range(3):
        center = mutate(root, 0.05)
        for k in range(8):
            s = mutate(center, 0.012)
            n_rate = rnd.choice([0, 0, 0.005, 0.02, 0.05])
            amb_rate = rnd.choice([0, 0.005, 0.01, 0.02, 0.04])
            for i in range(length):
                r = rnd.random()
                if r < n_rate:
                    s[i] = rnd.choice("NNN?")
                elif r < n_rate + amb_rate:
                    s[i] = rnd.choice("RYSWKMBDHV")
            for _ in range(rnd.choice([0, 1, 1, 2])):  # inner gap runs, some longer than 17
                start, run = rnd.randrange(length), rnd.choice([1, 3, 9, 20, 30])
                s[start:start + run] = "-" * len(s[start:start + run])
            lead, trail = rnd.choice([0, 0, 2, 15]), rnd.choice([0, 0, 1, 40])
            s[:lead] = "-" * lead
            s[length - trail:] = "-" * trail
            if rnd.random() < 0.25:
                s = [ch.lower() if rnd.random() < 0.3 else ch for ch in s]
            out.append((f"c{c}s{k}", "".join(s)))
    # Unrelated enough to saturate some models, and one mostly gaps
    out.append(("far", "".join(mutate(root, 0.6))))
    out.append(("gappy", "-" * 100 + "".join(root[100:130]) + "-" * (length - 130)))
    with open(path, "w") as f:
        for name, s in out:
            assert len(s) == length
            f.write(f">{name}\n{s}\n")


def build_tn93(workdir):
    src = os.path.join(workdir, "tn93")
    run = lambda *cmd, cwd=None: subprocess.run(cmd, cwd=cwd, check=True, stdout=subprocess.DEVNULL)
    run("git", "clone", "-q", "https://github.com/veg/tn93.git", src)
    run("git", "checkout", "-q", COMMIT, cwd=src)
    run("sed", "-i", "s/%g/%.17g/g", "src/TN93.cpp", cwd=src)
    build = os.path.join(src, "build")
    os.mkdir(build)
    run("cmake", "-Wno-dev", "-DCMAKE_BUILD_TYPE=Release", "..", cwd=build)
    run("make", "-j4", "tn93", cwd=build)
    return os.path.join(build, "tn93")


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:]
                seqs[name] = ""
            elif line:
                seqs[name] += line
    return seqs


def run_tn93(tn93, path, mode, fraction):
    cmd = [tn93, "-q", "-t", "1000", "-l", "1", "-f", "csv", "-a", mode]
    if fraction is not None:
        cmd += ["-g", repr(fraction)]
    out = subprocess.run(cmd + [path], check=True, capture_output=True, text=True).stdout
    lines = out.split("\n")
    assert lines[0] == "ID1,ID2,Distance", lines[0]
    found = {}
    for line in lines[1:]:
        if not line:
            break  # the summary follows a blank line
        a, b, d = line.split(",")
        found[(a, b)] = d
    return found


def main():
    write_random(os.path.join(HERE, "random.fasta"))
    with tempfile.TemporaryDirectory() as tmp:
        tn93 = os.environ.get("TN93") or build_tn93(tmp)
        rows = []
        for name in FILES:
            path = os.path.join(HERE, name)
            ids = list(read_fasta(path))
            for mode, fraction in CONFIGS:
                found = run_tn93(tn93, path, mode, fraction)
                for i, a in enumerate(ids):
                    for b in ids[i + 1:]:
                        d = found.pop((a, b), "NA")
                        rows.append((name, mode, "-" if fraction is None else repr(fraction), a, b, d))
                assert not found, f"unexpected pairs {found}"
    with open(OUT, "w") as f:
        f.write(f"# tn93 {COMMIT} (veg/tn93), written by tests/reference/tn93/generate.py\n")
        f.write("file\tmode\tfraction\tid1\tid2\ttn93\n")
        for row in rows:
            f.write("\t".join(row) + "\n")
    print(f"{len(rows)} rows -> {os.path.relpath(OUT)}", file=sys.stderr)


if __name__ == "__main__":
    main()
