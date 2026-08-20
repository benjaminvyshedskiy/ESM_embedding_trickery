#!/usr/bin/env python3
"""Map the Jager benchmark baits onto vORF library fragment IDs with jackhmmer.

Reads QC/Well_Known_interactions/jager_2011_HIV.csv (written by
extract_jager_benchmark.py), searches each unique viral bait sequence against the
12,598-fragment viral ORFeome in data/mmc2.xlsx, and writes back a
vorf_fragment_ids column plus a companion per-hit TSV.

Search: jackhmmer -E 1e-5 -N 3 (plus output-only flags).

Which hits reach the CSV: target coverage >= 0.5 and an HIV-1/HIV-2 organism.
The target-coverage cut is the important one. The library stores each isolate's
protein both as a standalone mature-chain fragment (P04585.1.1, 131 aa MA,
target coverage 1.00) and inside the gag / gag-pol fragments that contain it
(target coverage 0.23-0.36), so a bare E-value cut returns ~3 entries per isolate.
Requiring the hit to span the fragment keeps "this fragment IS the protein".
It also drops the HTLV gag and Trypanosoma cruzi hits that the NC query picks up.

Query coverage is deliberately NOT filtered on: the library truncates env at 570 aa
(the entry is 824 aa), so Gp160 only covers 0.67 of its query. The same truncation
means Gp41 aligns to just the tail of the env fragments (target coverage 0.12-0.14)
and therefore ends up with no fragment in the CSV -- see the TSV, where those three
env hits are kept with kept_in_csv=False.

Requires pandas + openpyxl and jackhmmer on PATH. On this machine:
/usr/local/bin/python3 (3.10.7) and HMMER 3.4 at /opt/homebrew/bin/jackhmmer.
"""

import csv
import gzip
import shutil
import subprocess
import tempfile
from pathlib import Path

import pandas as pd

REPO = Path(__file__).resolve().parent.parent
CSV = REPO / "QC/Well_Known_interactions/jager_2011_HIV.csv"
TSV = REPO / "QC/Well_Known_interactions/jager_2011_HIV_jackhmmer_hits.tsv"
LIBRARY = REPO / "data/mmc2.xlsx"
SCREEN = REPO / "data/bait_prey_binary_all13.tsv.gz"

EVALUE = "1e-5"
ITERATIONS = "3"
MIN_TARGET_COVERAGE = 0.5
KEEP_ORGANISMS = ("HIV-1", "HIV-2")


def union_length(intervals):
    """Total length covered by a set of closed 1-based intervals, counting overlaps once."""
    total, start, end = 0, None, None
    for a, b in sorted(intervals):
        if start is None:
            start, end = a, b
        elif a <= end + 1:
            end = max(end, b)
        else:
            total += end - start + 1
            start, end = a, b
    return total if start is None else total + end - start + 1


def parse_domtbl(path):
    """Collapse a jackhmmer --domtblout into one record per (query, target) pair."""
    pairs = {}
    for line in open(path):
        if line.startswith("#"):
            continue
        f = line.split()
        target, tlen, query, qlen = f[0], int(f[2]), f[3], int(f[5])
        rec = pairs.setdefault(
            (query, target),
            dict(
                query=query,
                target=target,
                qlen=qlen,
                tlen=tlen,
                evalue=float(f[6]),
                score=float(f[7]),
                hmm=[],
                env=[],
            ),
        )
        rec["hmm"].append((int(f[15]), int(f[16])))  # coords on the query profile
        rec["env"].append((int(f[19]), int(f[20])))  # envelope coords on the fragment
    for rec in pairs.values():
        rec["query_coverage"] = union_length(rec["hmm"]) / rec["qlen"]
        rec["target_coverage"] = union_length(rec["env"]) / rec["tlen"]
    return list(pairs.values())


def main():
    if shutil.which("jackhmmer") is None:
        raise SystemExit("jackhmmer not found on PATH")

    rows = list(csv.DictReader(open(CSV)))
    baits = {}
    for row in rows:
        baits.setdefault(row["viral_protein"], row["viral_sequence"])

    library = pd.read_excel(LIBRARY, sheet_name=0, header=1)
    library = library.set_index("Fragment ID")

    with gzip.open(SCREEN, "rt") as fh:
        screened = set(fh.readline().rstrip("\n").split("\t")[1:])

    work = Path(tempfile.mkdtemp(prefix="jager_jackhmmer_"))
    query_fa, db_fa, domtbl = work / "queries.fa", work / "vorf.fa", work / "hits.dom"
    with open(query_fa, "w") as fh:
        for name, seq in baits.items():
            fh.write(f">{name}\n{seq}\n")
    with open(db_fa, "w") as fh:
        for fid, seq in zip(library.index, library["Fragment sequence"]):
            fh.write(f">{fid}\n{seq}\n")

    subprocess.run(
        ["jackhmmer", "-E", EVALUE, "-N", ITERATIONS,
         "--noali", "--cpu", "4", "--domtblout", str(domtbl), "-o", "/dev/null",
         str(query_fa), str(db_fa)],
        check=True,
    )

    hits = parse_domtbl(domtbl)
    for hit in hits:
        entry = library.loc[hit["target"]]
        hit["organism"] = str(entry["Organism"])
        hit["annotation"] = str(entry["Annotation"])
        hit["fragment_length"] = int(entry["Fragment length"])
        hit["in_screen"] = hit["target"] in screened
        hit["kept"] = (
            hit["target_coverage"] >= MIN_TARGET_COVERAGE
            and any(org in hit["organism"] for org in KEEP_ORGANISMS)
        )
    hits.sort(key=lambda h: (list(baits).index(h["query"]), h["evalue"]))

    kept = {}
    for hit in hits:
        if hit["kept"]:
            kept.setdefault(hit["query"], []).append(hit["target"])

    # tolerate a re-run: drop any existing column before re-appending it
    fields = [f for f in rows[0] if f != "vorf_fragment_ids"] + ["vorf_fragment_ids"]
    with open(CSV, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            row["vorf_fragment_ids"] = ";".join(kept.get(row["viral_protein"], []))
            writer.writerow(row)

    with open(TSV, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(
            ["viral_protein", "vorf_fragment_id", "evalue", "bitscore",
             "query_coverage", "target_coverage", "organism", "annotation",
             "fragment_length", "in_screen", "kept_in_csv"]
        )
        for hit in hits:
            writer.writerow(
                [hit["query"], hit["target"], f"{hit['evalue']:.2g}", f"{hit['score']:.1f}",
                 f"{hit['query_coverage']:.3f}", f"{hit['target_coverage']:.3f}",
                 hit["organism"], hit["annotation"], hit["fragment_length"],
                 hit["in_screen"], hit["kept"]]
            )
    shutil.rmtree(work)

    print(f"{len(hits)} raw hits -> {TSV.relative_to(REPO)}")
    print(f"{sum(len(v) for v in kept.values())} fragments assigned in {CSV.relative_to(REPO)}")
    for name in baits:
        ids = kept.get(name, [])
        print(f"  {name:6s} {len(ids):2d} fragments, {sum(i in screened for i in ids)} screened")


if __name__ == "__main__":
    main()
