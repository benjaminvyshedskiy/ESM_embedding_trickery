#!/usr/bin/env python3
"""Build the Jager benchmark interaction table from the Nature 2012 supplement PDF.

Source: 41586_2012_BFnature10719_MOESM288_ESM.pdf (Jager et al., Nature 2012,
doi:10.1038/nature10719).

  page 25     Supplementary Table 3: 39 well-characterized HIV-human protein pairs
  pages 40-43 Amino acid sequences of the HIV-SF constructs used for AP-MS

Every construct on pp. 40-43 is an SF fusion. Per the legend, underlined residues are
the 2xStrepTagII-TEV-3xFLAG (SF) tag and, for gp41, the influenza HA signal peptide.
They have to go, or jackhmmer will hit on tag sequence shared by all 17 baits. The tag
boundary is not a fixed motif -- the underline starts at GAAAG..., GLEGGGG... or
LEGGGG... depending on the construct -- so we take the boundary from the underline
rectangles the PDF actually draws rather than from a string match. Bold residues (the
site-directed mutations) are kept: they are what was in the purified construct.

Requires PyMuPDF. On this machine: /usr/local/bin/python3 (3.10.7).
"""

import csv
import re
from pathlib import Path

import fitz  # PyMuPDF

REPO = Path(__file__).resolve().parent.parent
PDF = REPO / "QC/Well_Known_interactions/41586_2012_BFnature10719_MOESM288_ESM.pdf"
OUT = REPO / "QC/Well_Known_interactions/jager_2011_HIV.csv"

TABLE_PAGE = 24  # 0-based; PDF page 25
SEQ_PAGES = range(39, 43)  # 0-based; PDF pages 40-43

# Table 3 bait name -> FASTA record name on pp. 40-43.
BAIT_TO_RECORD = {
    "ma": "MA",
    "nc": "NC",
    "vif": "nVif",
    "vpr": "hVpr",
    "tat": "Tat",
    "nef": "Nef",
    "gp160": "gp160",
    "gp120": "gp120",
    "gp41": "SP-gp41",
}

# Vpu is absent from the supplement's sequence listing (17 records, no Vpu) even though
# Table 3 scores 4 Vpu interactions. Filled from the pNL4-3 vpu CDS, GenBank AAK08488.1
# (AF324493), 81 aa -- the right strain, since the paper's RT/IN constructs match pNL4-3
# pol (AAK08484.1) exactly. UniProt has no NL4-3 Vpu entry; the nearest is P05919 (HXB2,
# 82 aa), which differs by an indel and several substitutions. Not PDF-derived.
VPU_NL43 = (
    "MQPIIVAIVALVVAIIIAIVVWSIVIIEYRKILRQRKIDRLIDRLIERAEDSGNESEGEV"
    "SALVEMGVEMGHHAPWDIDDL"
)

FASTA_HEADER = re.compile(r"^>?([A-Za-z0-9()\-]+)-SF$")
SEQ_LINE = re.compile(r"^[A-Za-z*]{5,80}$")
AA = set("ACDEFGHIKLMNPQRSTVWY")


def parse_table(page):
    """Return [(bait, uniprot_id, human_name)] from Supplementary Table 3."""
    lines = [ln.strip() for ln in page.get_text().split("\n")]
    start = lines.index("Protein name") + 1
    rows, buf = [], []
    for line in lines[start:]:
        if not line:
            continue
        if line.startswith(("SUPPLEMENTARY", "RESEARCH", "WWW.NATURE", "doi:")):
            break
        buf.append(line)
        if len(buf) == 3:
            rows.append(tuple(buf))
            buf = []
    if buf:
        raise ValueError(f"table ended mid-row: {buf}")
    return rows


def underline_rects(page):
    """Thin filled rectangles -- the underlines drawn under tag/signal-peptide residues."""
    return [
        d["rect"]
        for d in page.get_drawings()
        if d["rect"].height < 1.5 and d["rect"].width > 3
    ]


def is_underlined(char_bbox, rects):
    x0, _, x1, y1 = char_bbox
    return any(
        y1 - 1.0 <= r.y0 <= y1 + 3.5 and r.x0 <= x0 + 1.0 and r.x1 >= x1 - 1.0
        for r in rects
    )


def parse_sequences(doc):
    """Return {record name: mature sequence} with underlined residues stripped."""
    records, current = {}, None
    for page_no in SEQ_PAGES:
        page = doc[page_no]
        rects = underline_rects(page)
        for block in page.get_text("rawdict")["blocks"]:
            if block["type"] != 0:
                continue
            for line in block["lines"]:
                chars = [c for span in line["spans"] for c in span["chars"]]
                text = "".join(c["c"] for c in chars).strip()
                header = FASTA_HEADER.match(text)
                if header:
                    current = header.group(1)
                    records[current] = []
                    continue
                if not SEQ_LINE.fullmatch(text):
                    if text:  # prose or page footer -- the record is over
                        current = None
                    continue
                if current is None:
                    continue
                for c in chars:
                    if c["c"] in " *":
                        continue
                    if not is_underlined(c["bbox"], rects):
                        records[current].append(c["c"])
    return {name: "".join(seq) for name, seq in records.items()}


def main():
    doc = fitz.open(PDF)
    interactions = parse_table(doc[TABLE_PAGE])
    sequences = parse_sequences(doc)

    rows = []
    for bait, uniprot, human_name in interactions:
        key = bait.lower()
        if key == "vpu":
            seq = VPU_NL43
        else:
            seq = sequences[BAIT_TO_RECORD[key]]
        if set(seq) - AA:
            raise ValueError(f"{bait}: non-standard residues {set(seq) - AA}")
        rows.append([bait, seq, uniprot, human_name])

    with open(OUT, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(
            ["viral_protein", "viral_sequence", "human_uniprot_id", "human_protein_name"]
        )
        writer.writerows(rows)

    print(f"wrote {len(rows)} interactions to {OUT.relative_to(REPO)}")
    seen = {}
    for bait, seq, _, _ in rows:
        seen.setdefault(bait, [0, len(seq)])[0] += 1
    for bait, (n, length) in seen.items():
        print(f"  {bait:6s} {n} pairs  {length:4d} aa")


if __name__ == "__main__":
    main()
