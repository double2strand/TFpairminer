#!/usr/bin/env python3
"""Export Gene_Register and Interaction_Priors sheets of the TF compendium xlsx to TSV."""
import argparse, csv
from pathlib import Path
import openpyxl

SHEETS = {"Gene_Register": "gene_register.tsv", "Interaction_Priors": "interaction_priors.tsv"}

def cell(v):
    if v is None:
        return ""
    if isinstance(v, float) and v.is_integer():
        v = int(v)
    return str(v).replace("\t", " ").replace("\r", " ").replace("\n", " ")

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("xlsx")
    p.add_argument("--outdir", default="compendium")
    a = p.parse_args()
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    wb = openpyxl.load_workbook(a.xlsx, read_only=True, data_only=True)
    for sheet, name in SHEETS.items():
        rows = wb[sheet].iter_rows(values_only=True)
        header = [cell(h) for h in next(rows)]
        while header and header[-1] == "":
            header.pop()
        n = 0
        with open(out / name, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(header)
            for r in rows:
                vals = [cell(v) for v in r[:len(header)]]
                if not any(vals):
                    continue
                w.writerow(vals); n += 1
        print(f"{sheet}: {n} rows -> {out/name}")

if __name__ == "__main__":
    main()
