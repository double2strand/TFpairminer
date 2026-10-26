"""tfpm v2 CLI.  python -m motif_cooccurrence.cli <cmd> ...   (run from mining/)

Manifest (TSV): sample_id dataset context kind inputs [celltype_col]
  kind=fragments : inputs = comma-separated replicate *.fragments.tsv.gz (paths or URLs)
  kind=corpus    : inputs = one Human-scATAC-Corpus cell x cCRE .h5ad (path or URL); each
                   cell type with >= --min-cells cells becomes its own sub-sample
"""
import argparse, json, os, shutil, subprocess, sys, time, urllib.request
from pathlib import Path

from . import engine as E

HERE = Path(__file__).resolve().parent


def read_manifest(path):
    rows = E.read_tsv(path)
    for r in rows:
        r["inputs"] = [x for x in r["inputs"].split(",") if x]
    return rows


def fetch(src, dest_dir):
    dest_dir = Path(dest_dir); dest_dir.mkdir(parents=True, exist_ok=True)
    if "://" not in src:
        return src
    out = dest_dir / src.rstrip("/").split("/")[-1]
    if not out.exists() or out.stat().st_size == 0:
        tmp = str(out) + ".part"
        subprocess.run(["curl", "-fL", "--retry", "5", "-sS", "-o", tmp, src], check=True)
        os.replace(tmp, out)
    return str(out)


def run_sample(a):
    man = {r["sample_id"]: r for r in read_manifest(a.manifest)}
    r = man[a.sample]; out = Path(a.out) / "samples" / r["sample_id"]; out.mkdir(parents=True, exist_ok=True)
    if (out / "DONE").exists() and not a.force:
        print("done already", out); return
    inputs = [fetch(x, Path(a.out) / "inputs") for x in r["inputs"]]
    chroms = a.chroms.split(",") if a.chroms else None
    panel = [x for x in a.panel.split(",") if x] if a.panel else None
    t0 = time.time(); units = []
    reuse = not a.redo_windows
    if r["kind"] == "fragments" and reuse and (out / "windows.tsv").exists() and (out / "windows_meta.json").exists():
        meta = json.loads((out / "windows_meta.json").read_text())
        units.append((r["sample_id"], r["context"], out / "windows.tsv", out))
    elif r["kind"] != "fragments" and reuse and (out / "celltypes" / "corpus_windows.json").exists():
        meta = json.loads((out / "celltypes" / "corpus_windows.json").read_text())
        for ct in meta["cell_types"]:
            units.append((f"{r['sample_id']}__{ct['slug']}", f"{r['context']} | {ct['cell_type']}",
                          out / "celltypes" / f"{ct['slug']}.windows.tsv", out / "celltypes" / ct["slug"]))
    elif r["kind"] == "fragments":
        meta = E.windows_from_fragments(inputs, a.genome, a.blacklist, out / "windows.tsv", a.window_bp,
                                        a.max_windows, a.min_cuts, a.min_rep_cuts, chroms, r["sample_id"])
        units.append((r["sample_id"], r["context"], out / "windows.tsv", out))
    else:
        meta = E.windows_from_corpus(inputs[0], a.genome, a.blacklist, out / "celltypes",
                                     r.get("celltype_col") or "Cell_type (HSC)", a.window_bp, a.max_windows,
                                     a.min_cells, label=r["sample_id"])
        for ct in meta["cell_types"]:
            units.append((f"{r['sample_id']}__{ct['slug']}", f"{r['context']} | {ct['cell_type']}",
                          out / "celltypes" / f"{ct['slug']}.windows.tsv", out / "celltypes" / ct["slug"]))
    (out / "windows_meta.json").write_text(json.dumps(meta, indent=1, default=str))
    sums = []
    for label, context, wt, ud in units:
        ud.mkdir(parents=True, exist_ok=True)
        if a.redo_windows or not (ud / "hits.npz").exists():
            E.scan_windows(wt, a.genome, a.meme, ud / "hits.npz", a.motif_pvalue, panel)

        s = E.analyse(wt, ud / "hits.npz", a.meme, ud, label, r["dataset"], context,
                      permutations=a.permutations, q_cut=a.q, min_lift=a.min_lift)
        s["kind"] = r["kind"]; s["parent_sample"] = r["sample_id"]
        (ud / "summary.json").write_text(json.dumps(s, indent=1)); sums.append(s)
        if not a.keep_hits:
            pass
    (out / "DONE").write_text(json.dumps({"seconds": round(time.time() - t0), "units": [s["sample"] for s in sums]}))
    print(json.dumps({"sample": r["sample_id"], "units": len(sums), "seconds": round(time.time() - t0)}))


def plan(a):
    """Print one command per manifest row (consumed by tfpm v2 run --runner wr|local)."""
    for r in read_manifest(a.manifest):
        extra = " ".join(a.passthrough or [])
        print(f"cd {HERE.parent} && {sys.executable} -m motif_cooccurrence.cli sample --manifest {Path(a.manifest).resolve()} "
              f"--sample {r['sample_id']} --out {Path(a.out).resolve()} --genome {a.genome} --meme {a.meme} "
              f"--blacklist {a.blacklist} {extra}".strip())


def status(a):
    man = read_manifest(a.manifest); done = 0
    for r in man:
        d = Path(a.out) / "samples" / r["sample_id"]
        st = "done" if (d / "DONE").exists() else ("started" if d.exists() else "pending")
        done += st == "done"; print(f"{r['sample_id']}\t{r['dataset']}\t{st}")
    print(f"{done}/{len(man)} samples done")


def main(argv=None):
    p = argparse.ArgumentParser(prog="tfpm v2"); sp = p.add_subparsers(dest="cmd", required=True)
    g = sp.add_parser("prep-genome"); g.add_argument("--fasta", required=True); g.add_argument("--out", required=True)
    for name in ("sample", "plan"):
        s = sp.add_parser(name)
        s.add_argument("--manifest", required=True); s.add_argument("--out", required=True)
        s.add_argument("--genome", required=True, help="dir from prep-genome")
        s.add_argument("--meme", required=True); s.add_argument("--blacklist", default="")
        if name == "sample":
            s.add_argument("--sample", required=True)
            s.add_argument("--window-bp", type=int, default=300); s.add_argument("--max-windows", type=int, default=20000)
            s.add_argument("--min-cuts", type=int, default=6); s.add_argument("--min-rep-cuts", type=int, default=2)
            s.add_argument("--min-cells", type=int, default=100)
            s.add_argument("--chroms", default=""); s.add_argument("--panel", default="", help="comma list of JASPAR ids/names; default all")
            s.add_argument("--motif-pvalue", type=float, default=1e-4); s.add_argument("--permutations", type=int, default=999)
            s.add_argument("--q", type=float, default=0.05); s.add_argument("--min-lift", type=float, default=1.25)
            s.add_argument("--keep-hits", action="store_true"); s.add_argument("--force", action="store_true", help="re-run analysis (reuses windows + scan)")
            s.add_argument("--redo-windows", action="store_true", help="recompute windows and motif scan too")
        else:
            s.add_argument("passthrough", nargs="*")
    st = sp.add_parser("status"); st.add_argument("--manifest", required=True); st.add_argument("--out", required=True)
    ag = sp.add_parser("aggregate"); ag.add_argument("--out", required=True); ag.add_argument("--site", required=True)
    ag.add_argument("--meme", required=True)
    a = p.parse_args(argv)
    if a.cmd == "prep-genome":
        E.prep_genome(a.fasta, a.out)
    elif a.cmd == "sample":
        run_sample(a)
    elif a.cmd == "plan":
        plan(a)
    elif a.cmd == "status":
        status(a)
    elif a.cmd == "aggregate":
        from .aggregate import aggregate
        aggregate(a.out, a.site, a.meme)


if __name__ == "__main__":
    main()
