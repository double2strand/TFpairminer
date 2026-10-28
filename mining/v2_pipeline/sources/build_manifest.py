#!/usr/bin/env python3
"""Enumerate every available input for TFpairminer v2 from the live sources and write samples.tsv.

  CCLA (scatlas HCL cancer cell lines, CNGB SCDS0000040): every *.fragments.tsv.gz (1 or 2 replicates per line)
  Human scATAC Corpus (Tsinghua): every sample (dataset x tissue) - cell x cCRE h5ad + standardized fragments
  ENCODE ChromBPNet ATAC models (HF kundajelab/encode-chrombpnet-ATAC-*): fold_0 bias-corrected model per CCLA line
Sizes come from HTTP Content-Length. MD5s are the publisher's where they exist (Corpus: per-dataset tarball only;
CNGB: none). HF models carry their LFS sha256. Stage 01 verifies size (+md5/sha256 when given) and writes a
sha256 lock (inputs.lock.tsv) so later re-downloads must match byte-for-byte.

Usage: python build_manifest.py --out samples.tsv [--no-head]   (network; ~3 min)
"""
import argparse, csv, json, re, subprocess, sys, urllib.request
from concurrent.futures import ThreadPoolExecutor

CNGB = "https://ftp.cngb.org/pub/SciRAID/cdcp/SCDS0000040/download"
HSC = "https://health.tsinghua.edu.cn/human-scatac-corpus"
UA = {"User-Agent": "Mozilla/5.0 TFpairminer-v2-pipeline"}
FIELDS = ["profile_id", "atlas", "source", "kind", "context", "tissue", "inputs", "sizes", "md5", "checksum_note",
          "fragments", "fragment_sizes", "celltype_col", "n_cells", "n_cell_types", "accession", "pmid",
          "cbp_pretrained", "cbp_pretrained_sha256", "rna_context", "mapping_level", "enabled", "disabled_reason"]

# ENCODE experiment to prefer when several ChromBPNet models exist (matches the published v2 HepG2/K562 runs).
PREFERRED = {"HepG2": "ENCSR291GJU", "K562": "ENCSR483RKN"}
# Cell line -> RNA context used by the targeted proxy screen (cancer line never 'exact').
LINE_CONTEXT = {"HepG2": "hepatocyte", "Huh7": "hepatocyte", "K562": "erythroid_myeloid_progenitor"}
LINE_DESC = {"K562": "chronic myelogenous leukemia", "HepG2": "hepatoblastoma/liver cancer", "Huh7": "hepatocellular carcinoma",
             "A549": "lung adenocarcinoma", "MCF7": "breast adenocarcinoma (ER+)", "HCT116": "colorectal carcinoma",
             "786-O": "clear-cell renal carcinoma", "HT-29": "colorectal adenocarcinoma", "SW480": "colorectal adenocarcinoma (primary)",
             "SW620": "colorectal adenocarcinoma (metastasis)", "DLD-1": "colorectal adenocarcinoma", "MDA-MB-468": "triple-negative breast cancer",
             "MDA-MB-231": "triple-negative breast cancer", "T-47D": "breast ductal carcinoma", "RPMI8226": "multiple myeloma",
             "Hap1": "near-haploid CML-derived", "RPE-1": "hTERT retinal pigment epithelium", "HK-2": "proximal tubule (HPV-immortalised)"}


def get(url, data=None):
    req = urllib.request.Request(url, data=data, headers=UA)
    with urllib.request.urlopen(req, timeout=120) as r:
        return r.read().decode()


def head_status(url):
    """Final HTTP status after redirects (HEAD); '000' on network failure."""
    return subprocess.run(["curl", "-sIL", "-m", "60", "-o", "/dev/null", "-w", "%{http_code}", "-A", UA["User-Agent"], url],
                          capture_output=True, text=True).stdout.strip() or "000"


def head_size(url):
    for _ in range(3):
        try:
            out = subprocess.run(["curl", "-sIL", "-m", "60", "-A", UA["User-Agent"], url], capture_output=True, text=True).stdout
            v = [l.split(":", 1)[1].strip() for l in out.splitlines() if l.lower().startswith("content-length")]
            if v:
                return v[-1]
        except Exception:
            pass
    return ""


def hf_models():
    out, url = [], "https://huggingface.co/api/models?author=kundajelab&search=encode-chrombpnet-ATAC&limit=1000"
    out = json.loads(get(url))
    return [m["id"] for m in out]


def hf_nobias(repo):
    try:
        tree = json.loads(get(f"https://huggingface.co/api/models/{repo}/tree/main/fold_0"))
    except Exception:
        return "", ""
    for x in tree:
        if re.search(r"chrombpnet_nobias\.fold_0\..*\.h5$", x["path"]):
            return f"https://huggingface.co/{repo}/resolve/main/{x['path']}", (x.get("lfs") or {}).get("oid", "")
    return "", ""


def norm(s):
    return re.sub(r"[^a-z0-9]", "", s.lower())


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--no-head", action="store_true")
    a = ap.parse_args()
    rows = []
    # ---- CCLA
    html = get(CNGB + "/")
    frags = sorted(set(re.findall(r'href="([^"/]+\.fragments\.tsv\.gz)"', html)))
    lines = {}
    for f in frags:
        lines.setdefault(re.sub(r"(_[12])?\.fragments\.tsv\.gz$", "", f), []).append(f)
    models = hf_models(); print(f"{len(lines)} CCLA lines, {len(models)} ENCODE ChromBPNet ATAC models", file=sys.stderr)
    for line, fs in sorted(lines.items()):
        urls = [f"{CNGB}/{f}" for f in sorted(fs)]
        pre = [m for m in models if norm(m.split("encode-chrombpnet-ATAC-")[1].rsplit("-ENCSR", 2)[0]) == norm(line)]
        pre = sorted(pre, key=lambda m: (PREFERRED.get(line, "~") not in m, m))
        purl, psha = hf_nobias(pre[0]) if pre else ("", "")
        rows.append({"profile_id": f"CCLA_{line}", "atlas": "CCLA", "source": "CNGB SCDS0000040", "kind": "fragments",
                     "context": f"{line} {LINE_DESC.get(line, 'cancer cell line')}", "tissue": "", "inputs": ",".join(urls),
                     "sizes": "", "md5": "", "checksum_note": "publisher gives no checksum; size checked + sha256 pinned on first fetch",
                     "fragments": ",".join(urls), "fragment_sizes": "", "celltype_col": "", "n_cells": "", "n_cell_types": "1",
                     "accession": "CNGB:SCDS0000040", "pmid": "", "cbp_pretrained": purl, "cbp_pretrained_sha256": psha,
                     "rna_context": LINE_CONTEXT.get(line, ""), "mapping_level": "proxy" if line in LINE_CONTEXT else "unmapped"})
    # ---- Corpus
    dl = get(HSC + "/download.php")
    ds_md5 = {}
    for d in sorted(set(re.findall(r"categorydata/Dataset/([^/]+)/[^\"]+-cell_by_cCRE\.tar\.gz\"", dl))):
        try:
            ds_md5[d] = get(f"{HSC}/data/categorydata/Dataset/{d}/{d}-cell_by_cCRE.tar.gz.md5").split()[0]
        except Exception:
            ds_md5[d] = ""
    js = json.loads(get(HSC + "/getResult_browse.php", b"filter=&item=&draw=1&start=0&length=-1"))["Data"]
    strip = lambda s: re.sub(r"<[^>]+>", "", s or "").replace("\xa0", " ").strip()
    for s in js:
        h5 = re.search(r"href='\./([^']+cell_by_cCRE\.h5ad)'", s["cell_by_cre_matrix"] or "")
        fr = re.search(r"href='\./([^']+fragments_standardized\.tar\.gz)'", s["fragment_file"] or "")
        if not h5:
            continue
        rows.append({"profile_id": "Corpus_" + re.sub(r"[^A-Za-z0-9]+", "_", s["sample"]).strip("_"), "atlas": "Corpus",
                     "source": "Human scATAC Corpus", "kind": "corpus", "context": s["sample"], "tissue": s["tissue"],
                     "inputs": f"{HSC}/{h5.group(1)}", "sizes": "", "md5": "",
                     "checksum_note": f"publisher md5 only for dataset tarball {s['dataset']}-cell_by_cCRE.tar.gz={ds_md5.get(s['dataset'], '')}; size checked + sha256 pinned",
                     "fragments": f"{HSC}/{fr.group(1)}" if fr else "", "fragment_sizes": "", "celltype_col": "Cell_type (HSC)",
                     "n_cells": s.get("number_of_cells", ""), "n_cell_types": s.get("number_of_cell_types", ""),
                     "accession": strip(s.get("accession")), "pmid": s.get("pmid", ""), "cbp_pretrained": "", "cbp_pretrained_sha256": "",
                     "rna_context": "auto", "mapping_level": "auto"})
    for r in rows:
        r.setdefault("enabled", "1"); r.setdefault("disabled_reason", "")
    if not a.no_head:
        def fill(r):
            urls = [u for c in ("inputs", "fragments", "cbp_pretrained") for u in r.get(c, "").split(",") if u]
            bad = [f"{u.rsplit('/', 1)[-1]}=HTTP {c}" for u in urls for c in [head_status(u)] if c != "200"]
            if bad:   # pipeline skips disabled rows instead of failing in stage 01
                r["enabled"], r["disabled_reason"] = "0", "unavailable at build time: " + "; ".join(bad)
            r["sizes"] = ",".join(head_size(u) for u in r["inputs"].split(","))
            r["fragment_sizes"] = ",".join(head_size(u) for u in r["fragments"].split(",") if u)
        with ThreadPoolExecutor(8) as ex:
            list(ex.map(fill, rows))
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, FIELDS, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
    print(f"wrote {len(rows)} rows -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
