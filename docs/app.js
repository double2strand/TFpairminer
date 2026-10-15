"use strict";
/* TF Pair Miner site — reads docs/data/runs/index.json (schema_version 2) and sharded per-run JSON. */
const D = {cache: {}}; const S = {mode: "visual", q: "", set: null, setLabel: "", sort: {k: "n_collections", dir: -1},
  f: {r: "", q: "", mq: "", lift: "", coll: "", prior: false, ct: ""}};
const $ = (s, el = document) => el.querySelector(s);
const esc = s => String(s ?? "").replace(/[&<>"']/g, c => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"}[c]));
const fmt = (v, d = 3) => v == null || v === "" ? "–" : (typeof v === "number" && !Number.isInteger(v) ? (Math.abs(v) < 1e-3 && v !== 0 ? v.toExponential(1) : v.toFixed(d)) : v);
const pct = v => v == null ? "–" : (v * 100).toFixed(1) + "%";
const base = () => `data/runs/${D.runId}/`;
async function J(path) { if (!D.cache[path]) D.cache[path] = fetch(path, {cache: "no-cache"}).then(r => { if (!r.ok) throw new Error(path + " " + r.status); return r.json(); }); return D.cache[path]; }
const rowsOf = o => o.rows.map(r => Object.fromEntries(o.cols.map((c, i) => [c, r[i]])));
function enrich(P) { P.forEach(p => { p.st = D.ST[p.s]; p.ct = p.st.cell_type; p.tissue = p.st.tissue; p.disease = p.st.disease; p.assay = p.st.assay; p.key = p.tf1 + ":" + p.tf2; }); return P; }

async function boot() {
  try { D.index = await J("data/runs/index.json"); } catch (e) { $("#app").innerHTML = `<div class="card">No published runs yet. Run <code>./tfpm publish</code>.</div>`; return; }
  const want = new URLSearchParams(location.search).get("run");
  const sel = $("#run"); sel.innerHTML = D.index.runs.map(r => `<option value="${esc(r.run_id)}">${esc(r.run_id === D.index.latest ? "★ " : "")}${esc(r.label || r.run_id)} · ${esc(r.exported_at.slice(0, 10))}${r.pilot ? " · pilot" : ""}</option>`).join("");
  sel.onchange = () => { const u = new URL(location); u.searchParams.set("run", sel.value); location.href = u; };
  D.runId = D.index.runs.some(r => r.run_id === want) ? want : D.index.latest; sel.value = D.runId;
  [D.run, D.strata, D.tfs, D.top, D.heat] = await Promise.all(["run.json", "strata.json", "tfs.json", "top_pairs.json", "heatmap.json"].map(f => J(base() + f)));
  D.ST = rowsOf(D.strata); D.TF = Object.fromEntries(rowsOf(D.tfs).map(t => [t.tf.toUpperCase(), t]));
  D.TOP = enrich(rowsOf(D.top)); S.set = D.TOP; S.setLabel = `top ${D.TOP.length} of ${D.top.total} candidate rows`;
  const sm = D.run.summary, b = $("#banner");
  const notLatest = D.runId !== D.index.latest ? " You are viewing an older run." : "";
  if (D.run.pilot) { b.hidden = false; b.innerHTML = `<b>Pilot data — not an atlas-wide result.</b> ${sm.n_completed_datasets}/${sm.n_planned_datasets} planned dataset(s) (selection “${esc(sm.dataset_selection)}”), max ${sm.run_config.max_strata || "∞"} strata per dataset. Replace by publishing the full run.${notLatest}`; }
  else if (notLatest) { b.hidden = false; b.textContent = notLatest; }
  window.addEventListener("hashchange", route); route();
}
const stratumLabel = st => `${st.cell_type} · ${st.tissue} · ${st.disease} · ${st.assay}`;

/* ---------- search: resolve query to a (sharded) row set ---------- */
async function runSearch() {
  const q = S.q.trim(); if (!q) { S.set = D.TOP; S.setLabel = `top ${D.TOP.length} of ${D.top.total} candidate rows`; return; }
  const parts = q.split(/[:\s,/+]+|-(?=[A-Za-z])/).filter(Boolean).map(x => x.toUpperCase());
  if (parts.length === 2 && D.TF[parts[0]] && D.TF[parts[1]]) {
    const t = await J(base() + `tf/${D.TF[parts[0]].file}.json`);
    S.set = enrich(rowsOf(t)).filter(p => [p.tf1.toUpperCase(), p.tf2.toUpperCase()].sort().join() === parts.slice().sort().join());
    S.setLabel = `pair ${parts.join(":")}`; return; }
  if (D.TF[q.toUpperCase()]) { const t = await J(base() + `tf/${D.TF[q.toUpperCase()].file}.json`); S.set = enrich(rowsOf(t)); S.setLabel = `all rows for ${q.toUpperCase()}`; return; }
  const ql = q.toLowerCase();
  const sts = D.ST.map((st, i) => [i, st]).filter(([i, st]) => stratumLabel(st).toLowerCase().includes(ql) || st.cell_type_id.toLowerCase() === ql).slice(0, 30);
  if (sts.length) { const shards = await Promise.all(sts.map(([i]) => J(base() + `ct/${i}.json`))); S.set = enrich(shards.flatMap(rowsOf));
    S.setLabel = `${sts.length} cell context(s) matching “${q}”`; S.ctProfiles = shards; return; }
  const tfs = Object.keys(D.TF).filter(t => t.startsWith(q.toUpperCase())).slice(0, 8);
  S.set = D.TOP.filter(p => p.tf1.toUpperCase().startsWith(q.toUpperCase()) || p.tf2.toUpperCase().startsWith(q.toUpperCase()));
  S.setLabel = `top rows with TF starting “${q}”` + (tfs.length ? ` — TFs: ${tfs.join(", ")}` : " — no exact TF / cell type match");
}
function filtered() {
  const f = S.f;
  return S.set.filter(p => (f.r === "" || (p.primary_r ?? p.r) >= +f.r) && (f.q === "" || (p.min_q != null && p.min_q <= +f.q)) && (f.mq === "" || (p.meta_q != null && p.meta_q <= +f.mq))
    && (f.lift === "" || p.lift >= +f.lift) && (f.coll === "" || p.n_collections >= +f.coll) && (!f.prior || p.prior) && (f.ct === "" || String(p.s) === f.ct));
}
async function route() {
  const h = decodeURIComponent(location.hash.slice(1)); const [m, arg] = h.split("/");
  if (["visual", "science", "info"].includes(m)) S.mode = m;
  document.querySelectorAll("nav button").forEach(b => b.classList.toggle("on", b.dataset.mode === S.mode && m !== "tf"));
  if (m === "tf" && arg) return tfPage(arg);
  await runSearch(); ({visual, science, info}[S.mode])();
}
function go(h) { if (location.hash === "#" + h) route(); else location.hash = h; }
function stat(v, l) { return `<div class="stat"><b>${esc(v)}</b>${esc(l)}</div>`; }
function sortBy(a, k, dir) { return a.sort((x, y) => { const u = x[k] ?? -Infinity, v = y[k] ?? -Infinity; return u > v ? dir : u < v ? -dir : ((y.primary_r ?? y.r) - (x.primary_r ?? x.r)); }); }

/* ---------- visual ---------- */
function visual() {
  const P = filtered(), c = D.run.counts;
  $("#app").innerHTML = `<div class="card"><div class="stats">${stat(c.pair_rows, "pair × context rows (run)")}${stat(c.unique_pairs, "unique TF pairs")}${stat(c.strata, "cell contexts")}${stat(c.pairs_with_prior, "rows with literature prior")}${stat(D.run.manifest.resolved_census_version, "Census release")}</div>
   <p class="muted">Showing <b>${esc(S.setLabel)}</b> (${P.length} rows). ${D.TF[S.q.trim().toUpperCase()] ? `<a onclick="go('tf/${esc(S.q.trim().toUpperCase())}')">Open ${esc(S.q.trim().toUpperCase())} page →</a>` : "Search a TF, pair (A:B) or cell type; click a node for its TF page."}</p></div>
   <div class="grid"><div class="card"><h2>Pair network</h2><p class="muted">Edges = candidate pairs (width ∝ primary r; green = literature prior). Top 150 edges by collections, then r.</p><div id="net"></div></div>
   <div class="card"><h2>TF × cell-context detection</h2><p class="muted">${S.q ? "TFs and contexts in the current result." : "Top TFs × contexts of the run."} Colour = fraction of cells detecting the TF.</p><div id="heat" class="heat"></div></div></div>
   ${controlsCard()}
   <div class="card"><h2>Top candidates</h2>${pairTable(sortBy(P.slice(), "n_collections", -1).slice(0, 25))}</div>`;
  network(P); heatmap(P);
}
function controlsCard() {
  const c = D.run.positive_controls || []; if (!c.length) return "";
  return `<div class="card"><h2>Positive controls</h2><p class="muted">Known partner pairs checked in this run (a sanity check of sensitivity; with pilot data most are expectedly absent).</p><div class="chips">${c.map(x => `<a onclick="S.q='${esc(x.pair.replace("-", ":"))}';$('#q').value=S.q;go('science')" title="contexts ${x.n_contexts_with_hit}, best r ${x.best_primary_r || "–"}, min q ${x.min_q || "–"}" style="${+x.n_contexts_with_hit ? "background:#dcfce7" : ""}">${esc(x.pair)} ${+x.n_contexts_with_hit ? "✓ " + x.n_contexts_with_hit : "–"}</a>`).join("")}</div></div>`;
}
function network(P) {
  const E = new Map();
  for (const p of P) { const r = p.primary_r ?? p.r; const e = E.get(p.key) || {a: p.tf1, b: p.tf2, r: -1, coll: 0, prior: 0, n: 0}; e.r = Math.max(e.r, r); e.coll = Math.max(e.coll, p.n_collections); e.prior |= p.prior; e.n++; E.set(p.key, e); }
  const edges = [...E.values()].sort((x, y) => y.coll - x.coll || y.r - x.r).slice(0, 150);
  const W = 640, H = 480; const ids = [...new Set(edges.flatMap(e => [e.a, e.b]))];
  if (!ids.length) { $("#net").innerHTML = '<p class="muted">No pairs match.</p>'; return; }
  const N = ids.map((id, i) => ({id, x: W / 2 + 200 * Math.cos(i * 2.4) * Math.sqrt((i + 1) / ids.length), y: H / 2 + 180 * Math.sin(i * 2.4) * Math.sqrt((i + 1) / ids.length), vx: 0, vy: 0, deg: 0}));
  const idx = Object.fromEntries(N.map((n, i) => [n.id, i])); edges.forEach(e => { N[idx[e.a]].deg++; N[idx[e.b]].deg++; });
  for (let it = 0; it < 300; it++) { const t = 1 - it / 300;
    for (let i = 0; i < N.length; i++) for (let j = i + 1; j < N.length; j++) { const a = N[i], b = N[j]; const dx = a.x - b.x, dy = a.y - b.y, d2 = dx * dx + dy * dy + .01, f = 900 / d2; a.vx += dx * f; a.vy += dy * f; b.vx -= dx * f; b.vy -= dy * f; }
    for (const e of edges) { const a = N[idx[e.a]], b = N[idx[e.b]]; const dx = b.x - a.x, dy = b.y - a.y, d = Math.hypot(dx, dy) + .01, f = (d - 60) * .02; a.vx += dx / d * f; a.vy += dy / d * f; b.vx -= dx / d * f; b.vy -= dy / d * f; }
    for (const n of N) { n.vx += (W / 2 - n.x) * .005; n.vy += (H / 2 - n.y) * .005; n.x += Math.max(-10, Math.min(10, n.vx)) * t; n.y += Math.max(-10, Math.min(10, n.vy)) * t; n.vx *= .5; n.vy *= .5; n.x = Math.max(20, Math.min(W - 40, n.x)); n.y = Math.max(12, Math.min(H - 12, n.y)); } }
  const q = S.q.trim().toUpperCase();
  $("#net").innerHTML = `<svg viewBox="0 0 ${W} ${H}" width="100%">${edges.map(e => { const a = N[idx[e.a]], b = N[idx[e.b]];
    return `<line x1="${a.x}" y1="${a.y}" x2="${b.x}" y2="${b.y}" stroke="${e.prior ? "#16a34a" : "#94a3b8"}" stroke-opacity=".7" stroke-width="${.5 + 4 * Math.max(0, e.r)}"><title>${esc(e.a)}–${esc(e.b)}: max r ${e.r.toFixed(3)}, ${e.n} context(s)${e.prior ? ", literature prior" : ""}</title></line>`; }).join("")}
    ${N.map(n => `<g class="node" onclick="go('tf/${esc(n.id)}')" style="cursor:pointer"><circle cx="${n.x}" cy="${n.y}" r="${3 + Math.sqrt(n.deg) * 2}" fill="${n.id.toUpperCase() === q ? "#dc2626" : "#2563eb"}"/><text x="${n.x + 6}" y="${n.y + 3}">${esc(n.id)}</text></g>`).join("")}</svg>`;
}
async function heatmap(P) {
  let tfs, ss, val = new Map();
  if (!S.q) { tfs = D.heat.tfs; ss = D.heat.strata; tfs.forEach((t, y) => ss.forEach((s, x) => val.set(t + "|" + s, D.heat.detection[y][x]))); }
  else { const deg = {}; P.forEach(p => { deg[p.tf1] = (deg[p.tf1] || 0) + 1; deg[p.tf2] = (deg[p.tf2] || 0) + 1; });
    tfs = Object.keys(deg).sort((a, b) => deg[b] - deg[a]).slice(0, 30); ss = [...new Set(P.map(p => p.s))].slice(0, 40);
    const shards = await Promise.all(tfs.filter(t => D.TF[t.toUpperCase()]).map(t => J(base() + `tf/${D.TF[t.toUpperCase()].file}.json`)));
    shards.forEach(sh => sh.profile.forEach(([s, d]) => val.set(sh.tf + "|" + s, d))); }
  const el = $("#heat"); if (!el) return;
  if (!tfs.length || !ss.length) { el.innerHTML = '<p class="muted">Nothing to show.</p>'; return; }
  const cw = 24, ch = 14, lw = 70, th = 150;
  const cells = []; tfs.forEach((tf, y) => ss.forEach((s, x) => { const d = val.get(tf + "|" + s);
    cells.push(`<rect x="${lw + x * cw}" y="${th + y * ch}" width="${cw - 1}" height="${ch - 1}" fill="${d == null ? "#f3f4f6" : `rgba(37,99,235,${.08 + .92 * Math.min(1, d)})`}"><title>${esc(tf)} in ${esc(stratumLabel(D.ST[s]))}: detection ${d == null ? "< " + D.run.profile_detection_floor + " or not measured" : pct(d)}</title></rect>`); }));
  el.innerHTML = `<svg width="${lw + ss.length * cw + 10}" height="${th + tfs.length * ch + 10}">${ss.map((s, x) => `<text transform="translate(${lw + x * cw + 13},${th - 4}) rotate(-60)">${esc(D.ST[s].cell_type.slice(0, 28))}</text>`).join("")}
   ${tfs.map((t, y) => `<text x="${lw - 4}" y="${th + y * ch + 10}" text-anchor="end" onclick="go('tf/${esc(t)}')" style="cursor:pointer">${esc(t)}</text>`).join("")}${cells.join("")}</svg>`;
}

/* ---------- TF page ---------- */
async function tfPage(tf) {
  const t = D.TF[tf.toUpperCase()]; if (!t) { $("#app").innerHTML = `<div class="card">${esc(tf)} is not in this run.</div>`; return; }
  const sh = await J(base() + `tf/${t.file}.json`); const P = enrich(rowsOf(sh)); tf = t.tf;
  const partners = {}; P.forEach(p => { const o = p.tf1 === tf ? p.tf2 : p.tf1; (partners[o] ||= []).push(p); });
  const prof = sh.profile.slice().sort((a, b) => b[1] - a[1]);
  $("#app").innerHTML = `<div class="card"><h2>${esc(tf)}</h2><div class="stats">${stat(Object.keys(partners).length, "candidate partners")}${stat(P.length, "pair × context rows")}${stat(prof.length, "contexts with detection ≥ " + D.run.profile_detection_floor)}</div>
   <h3>Partners</h3><div class="chips">${Object.entries(partners).sort((a, b) => b[1].length - a[1].length).map(([o, ps]) => `<a onclick="go('tf/${esc(o)}')" title="${ps.length} context(s)">${esc(o)} (${ps.length})${ps.some(p => p.prior) ? " ★" : ""}</a>`).join("") || '<span class="muted">none in this run</span>'}</div><p class="muted">★ = literature/physical prior in the compendium.</p></div>
   <div class="card"><h2>Pair rows involving ${esc(tf)}</h2>${pairTable(sortBy(P, "primary_r", -1))}</div>
   <div class="card"><h2>Expression across contexts</h2><div class="tablewrap"><table><tr><th>Cell context</th><th>Detection</th><th>Mean log1p CPTT</th><th>Cells</th><th></th></tr>
   ${prof.map(([s, d, mn]) => `<tr><td>${esc(stratumLabel(D.ST[s]))}</td><td class="n">${pct(d)}</td><td class="n">${mn.toFixed(3)}</td><td class="n">${D.ST[s].n_cells ?? "–"}</td><td><span class="bar" style="width:${Math.round(d * 120)}px"></span></td></tr>`).join("")}</table></div></div>`;
}

/* ---------- scientific ---------- */
const COLS = [["tf1", "TF1"], ["tf2", "TF2"], ["ct", "Cell type"], ["tissue", "Tissue"], ["disease", "Disease"], ["assay", "Assay"],
  ["primary_stat", "Primary stat"], ["primary_r", "Primary r"], ["min_q", "q (BH, stratum)"], ["r", "Cell r (raw)"], ["donor_r", "Pseudobulk r"], ["lift", "Lift"],
  ["tf1_det", "TF1 det."], ["tf2_det", "TF2 det."], ["co_detected", "Co-det. cells"], ["n_cells", "Cells"], ["n_datasets", "Datasets"], ["n_collections", "Collections"],
  ["k_meta", "Meta k"], ["meta_r", "Meta r (RE)"], ["meta_p", "Meta p"], ["meta_q", "Meta q"], ["i2", "I²"], ["prior", "Prior"], ["pair_id", "Pair_ID"], ["pmids", "PMIDs"]];
const NUM = {primary_r: 3, min_q: 3, r: 3, donor_r: 3, lift: 2, meta_r: 3, meta_p: 3, meta_q: 3, i2: 2};
function pairTable(P) {
  const head = COLS.map(([k, l]) => `<th data-k="${k}" class="${S.sort.k === k ? "sorted" : ""}" data-dir="${S.sort.dir > 0 ? " ▲" : " ▼"}">${l}</th>`).join("");
  const cell = (p, k) => { const v = p[k];
    if (k === "tf1" || k === "tf2") return `<td><a onclick="go('tf/${esc(v)}')">${esc(v)}</a></td>`;
    if (k === "ct") return `<td title="${esc(p.st.cell_type_id)}">${esc(v)}</td>`;
    if (k === "tf1_det" || k === "tf2_det") return `<td class="n">${pct(v)}</td>`;
    if (k in NUM) return `<td class="n">${fmt(v, NUM[k])}</td>`;
    if (k === "prior") return `<td>${v ? `<span class="pill" title="${esc(p.prior_evidence + " — " + p.prior_classes)}">${esc(p.prior_evidence || "yes")}</span>` : ""}</td>`;
    if (k === "pmids") { const ids = (v || "").split(" | ").filter(Boolean); return `<td>${ids.slice(0, 4).map(x => `<a href="https://pubmed.ncbi.nlm.nih.gov/${encodeURIComponent(x)}/" target="_blank" rel="noopener">${esc(x)}</a>`).join(" ")}${ids.length > 4 ? ` +${ids.length - 4}` : ""}</td>`; }
    return typeof v === "number" ? `<td class="n">${v}</td>` : `<td>${esc(v)}</td>`; };
  const body = P.slice(0, 2000).map(p => `<tr>${COLS.map(([k]) => cell(p, k)).join("")}</tr>`).join("");
  return `<div class="tablewrap"><table><tr>${head}</tr>${body}</table></div>${P.length > 2000 ? `<p class="muted">Showing 2,000 of ${P.length}; narrow the filters or download.</p>` : ""}`;
}
function science() {
  const P = sortBy(filtered(), S.sort.k, S.sort.dir), rc = D.run.summary.run_config, m = D.run.manifest;
  const ctOpts = [...new Set(S.set.map(p => p.s))].map(i => [i, stratumLabel(D.ST[i])]).sort((a, b) => a[1].localeCompare(b[1]));
  $("#app").innerHTML = `<div class="card"><p class="muted">Rows: <b>${esc(S.setLabel)}</b>. Full table: <a href="${base()}${esc(D.run.downloads["cell_type_pair_candidates.tsv"])}" download>download TSV</a>.</p><div class="filters">
   <label>Min primary r<input id="fr" type="number" step="0.05" value="${S.f.r}"></label><label>Max q (stratum)<input id="fq" type="number" step="0.01" value="${S.f.q}"></label>
   <label>Max meta q<input id="fmq" type="number" step="0.01" value="${S.f.mq}"></label><label>Min lift<input id="fl" type="number" step="0.1" value="${S.f.lift}"></label>
   <label>Min collections<input id="fc" type="number" min="1" value="${S.f.coll}"></label>
   <label>Cell context<select id="fct" style="width:260px"><option value="">all</option>${ctOpts.map(([i, l]) => `<option value="${i}" ${S.f.ct === String(i) ? "selected" : ""}>${esc(l)}</option>`).join("")}</select></label>
   <label><span>Prior only</span><input id="fp" type="checkbox" ${S.f.prior ? "checked" : ""} style="width:auto"></label>
   <span class="muted">${P.length} rows · <a id="dlf">download filtered TSV</a></span></div>
   <p class="muted">Screening (per dataset × stratum): primary statistic <b>${esc(rc.primary_stat || "cell_raw (v1)")}</b>${rc.primary_stat && rc.primary_stat !== "cell_raw" ? ` (metacells of ${rc.metacell_size} cells within donor; pseudobulk donors ≥ ${rc.min_donor_cells} cells; residualised on donor + log library size)` : " (unadjusted cell-level Pearson; q shown for reference only)"}; detection ≥ ${rc.min_detection}; primary r ≥ ${rc.min_rho}; ${rc.max_q != null && rc.primary_stat !== "cell_raw" ? `BH q ≤ ${rc.max_q} within stratum; ` : ""}lift ≥ ${rc.min_lift}; co-detected ≥ ${rc.min_codetected}; min cells ${rc.min_cells}${rc.max_strata ? `; max strata/dataset ${rc.max_strata} (pilot cap)` : ""}. Cross-dataset: Fisher-z fixed & DerSimonian–Laird random-effects meta-analysis, BH q within stratum. Census ${esc(m.resolved_census_version)}, TF set ${esc(m.tf_set)}.</p>
   ${pairTable(P)}</div>`;
  const upd = () => { S.f = {r: $("#fr").value, q: $("#fq").value, mq: $("#fmq").value, lift: $("#fl").value, coll: $("#fc").value, prior: $("#fp").checked, ct: $("#fct").value}; science(); };
  ["#fr", "#fq", "#fmq", "#fl", "#fc", "#fct", "#fp"].forEach(s => $(s).addEventListener("change", upd));
  document.querySelectorAll("th[data-k]").forEach(th => th.onclick = () => { const k = th.dataset.k; S.sort = {k, dir: S.sort.k === k ? -S.sort.dir : -1}; science(); });
  $("#dlf").onclick = () => { const cols = COLS.map(c => c[0]).concat(["cell_type_id", "tissue_id", "disease_id", "assay_id", "co_detected", "tau2", "prior_evidence", "prior_classes"]);
    const txt = [cols.join("\t"), ...P.map(p => cols.map(c => c.endsWith("_id") && c !== "pair_id" ? p.st[c] : (p[c] ?? "")).join("\t"))].join("\n");
    const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([txt], {type: "text/tab-separated-values"})); a.download = `tfpairs_${D.runId}_filtered.tsv`; a.click(); };
}

/* ---------- run info ---------- */
function info() {
  const r = D.run, m = r.manifest, s = r.summary;
  const row = (k, v) => `<tr><th>${esc(k)}</th><td style="white-space:normal">${v}</td></tr>`;
  $("#app").innerHTML = `<div class="card"><h2>Run information &amp; provenance</h2><table>
   ${row("Run ID", `<code>${esc(r.run_id)}</code> (schema v${r.schema_version})`)}${row("Status", r.pilot ? "<b>Pilot / partial</b> — not atlas-wide" : "<b>Atlas-wide</b> full run")}${row("Label", esc(r.label || "–"))}
   ${row("Census release", `${esc(m.resolved_census_version)} (requested ${esc(m.requested_census_version)})`)}${row("TF set", esc(m.tf_set))}
   ${row("Compendium", `${m.n_mapped_genes} of ${m.n_register_genes} genes mapped (${m.n_tf_mapped} TFs); register sha256 <code>${esc(m.register_sha256)}</code>`)}
   ${row("Datasets", `${s.n_completed_datasets} completed of ${s.n_planned_datasets} planned (selection ${esc(s.dataset_selection)}; complete ${s.complete}; atlas_wide ${s.atlas_wide})`)}
   ${row("Normalization", esc(m.normalization))}${row("Filters", esc(m.filters))}
   ${row("Run configuration", `<pre>${esc(JSON.stringify(s.run_config, null, 1))}</pre>`)}
   ${row("Positive controls", `${s.positive_controls_recovered ?? "–"} of ${s.positive_controls_tested ?? "–"} recovered`)}
   ${row("Published", `${r.counts.pair_rows} pair×context rows, ${r.counts.unique_pairs} unique pairs, ${r.counts.strata} contexts, ${r.counts.tfs} TFs, ${r.counts.pairs_with_prior} rows with prior`)}
   ${row("Warnings", esc((r.warnings || []).join("; ") || "none"))}
   ${row("Exported", `${esc(r.exported_at)} (code ${esc(r.code_commit_at_export || "–")})`)}
   ${row("Downloads", Object.entries(r.downloads).map(([k, v]) => `<a href="${base()}${esc(v)}" download>${esc(k)}</a>`).join(" · "))}
   ${row("All runs", D.index.runs.map(x => `<a href="?run=${encodeURIComponent(x.run_id)}">${esc(x.run_id)}</a>${x.run_id === D.index.latest ? " (latest)" : ""}`).join("<br>"))}</table></div>
   <div class="card"><h2>How to read this</h2><p>Each row is a TF pair that passed screening inside one cell-type × tissue × disease × assay stratum in at least one Census dataset.
   <b>Primary r</b>: partial correlation used for testing — by default over within-donor <i>metacells</i>, after removing donor means and log library size (falls back to residualised cell-level, then raw cell-level, when too few units). <b>q (BH)</b>: Benjamini–Hochberg over all TF pairs tested in that stratum (both TFs above the detection floor). <b>Cell r (raw)</b>: original v1 cell-level Pearson. <b>Pseudobulk r</b>: partial r over donor means (≥4 donors). <b>Lift</b>: observed/expected co-detection.
   <b>Meta</b>: Fisher-z meta-analysis over datasets that stored the pair (hits plus tested pairs with p ≤ store_p) — k, random-effects r and p, I²; datasets where the pair was not stored are not included, so meta-statistics are optimistic. <b>Collections</b>: independent source collections where the pair passed — the main replication signal. <b>Prior</b>: physical/association evidence in the compendium (Pair_ID).
   These are candidates for SCEA / Tabula Sapiens replication and scATAC evidence levels, not verified cooperative pairs.</p></div>`;
}

$("#q").addEventListener("input", e => { S.q = e.target.value; clearTimeout(S.t); S.t = setTimeout(() => { if (location.hash.startsWith("#tf/")) go(S.mode); else route(); }, 250); });
$("#q").addEventListener("keydown", e => { const t = S.q.trim().toUpperCase(); if (e.key === "Enter" && D.TF && D.TF[t]) go("tf/" + t); });
document.querySelectorAll("nav button").forEach(b => b.onclick = () => go(b.dataset.mode));
boot();
