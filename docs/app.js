"use strict";
/* TF Pair Miner — results explorer. Data: data/runs/index.json + per-run sharded JSON (schema_version 2). */
const D = {cache: {}}; const S = {sort: {k: "min_q", dir: 1}, f: {r: "", q: "", mq: "", lift: "", coll: "", prior: false, rep: false, atac: ""}};
const $ = (s, el = document) => el.querySelector(s);
const esc = s => String(s ?? "").replace(/[&<>"']/g, c => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;"}[c]));
const fmt = (v, d = 3) => v == null || v === "" ? "–" : (typeof v === "number" && !Number.isInteger(v) ? (Math.abs(v) < 1e-3 && v !== 0 ? v.toExponential(1) : v.toFixed(d)) : v);
const pct = v => v == null ? "–" : (v * 100).toFixed(1) + "%";
const base = () => `data/runs/${D.runId}/`;
const J = p => D.cache[p] ||= fetch(p, {cache: "no-cache"}).then(r => { if (!r.ok) throw new Error(p + " " + r.status); return r.json(); });
const rowsOf = o => o.rows.map(r => Object.fromEntries(o.cols.map((c, i) => [c, r[i]])));
const PL = {displaylogo: false, responsive: true};
const sig = p => p.min_q ?? p.meta_q ?? 1;
const rank = (a, b) => sig(a) - sig(b) || (b.primary_r ?? b.r) - (a.primary_r ?? a.r);
const label = st => `${st.cell_type} · ${st.tissue} · ${st.disease} · ${st.assay}`;
const short = st => `${st.cell_type} (${st.tissue}, ${st.assay})`;
function enrich(P) { P.forEach(p => { p.st = D.ST[p.s]; p.ct = p.st.cell_type; p.tissue = p.st.tissue; p.disease = p.st.disease; p.assay = p.st.assay; }); return P; }
const tfFile = t => D.TF[t.toUpperCase()] && base() + `tf/${D.TF[t.toUpperCase()].file}.json`;

async function boot() {
  try { D.index = await J("data/runs/index.json"); } catch (e) { $("#app").innerHTML = `<div class="card">No published runs yet.</div>`; return; }
  const want = new URLSearchParams(location.search).get("run");
  const sel = $("#run"); sel.innerHTML = D.index.runs.map(r => `<option value="${esc(r.run_id)}">${r.run_id === D.index.latest ? "★ " : ""}${esc(r.label || r.run_id)}</option>`).join("");
  sel.onchange = () => { const u = new URL(location); u.searchParams.set("run", sel.value); u.hash = "home"; location.href = u; };
  D.runId = D.index.runs.some(r => r.run_id === want) ? want : D.index.latest; sel.value = D.runId;
  [D.run, D.strata, D.tfs, D.top] = await Promise.all(["run.json", "strata.json", "tfs.json", "top_pairs.json"].map(f => J(base() + f)));
  D.ST = rowsOf(D.strata); D.TF = Object.fromEntries(rowsOf(D.tfs).map(t => [t.tf.toUpperCase(), t]));
  D.TOP = enrich(rowsOf(D.top));
  const sm = D.run.summary;
  if (D.run.pilot) { const b = $("#banner"); b.hidden = false; b.innerHTML = `<b>Pilot data, not an atlas-wide result.</b> ${esc(D.run.label || "")}: ${sm.n_completed_datasets}/${sm.n_planned_datasets} dataset(s), up to ${sm.run_config.max_strata || "∞"} strata per dataset. Use the run selector to switch runs. The full run replaces this one when it's published.`; }
  setupSearch(); window.addEventListener("hashchange", route); route();
}
function layerBadges(compact) {
  const L = D.run.layers || {}; const emb = D.run.embedding || {};
  const b = (st, t) => `<span class="badge ${st === "done" ? "b-done" : "b-todo"}">${esc(t)}: ${st === "done" ? "done" : "not yet run"}</span>`;
  return `<div class="layers">${b((L.discovery_cellxgene || {}).status || "done", "CELLxGENE discovery")}${b(emb.status || "not_run", "Per-cell UMAP")}${b((L.replication || {}).status || "not_run", "Replication (SCEA / Tabula Sapiens)")}${b((L.scatac || {}).status || "not_run", (L.scatac || {}).status === "done" ? "scATAC level B (A, C, D not run)" : "scATAC evidence A–D")}</div>`;
}
const REP = {replicated_donor_supported: ["replicated", "b-rep"], cell_level_support_only: ["cell-level only", "b-cell"], not_replicated: ["not replicated", "b-no"],
  inconclusive_low_expression: ["inconclusive", "b-inc"], insufficient_cells: ["insufficient cells", "b-inc"], gene_unavailable: ["gene unavailable", "b-inc"],
  no_cell_type_match: ["no TS cell type", "b-todo"], no_matching_context: ["no TS context", "b-todo"], no_allowed_assay: ["no TS assay", "b-todo"]};
const repBadge = (st, title) => { const [l, c] = REP[st] || [st, "b-todo"]; return `<span class="badge ${c}" title="${esc(title || st)}">${esc(l)}</span>`; };
const repLayerDone = () => ((D.run.layers || {}).replication || {}).status === "done";
const repCell = p => p.replication == null ? (repLayerDone() ? `<span class="badge b-todo">not tested</span>` : `<span class="badge b-todo">not yet run</span>`)
  : repBadge(p.replication.status, `${p.replication.assay || ""} ${p.replication.reason || ""}`);
const isRep = p => p.replication && p.replication.status === "replicated_donor_supported";
const ATACC = {catlas_coenriched: "b-rep", not_coenriched: "b-no", provisional_coenriched: "b-cell", provisional_not_coenriched: "b-no",
  proxy_coenriched: "b-proxy", proxy_not_coenriched: "b-proxyno", no_matching_context: "b-todo", no_exact_celltype: "b-todo", not_tested: "b-todo"};
const atacDone = () => ((D.run.layers || {}).scatac || {}).status === "done";
const atacBadge = (st, label, tip) => `<span class="badge ${ATACC[st] || "b-inc"}" title="${esc(tip || st)}">${esc(label || st)}</span>`;
const atacCell = p => { const a = p.atac_level;
  if (a == null) return `<span class="badge b-todo">${atacDone() ? "not tested" : "not yet run"}</span>`;
  return `<span class="muted" style="font-size:10px">B</span> ${atacBadge(a.status, a.label, "Level B: same-cCRE motif co-enrichment (not co-binding)")}${a.proxy_status ? " " + atacBadge(a.proxy_status, "fetal proxy: " + (a.proxy_label || "").replace(" in fetal proxy", ""), "Exploratory FETAL-proxy result, not adult validation") : ""}`; };
const atacKey = p => !p.atac_level ? "none" : p.atac_level.status;

/* ---------- search with autocomplete ---------- */
function setupSearch() {
  const q = $("#q"), ac = $("#ac"); let items = [], sel = -1;
  const ctNames = [...new Set(D.ST.map(s => s.cell_type))];
  const render = () => { ac.innerHTML = items.map((it, i) => `<div class="${i === sel ? "sel" : ""}" data-i="${i}"><span>${esc(it.t)}</span><small>${esc(it.k)}</small></div>`).join(""); ac.hidden = !items.length;
    ac.querySelectorAll("div").forEach(d => d.onmousedown = e => { e.preventDefault(); pick(items[+d.dataset.i]); }); };
  const pick = it => { ac.hidden = true; q.value = it.t; go(it.h); };
  q.addEventListener("input", () => {
    const v = q.value.trim(); sel = -1; if (!v) { items = []; return render(); }
    const parts = v.split(/[:\s,/+]+/).filter(Boolean).map(x => x.toUpperCase());
    items = [];
    if (parts.length === 2 && D.TF[parts[0]]) Object.keys(D.TF).filter(t => t.startsWith(parts[1]) && t !== parts[0]).slice(0, 6).forEach(t => items.push({t: `${D.TF[parts[0]].tf}:${D.TF[t].tf}`, k: "pair", h: `pair/${D.TF[parts[0]].tf}:${D.TF[t].tf}`}));
    const V = v.toUpperCase();
    Object.values(D.TF).filter(t => t.tf.toUpperCase().startsWith(V)).sort((a, b) => b.n_pair_rows - a.n_pair_rows || a.tf.localeCompare(b.tf)).slice(0, 8)
      .forEach(t => items.push({t: t.tf, k: `TF · ${t.n_partners} partners`, h: `tf/${t.tf}`}));
    D.ST.forEach((s, i) => { if (label(s).toLowerCase().includes(v.toLowerCase()) && items.length < 16) items.push({t: short(s), k: `cell type · ${s.n_pair_rows} pairs`, h: `ct/${i}`}); });
    render(); });
  q.addEventListener("keydown", e => { if (ac.hidden) return; if (e.key === "ArrowDown") { sel = Math.min(sel + 1, items.length - 1); render(); e.preventDefault(); }
    else if (e.key === "ArrowUp") { sel = Math.max(sel - 1, 0); render(); e.preventDefault(); } else if (e.key === "Enter" && items.length) pick(items[Math.max(sel, 0)]); else if (e.key === "Escape") ac.hidden = true; });
  q.addEventListener("blur", () => setTimeout(() => ac.hidden = true, 150));
}
function go(h) { if (location.hash === "#" + h) route(); else location.hash = h; }
async function route() {
  const h = decodeURIComponent(location.hash.slice(1) || "home"); const i = h.indexOf("/"); const m = i < 0 ? h : h.slice(0, i), arg = i < 0 ? "" : h.slice(i + 1);
  document.querySelectorAll("nav button").forEach(b => b.classList.toggle("on", b.dataset.mode === m));
  window.scrollTo(0, 0);
  try { await ({home, table, info, about, tf: () => tfPage(arg), ct: () => ctPage(+arg), pair: () => pairPage(arg)}[m === "table" && arg === "replicated" ? (S.f.rep = true, "table") : m] || home)(); }
  catch (e) { $("#app").innerHTML = `<div class="card">Error: ${esc(e.message)}</div>`; console.error(e); }
}
document.querySelectorAll("nav button").forEach(b => b.onclick = () => go(b.dataset.mode));
const stat = (v, l) => `<div class="stat"><b>${esc(v)}</b>${esc(l)}</div>`;

/* ---------- shared widgets ---------- */
function network(el, P, focus) {
  const E = new Map();
  for (const p of P) { const k = p.tf1 + ":" + p.tf2, r = p.primary_r ?? p.r; const e = E.get(k) || {a: p.tf1, b: p.tf2, r: -1, prior: 0, n: 0, q: 1}; e.r = Math.max(e.r, r); e.q = Math.min(e.q, sig(p)); e.prior |= p.prior; e.n++; E.set(k, e); }
  const edges = [...E.values()].sort((x, y) => x.q - y.q || y.r - x.r).slice(0, 120);
  const W = 640, H = 420, ids = [...new Set(edges.flatMap(e => [e.a, e.b]))];
  if (!ids.length) { el.innerHTML = '<p class="muted">No pairs.</p>'; return; }
  const N = ids.map((id, i) => ({id, x: W / 2 + 180 * Math.cos(i * 2.4) * Math.sqrt((i + 1) / ids.length), y: H / 2 + 160 * Math.sin(i * 2.4) * Math.sqrt((i + 1) / ids.length), vx: 0, vy: 0, deg: 0}));
  const ix = Object.fromEntries(N.map((n, i) => [n.id, i])); edges.forEach(e => { N[ix[e.a]].deg++; N[ix[e.b]].deg++; });
  for (let it = 0; it < 300; it++) { const t = 1 - it / 300;
    for (let i = 0; i < N.length; i++) for (let j = i + 1; j < N.length; j++) { const a = N[i], b = N[j], dx = a.x - b.x, dy = a.y - b.y, d2 = dx * dx + dy * dy + .01, f = 900 / d2; a.vx += dx * f; a.vy += dy * f; b.vx -= dx * f; b.vy -= dy * f; }
    for (const e of edges) { const a = N[ix[e.a]], b = N[ix[e.b]], dx = b.x - a.x, dy = b.y - a.y, d = Math.hypot(dx, dy) + .01, f = (d - 60) * .02; a.vx += dx / d * f; a.vy += dy / d * f; b.vx -= dx / d * f; b.vy -= dy / d * f; }
    for (const n of N) { n.vx += (W / 2 - n.x) * .006; n.vy += (H / 2 - n.y) * .006; n.x += Math.max(-10, Math.min(10, n.vx)) * t; n.y += Math.max(-10, Math.min(10, n.vy)) * t; n.vx *= .5; n.vy *= .5; n.x = Math.max(20, Math.min(W - 45, n.x)); n.y = Math.max(12, Math.min(H - 12, n.y)); } }
  el.innerHTML = `<svg viewBox="0 0 ${W} ${H}" width="100%">${edges.map(e => { const a = N[ix[e.a]], b = N[ix[e.b]];
    return `<line onclick="go('pair/${esc(e.a)}:${esc(e.b)}')" style="cursor:pointer" x1="${a.x}" y1="${a.y}" x2="${b.x}" y2="${b.y}" stroke="${e.prior ? "#16a34a" : "#94a3b8"}" stroke-opacity=".75" stroke-width="${1 + 5 * Math.max(0, e.r)}"><title>${esc(e.a)}–${esc(e.b)}: best q ${fmt(e.q)}, max r ${e.r.toFixed(3)}${e.prior ? ", literature prior" : ""} (click for pair page)</title></line>`; }).join("")}
    ${N.map(n => `<g class="node" onclick="go('tf/${esc(n.id)}')" style="cursor:pointer"><circle cx="${n.x}" cy="${n.y}" r="${3 + Math.sqrt(n.deg) * 2.2}" fill="${n.id === focus ? "#dc2626" : "#2563eb"}"/><text x="${n.x + 6}" y="${n.y + 3}">${esc(n.id)}</text></g>`).join("")}</svg>
    <p class="muted">Width = r; green = literature prior. Click an edge to open the pair, or a node to open the TF.</p>`;
}
const COLS = [["tf1", "TF1"], ["tf2", "TF2"], ["ct", "Cell type"], ["tissue", "Tissue"], ["disease", "Disease"], ["assay", "Assay"], ["primary_stat", "Primary stat"], ["primary_r", "Primary r"], ["min_q", "q (BH)"],
  ["r", "Cell r (raw)"], ["donor_r", "Pseudobulk r"], ["lift", "Lift"], ["tf1_det", "TF1 det."], ["tf2_det", "TF2 det."], ["co_detected", "Co-det. cells"], ["n_cells", "Cells"], ["n_datasets", "Datasets"], ["n_collections", "Collections"],
  ["k_meta", "Meta k"], ["meta_r", "Meta r"], ["meta_q", "Meta q"], ["i2", "I²"], ["prior", "Prior"], ["pair_id", "Pair_ID"], ["pmids", "PMIDs"], ["replication", "Replication"], ["atac_level", "scATAC"]];
const NUM = {primary_r: 3, min_q: 3, r: 3, donor_r: 3, lift: 2, meta_r: 3, meta_q: 3, i2: 2};
function pairTable(P, opts = {}) {
  const cols = opts.cols ? COLS.filter(c => opts.cols.includes(c[0])) : COLS;
  const head = cols.map(([k, l]) => `<th data-k="${k}" class="${S.sort.k === k ? "sorted" : ""}" data-dir="${S.sort.dir > 0 ? " ▲" : " ▼"}">${l}</th>`).join("");
  const cell = (p, k) => { const v = p[k];
    if (k === "tf1" || k === "tf2") return `<td><a onclick="go('tf/${esc(v)}')">${esc(v)}</a> <a class="muted" title="pair page" onclick="go('pair/${esc(p.tf1)}:${esc(p.tf2)}')">⇄</a></td>`;
    if (k === "ct") return `<td><a onclick="go('ct/${p.s}')" title="${esc(p.st.cell_type_id)}">${esc(v)}</a></td>`;
    if (k === "tf1_det" || k === "tf2_det") return `<td class="n">${pct(v)}</td>`;
    if (k in NUM) return `<td class="n">${fmt(v, NUM[k])}</td>`;
    if (k === "prior") return `<td>${v ? `<span class="pill" title="${esc(p.prior_evidence + " — " + p.prior_classes)}">${esc(p.prior_evidence || "yes")}</span>` : ""}</td>`;
    if (k === "pmids") { const ids = (v || "").split(" | ").filter(Boolean); return `<td>${ids.slice(0, 3).map(x => `<a href="https://pubmed.ncbi.nlm.nih.gov/${encodeURIComponent(x)}/" target="_blank" rel="noopener">${esc(x)}</a>`).join(" ")}${ids.length > 3 ? ` +${ids.length - 3}` : ""}</td>`; }
    if (k === "replication") return `<td>${repCell(p)}</td>`; if (k === "atac_level") return `<td>${atacCell(p)}</td>`;
    return typeof v === "number" ? `<td class="n">${v}</td>` : `<td>${esc(v)}</td>`; };
  const lim = opts.limit || 2000;
  return `<div class="tablewrap"><table><tr>${head}</tr>${P.slice(0, lim).map(p => `<tr>${cols.map(([k]) => cell(p, k)).join("")}</tr>`).join("")}</table></div>${P.length > lim ? `<p class="muted">Showing ${lim} of ${P.length}.</p>` : ""}`;
}
function sortBy(a, k, dir) { return a.sort((x, y) => { const u = x[k] ?? (dir > 0 ? Infinity : -Infinity), v = y[k] ?? (dir > 0 ? Infinity : -Infinity); return u > v ? dir : u < v ? -dir : rank(x, y); }); }

/* ---------- embedding helpers ---------- */
async function embedFor(genes, sPref) {
  const E = (D.run.embedding || {}).datasets || []; const out = [];
  for (const d of E) { if (!genes.every(g => d.genes.includes(g))) continue;
    const cells = await J(base() + `embed/${d.dataset_id}/cells.json`); const nS = sPref == null ? 0 : cells.rows.filter(r => r[3] === sPref).length; out.push({d, cells, nS}); }
  out.sort((a, b) => b.nS - a.nS); return out;
}
async function geneVec(did, g, n) { const j = await J(base() + `embed/${did}/genes/${D.TF[g.toUpperCase()] ? D.TF[g.toUpperCase()].file : g}.json`); const v = new Float32Array(n); j.i.forEach((i, k) => v[i] = j.v[k]); return v; }
function umap(el, cells, values, title, cs = "Viridis", hl) {
  const x = cells.rows.map(r => r[0]), y = cells.rows.map(r => r[1]);
  const text = cells.rows.map((r, i) => `${cells.cell_types[r[2]]}${values ? "<br>" + values[i].toFixed(2) : ""}`);
  const ord = [...x.keys()].sort((a, b) => (values ? values[a] - values[b] : 0));
  const tr = values ? {x: ord.map(i => x[i]), y: ord.map(i => y[i]), text: ord.map(i => text[i]), mode: "markers", type: "scattergl", hoverinfo: "text",
    marker: {size: 3, color: ord.map(i => values[i]), colorscale: cs, showscale: true, colorbar: {thickness: 8, len: .6}}} : null;
  const data = values ? [tr] : cells.cell_types.map((c, k) => { const idx = cells.rows.map((r, i) => r[2] === k ? i : -1).filter(i => i >= 0);
    return {name: c.length > 32 ? c.slice(0, 30) + "…" : c, x: idx.map(i => x[i]), y: idx.map(i => y[i]), mode: "markers", type: "scattergl", marker: {size: 3, opacity: hl == null ? .8 : .25}, hovertext: idx.map(i => c), hoverinfo: "text"}; });
  if (hl != null) { const idx = cells.rows.map((r, i) => r[3] === hl ? i : -1).filter(i => i >= 0); data.push({name: "this stratum", x: idx.map(i => x[i]), y: idx.map(i => y[i]), mode: "markers", type: "scattergl", marker: {size: 4, color: "#dc2626"}}); }
  Plotly.newPlot(el, data, {title: {text: title, font: {size: 13}}, margin: {l: 10, r: 10, t: 30, b: 10}, xaxis: {visible: false}, yaxis: {visible: false}, showlegend: !values, legend: {font: {size: 9}}, height: 340}, PL);
}
function embedNote(e) { return `<p class="muted">UMAP from <b>${esc(e.d.title || e.d.dataset_id)}</b>: ${e.cells.n_cells} of ${e.cells.n_primary_cells} primary cells, downsampled by cell type. ${esc(e.cells.method)}. Colour shows log1p CPTT.</p>`; }

/* ---------- pages ---------- */
async function home() {
  const c = D.run.counts, P = D.TOP.slice().sort(rank);
  const tfDeg = {}; P.forEach(p => { tfDeg[p.tf1] = (tfDeg[p.tf1] || 0) + 1; tfDeg[p.tf2] = (tfDeg[p.tf2] || 0) + 1; });
  const topTF = Object.keys(tfDeg).sort((a, b) => tfDeg[b] - tfDeg[a]).slice(0, 12);
  $("#app").innerHTML = `<div class="card"><div class="stats">${stat(c.unique_pairs, "unique TF pairs")}${stat(c.pair_rows, "pair × cell-type rows")}${stat(c.strata, "cell contexts")}${stat(c.pairs_with_prior, "with literature prior")}${stat(D.run.manifest.resolved_census_version, "Census release")}</div>${layerBadges()}
   <p>Search above, or start from a TF: <span class="chips">${topTF.map(t => `<a onclick="go('tf/${t}')">${t}</a>`).join("")}</span></p>
   <p>Or a cell context: <span class="chips">${D.ST.map((s, i) => [s, i]).filter(([s]) => s.n_pair_rows).map(([s, i]) => `<a onclick="go('ct/${i}')">${esc(short(s))} (${s.n_pair_rows})</a>`).join("") || '<span class="muted">no cell context has passing pairs in this run</span>'}</span></p></div>
   <div class="grid"><div class="card"><h2>Pair network (most significant)</h2><div id="net"></div></div>
   <div class="card"><h2>Positive controls</h2><p class="muted">Known partner pairs checked in this run.</p><div class="chips">${(D.run.positive_controls || []).map(x => `<a ${+x.n_contexts_with_hit ? `style="background:#dcfce7" onclick="go('pair/${x.pair.replace("-", ":")}')"` : ""} title="best r ${x.best_primary_r || "–"}, min q ${x.min_q || "–"}">${esc(x.pair)} ${+x.n_contexts_with_hit ? "✓" : "–"}</a>`).join("")}</div>
   <h2 style="margin-top:14px">Embeddings available</h2>${((D.run.embedding || {}).datasets || []).map(d => `<div class="muted">• ${esc(d.title || d.dataset_id)}: ${d.n_cells} cells, ${d.genes.length} TFs</div>`).join("") || '<p class="muted">None for this run.</p>'}</div></div>
   <div class="card"><h2>Top candidates <span class="muted">(ranked by q, then r)</span></h2>${pairTable(P, {limit: 30, cols: ["tf1", "tf2", "ct", "assay", "primary_r", "min_q", "lift", "tf1_det", "tf2_det", "co_detected", "n_collections", "prior", "replication", "atac_level"]})}</div>`;
  network($("#net"), P);
}

async function tfPage(tf) {
  const t = D.TF[tf.toUpperCase()]; if (!t) { $("#app").innerHTML = `<div class="card">${esc(tf)} isn't in this run.</div>`; return; }
  tf = t.tf; const sh = await J(tfFile(tf)); const P = enrich(rowsOf(sh)).sort(rank);
  const partners = {}; P.forEach(p => { const o = p.tf1 === tf ? p.tf2 : p.tf1; (partners[o] ||= []).push(p); });
  const plist = Object.entries(partners).map(([o, ps]) => ({o, ps, best: ps.slice().sort(rank)[0]})).sort((a, b) => rank(a.best, b.best));
  $("#app").innerHTML = `<div class="crumb"><a onclick="go('home')">Explore</a> › TF</div><div class="card"><h2>${esc(tf)}</h2><div class="stats">${stat(plist.length, "candidate partners")}${stat(P.length, "pair × cell-type rows")}${stat(sh.profile.length, "contexts expressed (≥" + D.run.profile_detection_floor * 100 + "%)")}</div>${layerBadges()}</div>
   <div class="grid"><div class="card"><h2>Partners ranked by significance</h2><div class="tablewrap" style="max-height:420px"><table><tr><th>#</th><th>Partner</th><th>Best q</th><th>Best r</th><th>Cell types</th><th>Prior</th><th>Replication</th><th>scATAC</th></tr>
   ${plist.map((x, i) => `<tr><td>${i + 1}</td><td><a onclick="go('pair/${esc(tf)}:${esc(x.o)}')">${esc(x.o)}</a></td><td class="n">${fmt(sig(x.best))}</td><td class="n">${fmt(x.best.primary_r ?? x.best.r)}</td><td>${[...new Set(x.ps.map(p => p.ct))].map(esc).join(", ")}</td><td>${x.ps.some(p => p.prior) ? '<span class="pill">yes</span>' : ""}</td><td>${repCell(x.best)}</td><td>${atacCell(x.best)}</td></tr>`).join("") || '<tr><td colspan="8" class="muted">No passing partners in this run</td></tr>'}</table></div></div>
   <div class="card"><h2>Partner × cell-type dot plot</h2><div id="dot" class="plot"></div></div></div>
   <div class="grid"><div class="card"><h2>${esc(tf)} detection by cell context</h2><div id="bar" class="plot"></div></div><div class="card"><h2>${esc(tf)} on UMAP</h2><div id="um" class="plot"></div><div id="umn"></div></div></div>
   <div class="card"><h2>All rows</h2>${pairTable(P)}</div>`;
  if (P.length) { const y = P.map(p => short(p.st)); Plotly.newPlot("dot", [{type: "scatter", mode: "markers", x: P.map(p => p.tf1 === tf ? p.tf2 : p.tf1), y,
      marker: {size: P.map(p => 6 + Math.min(24, -Math.log10(Math.max(sig(p), 1e-300)) * 1.5)), color: P.map(p => p.primary_r ?? p.r), colorscale: "Reds", showscale: true, colorbar: {title: "r", thickness: 8}, line: {width: P.map(p => p.prior ? 2 : 0), color: "#16a34a"}},
      text: P.map(p => `${p.tf1}:${p.tf2}<br>${short(p.st)}<br>r ${fmt(p.primary_r ?? p.r)}, q ${fmt(sig(p))}`), hoverinfo: "text"}],
      {margin: {l: 200, r: 10, t: 10, b: 80}, height: Math.max(300, 120 + 40 * new Set(y).size), xaxis: {tickangle: -45, title: "partner"}, yaxis: {automargin: true}}, PL);
    document.getElementById("dot").on("plotly_click", e => go(`pair/${tf}:${e.points[0].x}`)); }
  else $("#dot").innerHTML = '<p class="muted">No pairs.</p>';
  const prof = sh.profile.slice().sort((a, b) => b[1] - a[1]).slice(0, 40);
  Plotly.newPlot("bar", [{type: "bar", orientation: "h", y: prof.map(r => short(D.ST[r[0]])), x: prof.map(r => r[1] * 100), marker: {color: "#2563eb"}, text: prof.map(r => `mean log1p ${r[2].toFixed(2)}`), hovertemplate: "%{y}<br>%{x:.1f}% cells<br>%{text}<extra></extra>"}],
    {margin: {l: 220, r: 10, t: 10, b: 40}, height: Math.max(260, 60 + 24 * prof.length), xaxis: {title: "% cells detecting"}, yaxis: {autorange: "reversed", automargin: true}}, PL);
  const E = await embedFor([tf]); if (E.length) { const e = E[0]; umap("um", e.cells, await geneVec(e.d.dataset_id, tf, e.cells.rows.length), `${tf} (log1p CPTT)`); $("#umn").innerHTML = embedNote(e); }
  else $("#um").innerHTML = '<p class="muted">No per-cell embedding layer for this TF in this run.</p>';
  bindSort(P, "tfPage", tf);
}

async function ctPage(i) {
  const st = D.ST[i]; if (!st) { $("#app").innerHTML = '<div class="card">Unknown cell context.</div>'; return; }
  const sh = await J(base() + `ct/${i}.json`); const P = enrich(rowsOf(sh)).sort(rank);
  $("#app").innerHTML = `<div class="crumb"><a onclick="go('home')">Explore</a> › Cell type</div><div class="card"><h2>${esc(st.cell_type)}</h2><div class="kv"><span><b>Tissue</b> ${esc(st.tissue)}</span><span><b>Disease</b> ${esc(st.disease)}</span><span><b>Assay</b> ${esc(st.assay)}</span><span><b>Ontology</b> ${esc(st.cell_type_id)}</span><span><b>Cells</b> ${st.n_cells ?? "–"}</span></div>
   <div class="stats" style="margin-top:8px">${stat(P.length, "candidate pairs")}${stat(new Set(P.flatMap(p => [p.tf1, p.tf2])).size, "TFs involved")}${stat(P.filter(p => p.prior).length, "with literature prior")}</div>${layerBadges()}</div>
   <div class="grid"><div class="card"><h2>Pair network</h2><div id="net"></div></div><div class="card"><h2>Most detected TFs here</h2><div id="bar" class="plot"></div></div></div>
   <div class="card"><h2>Where these cells sit</h2><div id="um" class="plot"></div><div id="umn"></div></div>
   <div class="card"><h2>Pairs ranked (q, then r)</h2>${pairTable(P)}</div>`;
  network($("#net"), P);
  const prof = sh.profile.slice(0, 30);
  Plotly.newPlot("bar", [{type: "bar", orientation: "h", y: prof.map(r => r[0]), x: prof.map(r => r[1] * 100), marker: {color: "#7c3aed"}, hovertemplate: "%{y}: %{x:.1f}%<extra></extra>"}],
    {margin: {l: 80, r: 10, t: 10, b: 40}, height: 60 + 20 * prof.length, xaxis: {title: "% cells detecting"}, yaxis: {autorange: "reversed"}}, PL);
  document.getElementById("bar").on("plotly_click", e => go("tf/" + e.points[0].y));
  const E = await embedFor([], i); if (E.length && E[0].nS) { umap("um", E[0].cells, null, "Cell types (red = this stratum)", null, i); $("#umn").innerHTML = embedNote(E[0]); }
  else $("#um").innerHTML = '<p class="muted">No per-cell embedding contains this stratum.</p>';
  bindSort(P, "ctPage", i);
}

async function pairPage(arg) {
  const [a0, b0] = arg.split(":"); const A = D.TF[(a0 || "").toUpperCase()], B = D.TF[(b0 || "").toUpperCase()];
  if (!A || !B) { $("#app").innerHTML = `<div class="card">Pair ${esc(arg)} isn't in this run.</div>`; return; }
  const [sa, sb] = await Promise.all([J(tfFile(A.tf)), J(tfFile(B.tf))]);
  const P = enrich(rowsOf(sa)).filter(p => (p.tf1 === A.tf && p.tf2 === B.tf) || (p.tf1 === B.tf && p.tf2 === A.tf)).sort(rank);
  const best = P[0];
  $("#app").innerHTML = `<div class="crumb"><a onclick="go('home')">Explore</a> › <a onclick="go('tf/${A.tf}')">${A.tf}</a> › pair</div>
   <div class="card"><h2>${esc(A.tf)} : ${esc(B.tf)}</h2>${best ? `<div class="stats">${stat(P.length, "cell contexts with a hit")}${stat(fmt(sig(best)), "best q")}${stat(fmt(best.primary_r ?? best.r), "best r")}${stat(Math.max(...P.map(p => p.n_collections)), "max collections")}${stat(best.prior ? "yes" : "no", "literature prior")}</div>` : '<p class="muted">This pair did not pass screening in this run. The expression panels below are still real data.</p>'}
   <h3>Evidence layers</h3><table><tr><th>Layer</th><th>Status</th><th>Evidence</th></tr>
   <tr><td>CELLxGENE discovery</td><td><span class="badge b-done">done</span></td><td>${best ? `${P.length} context(s); best ${esc(best.primary_stat)} r = ${fmt(best.primary_r)}, q = ${fmt(best.min_q)}; meta k = ${best.k_meta ?? 0}` : "no passing context"}</td></tr>
   <tr><td>Literature / physical prior</td><td><span class="badge ${best && best.prior ? "b-done" : "b-todo"}">${best && best.prior ? "recorded" : "none recorded"}</span></td><td>${best && best.prior ? `${esc(best.prior_evidence)}; ${esc(best.prior_classes)}; Pair_ID ${esc(best.pair_id || "–")}; PMIDs ${(best.pmids || "").split(" | ").filter(Boolean).slice(0, 8).map(x => `<a href="https://pubmed.ncbi.nlm.nih.gov/${encodeURIComponent(x)}/" target="_blank" rel="noopener">${esc(x)}</a>`).join(" ")}` : "Not recorded in the compendium. That doesn't mean the pair doesn't interact."}</td></tr>
   ${repRow(P)}
   ${atacRow(P)}</table></div>
   ${atacPanel(P)}
   ${repTable(P)}
   <div class="card"><h2>All evidence rows</h2>${P.length ? pairTable(P) : '<p class="muted">none</p>'}</div>
   <div class="card"><h2>Detection by cell context</h2><div id="det" class="plot"></div></div>
   <div class="card"><h2>Per-cell expression</h2><div id="embsel"></div><div class="grid3"><div id="uA" class="plot"></div><div id="uB" class="plot"></div><div id="uC" class="plot"></div></div><div id="umn"></div>
   <div class="grid"><div id="vio" class="plot"></div><div id="sc" class="plot"></div></div></div>`;
  const pa = new Map(sa.profile.map(r => [r[0], r[1]])), pb = new Map(sb.profile.map(r => [r[0], r[1]]));
  const ss = [...new Set([...pa.keys(), ...pb.keys()])].sort((x, y) => ((pa.get(y) || 0) + (pb.get(y) || 0)) - ((pa.get(x) || 0) + (pb.get(x) || 0))).slice(0, 30);
  const hit = new Set(P.map(p => p.s)); const yl = ss.map(s => (hit.has(s) ? "★ " : "") + short(D.ST[s]));
  Plotly.newPlot("det", [{type: "bar", orientation: "h", name: A.tf, y: yl, x: ss.map(s => (pa.get(s) || 0) * 100)}, {type: "bar", orientation: "h", name: B.tf, y: yl, x: ss.map(s => (pb.get(s) || 0) * 100)}],
    {barmode: "group", margin: {l: 240, r: 10, t: 10, b: 40}, height: 80 + 34 * ss.length, xaxis: {title: `% cells detecting (≥${D.run.profile_detection_floor * 100}% shown; ★ = pair hit)`}, yaxis: {autorange: "reversed", automargin: true}}, PL);
  const E = await embedFor([A.tf, B.tf], best ? best.s : null);
  if (!E.length) { $("#uA").innerHTML = '<p class="muted">No per-cell embedding with both TFs in this run.</p>'; return; }
  $("#embsel").innerHTML = E.length > 1 ? `<label class="muted">Dataset <select id="es">${E.map((e, k) => `<option value="${k}">${esc(e.d.title || e.d.dataset_id)} (${e.nS} cells in hit stratum)</option>`).join("")}</select></label>` : "";
  const draw = async k => { const e = E[k], n = e.cells.rows.length; const [va, vb] = await Promise.all([geneVec(e.d.dataset_id, A.tf, n), geneVec(e.d.dataset_id, B.tf, n)]);
    const co = va.map((v, i) => Math.sqrt(v * vb[i]));
    umap("uA", e.cells, va, `${A.tf}`); umap("uB", e.cells, vb, `${B.tf}`); umap("uC", e.cells, co, `co-expression √(${A.tf}·${B.tf})`, "Hot");
    $("#umn").innerHTML = embedNote(e);
    const cts = e.cells.cell_types; const ctOf = e.cells.rows.map(r => cts[r[2]]); const top = cts.map((c, k2) => [c, ctOf.filter(x => x === c).length]).sort((x, y) => y[1] - x[1]).slice(0, 12).map(x => x[0]);
    const vt = (v, nm, side) => ({type: "violin", name: nm, legendgroup: nm, scalegroup: nm, side, x: ctOf.filter(c => top.includes(c)), y: [...v].filter((_, i) => top.includes(ctOf[i])), points: false, spanmode: "hard", meanline: {visible: true}});
    Plotly.newPlot("vio", [vt(va, A.tf, "negative"), vt(vb, B.tf, "positive")], {title: {text: "Expression by cell type (log1p CPTT)", font: {size: 13}}, violinmode: "overlay", violingap: 0, margin: {l: 40, r: 10, t: 30, b: 140}, height: 420, xaxis: {tickangle: -40}}, PL);
    const inS = e.cells.rows.map(r => best && r[3] === best.s);
    const jit = () => (Math.random() - .5) * .08;
    Plotly.newPlot("sc", [{name: "other cells", x: [...va].filter((_, i) => !inS[i]).map(v => v + jit()), y: [...vb].filter((_, i) => !inS[i]).map(v => v + jit()), mode: "markers", type: "scattergl", marker: {size: 3, color: "#94a3b8", opacity: .5}},
      {name: best ? "hit stratum" : "", x: [...va].filter((_, i) => inS[i]).map(v => v + jit()), y: [...vb].filter((_, i) => inS[i]).map(v => v + jit()), mode: "markers", type: "scattergl", marker: {size: 4, color: "#dc2626"}}],
      {title: {text: `${A.tf} vs ${B.tf} per cell (jittered)`, font: {size: 13}}, xaxis: {title: A.tf}, yaxis: {title: B.tf}, margin: {l: 50, r: 10, t: 30, b: 40}, height: 420}, PL); };
  await draw(0); if ($("#es")) $("#es").onchange = e => draw(+e.target.value);
}

function repRow(P) {
  const L = (D.run.layers || {}).replication || {};
  if (L.status !== "done") return `<tr><td>Replication (SCEA / Tabula Sapiens)</td><td><span class="badge b-todo">not yet run</span></td><td class="muted">This layer will be filled in by a later pipeline stage</td></tr>`;
  const R = P.filter(p => p.replication); const best = R.sort((a, b) => (isRep(b) - isRep(a)))[0];
  return `<tr><td>Replication (${esc(L.dataset || "Tabula Sapiens")})</td><td>${best ? repBadge(best.replication.status) : '<span class="badge b-todo">not tested</span>'}</td><td>${best ? `${R.filter(isRep).length} of ${R.length} discovery context(s) replicated; best: ${esc(best.replication.assay || "")}, ${best.replication.n_cells ?? "–"} cells, ${best.replication.n_donors ?? "–"} donors, r = ${fmt(best.replication.r)}, lift = ${fmt(best.replication.lift)}, positive donors ${fmt(best.replication.pos_donor_frac)}` : "No candidate row for this pair was sent to replication"}</td></tr>`;
}
function repTable(P) {
  const L = (D.run.layers || {}).replication || {}; if (L.status !== "done") return "";
  const R = P.filter(p => p.replication); if (!R.length) return "";
  const t = L.thresholds || {};
  return `<div class="card"><h2>Replication in ${esc(L.dataset || "Tabula Sapiens")} by assay</h2><p class="muted">Same cell type × tissue × disease; independent donors and datasets (discovery datasets excluded). Thresholds: ≥${t.min_cells} cells, detection ≥${t.min_detection}, r ≥${t.min_rho}, lift ≥${t.min_lift}, co-detected ≥${t.min_codetected}; donor support = ≥${t.min_donors} donors with ≥${t.min_donor_cells} cells and ≥${t.min_positive_donor_fraction} of donors positive. ${esc(L.independence_note || "")}</p>
   <table><tr><th>Discovery context</th><th>TS assay</th><th>Status</th><th>Cells</th><th>Donors</th><th>r</th><th>Lift</th><th>Co-detected</th><th>Det A</th><th>Det B</th><th>Donor mean r</th><th>Donors tested</th><th>Positive donors</th><th>Reason</th></tr>
   ${R.flatMap(p => (p.replication.assays.length ? p.replication.assays : [p.replication]).map(a => `<tr><td>${esc(p.ct)} · ${esc(p.assay)}</td><td>${esc(a.assay || "–")}</td><td>${repBadge(a.status)}</td><td class="n">${a.n_cells ?? "–"}</td><td class="n">${a.n_donors ?? "–"}</td><td class="n">${fmt(a.r)}</td><td class="n">${fmt(a.lift)}</td><td class="n">${a.codetected ?? "–"}</td><td class="n">${fmt(a.tf1_det)}</td><td class="n">${fmt(a.tf2_det)}</td><td class="n">${fmt(a.donor_mean_r)}</td><td class="n">${a.donors_tested ?? "–"}</td><td class="n">${fmt(a.pos_donor_frac)}</td><td class="muted">${esc(a.reason || "")}</td></tr>`)).join("")}</table></div>`;
}

function atacOk(p, v) { const a = p.atac_level; if (!v) return true; if (!a) return false;
  if (v === "tested") return ["catlas_coenriched", "not_coenriched"].includes(a.status);
  if (v === "gap") return ["no_matching_context", "no_exact_celltype", "ontology_present_cluster_unresolved"].includes(a.status);
  if (v === "proxy") return !!a.proxy_status; if (v.startsWith("proxy_")) return a.proxy_status === v; return a.status === v; }
function atacRow(P) {
  const L = (D.run.layers || {}).scatac || {};
  if (L.status !== "done") return `<tr><td>scATAC evidence (A–D)</td><td><span class="badge b-todo">not yet run</span></td><td class="muted">Levels A–D: motif co-accessibility, with motif ambiguity handled</td></tr>`;
  const R = P.filter(p => p.atac_level && p.atac_level.status !== "not_tested"), a = R[0] && R[0].atac_level;
  return `<tr><td>scATAC level B (CATlas)</td><td>${a ? atacCell(R[0]) : '<span class="badge b-todo">not tested</span>'}</td><td>${a ? `${esc(a.label)}${a.fetal_clusters ? ` – CATlas has only <i>${esc(a.fetal_clusters)}</i>` : ""}${a.proxy_label ? `; exploratory fetal proxy: ${esc(a.proxy_label)}` : ""}. Levels A, C, D not run.` : "Only Tabula Sapiens–replicated rows are sent to CATlas"}</td></tr>`;
}
function atacPanel(P) {
  const L = (D.run.layers || {}).scatac || {}; if (L.status !== "done") return "";
  const R = P.filter(p => p.atac_level && p.atac_level.status !== "not_tested");
  const st = x => x == null ? "–" : fmt(x);
  const rowsOf2 = (p, list, kind) => list.map(m => `<tr><td>${esc(p.ct)} · ${esc(p.assay)}</td><td>${kind === "proxy" ? '<span class="badge b-proxy">fetal proxy</span>' : "adult/exact"}</td><td>${esc(m.cluster || "—")}</td><td>${esc(m.stage || "–")}</td><td>${esc(m.tier || "–")}${m.reviewed === "Yes" ? " (reviewed)" : ""}</td><td>${atacBadge(m.result, m.label)}</td><td class="n">${m.n_matched_foreground ?? "–"} / ${m.n_matched_background ?? "–"}</td><td class="n">${m.tf1_fg_hits ?? "–"}/${m.tf1_bg_hits ?? "–"} · OR ${st(m.tf1_odds_ratio)} · q ${st(m.tf1_q)}</td><td class="n">${m.tf2_fg_hits ?? "–"}/${m.tf2_bg_hits ?? "–"} · OR ${st(m.tf2_odds_ratio)} · q ${st(m.tf2_q)}</td><td class="n">${m.joint_fg_hits ?? "–"}/${m.joint_bg_hits ?? "–"} · OR ${st(m.joint_odds_ratio)} · q ${st(m.joint_q)}</td><td class="n">${st(m.joint_lift_fg)} vs ${st(m.joint_lift_bg)}</td><td class="muted">${esc(m.tf1_motifs)} | ${esc(m.tf2_motifs)}</td></tr>`).join("");
  const body = R.flatMap(p => [rowsOf2(p, p.atac_level.exact || [], "exact"), rowsOf2(p, p.atac_level.proxy_fetal || [], "proxy")]).join("");
  const rv = ((L.exact || L.proxy_fetal || {}).mapping_review || []);
  const cw = ((L.proxy_fetal || {}).crosswalk || []);
  return `<div class="card"><h2>scATAC evidence – level B: CATlas same-cCRE motif co-enrichment</h2>
   <p class="caveat"><b>Caveat.</b> ${esc(L.caveat || "")}</p>
   ${L.proxy_fetal ? `<p class="caveat proxy"><b>Fetal proxy.</b> ${esc(L.proxy_fetal.note)}</p>` : ""}
   <h3>Cell-type mapping (RNA → CATlas)</h3><table><tr><th>RNA cell type</th><th>Tissue</th><th>Exact CATlas cluster (ontology)</th><th>Nearest names</th><th>Fetal-only counterpart</th><th>Proxy crosswalk</th></tr>
   ${rv.map(r => `<tr><td>${esc(r.cell_type)} <span class="muted">${esc(r.cell_type_id)}</span></td><td>${esc(r.tissue)}</td><td>${esc(r.exact_ontology_clusters || "none")}</td><td class="muted">${esc(r.name_suggestions)}</td><td>${esc(r.fetal_only_clusters || "–")}</td><td>${cw.filter(c => c.cell_type_id === r.cell_type_id).map(c => `${esc(c.catlas_cluster)} (${esc(c.mapping_tier)})`).join(", ") || "–"}</td></tr>`).join("")}</table>
   <h3>Enrichment (accessible vs GC/width-matched inaccessible distal cCREs; one-sided Fisher, BH q)</h3>
   ${R.length ? `<table><tr><th>Discovery context</th><th>Run</th><th>CATlas cluster</th><th>Stage</th><th>Mapping tier</th><th>Result</th><th>Matched fg / bg cCREs</th><th>TF A hits fg/bg</th><th>TF B hits fg/bg</th><th>Joint (same cCRE) fg/bg</th><th>Joint lift fg vs bg</th><th>JASPAR motifs</th></tr>${body}</table>` : '<p class="muted">No replicated row for this pair was sent to CATlas.</p>'}
   <p class="muted">CATlas mirror ${esc(((L.exact || L.proxy_fetal || {}).manifest || {}).catlas_mirror_repo || "")}@${esc((((L.exact || L.proxy_fetal || {}).manifest || {}).catlas_mirror_revision || "").slice(0, 10))}; JASPAR 2026 CORE vertebrates; hg38. Thresholds: q ≤ ${(((L.exact || {}).summary) || {}).q_threshold ?? 0.05}, joint accessible cCREs ≥ ${(((L.exact || {}).summary) || {}).min_joint_peaks ?? 20}. Levels A, C and D are not run.</p></div>`;
}

function bindSort(P, fn, arg) { document.querySelectorAll("th[data-k]").forEach(th => th.onclick = () => { const k = th.dataset.k; S.sort = {k, dir: S.sort.k === k ? -S.sort.dir : (k === "min_q" || k === "meta_q" ? 1 : -1)}; route(); }); }

async function table() {
  const f = S.f;
  const P = sortBy(D.TOP.filter(p => (f.r === "" || (p.primary_r ?? p.r) >= +f.r) && (f.q === "" || (p.min_q != null && p.min_q <= +f.q)) && (f.mq === "" || (p.meta_q != null && p.meta_q <= +f.mq))
    && (f.lift === "" || p.lift >= +f.lift) && (f.coll === "" || p.n_collections >= +f.coll) && (!f.prior || p.prior) && (!f.rep || isRep(p)) && atacOk(p, f.atac)), S.sort.k, S.sort.dir);
  const rc = D.run.summary.run_config, m = D.run.manifest;
  $("#app").innerHTML = `<div class="card"><h2>Scientific table</h2><p class="muted">Top ${D.TOP.length} of ${D.top.total} rows. <a href="${base()}${esc(D.run.downloads["cell_type_pair_candidates.tsv"])}" download>Download the full TSV</a>.</p><div class="filters">
   <label>Min primary r<input id="fr" type="number" step="0.05" value="${f.r}"></label><label>Max q<input id="fq" type="number" step="0.01" value="${f.q}"></label><label>Max meta q<input id="fmq" type="number" step="0.01" value="${f.mq}"></label>
   <label>Min lift<input id="fl" type="number" step="0.1" value="${f.lift}"></label><label>Min collections<input id="fc" type="number" min="1" value="${f.coll}"></label><label><span>Prior only</span><input id="fp" type="checkbox" ${f.prior ? "checked" : ""} style="width:auto"></label><label><span>Replicated only</span><input id="frep" type="checkbox" ${f.rep ? "checked" : ""} ${repLayerDone() ? "" : "disabled"} style="width:auto"></label><label>ATAC (level B)<select id="fatac" ${atacDone() ? "" : "disabled"}>${[["", "any"], ["tested", "tested (adult)"], ["catlas_coenriched", "co-enriched (adult)"], ["gap", "no CATlas match"], ["proxy", "fetal proxy tested"], ["proxy_coenriched", "fetal proxy co-enriched"], ["proxy_not_coenriched", "fetal proxy not co-enriched"]].map(([v, l]) => `<option value="${v}" ${f.atac === v ? "selected" : ""}>${l}</option>`).join("")}</select></label><span class="muted">${P.length} rows</span></div>
   <p class="muted">Primary statistic <b>${esc(rc.primary_stat || "cell_raw")}</b>; detection ≥ ${rc.min_detection}; r ≥ ${rc.min_rho}; ${rc.primary_stat && rc.primary_stat !== "cell_raw" ? `BH q ≤ ${rc.max_q}; ` : ""}lift ≥ ${rc.min_lift}; co-detected ≥ ${rc.min_codetected}; min cells ${rc.min_cells}. Census ${esc(m.resolved_census_version)}.</p>${pairTable(P)}</div>`;
  const upd = () => { S.f = {r: $("#fr").value, q: $("#fq").value, mq: $("#fmq").value, lift: $("#fl").value, coll: $("#fc").value, prior: $("#fp").checked, rep: $("#frep").checked, atac: $("#fatac").value}; table(); };
  ["#fr", "#fq", "#fmq", "#fl", "#fc", "#fp", "#frep", "#fatac"].forEach(s => $(s).addEventListener("change", upd)); bindSort();
}

function info() {
  const r = D.run, m = r.manifest, s = r.summary; const row = (k, v) => `<tr><th>${esc(k)}</th><td style="white-space:normal">${v}</td></tr>`;
  $("#app").innerHTML = `<div class="card"><h2>Run information &amp; provenance</h2><table>
   ${row("Run", `<code>${esc(r.run_id)}</code>: ${esc(r.label || "")} (schema v${r.schema_version})`)}${row("Status", r.pilot ? "<b>Pilot / partial</b>, not atlas-wide" : "<b>Atlas-wide</b>")}
   ${row("Census", esc(m.resolved_census_version))}${row("TF set", esc(m.tf_set))}${row("Compendium", `${m.n_mapped_genes}/${m.n_register_genes} genes mapped (${m.n_tf_mapped} TFs); sha256 <code>${esc(m.register_sha256)}</code>`)}
   ${row("Datasets", `${s.n_completed_datasets}/${s.n_planned_datasets} (${esc(s.dataset_selection)})`)}${row("Normalization", esc(m.normalization))}
   ${row("Run configuration", `<pre>${esc(JSON.stringify(s.run_config, null, 1))}</pre>`)}${row("Positive controls", `${s.positive_controls_recovered ?? "–"} of ${s.positive_controls_tested ?? "–"} recovered`)}
   ${row("Layers", layerBadges())}${row("Embedding", ((r.embedding || {}).datasets || []).map(d => `${esc(d.title)}: ${d.n_cells}/${d.n_primary_cells} cells, ${d.genes.length} TFs`).join("<br>") || "not run")}
   ${row("Exported", `${esc(r.exported_at)} (code ${esc(r.code_commit_at_export || "–")})`)}${row("Downloads", Object.entries(r.downloads).map(([k, v]) => `<a href="${base()}${esc(v)}" download>${esc(k)}</a>`).join(" · "))}
   ${row("All runs", D.index.runs.map(x => `<a href="?run=${encodeURIComponent(x.run_id)}#home">${esc(x.label || x.run_id)}</a>${x.run_id === D.index.latest ? " (default)" : ""}`).join("<br>"))}</table></div>`;
}
function about() {
  $("#app").innerHTML = `<div class="card"><h2>About</h2><p>This tool screens human single-cell RNA data in the CELLxGENE Census for transcription-factor pairs that are co-expressed within the same cell type × tissue × disease × assay context. For each pair it reports correlation (raw cell-level, or donor- and library-adjusted metacell/pseudobulk), co-detection lift, Benjamini–Hochberg q within each context, a cross-dataset meta-analysis, and literature priors from the TF compendium. The results are candidates for wet-lab testing, not proof that the TFs cooperate. Replication in Tabula Sapiens (independent donors, same cell type × tissue × disease; discovery datasets excluded) is shown as badges and a per-assay table when a run has it; scATAC level B (CATlas same-cCRE motif co-enrichment for replicated pairs; fetal-proxy results are labelled separately and are never adult validation) is shown as a badge and a pair-page panel; SCEA replication and scATAC levels A, C, D are later layers.</p></div>
   <div class="card"><h2>Run it</h2><pre>git clone https://github.com/double2strand/TFpairminer.git && cd TFpairminer
./tfpm setup                    # venv, deps, tests, Census plan
./tfpm mine --runner wr         # cluster (wr manager must be running); or --runner local --parallel 6
./tfpm status
./tfpm publish --embed          # validate, export versioned JSON (+ optional UMAP layer), commit, push</pre>
   <p class="muted">Full details: <a href="https://github.com/double2strand/TFpairminer#readme">README</a> · <a href="https://github.com/double2strand/TFpairminer/blob/main/mining/METHOD.md">METHOD</a></p></div>`;
}
boot();
