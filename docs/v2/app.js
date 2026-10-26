"use strict";
const $ = s => document.querySelector(s), main = () => $("#main");
const D = {index: null, tf: {}, sample: {}};
const RL = {strong: "Strong hypothesis", moderate: "Moderate hypothesis", preliminary: "Preliminary hypothesis", family: "Family-level only", weak: "Weak / single-replicate"};
const RT = {strong: "Distinct motifs, enriched with replicate support in ≥3 profiles from ≥2 independent samples",
  moderate: "Distinct, replicate-supported in 2 profiles (or 1 with permutation p ≤ 0.01)", preliminary: "Distinct, replicate-supported in a single profile",
  family: "Enriched only where motifs are same-family / similar PWM / overlapping — may be one factor family", weak: "Enriched but not replicate-supported"};
const esc = s => String(s ?? "").replace(/[&<>"]/g, c => ({"&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;"}[c]));
const fmt = (x, d = 3) => x == null || x === "" ? "–" : (Math.abs(x) < 1e-3 && x != 0 ? Number(x).toExponential(1) : (+x).toFixed(d).replace(/\.?0+$/, ""));
const slug = s => String(s).replace(/[^A-Za-z0-9]+/g, "_").replace(/^_|_$/g, "");
const J = async p => { const r = await fetch(p); if (!r.ok) throw new Error(p + " " + r.status); return r.json(); };
const tfLink = t => `<a href="#/tf/${encodeURIComponent(t)}">${esc(t)}</a>`;
const pairLink = (a, b, txt) => `<a href="#/pair/${encodeURIComponent(a)}/${encodeURIComponent(b)}">${txt || esc(a) + " + " + esc(b)}</a>`;
const unitById = id => D.index.units.find(u => u.id === id) || {id, context: id};
const unitLink = id => `<a href="#/sample/${encodeURIComponent(id)}">${esc(shortU(id))}</a>`;
const shortU = id => id.replace(/^Corpus_/, "").replace(/^CCLA_/, "").replace("__", " · ").replace(/_/g, " ");
const rbadge = r => `<span class="badge r-${r}" title="${esc(RT[r])}">${RL[r]}</span>`;
const flagBadges = f => (f || []).filter(Boolean).map(x => `<span class="badge ${x === "repeat_rich" ? "fl-rep" : x === "broad" ? "fl-broad" : "fl"}" title="${esc(FT[x] || x)}">${esc(x.replace(/_/g, " "))}</span>`).join("");
const FT = {same_family: "Both TFs are in the same TF family — motif cannot distinguish them", similar_pwm: "PWMs are ≥0.75 similar", overlapping_sites: "≥50% of joint windows have overlapping sites",
  dimer_motif: "One motif is a heterodimer (::) motif", repeat_rich: "≥50% of joint windows are repeat-rich", broad: "Enriched in ≥75% of profiles — may reflect sequence composition rather than context-specific cooperation"};
const CB = {supported: ["cbp-yes", "ChromBPNet: contribution-supported"], not_supported: ["cbp-no", "ChromBPNet: not supported"],
  untestable: ["cbp-na", "ChromBPNet: too few joint windows"], not_run: ["cbp-nr", "Motif-only · ChromBPNet not yet run"]};
const cbadge = (k, short) => { const c = CB[k || "not_run"]; return `<span class="badge ${c[0]}" title="${esc(c[1])}">${short ? ({supported: "contrib ✓", not_supported: "contrib ✗", untestable: "untestable", not_run: "not yet run"})[k || "not_run"] : c[1]}</span>`; };
const COOP = {"synergistic": "synergistic (joint ablation > sum of singles)", "sub-additive": "sub-additive (joint < sum)", "additive/ns": "additive / not significant"};
function trackSVG(t) {
  const W = 760, H = 110, n = t.contrib.length, x = i => 30 + i * (W - 40) / n, mx = Math.max(1e-6, ...t.contrib.map(Math.abs)), y0 = 70, sc = v => v / mx * 55;
  const colsC = ["#16a34a", "#2563eb", "#f59e0b", "#dc2626"];
  let s = `<svg viewBox="0 0 ${W} ${H}" class="trk"><line x1="30" x2="${W - 10}" y1="${y0}" y2="${y0}" stroke="#cbd5e1"/>`;
  t.sites.forEach((st, k) => s += `<rect x="${x(st.pos)}" y="4" width="${Math.max(2, x(st.pos + st.len) - x(st.pos))}" height="${H - 24}" fill="${["#fde68a", "#bfdbfe", "#fbcfe8"][k % 3]}" opacity=".7"/><text x="${x(st.pos)}" y="${H - 6}" font-size="10">${esc(st.tf)}</text>`);
  t.contrib.forEach((v, i) => { const h = sc(v); s += `<rect x="${x(i)}" y="${h > 0 ? y0 - h : y0}" width="${(W - 40) / n}" height="${Math.abs(h) || .3}" fill="${colsC["ACGT".indexOf(t.seq[i])] || "#999"}"/>`; });
  return s + `</svg><div class="tag">${esc(t.chrom)}:${t.start.toLocaleString()}–${(t.start + n).toLocaleString()} · ${esc(shortU(t.u))} · bar height = projected contribution to predicted log-counts (A green, C blue, G amber, T red) · shaded = motif sites</div>`;
}
function cbpPanel(per, tracks, kind) {
  const ent = Object.entries(per || {}), tested = ent.filter(([, e]) => e.st === "tested"), ran = Object.keys(D.index.cbp_units || {});
  const notRun = D.index.units.filter(u => !ran.includes(u.id) && !(per || {})[u.id]).map(u => shortU(u.id));
  const rows = ent.map(([u, e]) => e.st !== "tested" ? `<tr><td>${unitLink(u)}</td><td colspan="6" class="tag">${e.st === "too_few_windows" ? `only ${e.n} windows with all sites (≥8 needed)` : esc(e.st)}</td></tr>` :
    `<tr><td>${unitLink(u)}</td><td>${e.sup === "yes" ? '<b class="ok">yes</b>' : "no"}</td><td class="n">${(e.pct || []).map(v => fmt(v, 2)).join(" / ")}</td><td class="n">${fmt(e.syn, 3)}</td><td class="n">${fmt(e.q)}</td><td>${esc(COOP[e.coop] || e.coop)}</td><td class="n">${e.n}</td></tr>`).join("");
  return `<div class="card cbp"><h2>ChromBPNet evidence layer <span class="tag">separate from motif co-occurrence</span></h2>
   ${ent.length ? `<div class="tablewrap"><table><thead><tr><th>Profile</th><th>Contribution-supported</th><th>Median site contribution percentile</th><th>Synergy (Δjoint − ΣΔsingle, log-counts)</th><th>q</th><th>Ablation call</th><th>Joint windows</th></tr></thead><tbody>${rows}</tbody></table></div>` : `<p>${cbadge("not_run")} This ${kind} was not among the candidates scored with ChromBPNet in any profile.</p>`}
   ${(tracks || []).length ? `<h3>Example contribution tracks</h3>${tracks.slice(0, 3).map(trackSVG).join("")}` : ""}
   <p class="tag"><b>Not yet run:</b> ${notRun.length ? esc(notRun.join(", ")) : "none"}. Contribution-supported = the median contribution percentile of every motif site ≥ 0.75 relative to random same-length positions in the same windows. Synergy is from in-silico ablation (site bases shuffled, 2 seeds) on the model's count head; Wilcoxon over windows, BH across candidates. Models: ENCODE pretrained bias-corrected ChromBPNet (bulk ATAC of the same cell line, fold 0), not trained on CCLA fragments.</p></div>`;
}
const status = {D: "distinct", F: "family-level", E: "exploratory"};
function color(l) { if (l == null) return "#f3f4f6"; const v = Math.max(-1, Math.min(1, Math.log2(l) / 1.5)); return v >= 0 ? `rgba(220,38,38,${0.12 + 0.88 * v})` : `rgba(37,99,235,${0.12 - 0.88 * v})`; }

async function boot() {
  D.index = await J("data/index.json"); $("#gen").textContent = "data generated " + D.index.generated;
  D.pf = D.index.pair_fields; D.pairs = D.index.pairs.map(p => Object.fromEntries(D.pf.map((k, i) => [k, p[i]])));
  D.tfs = Object.keys(D.index.tf_slugs).sort();
  document.querySelectorAll("nav button").forEach(b => b.onclick = () => location.hash = b.dataset.r);
  setupSearch(); window.onhashchange = route; route();
}
async function tfData(t) { if (!D.tf[t]) D.tf[t] = await J(`data/tf/${D.index.tf_slugs[t] || slug(t)}.json`); return D.tf[t]; }
async function sampleData(u) { if (!D.sample[u]) D.sample[u] = await J(`data/samples/${slug(u)}.json`); return D.sample[u]; }

function setupSearch() {
  const q = $("#q"), ac = $("#ac"); let sel = 0, items = [];
  const go = it => { ac.hidden = true; q.value = ""; location.hash = it.h; };
  q.oninput = () => {
    const v = q.value.trim().toLowerCase(); if (!v) { ac.hidden = true; return; }
    const toks = v.split(/[\s+,/]+/).filter(Boolean); items = [];
    if (toks.length >= 2) D.pairs.filter(p => toks.every(t => p.tf_a.toLowerCase().startsWith(t) || p.tf_b.toLowerCase().startsWith(t))).slice(0, 8)
      .forEach(p => items.push({h: `#/pair/${encodeURIComponent(p.tf_a)}/${encodeURIComponent(p.tf_b)}`, l: `${p.tf_a} + ${p.tf_b}`, s: RL[p.rating]}));
    D.tfs.filter(t => t.toLowerCase().startsWith(toks[0])).slice(0, 10).forEach(t => items.push({h: "#/tf/" + encodeURIComponent(t), l: t, s: "TF"}));
    D.tfs.filter(t => !t.toLowerCase().startsWith(toks[0]) && t.toLowerCase().includes(toks[0])).slice(0, 5).forEach(t => items.push({h: "#/tf/" + encodeURIComponent(t), l: t, s: "TF"}));
    D.index.units.filter(u => (u.id + " " + u.context).toLowerCase().includes(v)).slice(0, 8).forEach(u => items.push({h: "#/sample/" + encodeURIComponent(u.id), l: shortU(u.id), s: u.dataset}));
    sel = 0; ac.innerHTML = items.map((it, i) => `<div class="${i ? "" : "sel"}" data-i="${i}">${esc(it.l)}<small>${esc(it.s)}</small></div>`).join("") || "<div><small>no match</small></div>";
    ac.hidden = false; ac.querySelectorAll("div[data-i]").forEach(d => d.onmousedown = () => go(items[+d.dataset.i]));
  };
  q.onkeydown = e => { if (e.key === "Enter" && items[sel]) go(items[sel]); if (e.key === "ArrowDown" || e.key === "ArrowUp") { sel = Math.max(0, Math.min(items.length - 1, sel + (e.key === "ArrowDown" ? 1 : -1))); ac.querySelectorAll("div").forEach((d, i) => d.classList.toggle("sel", i === sel)); e.preventDefault(); } };
  q.onblur = () => setTimeout(() => ac.hidden = true, 150);
}
async function route() {
  const h = decodeURIComponent(location.hash || "#/"), parts = h.slice(2).split("/");
  document.querySelectorAll("nav button").forEach(b => b.classList.toggle("on", h === b.dataset.r || (b.dataset.r !== "#/" && h.startsWith(b.dataset.r))));
  try {
    if (parts[0] === "tf") await viewTF(parts[1]);
    else if (parts[0] === "pair") await viewPair(parts[1], parts[2]);
    else if (parts[0] === "sample") await viewSample(parts.slice(1).join("/"));
    else if (parts[0] === "samples") viewSamples();
    else if (parts[0] === "pairs") viewPairs();
    else if (parts[0] === "sets") viewSets();
    else if (parts[0] === "set") await viewSet(parts[1]);
    else if (parts[0] === "about") viewAbout();
    else viewHome();
  } catch (e) { main().innerHTML = `<div class="card">Could not load: ${esc(e.message)}</div>`; }
  window.scrollTo(0, 0);
}
function sortable(el) {
  el.querySelectorAll("table.sort").forEach(t => t.querySelectorAll("th").forEach((th, ci) => th.onclick = () => {
    const rows = [...t.tBodies[0].rows], dir = th.dataset.dir === "▲" ? -1 : 1;
    rows.sort((a, b) => { const x = a.cells[ci].dataset.v ?? a.cells[ci].textContent, y = b.cells[ci].dataset.v ?? b.cells[ci].textContent;
      const nx = parseFloat(x), ny = parseFloat(y); return (isNaN(nx) || isNaN(ny) ? x.localeCompare(y) : nx - ny) * dir; });
    t.querySelectorAll("th").forEach(o => { o.classList.remove("sorted"); o.dataset.dir = ""; }); th.classList.add("sorted"); th.dataset.dir = dir === 1 ? "▲" : "▼";
    rows.forEach(r => t.tBodies[0].appendChild(r)); }));
}
function pairsTable(list, opts = {}) {
  return `<div class="tablewrap"><table class="sort"><thead><tr><th>Pair</th><th>Rating</th><th>Profiles enriched</th><th>Distinct + replicated</th><th>Median lift</th><th>Best q</th><th>Spacing (bp)</th><th>Datasets</th><th>Flags</th><th>Evidence</th></tr></thead><tbody>${
    list.map(p => `<tr><td>${pairLink(p.tf_a, p.tf_b)}</td><td data-v="${["strong", "moderate", "preliminary", "family", "weak"].indexOf(p.rating)}">${rbadge(p.rating)}</td><td class="n">${p.n_sig}/${D.index.units.length}</td><td class="n">${p.n_distinct_rep}</td><td class="n">${fmt(p.median_lift ?? p.median_lift_sig, 2)}</td><td class="n" data-v="${p.best_q}">${fmt(p.best_q)}</td><td class="n">${p.median_gap == null ? "–" : Math.round(p.median_gap)}</td><td>${esc(p.datasets)}</td><td>${flagBadges(typeof p.flags === "string" ? p.flags.split("|") : p.family_flags)}</td><td data-v="${["supported", "not_supported", "untestable", "not_run"].indexOf(p.cbp || "not_run")}">${cbadge(p.cbp, 1)}${(p.cbp_coop || "").includes("synergistic") && !(p.cbp_coop || "").includes("sub") ? ' <span class="badge cbp-syn">synergy</span>' : ""}</td></tr>`).join("")}</tbody></table></div>`;
}
function viewHome() {
  const I = D.index, c = I.counts, us = I.units, ds = [...new Set(us.map(u => u.dataset))];
  const top = D.pairs.filter(p => p.rating === "strong" || p.rating === "moderate").slice(0, 25);
  main().innerHTML = `
  <div class="card"><h2>What v2 does</h2><p>TFpairminer v2 mines scATAC accessible windows for TF <b>motif pairs</b>, <b>triplet sets</b> and <b>modules</b> that co-occur more often than expected after matching GC, accessibility, motif density and repeat content. Two independent sources are analysed in parallel, separate from the v1 RNA co-expression layers:</p>
  <ul><li><b>CCLA / scatlas HCL cancer cell lines</b>: raw scATAC fragments, two replicate libraries per line.</li><li><b>Human scATAC Corpus</b>: per-cell-type pseudobulk cCREs, with random half-cell pseudoreplicates.</li></ul>
  <p class="muted">Ratings combine replicate support, cross-sample consistency and family ambiguity. They rank <i>hypotheses worth testing</i> and are not confidence in a physical interaction.</p>
  <div class="stats"><div class="stat"><b>${us.length}</b>profiles analysed</div><div class="stat"><b>${new Set(us.map(u => u.parent)).size}</b>samples / cell lines</div><div class="stat"><b>${ds.length}</b>datasets</div><div class="stat"><b>${c.pairs_listed.toLocaleString()}</b>enriched pairs listed</div><div class="stat"><b>${c.tfs}</b>TF motifs with partners</div><div class="stat"><b>${I.top_sets.length}</b>synergistic sets (top)</div></div>
  <h3>Pairs by rating</h3><div>${Object.entries(c.by_rating).map(([r, n]) => `${rbadge(r)} ${n.toLocaleString()}`).join(" &nbsp; ")}</div>
  <h3>Evidence layer</h3><div>${Object.entries(I.cbp_counts || {}).map(([k, n]) => `${cbadge(k)} ${n.toLocaleString()}`).join(" &nbsp; ")}</div>
  <p class="tag">ChromBPNet has been run for ${Object.keys(I.cbp_units || {}).map(shortU).join(", ") || "no profile"}. Ratings come from motif co-occurrence only; the evidence badge says whether a sequence-to-accessibility model also attributes accessibility to both motif sites.</p></div>
  <div class="grid"><div class="card"><h2>Search</h2><p>Use the search bar for a <b>TF</b> (ranked partners, sets, consistency heatmap), a <b>pair</b> (type both names) or a <b>sample, cell line or cell type</b> (its candidates).</p>
   <div class="chips">${["FOXA2", "HNF1A", "HNF4A", "CEBPB", "GATA1", "TAL1::TCF3", "FOS::JUN", "TEAD1", "KLF5", "GRHL2"].filter(t => D.index.tf_slugs[t]).map(tfLink).join("")}</div>
   <h3>Samples</h3><div class="chips">${us.map(u => unitLink(u.id)).join("")}</div></div>
   <div class="card"><h2>Profiles</h2>${unitTable(us)}</div></div>
  <div class="card"><h2>Top consistent pair hypotheses</h2>${pairsTable(top)}<p><a href="#/pairs">All pairs →</a></p></div>`;
  sortable(main());
}
function unitTable(us) {
  return `<div class="tablewrap" style="max-height:420px"><table class="sort"><thead><tr><th>Profile</th><th>Dataset</th><th>Windows</th><th>Distinct pairs</th><th>Family-level</th><th>Synergistic sets</th></tr></thead><tbody>${us.map(u => `<tr><td>${unitLink(u.id)}<div class="tag">${esc(u.context)}</div></td><td>${esc(u.dataset.split(" (")[0])}</td><td class="n">${(u.n_windows || {}).pooled ?? "–"}</td><td class="n">${u.n_priority_distinct}</td><td class="n">${u.n_family_level}</td><td class="n">${u.n_synergistic_sets}</td></tr>`).join("")}</tbody></table></div>`;
}
function viewPairs() {
  main().innerHTML = `<div class="card"><h2>All enriched pairs (${D.pairs.length.toLocaleString()})</h2><div class="filters">
   <label>Rating<select id="fr"><option value="">any</option>${Object.keys(RL).map(r => `<option value="${r}">${RL[r]}</option>`).join("")}</select></label>
   <label>TF contains<input id="ft"></label><label>Min profiles<input id="fn" type="number" value="1" min="1"></label>
   <label>Dataset<select id="fd"><option value="">any</option>${[...new Set(D.index.units.map(u => u.dataset))].map(d => `<option>${esc(d)}</option>`).join("")}</select></label>
   <label>Hide broad/family<select id="fh"><option value="">no</option><option value="1">yes</option></select></label>
   <label>Evidence<select id="fe"><option value="">any</option><option value="supported">ChromBPNet contribution-supported</option><option value="tested">ChromBPNet tested (any result)</option><option value="not_supported">tested, not supported</option><option value="synergy">ablation synergy</option><option value="not_run">motif-only (not yet run)</option></select></label></div><div id="pt"></div></div>`;
  const upd = () => { const e = $("#fe").value, r = $("#fr").value, t = $("#ft").value.toLowerCase(), n = +$("#fn").value || 1, d = $("#fd").value, h = $("#fh").value;
    const l = D.pairs.filter(p => (!r || p.rating === r) && (!t || (p.tf_a + " " + p.tf_b).toLowerCase().includes(t)) && p.n_sig >= n && (!d || p.datasets.includes(d)) && (!h || (p.rating !== "family" && !p.flags.includes("broad"))) && (!e || (e === "tested" ? ["supported", "not_supported"].includes(p.cbp) : e === "synergy" ? (p.cbp_coop || "").split(",").includes("synergistic") : (p.cbp || "not_run") === e)));
    $("#pt").innerHTML = `<p class="muted">${l.length.toLocaleString()} pairs${l.length > 1500 ? " (first 1500 shown)" : ""}</p>` + pairsTable(l.slice(0, 1500)); sortable($("#pt")); };
  ["#fr", "#ft", "#fn", "#fd", "#fh", "#fe"].forEach(s => $(s).oninput = upd); upd();
}
function heatmap(rows, units, getter, rowLabel) {
  const legend = `<div class="legend tag">log2 lift: ${[0.5, 0.75, 1, 1.5, 2, 3].map(l => `<span style="background:${color(l)}"></span>${l}×`).join(" ")} &nbsp; ● = significant distinct + replicated · ○ = significant (family-level or single-replicate) · blank = not enriched / not tested</div>`;
  return legend + `<div class="heat"><table class="hm"><thead><tr><th></th>${units.map(u => `<th class="v" title="${esc(u.context)}">${esc(shortU(u.id))}</th>`).join("")}</tr></thead><tbody>${
    rows.map(r => `<tr><th style="text-align:right">${rowLabel(r)}</th>${units.map(u => { const v = getter(r, u.id); if (!v) return `<td style="background:#f9fafb"></td>`;
      const mark = (v.s === "D" && v.r === "y") ? "●" : (v.s === "D" || v.s === "F") ? "○" : "";
      return `<td style="background:${color(v.l)}" title="${esc(shortU(u.id))}: lift ${fmt(v.l, 2)}, q ${fmt(v.q)}, ${status[v.s]}, replicate ${v.r === "y" ? "yes" : "no"}">${mark}</td>`; }).join("")}</tr>`).join("")}</tbody></table></div>`;
}
async function viewTF(t) {
  const d = await tfData(t), units = D.index.units;
  const ps = d.pairs.map(p => ({...p, partner: p.tf_a === t ? p.tf_b : p.tf_a, flags: p.family_flags.join("|")}));
  const best = ps.filter(p => p.rating !== "weak").slice(0, 40);
  main().innerHTML = `<div class="crumb"><a href="#/">Overview</a> › TF</div><div class="card"><p class="big"><b>${esc(t)}</b> <span class="tag">motif partners from scATAC co-occurrence</span></p>
   <div class="kv"><span><b>${ps.length}</b> enriched partners</span><span><b>${ps.filter(p => p.rating === "strong").length}</b> strong</span><span><b>${ps.filter(p => p.rating === "moderate").length}</b> moderate</span><span><b>${d.sets.length}</b> sets</span><span><b>${d.modules.length}</b> module memberships</span></div></div>
  <div class="card"><h2>Ranked partners</h2>${pairsTable(ps.map(p => ({...p, median_lift: p.median_lift_sig})))}</div>
  <div class="card"><h2>Cross-sample consistency (top ${Math.min(40, best.length)} partners)</h2>${heatmap(best, units, (p, u) => p.per[u], p => pairLink(p.tf_a, p.tf_b, esc(p.partner)))}</div>
  <div class="grid"><div class="card"><h2>Higher-order sets containing ${esc(t)}</h2>${setsTable(d.sets.slice(0, 100))}</div>
  <div class="card"><h2>Modules (cliques of distinct, replicated pairs)</h2>${modulesTable(d.modules.slice(0, 60))}</div></div>`;
  sortable(main());
}
function setsTable(s) {
  if (!s.length) return `<p class="muted">No triplet passed the synergy test (q ≤ 0.05 against the best constituent pair).</p>`;
  return `<div class="tablewrap" style="max-height:420px"><table class="sort"><thead><tr><th>Set</th><th>Profile</th><th>Windows with all 3</th><th>Lift vs independent</th><th>Synergy lift vs best pair</th><th>q</th></tr></thead><tbody>${
    s.map(x => `<tr><td><a href="#/set/${encodeURIComponent(x.tfs.join("+"))}">${x.tfs.map(esc).join(" + ")}</a></td><td>${unitLink(x.u)}</td><td class="n">${x.n}</td><td class="n">${fmt(x.li, 2)}</td><td class="n">${fmt(x.sl, 2)}</td><td class="n" data-v="${x.q}">${fmt(x.q)}</td></tr>`).join("")}</tbody></table></div>`;
}
function modulesTable(m) {
  if (!m.length) return `<p class="muted">No modules.</p>`;
  return `<div class="tablewrap" style="max-height:420px"><table class="sort"><thead><tr><th>Module</th><th>Size</th><th>Profile</th><th>Windows with all motifs</th></tr></thead><tbody>${
    m.map(x => `<tr><td>${x.tfs.map(tfLink).join(" + ")}</td><td class="n">${x.tfs.length}</td><td>${unitLink(x.u)}</td><td class="n">${x.n}</td></tr>`).join("")}</tbody></table></div>`;
}
function liftPlot(p, units) {
  const W = 760, rowH = 18, L = 190, rows = units.filter(u => p.per[u.id]), H = rows.length * rowH + 40;
  const lifts = rows.map(u => p.per[u.id].l).concat([p.per[rows[0]?.id]?.r1, 2]).filter(Boolean);
  const mx = Math.max(2, ...rows.flatMap(u => [p.per[u.id].l, p.per[u.id].r1, p.per[u.id].r2]).filter(Boolean)) * 1.05;
  const x = v => L + (W - L - 120) * Math.log2(Math.max(v, 0.5) / 0.5) / Math.log2(mx / 0.5);
  let s = `<svg viewBox="0 0 ${W} ${H}" class="plot" style="min-height:${H}px"><line x1="${x(1)}" x2="${x(1)}" y1="10" y2="${H - 20}" stroke="#9ca3af" stroke-dasharray="3"/>`;
  [0.5, 1, 1.5, 2, 3, 4, 6].filter(v => v <= mx).forEach(v => s += `<text x="${x(v)}" y="${H - 6}" text-anchor="middle">${v}×</text>`);
  rows.forEach((u, i) => { const v = p.per[u.id], y = 20 + i * rowH, c = v.s === "D" && v.r === "y" ? "#166534" : v.s === "E" ? "#9ca3af" : "#d97706";
    s += `<text x="${L - 6}" y="${y + 4}" text-anchor="end">${esc(shortU(u.id).slice(0, 34))}</text>`;
    [v.r1, v.r2].forEach((r, k) => { if (r) s += `<circle cx="${x(r)}" cy="${y + (k ? 3 : -3)}" r="2.5" fill="none" stroke="#6b7280"><title>rep${k + 1} lift ${fmt(r, 2)}</title></circle>`; });
    s += `<line x1="${x(1)}" x2="${x(v.l)}" y1="${y}" y2="${y}" stroke="${c}" stroke-width="3"/><circle cx="${x(v.l)}" cy="${y}" r="5" fill="${c}"><title>lift ${fmt(v.l, 2)} q ${fmt(v.q)}</title></circle>
      <text x="${W - 112}" y="${y + 4}">q ${fmt(v.q)} ${v.s === "D" ? "" : "(" + status[v.s] + ")"}</text>`; });
  return s + `</svg><div class="tag">Filled dot = matched-background lift in pooled windows (green: distinct + replicated; amber: significant but family-level or single-replicate; grey: exploratory). Open circles = rep1/rep2 lifts. Dashed line = no enrichment.</div>`;
}
function volcano(list, onclickBase) {
  const W = 520, H = 300, pts = list.filter(p => p.l && p.q != null);
  if (!pts.length) return "";
  const mx = Math.max(...pts.map(p => p.l)) * 1.05, my = Math.max(...pts.map(p => -Math.log10(Math.max(p.q, 1e-300)))) * 1.05 || 1;
  const x = v => 40 + (W - 50) * (Math.log2(v) - Math.log2(0.9)) / (Math.log2(mx) - Math.log2(0.9)), y = v => H - 25 - (H - 40) * v / my;
  let s = `<svg viewBox="0 0 ${W} ${H}" class="plot" style="min-height:${H}px"><line x1="40" x2="${W}" y1="${H - 25}" y2="${H - 25}" stroke="#9ca3af"/><line x1="40" x2="40" y1="10" y2="${H - 25}" stroke="#9ca3af"/>
   <text x="${W / 2}" y="${H - 5}" text-anchor="middle">lift (log scale)</text><text x="12" y="${H / 2}" transform="rotate(-90 12 ${H / 2})" text-anchor="middle">−log10 q</text>`;
  [1, 1.5, 2, 3, 4, 6, 10].filter(v => v <= mx).forEach(v => s += `<text x="${x(v)}" y="${H - 14}" text-anchor="middle">${v}</text>`);
  pts.forEach(p => { const c = p.s === "D" && p.r === "y" ? "#166534" : p.s === "D" ? "#22c55e" : "#d97706";
    s += `<a href="#/pair/${encodeURIComponent(p.tf_a)}/${encodeURIComponent(p.tf_b)}"><circle cx="${x(p.l)}" cy="${y(-Math.log10(Math.max(p.q, 1e-300)))}" r="3" fill="${c}" fill-opacity=".6"><title>${esc(p.tf_a)} + ${esc(p.tf_b)}: lift ${fmt(p.l, 2)}, q ${fmt(p.q)}</title></circle></a>`; });
  return s + `</svg><div class="tag">dark green = distinct + replicated · light green = distinct, single replicate · amber = family-level</div>`;
}
async function viewPair(a, b) {
  const d = await tfData(a), p = d.pairs.find(x => (x.tf_a === a && x.tf_b === b) || (x.tf_a === b && x.tf_b === a));
  if (!p) { main().innerHTML = `<div class="card">${esc(a)} + ${esc(b)} is not among the enriched pairs.</div>`; return; }
  const units = D.index.units, fam = [p.fam_a, p.fam_b], per = Object.entries(p.per);
  const anyv = {ma: p.ma, mb: p.mb}, famBadge = p.rating === "family" || p.family_flags.length ? `<span class="badge fl">family-level caution</span>` : `<span class="badge r-preliminary">motifs distinguishable</span>`;
  main().innerHTML = `<div class="crumb"><a href="#/">Overview</a> › ${tfLink(p.tf_a)} › pair</div>
  <div class="card"><p class="big"><b>${tfLink(p.tf_a)} + ${tfLink(p.tf_b)}</b> ${rbadge(p.rating)} ${famBadge} ${flagBadges(p.family_flags)} ${cbadge(p.cbp)}</p>
   <div class="kv tag"><span>${esc(p.tf_a)}: ${esc(anyv.ma)} · family ${esc(p.fam_a || "unassigned")}</span><span>${esc(p.tf_b)}: ${esc(anyv.mb)} · family ${esc(p.fam_b || "unassigned")}</span>
   <span><a href="https://jaspar.elixir.no/matrix/${esc(anyv.ma)}/" target="_blank">JASPAR ${esc(anyv.ma)}</a> · <a href="https://jaspar.elixir.no/matrix/${esc(anyv.mb)}/" target="_blank">JASPAR ${esc(anyv.mb)}</a></span></div>
   <p class="summary">${esc(p.summary)}</p>
   <div class="stats"><div class="stat"><b>${p.n_sig}/${units.length}</b>profiles enriched</div><div class="stat"><b>${p.n_distinct_rep}</b>distinct + replicated</div><div class="stat"><b>${fmt(p.median_lift_sig, 2)}×</b>median lift</div><div class="stat"><b>${fmt(p.best_q)}</b>best q</div><div class="stat"><b>${p.median_gap == null ? "–" : Math.round(p.median_gap) + " bp"}</b>median spacing</div><div class="stat"><b>${esc(p.datasets.join(", "))}</b>datasets</div></div>
   <p class="tag">Rating: ${esc(RT[p.rating])}.</p></div>
  <div class="card"><h2>Lift and q across profiles</h2>${liftPlot(p, units)}</div>
  <div class="card"><h2>Cross-sample consistency</h2>${heatmap([p], units, (r, u) => r.per[u], () => esc(p.tf_a + " + " + p.tf_b))}
   <div class="tablewrap" style="max-height:380px;margin-top:8px"><table class="sort"><thead><tr><th>Profile</th><th>Status</th><th>Joint windows</th><th>Expected</th><th>Lift</th><th>q</th><th>Perm p</th><th>Rep1 lift</th><th>Rep2 lift</th><th>Replicated</th><th>Spacing</th><th>Flags</th></tr></thead><tbody>${
    per.map(([u, v]) => `<tr><td>${unitLink(u)}</td><td>${status[v.s]}</td><td class="n">${v.n}</td><td class="n">${fmt(v.e, 1)}</td><td class="n">${fmt(v.l, 2)}</td><td class="n" data-v="${v.q}">${fmt(v.q)}</td><td class="n">${fmt(v.pp)}</td><td class="n">${fmt(v.r1, 2)}</td><td class="n">${fmt(v.r2, 2)}</td><td>${v.r === "y" ? "yes" : "no"}</td><td class="n">${v.g == null ? "–" : Math.round(v.g)}</td><td>${flagBadges((v.amb || "").split(";").map(f => f.split(":")[0]))}</td></tr>`).join("")}</tbody></table></div>
   <p class="tag">Profiles not listed had no enrichment (q > 0.25) or one motif occurred in fewer than 20 windows.</p></div>
  <div class="card"><h2>Site examples</h2>${p.examples.length ? `<table class="ex"><thead><tr><th>Profile</th><th>Window</th><th>${esc(p.tf_a)} site</th><th>${esc(p.tf_b)} site</th><th>Gap</th><th>Cuts (rep1/rep2)</th><th></th></tr></thead><tbody>${
     p.examples.map(e => `<tr><td>${unitLink(e.u)}</td><td>${esc(e.w)}</td><td>${e.chrom}:${e.a[0]}-${e.a[1]}</td><td>${e.chrom}:${e.b[0]}-${e.b[1]}</td><td>${e.gap} bp</td><td>${e.cuts} (${esc(e.rc)})</td><td><a target="_blank" href="https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&position=${encodeURIComponent(e.w)}">UCSC</a></td></tr>`).join("")}</tbody></table>` : `<p class="muted">Site examples are exported for the top 60 significant pairs per profile. This pair has none.</p>`}
   <p class="tag">Coordinates are hg38, 0-based starts, and the most accessible joint windows come first. A site example shows co-occurring motif matches, not bound factors.</p></div>
  ${cbpPanel(p.cbp_per, p.cbp_tracks, "pair")}
  <div class="card"><h2>How to read this</h2><p>The lift is observed / expected joint windows, with the expectation computed within strata of GC, accessibility, motif density and repeat content. q is BH-adjusted across all pairs in that profile. "Distinct" means the two motifs are not the same family, are not highly similar PWMs, and do not mostly overlap. "Replicated" means the enrichment holds independently in both replicate (or pseudoreplicate) window sets. <b>Next tests:</b> co-binding (ChIP/CUT&amp;RUN) at the example sites, motif mutagenesis or MPRA, and, where the ChromBPNet panel shows support, joint motif mutagenesis. The rating above is motif-only. ChromBPNet evidence is a separate layer and does not change it.</p></div>`;
  sortable(main());
}
function viewSamples() { main().innerHTML = `<div class="card"><h2>Samples, cell lines and cell types</h2>${unitTable(D.index.units)}</div>`; sortable(main()); }
async function viewSample(u) {
  const d = await sampleData(u), U = d.unit, s = U.summary, W = U.windows_meta || {};
  const ps = d.pairs, frag = (W.fragments || []).map(f => f.n_fragments ? `${(f.n_fragments / 1e6).toFixed(1)}M fragments` : "").filter(Boolean).join(" / ");
  main().innerHTML = `<div class="crumb"><a href="#/samples">Samples</a> › ${esc(shortU(u))}</div><div class="card"><p class="big"><b>${esc(shortU(u))}</b> <span class="tag">${esc(U.dataset)}</span></p><p>${esc(U.context)}</p>
  <div class="stats"><div class="stat"><b>${(s.n_windows || {}).pooled}</b>pooled windows</div><div class="stat"><b>${(s.n_windows || {}).rep1 ?? "–"}/${(s.n_windows || {}).rep2 ?? "–"}</b>rep1/rep2 windows</div><div class="stat"><b>${s.n_motifs_tested}</b>motifs tested</div><div class="stat"><b>${(s.n_pairs_tested || 0).toLocaleString()}</b>pairs tested</div><div class="stat"><b>${s.n_priority_distinct}</b>distinct pairs</div><div class="stat"><b>${s.n_family_level}</b>family-level</div><div class="stat"><b>${s.n_synergistic_sets}</b>synergistic sets</div><div class="stat"><b>${s.n_modules}</b>modules</div></div>
  <p class="tag">${frag ? "Input: " + frag + ". " : ""}${W.pseudoreplicates ? "Replicates: " + esc(W.pseudoreplicates) + ". " : ""}Motif p-value 1e-4, JASPAR 2024 CORE vertebrates non-redundant.</p></div>
  ${d.cbp ? `<div class="card cbp"><h2>ChromBPNet evidence layer</h2><div class="stats"><div class="stat"><b>${d.cbp.n_candidates}</b>candidates scored</div><div class="stat"><b>${d.cbp.n_tested}</b>testable</div><div class="stat"><b>${d.cbp.n_contribution_supported}</b>contribution-supported</div><div class="stat"><b>${d.cbp.n_synergistic}</b>synergistic</div><div class="stat"><b>${d.cbp.n_subadditive}</b>sub-additive</div></div><p class="tag">Model: ${esc(d.cbp.model)}. ${esc(d.cbp.method)}</p>${(d.cbp.modisco || []).length ? `<h3>TF-MoDISco patterns (de novo, from contribution scores)</h3><div class="tablewrap" style="max-height:300px"><table><thead><tr><th>Pattern</th><th>Seqlets</th><th>Consensus</th><th>Best JASPAR match</th></tr></thead><tbody>${d.cbp.modisco.map(m => `<tr><td>${esc(m.pattern)}</td><td class="n">${m.n_seqlets}</td><td><code>${esc(m.consensus)}</code></td><td>${esc(m.best_jaspar)} <span class="tag">${esc(m.best_jaspar_id)}</span></td></tr>`).join("")}</tbody></table></div><p class="tag">Matches for patterns shorter than about 6 bp are unreliable.</p>` : ""}</div>` : `<div class="card cbp"><p>${cbadge("not_run")} ChromBPNet has not yet been run for this profile. All calls here are motif co-occurrence only.</p></div>`}
  <div class="grid"><div class="card"><h2>Lift vs q (top ${ps.length} significant pairs)</h2>${volcano(ps)}</div>
  <div class="card"><h2>Candidates</h2><div class="filters"><label>TF contains<input id="sf"></label><label>Show<select id="ss"><option value="D">distinct only</option><option value="DF">incl. family-level</option></select></label><label>Replicated only<select id="sr"><option value="1">yes</option><option value="">no</option></select></label></div><div id="st"></div></div></div>
  <div class="grid"><div class="card"><h2>Synergistic / higher-order sets</h2>${setsTable(d.sets.slice(0, 150))}</div><div class="card"><h2>Modules</h2>${modulesTable(d.modules)}</div></div>`;
  const upd = () => { const f = $("#sf").value.toLowerCase(), s_ = $("#ss").value, r = $("#sr").value;
    const l = ps.filter(p => s_.includes(p.s) && (!r || p.r === "y") && (!f || (p.tf_a + " " + p.tf_b).toLowerCase().includes(f)));
    $("#st").innerHTML = `<div class="tablewrap" style="max-height:520px"><table class="sort"><thead><tr><th>Pair</th><th>Status</th><th>Joint</th><th>Lift</th><th>q</th><th>Rep1/Rep2 lift</th><th>Spacing</th><th>Flags</th><th>ChromBPNet</th></tr></thead><tbody>${
      l.slice(0, 600).map(p => `<tr><td>${pairLink(p.tf_a, p.tf_b)}</td><td>${status[p.s]}</td><td class="n">${p.n}</td><td class="n">${fmt(p.l, 2)}</td><td class="n" data-v="${p.q}">${fmt(p.q)}</td><td class="n">${fmt(p.r1, 2)} / ${fmt(p.r2, 2)}</td><td class="n">${p.g == null ? "–" : Math.round(p.g)}</td><td>${flagBadges((p.amb || "").split(";").map(x => x.split(":")[0]))}</td><td>${p.cbp_e ? (p.cbp_e.st === "tested" ? `${p.cbp_e.sup === "yes" ? "✓" : "✗"} syn ${fmt(p.cbp_e.syn, 3)}${p.cbp_e.coop === "synergistic" ? "*" : ""}` : "untestable") : '<span class="tag">–</span>'}</td></tr>`).join("")}</tbody></table></div>`; sortable($("#st")); };
  ["#sf", "#ss", "#sr"].forEach(x => $(x).oninput = upd); upd();
}
function viewSets() {
  const s = D.index.top_sets, cnt = {};
  s.forEach(x => { const k = x.tfs.join("+"); cnt[k] = (cnt[k] || 0) + 1; });
  const rec = Object.entries(cnt).filter(([, n]) => n > 1).sort((a, b) => b[1] - a[1]).slice(0, 30);
  main().innerHTML = `<div class="card"><h2>Higher-order motif sets</h2><p>A triplet is called <b>synergistic</b> when the windows containing all three motifs exceed the count expected from the <i>best constituent pair</i> × the third motif's stratum frequency (Poisson, BH q ≤ 0.05, synergy lift ≥ 1.2). Triplets are built from one motif per family or PWM cluster, so they are not stacks of homologous PWMs. <b>Modules</b> are cliques of distinct, replicated pairs; see each TF or sample page.</p>
  ${rec.length ? `<h3>Sets synergistic in more than one profile</h3><div class="chips">${rec.map(([k, n]) => `<a href="#/set/${encodeURIComponent(k)}">${esc(k.replace(/\+/g, " + "))} (${n})</a>`).join("")}</div>` : ""}
  <h3>Top synergistic sets</h3>${setsTable(s)}</div>`; sortable(main());
}
async function viewSet(k) {
  const tfs = k.split("+"), d = await tfData(tfs[0]), hits = d.sets.filter(x => x.tfs.join("+") === k);
  const pairs = [[0, 1], [0, 2], [1, 2]].map(([i, j]) => D.pairs.find(p => (p.tf_a === tfs[i] && p.tf_b === tfs[j]) || (p.tf_a === tfs[j] && p.tf_b === tfs[i])) || {tf_a: tfs[i], tf_b: tfs[j], rating: "weak", missing: 1});
  main().innerHTML = `<div class="crumb"><a href="#/sets">Sets</a></div><div class="card"><p class="big"><b>${tfs.map(tfLink).join(" + ")}</b></p>
   <p class="summary">These three motifs occur together in accessible windows more often than expected from their strongest pairwise co-occurrence, in ${hits.length} profile(s). That makes this a candidate for a higher-order regulatory grammar (a "mixture") to test, for example by combinatorial motif mutagenesis. It is not evidence of a protein complex.</p>
   ${cbadge((hits.find(h => h.cbp && h.cbp !== "not_run") || hits[0] || {}).cbp)}
   ${setsTable(hits)}<h3>Constituent pairs</h3><ul>${pairs.map(p => `<li>${pairLink(p.tf_a, p.tf_b)} ${p.missing ? '<span class="tag">not individually enriched</span>' : rbadge(p.rating) + " " + cbadge(p.cbp, 1)}</li>`).join("")}</ul></div>
   ${cbpPanel(Object.fromEntries(hits.filter(h => h.cbp_e).map(h => [h.u, h.cbp_e])), (d.set_tracks || {})[[...tfs].sort().join("+")], "set")}`;
}
function viewAbout() {
  main().innerHTML = `<div class="card"><h2>Method (summary)</h2><ol>
  <li><b>Windows.</b> <i>CCLA fragments:</i> Tn5 insertions (+4/−5) are binned at 50 bp and summits of smoothed hotspots are taken as 300 bp windows. Blacklist and N filters are applied. Three window sets are built: pooled (top 20k, requiring ≥2 cuts in each replicate), rep1 and rep2. <i>Corpus:</i> per-cell-type pseudobulk cCREs (top 20k by CPM, ≥100 cells), with two random half-cell pseudoreplicates.</li>
  <li><b>Motifs.</b> All 879 JASPAR 2024 CORE vertebrate non-redundant PWMs are scanned on both strands, using exact p = 1e-4 thresholds.</li>
  <li><b>Pairs.</b> A stratified hypergeometric z-test is run with matched background strata: GC quintile × accessibility tertile × motif-density tertile × repeat-rich (soft-masked hg38). BH correction gives q, and the top 300 pairs also get stratum-preserving permutations (999). Replicate support requires lift > 1 and p < 0.05 in both rep1 and rep2. <i>Distinct</i> requires q ≤ 0.05, matched lift ≥ 1.05, pilot-style lift ≥ 1.25, ≥10 joint windows, ≥5 non-overlapping sites ≤150 bp apart, and no family, PWM-similarity (≥0.75), overlap or repeat flag.</li>
  <li><b>Sets & modules.</b> Triplets are tested against the best constituent pair (Poisson + BH). Modules are maximal cliques of replicated distinct pairs.</li>
  <li><b>Consistency & ratings.</b> Pairs are aggregated over all profiles. "broad" marks enrichment in ≥75% of profiles, which may be compositional.</li>
  <li><b>ChromBPNet evidence layer (separate).</b> For profiles with a model, candidates (top distinct pairs, cross-profile strong/moderate pairs, broad controls and synergistic sets) are scored as follows. Contribution scores (integrated gradients from dinucleotide-shuffled references on the box; DeepSHAP via <code>chrombpnet contribs_bw</code> on the cluster) are computed over accessible windows. A pair is <i>contribution-supported</i> when every motif site's median contribution percentile is ≥ 0.75. Joint vs single in-silico ablation of the sites gives a synergy score, Δjoint − ΣΔsingle on predicted log-counts. TF-MoDISco patterns are mapped to JASPAR as a de novo check.</li></ol>
  <h3>Pilot</h3><p>The earlier HepG2 chr22 pilot (500 windows; GC × insertion strata; 999 permutations) was re-run on this build box: FOXA2–HNF1A had 144 joint windows, lift 1.63 and permutation q 0.018 (pilot priority). FOS–JUN (lift 2.02), CEBPA–CEBPB (2.12) and MAX–MYC (2.42) were family-level. Genome-wide v2 runs use the stricter matched background described above, so lifts are lower.</p>
  <h3>Toolkit</h3><pre>./tfpm v2 setup && ./tfpm v2 discover
./tfpm v2 run --runner wr      # or --runner local --parallel 1
./tfpm v2 publish</pre><p><a href="https://github.com/double2strand/TFpairminer/tree/main/mining/motif_cooccurrence">mining/motif_cooccurrence/README.md</a> has the cluster and workstation steps.</p>
  <h3>Not done / caveats</h3><ul><li>ChromBPNet has been run only for the profiles listed on the overview, using ENCODE pretrained bulk-ATAC models (fold 0) rather than models trained on the CCLA fragments. Every other profile shows "not yet run".</li><li>Lifts are window-level co-occurrence and say nothing about binding.</li><li>Several CCLA lines (HCT116, MCF7, HeLa, GM12878) were unavailable (HTTP 404).</li><li>The Corpus provides only cell × cCRE matrices, not fragments, so its windows are cCREs and its replicates are pseudoreplicates.</li></ul></div>`;
}
boot();
