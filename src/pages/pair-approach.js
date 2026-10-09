import { initPage, h, field, tex, segmented, slider, fmt } from "../ui/page.js";
import { createMatrixView } from "../ui/matrixView.js";
import { createPlot } from "../ui/plot.js";
import { joPolynomials, joXA } from "../physics/jo.js";
import { axialCollocation } from "../physics/collocation.js";
import { resistanceCoefficients, pairMatrix } from "../models/pairA.js";
import { approachTrajectory } from "../models/approach.js";

initPage("p9");

const a = 1, mu = 1, F = 1;
const unit = 6 * Math.PI * mu * a;
const Tmax = () => 120; // simulated time window
const polys60 = joPolynomials(60);

const state = { h0: 2, hc: 0.2, K: 40, L: 32, matrixModel: "lub", tFrac: 0, colloc: [] };

const MODELS = {
  rpy: { name: "RPY の逆行列", color: "var(--c1)" },
  lub: { name: "RPY ＋ 教材潤滑", color: "var(--c2)" },
  jo: { name: "JO 級数（軸方向）", color: "var(--c3)" },
};

function coeffs(model, h) {
  const r = 2 * a + h;
  const opts = { a, mu };
  if (model === "lub") Object.assign(opts, { lubrication: true, hc: state.hc * a });
  if (model === "jo") Object.assign(opts, { axial: "jo", polys: polys60.slice(0, state.K + 1) });
  return resistanceCoefficients(r, opts);
}
// X11 - X12 (normalised) and whether the value can be trusted
function relative(model, h) {
  const c = coeffs(model, h);
  return { value: (c.R.par.self - c.R.par.cross) / unit, ok: model !== "jo" || c.jo.converged, c };
}

// Gap trajectories, cached per model; a control invalidates only the models it affects.
function trajectory(model) {
  return approachTrajectory((hh) => relative(model, hh), {
    h0: state.h0, Tmax: Tmax(), F, a, mu,
    // RPY alone reaches contact in finite time; the lubricated gap only decays.
    contactBelow: model === "rpy" ? 1e-9 : undefined,
  });
}

// ---------- controls ----------
const controls = document.getElementById("controls");
controls.append(
  slider({ label: "潤滑の切替え h_c/a", min: 0.05, max: 0.5, step: 0.01, value: state.hc, format: (v) => v.toFixed(2), onInput: (v) => { state.hc = v; update(["lub"]); } }),
  field("JO の次数 K", segmented([[10, "10"], [20, "20"], [40, "40"], [60, "60"]], state.K, (v) => { state.K = v; update(["jo"]); renderColloc(); })));

const cc = document.getElementById("colloc-controls");
const collocStatus = h("span", { class: "caption" });
const collocButton = h("button", { class: "action", type: "button", onclick: () => runCollocation() }, "計算する");
cc.append(
  field("次数 L", segmented([[16, "16"], [32, "32"], [64, "64"]], state.L, (v) => { state.L = v; })),
  collocButton, collocStatus);

const mc = document.getElementById("motion-controls");
const tSlider = slider({ label: "時刻", min: 0, max: 1, step: 0.001, value: 0, format: (v) => `t = ${(v * Tmax()).toFixed(1)}`, onInput: (v) => { state.tFrac = v; renderMatrix(); renderGap(); } });
mc.append(
  slider({ label: "最初の隙間 h₀/a", min: 0.3, max: 4, step: 0.1, value: state.h0, format: (v) => v.toFixed(1), onInput: (v) => { state.h0 = v; renderMotion(); } }),
  tSlider);
document.getElementById("matrix-controls").append(
  field("モデル", segmented(Object.entries(MODELS).map(([k, m]) => [k, m.name]), state.matrixModel, (v) => { state.matrixModel = v; renderMatrix(); renderGap(); })));

const resPlot = createPlot(document.getElementById("resplot"), { height: 320, xlog: true, ylog: true, xlabel: "隙間 h/a", ylabel: "X₁₁ − X₁₂" });
const gapPlot = createPlot(document.getElementById("gapplot"), { height: 280, ylog: true, xlabel: "時間 t（μa²/F）", ylabel: "隙間 h/a" });
const diffPlot = createPlot(document.getElementById("diffplot"), { height: 260, xlog: true, ylog: true, xlabel: "中心間距離 s = r/a", ylabel: "|X₁₁(RPY⁻¹) − X^A₁₁(JO)|" });
const LAB = ["1x", "1y", "1z", "2x", "2y", "2z"];
const matrixView = createMatrixView(document.getElementById("matrix"), { labels: LAB, explain: explainEntry });
matrixView.select(0, 3);

// ---------- resistance vs gap ----------
const hgrid = Array.from({ length: 160 }, (_, i) => 10 ** (-3 + (i * 4) / 159));
// One cached curve per model, replaced when its parameter (h_c or K) changes.
const curveCache = new Map();
function curve(model) {
  const key = model === "lub" ? state.hc : model === "jo" ? state.K : "";
  const hit = curveCache.get(model);
  if (hit && hit.key === key) return hit.data;
  const data = hgrid.map((hh) => ({ hh, ...relative(model, hh) }));
  curveCache.set(model, { key, data });
  return data;
}
function renderResistance() {
  const series = [];
  for (const model of ["rpy", "lub"]) series.push({ name: MODELS[model].name, color: MODELS[model].color, points: curve(model).map((p) => [p.hh, p.value]) });
  const jo = curve("jo");
  series.push({ name: `JO 級数 K = ${state.K}（収束範囲）`, color: MODELS.jo.color, width: 2.5, points: jo.map((p) => [p.hh, p.ok ? p.value : NaN]) });
  series.push({ name: "同（未収束）", color: MODELS.jo.color, width: 1, dash: "3 3", points: jo.map((p) => [p.hh, p.ok ? NaN : p.value]) });
  series.push({ name: "潤滑の主要項 a/(2h)", color: "var(--muted)", dash: "6 4", width: 1.5, points: hgrid.map((hh) => [hh, a / (2 * hh)]) });
  series.push({ name: `境界条件解法（L = ${state.colloc[0]?.L ?? state.L}）`, color: "var(--ink)", marker: true, r: 4, points: state.colloc.map((c) => [c.h, c.value]) });
  resPlot.update({ series, ydomain: [0.5, 2e3], xdomain: [1e-3, 10], vlines: [{ x: state.hc, color: MODELS.lub.color }] });
}

let collocBusy = false;
function runCollocation() {
  if (collocBusy) return; // one solve at a time; a long L = 64 run must not be queued twice
  collocBusy = true;
  collocButton.disabled = true;
  collocStatus.textContent = "計算中…";
  setTimeout(() => {
    const L = state.L, rows = [];
    const t0 = performance.now();
    for (const hh of [2, 1, 0.5, 0.2, 0.1, 0.05]) {
      const res = axialCollocation([0, 2 * a + hh], [a, a], { L, withCondition: false });
      rows.push({ h: hh, L, value: (res.R[0][0] - res.R[0][1]) / unit, x11: res.R[0][0] / unit, x12: res.R[0][1] / unit, residual: res.boundaryError });
    }
    state.colloc = rows;
    collocBusy = false;
    collocButton.disabled = false;
    collocStatus.textContent = `完了（${((performance.now() - t0) / 1000).toFixed(1)} 秒）`;
    renderColloc();
    renderResistance();
  }, 30);
}
function renderColloc() {
  const box = document.getElementById("colloc");
  if (!state.colloc.length) { box.replaceChildren(h("p", { class: "caption" }, "「計算する」を押すと、隙間 h/a = 2, 1, 0.5, 0.2, 0.1, 0.05 で二球の境界値問題を解く（L = 64 では数秒かかる）。")); return; }
  const table = h("table", { class: "data" }, h("tr", {}, h("th", {}, "h/a"), h("th", {}, "X₁₁"), h("th", {}, "X₁₂"), h("th", {}, "JO（K）との差"), h("th", {}, "境界残差")));
  for (const c of state.colloc) {
    const jo = relative("jo", c.h);
    table.append(h("tr", {}, h("td", {}, String(c.h)), h("td", {}, fmt(c.x11, 7)), h("td", {}, fmt(c.x12, 7)),
      h("td", {}, jo.ok ? fmt(Math.abs(jo.value - c.value) / c.value, 2) : "未収束"),
      h("td", { style: { color: c.residual > 1e-3 ? "var(--warn)" : "" } }, fmt(c.residual, 2))));
  }
  box.replaceChildren(table);
}

// ---------- motion ----------
const traj = {};
function renderMotion(models = Object.keys(MODELS)) {
  for (const m of models) traj[m] = trajectory(m);
  renderGap();
  renderMatrix();
}
function hAt(model, t) {
  const { ts, hs } = traj[model];
  if (t >= ts[ts.length - 1]) return { h: hs[hs.length - 1], past: true };
  let i = 1;
  while (ts[i] < t) i++;
  const w = (t - ts[i - 1]) / (ts[i] - ts[i - 1]);
  return { h: hs[i - 1] + w * (hs[i] - hs[i - 1]), past: false };
}
function renderGap() {
  const t = state.tFrac * Tmax();
  const series = Object.entries(MODELS).map(([k, m]) => ({ name: m.name, color: m.color, points: traj[k].ts.map((tt, i) => [tt, traj[k].hs[i]]) }));
  const K = (2 * F) / (3 * Math.PI * mu * a * a);
  series.push({ name: "h_c 以下での主要項 h ∝ e^{−Kt}", color: "var(--muted)", dash: "6 4", width: 1.2, points: lubTail(K) });
  const marks = Object.entries(MODELS).map(([k, m]) => ({ x: Math.min(t, traj[k].ts.at(-1)), y: hAt(k, t).h, color: m.color }));
  // lower bound: the decade below the smallest gap reached by any model (at most 1e-9)
  const hmin = Math.min(...Object.values(traj).flatMap((tr) => tr.hs));
  const floor = Math.min(1e-9, 10 ** Math.floor(Math.log10(hmin)));
  gapPlot.update({ series, points: marks, ydomain: [floor, Math.max(5, state.h0 * 1.5)], xdomain: [0, Tmax()], vlines: [{ x: t }] });
  const notes = [];
  if (traj.rpy.stop === "contact") notes.push(`RPY だけでは t ≈ ${traj.rpy.ts.at(-1).toFixed(1)} で接触する（接近速度が h → 0 でも有限）。`);
  if (traj.jo.stop === "series") notes.push(`JO 級数（K = ${state.K}）は h ≈ ${traj.jo.hs.at(-1).toPrecision(2)}a で収束しなくなるので、そこで止めた。`);
  notes.push("教材潤滑を加えると、h_c より近くで隙間は指数的に減り、有限時間では接触しない。");
  document.getElementById("gap-caption").textContent = notes.join(" ");
}
// reference line: starting where the lubricated trajectory crosses h_c
function lubTail(K) {
  const { ts, hs } = traj.lub;
  const i = hs.findIndex((x) => x < state.hc);
  if (i < 0) return [];
  const t0 = ts[i], h0 = hs[i];
  return Array.from({ length: 50 }, (_, j) => { const tt = t0 + (j * (Tmax() - t0)) / 49; return [tt, h0 * Math.exp(-K * (tt - t0))]; });
}

function renderMatrix() {
  const t = state.tFrac * Tmax();
  const m = state.matrixModel;
  const { h: hh, past } = hAt(m, t);
  const c = coeffs(m, hh);
  const norm = { par: { self: c.R.par.self / unit, cross: c.R.par.cross / unit }, perp: { self: c.R.perp.self / unit, cross: c.R.perp.cross / unit } };
  current = { c, hh, past, m };
  document.getElementById("matrix-title").textContent = `抵抗行列 R/(6πμa)、h = ${fmt(hh, 3)}a${past ? "（計算を止めた後）" : ""}`;
  matrixView.update(pairMatrix(norm, [1, 0, 0]), { labels: LAB });
}
let current = null;

function explainEntry(i, j, v) {
  const { c, hh, m } = current;
  const k = i % 3, l = j % 3, same = Math.floor(i / 3) === Math.floor(j / 3);
  const wrap = h("div", {}, h("div", {}, h("strong", {}, `${LAB[i]}, ${LAB[j]} 成分 = ${fmt(v)}`), `（h = ${fmt(hh, 3)}a）`));
  if (k !== l) { wrap.append(h("div", {}, "中心線（x 軸）方向と垂直方向は結合しない。0。")); return wrap; }
  if (k === 0) {
    wrap.append(tex(same ? "X_{11}" : "X_{12}", { display: true }), h("div", { class: "caption" }, `出典：${c.source.par}`));
    if (m === "lub" && c.zeta > 0) wrap.append(h("div", { class: "val" }, `教材潤滑 ζₙ/(6πμa) = ${fmt(c.zeta / unit)} を ${same ? "加えた" : "引いた"}。`));
    if (m === "jo") wrap.append(h("div", { class: "caption" }, c.jo.converged ? `K = ${c.jo.K} 次で収束（${c.jo.Kc} 次との相対差 ${fmt(c.jo.change, 2)}）。` : "級数が収束していない。この値は使えない。"));
  } else {
    wrap.append(tex(same ? "Y_{11}" : "Y_{12}", { display: true }),
      h("div", { class: "caption" }, `出典：${c.source.perp}。案Aでは垂直方向の厳密な関数がないので、どのモデルでも RPY の値を使っている。`));
  }
  return wrap;
}

// ---------- RPY vs JO difference ----------
function renderDiff() {
  const ss = Array.from({ length: 80 }, (_, i) => 10 ** (Math.log10(2.5) + (i * (Math.log10(60) - Math.log10(2.5))) / 79));
  const pts = ss.map((s) => {
    const rpy = resistanceCoefficients(s * a, { a, mu }).R.par.self / unit;
    return [s, Math.abs(rpy - joXA(s, 1, polys60).x11)];
  });
  const ref = pts.find((p) => p[0] > 20);
  const C = ref[1] * ref[0] ** 4;
  diffPlot.update({
    series: [
      { name: "差", color: "var(--c3)", points: pts },
      { name: "傾き −4 の線", color: "var(--muted)", dash: "6 4", width: 1.2, points: ss.map((s) => [s, C / s ** 4]) },
    ],
  });
}

// Recompute only what depends on the changed models.
function update(models) { renderResistance(); renderMotion(models); }
renderResistance();
renderColloc();
renderMotion();
renderDiff();
