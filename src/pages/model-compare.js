import { initPage, h, field, segmented, slider, fmt, frameThrottle } from "../ui/page.js";
import { createPlot } from "../ui/plot.js";
import { createMatrixView } from "../ui/matrixView.js";
import { joScalarAtGap, nearContactAt } from "../physics/joFull.js";
import { pairScalarsUnequal, resistance12, torqueFree, NAMES16 } from "../models/pairB.js";

initPage("p13");

// sphere 1 has radius 1; sphere 2 has radius lambda (the size ratio a2/a1)
const RATIOS = [[1, "1 : 1"], [0.5, "1 : 1/2"], [0.25, "1 : 1/4"]];
const FAMILIES = {
  XA: { label: "X^A（中心線方向の並進）", names: ["XA11", "XA12", "XA22"], ydomain: [-30, 30] },
  YA: { label: "Y^A（垂直方向の並進）", names: ["YA11", "YA12", "YA22"] },
  YB: { label: "Y^B（並進と回転の結合）", names: ["YB11", "YB12", "YB21", "YB22"] },
  XC: { label: "X^C（中心線まわりの回転）", names: ["XC11", "XC12", "XC22"] },
  YC: { label: "Y^C（垂直な軸まわりの回転）", names: ["YC11", "YC12", "YC22"] },
};
const COLOR = { 11: "var(--c1)", 12: "var(--c2)", 21: "var(--c4)", 22: "var(--c3)" };
const st = { lambda: 1, fam: "XA", eval: "full", gap: 0.05, theta: 30, show: "B" };
const gaps = Array.from({ length: 121 }, (_, k) => 10 ** (-3 + (4 * k) / 120));
const radii = () => ({ a1: 1, a2: st.lambda });
const shown = (names) => (st.lambda === 1 ? names.filter((n) => !n.endsWith("22") && !n.endsWith("21")) : names);

// model-B value of any of the 16 scalars, plain series or full; "22"/"21" come from 1/lambda (JO 1.9)
function plainValue(n, xi) {
  const l = st.lambda, k = n.slice(2), base = n.slice(0, 2);
  if (k === "11" || k === "12") return joScalarAtGap(n, xi, { plain: true, lambda: l });
  const v = joScalarAtGap(base + (k === "22" ? "11" : "12"), xi, { plain: true, lambda: 1 / l });
  return base === "YB" ? -v : v;
}
function leadValue(n, xi) {
  const l = st.lambda, k = n.slice(2), base = n.slice(0, 2);
  if (base === "XC") return NaN;
  if (k === "11" || k === "12") return nearContactAt(n, xi, l);
  const v = nearContactAt(base + (k === "22" ? "11" : "12"), xi, 1 / l);
  return base === "YB" ? -v : v;
}

// cached curves for the current ratio
let grid = null;
function curves() {
  if (grid && grid.lambda === st.lambda) return grid;
  const A = gaps.map((g) => pairScalarsUnequal(g, "A", radii())), B = gaps.map((g) => pairScalarsUnequal(g, "B", radii()));
  grid = { lambda: st.lambda, A, B, plain: {} };
  return grid;
}
const plainCurve = (n) => (grid.plain[n] ??= gaps.map((g) => plainValue(n, g)));

// scalars at the gap chosen with the slider, computed once per gap and ratio
let atGap = null;
function scalarsAtGap() {
  if (!atGap || atGap.gap !== st.gap || atGap.lambda !== st.lambda)
    atGap = { gap: st.gap, lambda: st.lambda, A: pairScalarsUnequal(st.gap, "A", radii()), B: pairScalarsUnequal(st.gap, "B", radii()), plain: null };
  if (st.eval === "plain" && !atGap.plain) atGap.plain = Object.fromEntries(NAMES16.map((n) => [n, plainValue(n, st.gap)]));
  return atGap;
}

// ---------------------------------------------------------------------
// functions of the gap
// ---------------------------------------------------------------------
const famPlot = createPlot(document.getElementById("fam-plot"), { width: 960, height: 320, xlog: true, xlabel: "すき間 ξ = 2h/(a₁+a₂)", ylabel: "無次元の抵抗関数" });
const scheduleFam = frameThrottle(() => { renderFam(); renderTable(); });
document.getElementById("ratio-controls").append(
  field("半径の比 a₁ : a₂", segmented(RATIOS, st.lambda, (v) => { st.lambda = v; renderAll(); })));
document.getElementById("fam-controls").append(
  field("関数", segmented(Object.entries(FAMILIES).map(([k, f]) => [k, f.label]), st.fam, (v) => { st.fam = v; renderFam(); })),
  field("案B の計算", segmented([["full", "特異項を分けて足す"], ["plain", "級数だけ（200 項）"]], st.eval, (v) => { st.eval = v; renderFam(); renderTable(); })));
function renderFam() {
  const c = curves(), f = FAMILIES[st.fam];
  const series = [];
  shown(f.names).forEach((n) => {
    const tag = n.slice(2), color = COLOR[tag];
    series.push({ name: `案B ${tag}`, color, width: 3, points: gaps.map((g, k) => [g, st.eval === "plain" ? plainCurve(n)[k] : c.B[k][n]]) });
    series.push({ name: `案A ${tag}`, color, dash: "6 4", points: gaps.map((g, k) => [g, c.A[k][n]]) });
    if (st.fam !== "XC") series.push({ color: "var(--muted)", dash: "2 3", points: gaps.filter((g) => g < 0.3).map((g) => [g, leadValue(n, g)]) });
  });
  series.push({ name: "近接の主要項", color: "var(--muted)", dash: "2 3", points: [] });
  famPlot.update({ series, hlines: [{ y: 0 }], vlines: [{ x: st.gap }], ydomain: f.ydomain });
}

const gapSlider = slider({ label: "すき間 ξ", min: -3, max: 1, step: 0.01, value: Math.log10(st.gap), format: (v) => fmt(10 ** v, 3), onInput: (v) => { st.gap = 10 ** v; scheduleFam(); scheduleMat(); } });
document.getElementById("gap-controls").append(gapSlider);
function renderTable() {
  const { A, B, plain } = scalarsAtGap();
  const rows = shown(NAMES16).map((n) => {
    const bv = st.eval === "plain" ? plain[n] : B[n];
    return h("tr", {}, h("td", {}, n), h("td", {}, fmt(A[n], 4)), h("td", {}, fmt(bv, 4)), h("td", {}, Math.abs(A[n]) > 1e-12 ? fmt(bv / A[n], 3) : "—"));
  });
  document.getElementById("fam-table").replaceChildren(h("table", { class: "data" },
    h("tr", {}, h("th", {}, `ξ = ${fmt(st.gap, 3)}`), h("th", {}, "案A"), h("th", {}, st.eval === "plain" ? "案B（級数だけ）" : "案B"), h("th", {}, "B ÷ A")), ...rows));
}

// ---------------------------------------------------------------------
// the 12 x 12 matrix
// ---------------------------------------------------------------------
const LABELS = ["U₁x", "U₁y", "U₁z", "U₂x", "U₂y", "U₂z", "Ω₁x", "Ω₁y", "Ω₁z", "Ω₂x", "Ω₂y", "Ω₂z"];
const ROWNAME = ["F₁", "F₁", "F₁", "F₂", "F₂", "F₂", "L₁", "L₁", "L₁", "L₂", "L₂", "L₂"];
let lastMats = null;
const matView = createMatrixView(document.getElementById("mat"), {
  labels: LABELS, blocks: [3, 3, 3, 3],
  format: (v) => (Math.abs(v) < 5e-4 ? "0" : Math.abs(v) >= 100 ? v.toFixed(0) : v.toFixed(2)),
  explain: (i, j, v) => explainEntry(i, j, v),
});
const scheduleMat = frameThrottle(() => renderMat());
document.getElementById("mat-controls").append(
  slider({ label: "対の向き θ（z 軸から）", min: 0, max: 90, step: 1, value: st.theta, format: (v) => `${v}°`, onInput: (v) => { st.theta = v; scheduleMat(); } }),
  field("表示", segmented([["A", "案A"], ["B", "案B"], ["D", "差 B − A"]], st.show, (v) => { st.show = v; renderMat(); })));

// each 3 x 3 block divided by its JO unit (3 pi mu (a_a + a_b), pi mu (a_a + a_b)^2, pi mu (a_a + a_b)^3),
// so that the entries are the scalars times direction factors
function normalised(R) {
  const rad = [1, st.lambda];
  const sphere = (k) => Math.floor(k / 3) % 2;           // 0 or 1
  const power = (i, j) => (i < 6 && j < 6 ? 1 : i >= 6 && j >= 6 ? 3 : 2);
  return R.map((row, i) => row.map((v, j) => {
    const s = rad[sphere(i)] + rad[sphere(j)], p = power(i, j);
    return v / ((p === 1 ? 3 : 1) * Math.PI * s ** p);
  }));
}
function renderMat() {
  const th = (st.theta * Math.PI) / 180, e = [Math.sin(th), 0, Math.cos(th)];
  const { A, B } = scalarsAtGap();
  const RA = normalised(resistance12(A, e, radii())), RB = normalised(resistance12(B, e, radii()));
  const RD = RB.map((row, i) => row.map((v, j) => v - RA[i][j]));
  lastMats = { RA, RB, RD, e };
  const M = st.show === "A" ? RA : st.show === "B" ? RB : RD;
  matView.update(M, { vmax: Math.max(...M.flat().map(Math.abs), 1e-12) });
}
function explainEntry(i, j, v) {
  const { RA, RB } = lastMats;
  const kind = i < 6 && j < 6 ? "A（力と並進速度）" : i >= 6 && j >= 6 ? "C（トルクと角速度）" : "B（トルクと並進、または力と回転）";
  const pa = Math.floor(i / 3) % 2 + 1, pb = Math.floor(j / 3) % 2 + 1;
  return h("div", {}, h("strong", {}, `${ROWNAME[i]}${"xyz"[i % 3]} ← ${LABELS[j]}：${fmt(v, 4)}`),
    h("div", { class: "caption" }, `ブロック ${kind}、粒子 ${pa} と ${pb}。案A ${fmt(RA[i][j], 4)}、案B ${fmt(RB[i][j], 4)}。` +
      (kind.startsWith("B") ? "この成分は ε·e の形で Y^B に比例する。" : "X（中心線方向）と Y（垂直方向）の重ね合わせ。")));
}

// ---------------------------------------------------------------------
// where the models part: |B - A| versus s
// ---------------------------------------------------------------------
const farPlot = createPlot(document.getElementById("far-plot"), { width: 960, height: 320, xlog: true, ylog: true, xlabel: "中心間距離 s = 2r/(a₁+a₂)", ylabel: "|案B − 案A|" });
function renderFar() {
  const ss = Array.from({ length: 41 }, (_, k) => 10 ** (Math.log10(2.2) + ((Math.log10(40) - Math.log10(2.2)) * k) / 40));
  const pick = [["XA11", "var(--c1)", "r⁻⁴"], ["XA12", "var(--c1)", "r⁻⁵"], ["YA11", "var(--c2)", "r⁻⁶"], ["YC11", "var(--c3)", "r⁻⁶"], ["YB11", "var(--ink)", "r⁻⁷"]];
  if (st.lambda !== 1) pick.push(["XA22", "var(--c4)", "r⁻⁴"]);
  const sA = ss.map((s) => pairScalarsUnequal(s - 2, "A", radii())), sB = ss.map((s) => pairScalarsUnequal(s - 2, "B", radii()));
  farPlot.update({
    series: pick.map(([n, color, slope], i) => ({ name: `${n}（${slope}）`, color, dash: i === 1 ? "6 4" : null, points: ss.map((s, k) => [s, Math.abs(sB[k][n] - sA[k][n])]) })),
    ydomain: [1e-12, 1e2],
  });
}

// ---------------------------------------------------------------------
// rotation fixed or torque free
// ---------------------------------------------------------------------
const tfPlot = createPlot(document.getElementById("tf-plot"), { width: 960, height: 300, xlog: true, xlabel: "すき間 ξ", ylabel: "垂直方向の自己抵抗 ÷ 6πμa_α" });
function renderTf() {
  const c = curves(), e = [0, 0, 1];
  const tf = (sc) => torqueFree(resistance12(sc, e, radii()));
  const series = [];
  const spheres = st.lambda === 1 ? [0] : [0, 1];
  for (const p of spheres) {
    const k = 3 * p, unit = 6 * Math.PI * (p === 0 ? 1 : st.lambda), tag = p === 0 ? "球 1" : "球 2（小）";
    series.push({ name: `案B 回転を固定（${tag}）`, color: p === 0 ? "var(--c1)" : "var(--c3)", width: 3, points: gaps.map((g, i) => [g, p === 0 ? c.B[i].YA11 : c.B[i].YA22]) });
    series.push({ name: `案B トルクをゼロ（${tag}）`, color: p === 0 ? "var(--c2)" : "var(--c4)", width: 3, points: gaps.map((g, i) => [g, tf(c.B[i])[k][k] / unit]) });
    series.push({ name: `案A トルクをゼロ（${tag}）`, color: p === 0 ? "var(--c2)" : "var(--c4)", dash: "6 4", points: gaps.map((g, i) => [g, tf(c.A[i])[k][k] / unit]) });
  }
  tfPlot.update({ series, hlines: [{ y: 1 }] });
}

function renderAll() {
  atGap = null;
  renderFam(); renderTable(); renderMat(); renderFar(); renderTf();
}
renderAll();
