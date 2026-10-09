import { initPage, h, field, segmented, slider, fmt, frameThrottle } from "../ui/page.js";
import { createPlot } from "../ui/plot.js";
import { createMatrixView } from "../ui/matrixView.js";
import { joScalarAtGap, nearContact } from "../physics/joFull.js";
import { pairScalarsAtGap, resistance12, torqueFree, NAMES } from "../models/pairB.js";

initPage("p13");

const FAMILIES = {
  XA: { label: "X^A（中心線方向の並進）", names: ["XA11", "XA12"], ydomain: [-30, 30] },
  YA: { label: "Y^A（垂直方向の並進）", names: ["YA11", "YA12"] },
  YB: { label: "Y^B（並進と回転の結合）", names: ["YB11", "YB12"] },
  XC: { label: "X^C（中心線まわりの回転）", names: ["XC11", "XC12"] },
  YC: { label: "Y^C（垂直な軸まわりの回転）", names: ["YC11", "YC12"] },
};
const st = { fam: "XA", eval: "full", gap: 0.05, theta: 30, show: "B" };
const gaps = Array.from({ length: 121 }, (_, k) => 10 ** (-3 + (4 * k) / 120));

// cache both models on the plotting grid (they do not depend on the controls)
const gridA = gaps.map((g) => pairScalarsAtGap(g, "A")), gridB = gaps.map((g) => pairScalarsAtGap(g, "B"));
const gridPlain = Object.fromEntries(NAMES.map((n) => [n, gaps.map((g) => joScalarAtGap(n, g, { plain: true }))]));

// ---------------------------------------------------------------------
// functions of the gap
// ---------------------------------------------------------------------
const famPlot = createPlot(document.getElementById("fam-plot"), { width: 960, height: 320, xlog: true, xlabel: "すき間 h/a", ylabel: "無次元の抵抗関数" });
const scheduleFam = frameThrottle(() => { renderFam(); renderTable(); });
document.getElementById("fam-controls").append(
  field("関数", segmented(Object.entries(FAMILIES).map(([k, f]) => [k, f.label]), st.fam, (v) => { st.fam = v; renderFam(); })),
  field("案B の計算", segmented([["full", "特異項を分けて足す"], ["plain", "級数だけ（200 項）"]], st.eval, (v) => { st.eval = v; renderFam(); renderTable(); })));
function renderFam() {
  const f = FAMILIES[st.fam];
  const series = [];
  f.names.forEach((n, i) => {
    const color = i === 0 ? "var(--c1)" : "var(--c2)", tag = n.slice(2);
    series.push({ name: `案B ${tag}`, color, width: 3, points: gaps.map((g, k) => [g, st.eval === "plain" ? gridPlain[n][k] : gridB[k][n]]) });
    series.push({ name: `案A ${tag}`, color, dash: "6 4", points: gaps.map((g, k) => [g, gridA[k][n]]) });
    if (nearContact[n]) series.push({ name: `近接の主要項 ${tag}`, color: "var(--muted)", dash: "2 3", points: gaps.filter((g) => g < 0.3).map((g) => [g, nearContact[n](g)]) });
  });
  famPlot.update({ series, hlines: [{ y: 0 }], vlines: [{ x: st.gap }], ydomain: f.ydomain });
}

const gapSlider = slider({ label: "すき間 h/a", min: -3, max: 1, step: 0.01, value: Math.log10(st.gap), format: (v) => fmt(10 ** v, 3), onInput: (v) => { st.gap = 10 ** v; scheduleFam(); scheduleMat(); } });
document.getElementById("gap-controls").append(gapSlider);
function renderTable() {
  const A = pairScalarsAtGap(st.gap, "A"), B = pairScalarsAtGap(st.gap, "B");
  const rows = NAMES.map((n) => {
    const bv = st.eval === "plain" ? joScalarAtGap(n, st.gap, { plain: true }) : B[n];
    return h("tr", {}, h("td", {}, n), h("td", {}, fmt(A[n], 4)), h("td", {}, fmt(bv, 4)), h("td", {}, Math.abs(A[n]) > 1e-12 ? fmt(bv / A[n], 3) : "—"));
  });
  document.getElementById("fam-table").replaceChildren(h("table", { class: "data" },
    h("tr", {}, h("th", {}, `h/a = ${fmt(st.gap, 3)}`), h("th", {}, "案A"), h("th", {}, st.eval === "plain" ? "案B（級数だけ）" : "案B"), h("th", {}, "B ÷ A")), ...rows));
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

// each 3 x 3 block divided by its own unit, so that the entries are the JO scalars times direction factors
function normalised(R) {
  const u = [6 * Math.PI, 4 * Math.PI, 4 * Math.PI, 8 * Math.PI]; // A, B~, B, C blocks (a = mu = 1)
  return R.map((row, i) => row.map((v, j) => v / (i < 6 && j < 6 ? u[0] : i >= 6 && j >= 6 ? u[3] : u[1])));
}
function renderMat() {
  const th = (st.theta * Math.PI) / 180, e = [Math.sin(th), 0, Math.cos(th)];
  const RA = normalised(resistance12(pairScalarsAtGap(st.gap, "A"), e)), RB = normalised(resistance12(pairScalarsAtGap(st.gap, "B"), e));
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
// where the models part: |B - A| versus r
// ---------------------------------------------------------------------
const farPlot = createPlot(document.getElementById("far-plot"), { width: 960, height: 320, xlog: true, ylog: true, xlabel: "中心間距離 r/a", ylabel: "|案B − 案A|" });
{
  const rs = Array.from({ length: 41 }, (_, k) => 10 ** (Math.log10(2.2) + ((Math.log10(40) - Math.log10(2.2)) * k) / 40));
  const pick = [["XA11", "var(--c1)", "r⁻⁴"], ["XA12", "var(--c1)", "r⁻⁵"], ["YA11", "var(--c2)", "r⁻⁶"], ["YC11", "var(--c3)", "r⁻⁶"], ["YB11", "var(--ink)", "r⁻⁷"]];
  const sA = rs.map((r) => pairScalarsAtGap(r - 2, "A")), sB = rs.map((r) => pairScalarsAtGap(r - 2, "B"));
  farPlot.update({
    series: pick.map(([n, color, slope], i) => ({ name: `${n}（${slope}）`, color, dash: i === 1 ? "6 4" : null, points: rs.map((r, k) => [r, Math.abs(sB[k][n] - sA[k][n])]) })),
    ydomain: [1e-12, 1e2],
  });
}

// ---------------------------------------------------------------------
// rotation fixed or torque free
// ---------------------------------------------------------------------
const tfPlot = createPlot(document.getElementById("tf-plot"), { width: 960, height: 300, xlog: true, xlabel: "すき間 h/a", ylabel: "垂直方向の自己抵抗 ÷ 6πμa" });
{
  const u = 6 * Math.PI;
  const tfB = gridB.map((sc) => torqueFree(resistance12(sc, [0, 0, 1]))[0][0] / u);
  const tfA = gridA.map((sc) => torqueFree(resistance12(sc, [0, 0, 1]))[0][0] / u);
  tfPlot.update({
    series: [
      { name: "案B 回転を固定（Y^A₁₁）", color: "var(--c1)", width: 3, points: gaps.map((g, k) => [g, gridB[k].YA11]) },
      { name: "案B トルクをゼロ（自由に回る）", color: "var(--c2)", width: 3, points: gaps.map((g, k) => [g, tfB[k]]) },
      { name: "案A トルクをゼロ（ページ8〜11の RPY 逆行列）", color: "var(--c2)", dash: "6 4", points: gaps.map((g, k) => [g, tfA[k]]) },
    ],
    hlines: [{ y: 1 }],
  });
}

renderFam();
renderTable();
renderMat();
