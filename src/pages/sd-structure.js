import { initPage, h, slider, segmented, fmt } from "../ui/page.js";
import { createMatrixView } from "../ui/matrixView.js";
import { createPlot } from "../ui/plot.js";
import { matmul, transpose, solve } from "../core/linalg.js";
import { mulberry32, normalSampler } from "../core/random.js";
import { stf, stf5, ddot, naive5 } from "../physics/stf.js";
import { schurReduce, fullSolve, sdResistance, isolatedResistance11 } from "../models/sd.js";
import { axialManyBody, resistance, pairwiseSumResistance } from "../models/manyBody.js";

initPage("p11");

const a = 1, mu = 1, unit = 6 * Math.PI * mu * a;

// ---------------------------------------------------------------------
// Isolated sphere, 11 x 11
// ---------------------------------------------------------------------
const L11 = ["Ux", "Uy", "Uz", "Ωx", "Ωy", "Ωz", "e₁", "e₂", "e₃", "e₄", "e₅"];
const r0View = createMatrixView(document.getElementById("r0"), {
  labels: L11, blocks: [3, 3, 5],
  format: (v) => (v === 0 ? "0" : v.toFixed(1)),
  explain: (i, j, v) => {
    if (i !== j) return h("div", {}, h("strong", {}, `${L11[i]}, ${L11[j]}`), "：孤立球では異なる成分は結合しない（球の対称性）。0。");
    const [name, formula, page] = i < 3 ? ["並進の抵抗", "6πμa", "第4章"] : i < 6 ? ["回転の抵抗", "8πμa³", "第4章"] : ["ひずみに対するストレスレットの係数", "20πμa³/3", "第5章"];
    return h("div", {}, h("strong", {}, `${L11[i]}, ${L11[j]}：${name}`), h("div", { class: "val" }, `${formula} = ${v.toFixed(4)}（a = μ = 1）`),
      h("div", { class: "caption" }, `一球の表面力を積分して得る係数（教科書${page}）。二粒子以上では、この対角ブロックに相手の影響が加わり、非対角ブロックが現れる。`));
  },
});
r0View.update(isolatedResistance11({ a, mu }));
r0View.select(6, 6);

// ---------------------------------------------------------------------
// STF 5-component representation
// ---------------------------------------------------------------------
let stfRng = mulberry32(5);
function renderStf() {
  const randn = normalSampler(stfRng);
  const rand3 = () => [0, 1, 2].map(() => [0, 1, 2].map(() => Math.round(randn() * 10) / 10));
  const S = stf(rand3()), E = stf(rand3());
  const s = stf5(S), e = stf5(E), ns = naive5(S), ne = naive5(E);
  const dot = (x, y) => x.reduce((acc, v, k) => acc + v * y[k], 0);
  const row = (label, v) => h("tr", {}, h("td", {}, label), ...v.map((x) => h("td", {}, fmt(x, 3))));
  document.getElementById("stf").replaceChildren(
    h("table", { class: "data" },
      h("tr", {}, h("th", {}, ""), ...[1, 2, 3, 4, 5].map((k) => h("th", {}, `${k}`))),
      row("s（正規直交基底）", s), row("e（正規直交基底）", e),
      row("S の 5 成分を抜き出しただけ", ns), row("E の 5 成分を抜き出しただけ", ne)),
    h("table", { class: "data", style: { marginTop: "8px" } },
      h("tr", {}, h("td", {}, "S:E（二重縮約）"), h("td", { class: "val" }, fmt(ddot(S, E), 6))),
      h("tr", {}, h("td", {}, "s·e"), h("td", { class: "val" }, fmt(dot(s, e), 6))),
      h("tr", {}, h("td", {}, "抜き出した 5 成分どうしの内積"), h("td", { class: "val", style: { color: "var(--warn)" } }, fmt(dot(ns, ne), 6)))));
}
document.getElementById("stf-new").addEventListener("click", renderStf);
renderStf();

// ---------------------------------------------------------------------
// Schur complement on a generic SPD generalised mobility (algebra check)
// ---------------------------------------------------------------------
const NV = 6, NS = 5;
const base = (() => {
  const randn = normalSampler(mulberry32(14));
  const B = Array.from({ length: NV + NS }, () => Array.from({ length: NV + NS }, randn));
  return matmul(B, transpose(B)).map((r, i) => r.map((v, j) => (v + (i === j ? 2 : 0)) / 10));
})();
const gv = [0, 0, -1, 0, 0, 0]; // a force along -z, no torque
const schur = { eps: 1, einf: false };
// Scaling the off-diagonal blocks by eps in [0, 1] keeps the matrix positive definite.
const coupled = () => base.map((r, i) => r.map((v, j) => ((i < NV) !== (j < NV) ? schur.eps * v : v)));
const schurView = createMatrixView(document.getElementById("schur-matrix"), {
  labels: ["Ux", "Uy", "Uz", "Ωx", "Ωy", "Ωz", "e₁", "e₂", "e₃", "e₄", "e₅"], blocks: [6, 5],
  explain: (i, j, v) => {
    const part = (k) => (k < NV ? "v" : "s");
    return h("div", {}, h("strong", {}, `M_${part(i)}${part(j)} のブロック`), h("div", { class: "val" }, fmt(v)),
      h("div", { class: "caption" }, part(i) !== part(j) ? `運動とひずみを結ぶ成分。結合の強さ ε = ${schur.eps.toFixed(2)} を掛けている。ε = 0 なら誘起ストレスレットは生じない。` : "例として作った対称正定値行列の成分（物理的な係数ではない）。"));
  },
});
document.getElementById("schur-controls").append(
  slider({ label: "結合の強さ ε（M_vs, M_sv に掛ける）", min: 0, max: 1, step: 0.01, value: 1, format: (v) => v.toFixed(2), onInput: (v) => { schur.eps = v; renderSchur(); } }),
  h("label", {}, "背景のひずみ e∞", segmented([[false, "なし"], [true, "あり"]], schur.einf, (v) => { schur.einf = v; renderSchur(); })));
function renderSchur() {
  const M = coupled();
  const einf = schur.einf ? [0.5, 0, 0.3, 0, 0] : [0, 0, 0, 0, 0];
  const red = schurReduce(M, NV, { gv, einf }), full = fullSolve(M, NV, { gv, einf });
  schurView.update(M);
  const diff = Math.max(...red.qv.map((x, i) => Math.abs(x - full.qv[i])));
  const gsNorm = Math.hypot(...red.gs);
  document.getElementById("schur-out").replaceChildren(
    h("table", { class: "data" },
      h("tr", {}, h("th", {}, ""), ...[1, 2, 3, 4, 5].map((k) => h("th", {}, `${k}`))),
      h("tr", {}, h("td", {}, "g_s"), ...red.gs.map((x) => h("td", {}, fmt(x, 3))))),
    h("p", {}, `|g_s| = ${fmt(gsNorm, 3)}。`, gsNorm < 1e-12 ? "結合がなく背景のひずみもないので、ストレスレットは誘起されない。" : "加えたのは −z 方向の力だけだが、球が変形しないための表面力がストレスレットとして現れる。"),
    h("p", { class: "caption" }, `Schur 補行列で消去した結果と、元の連立方程式を直接解いた結果の差：最大 ${fmt(diff, 2)}（丸め誤差の範囲）。運動の応答 q_v の z 成分は ${fmt(red.qv[2], 4)}（結合なしなら ${fmt(-M[2][2], 4)}）。`));
}
renderSchur();

// ---------------------------------------------------------------------
// SD-type resistance for three coaxial spheres
// ---------------------------------------------------------------------
const sd = { g12: 1, g23: 1, L: 24 };
const pairCache = new Map();
// axial (zz) components of a 3N x 3N translational matrix
const zz = (R) => { const idx = Array.from({ length: R.length / 3 }, (_, k) => 3 * k + 2); return idx.map((i) => idx.map((j) => R[i][j])); };
function methodsFor(z, L) {
  if (pairCache.size > 500) pairCache.clear();
  const ex = axialManyBody(z, z.map(() => a), { L, mu, pairCache });
  const X = z.map((v) => [0, 0, v]);
  const Rr = zz(resistance(X, { a, mu })), R2r = zz(pairwiseSumResistance(X, { a, mu }));
  return {
    residual: ex.boundaryError,
    list: [
      { name: "厳密（境界条件解法）", R: ex.R },
      { name: "RPY の逆行列（遠方だけ）", R: Rr },
      { name: "二体の和", R: ex.R2B },
      { name: "SD 型", R: sdResistance(Rr, ex.R2B, R2r) },
    ],
  };
}
const maxErr = (R, Rex) => Math.max(...R.flatMap((row, i) => row.map((v, j) => Math.abs(v - Rex[i][j])))) / Math.abs(Rex[0][0]);

let sdPending = false;
const scheduleSd = () => { if (!sdPending) { sdPending = true; requestAnimationFrame(() => { sdPending = false; renderSd(); }); } };
document.getElementById("sd-controls").append(
  slider({ label: "隙間 h₁₂/a", min: 0.1, max: 4, step: 0.05, value: sd.g12, format: (v) => v.toFixed(2), onInput: (v) => { sd.g12 = v; scheduleSd(); } }),
  slider({ label: "隙間 h₂₃/a", min: 0.1, max: 4, step: 0.05, value: sd.g23, format: (v) => v.toFixed(2), onInput: (v) => { sd.g23 = v; scheduleSd(); } }),
  h("label", {}, "次数 L", segmented([[16, "16"], [24, "24"], [32, "32"]], sd.L, (v) => { sd.L = v; scheduleSd(); })));

function renderSd() {
  const z = [0, 2 * a + sd.g12, 4 * a + sd.g12 + sd.g23];
  const { residual, list } = methodsFor(z, sd.L);
  const Rex = list[0].R;
  const entries = [["R₁₁", 0, 0], ["R₂₂", 1, 1], ["R₁₂", 0, 1], ["R₂₃", 1, 2], ["R₁₃", 0, 2]];
  const speed = (R) => solve(R, [-1, -1, -1]).map((u) => -u * unit);
  const table = h("table", { class: "data" },
    h("tr", {}, h("th", {}, ""), ...entries.map(([lab]) => h("th", {}, lab)), h("th", {}, "最大誤差"), h("th", {}, "球 2 の沈降速度")));
  for (const m of list) {
    const isExact = m === list[0];
    table.append(h("tr", {}, h("td", { style: { whiteSpace: "nowrap" } }, m.name),
      ...entries.map(([, i, j]) => {
        const v = m.R[i][j] / unit, wrongSign = !isExact && Math.sign(m.R[i][j]) !== Math.sign(Rex[i][j]);
        return h("td", { style: { color: wrongSign ? "var(--warn)" : "" } }, fmt(v));
      }),
      h("td", {}, isExact ? "—" : `${(100 * maxErr(m.R, Rex)).toFixed(1)}%`),
      h("td", {}, fmt(speed(m.R)[1]))));
  }
  document.getElementById("sd-table").replaceChildren(h("div", { style: { overflowX: "auto" } }, table));
  const sdR = list[3].R;
  const notes = [`境界残差の最大値 ${fmt(residual, 2)}（L = ${sd.L}）。沈降速度は孤立球の値 m₀F を単位とする。`];
  if (Math.sign(sdR[0][2]) !== Math.sign(Rex[0][2]))
    notes.push("SD 型の R₁₃ は符号が厳密解と逆になっている（橙色）。遠方の部分を並進の力だけで作り、ストレスレットを含めていないためである。球 1 と球 3 の結合は、間の球を経由する遠方の反射で決まるので、近接の二体補正では直らない。");
  document.getElementById("sd-note").textContent = notes.join(" ");
}
renderSd();

const sdPlot = createPlot(document.getElementById("sd-plot"), { height: 280, xlog: true, ylog: true, xlabel: "隙間 h/a（二つとも同じ）", ylabel: "最大相対誤差" });
document.getElementById("sd-scan").addEventListener("click", (ev) => {
  const btn = ev.currentTarget, status = document.getElementById("sd-scan-status");
  btn.disabled = true;
  status.textContent = "計算中…";
  setTimeout(() => {
    try {
      const gaps = [4, 2.5, 1.5, 1, 0.6, 0.4, 0.25, 0.15, 0.1];
      const rows = gaps.map((g) => {
        const { list } = methodsFor([0, 2 + g, 4 + 2 * g], 32);
        return { g, errs: list.slice(1).map((m) => maxErr(m.R, list[0].R)) };
      });
      const colors = ["var(--c1)", "var(--c2)", "var(--c3)"];
      const names = ["RPY の逆行列（遠方だけ）", "二体の和", "SD 型"];
      sdPlot.update({ series: names.map((name, k) => ({ name, color: colors[k], points: rows.map((r) => [r.g, r.errs[k]]) })) });
      status.textContent = "完了（L = 32）";
    } finally {
      btn.disabled = false;
    }
  }, 30);
});
