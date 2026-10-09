import { initPage, h, field, segmented, slider, fmt } from "../ui/page.js";
import { createMatrixView } from "../ui/matrixView.js";
import { createPlot } from "../ui/plot.js";
import { createScene } from "../ui/scene2d.js";
import { rotation2, transformVector, transformTensor, bilinear, trace, doubleContraction, symPart, eigSym2 } from "../physics/tensor.js";

initPage("p1");

const PRESETS = {
  b16: { name: "問1.6 の B", T: [[2, 1], [-3, 4]] },
  sym: { name: "対称な例", T: [[3, 1], [1, 1]] },
  iso: { name: "等方 2I", T: [[2, 0], [0, 2]] },
};
const st = { preset: "b16", T: PRESETS.b16.T.map((r) => r.slice()), phi: 30, aDeg: 20, bDeg: 110 };
const deg = (d) => (d * Math.PI) / 180;
const unit = (d) => [Math.cos(deg(d)), Math.sin(deg(d))];

let presetSeg = null;
const entrySliders = [];
const controls = document.getElementById("controls");
controls.append(
  field("テンソル", presetSeg = segmented(Object.entries(PRESETS).map(([k, p]) => [k, p.name]), st.preset, (v) => {
    st.preset = v;
    st.T = PRESETS[v].T.map((r) => r.slice());
    entrySliders.forEach((s, k) => s.set(st.T[k >> 1][k & 1]));
    render();
  })),
  ...[[0, 0], [0, 1], [1, 0], [1, 1]].map(([i, j]) => {
    const s = slider({ label: `T${i + 1}${j + 1}`, min: -4, max: 4, step: 0.1, value: st.T[i][j], format: (v) => v.toFixed(1), onInput: (v) => { st.T[i][j] = v; presetSeg.set(null); render(); } });
    entrySliders.push(s);
    return s;
  }),
  slider({ label: "基底の回転 φ", min: 0, max: 180, step: 1, value: st.phi, format: (v) => `${v}°`, onInput: (v) => { st.phi = v; render(); } }),
  slider({ label: "a の向き", min: 0, max: 360, step: 1, value: st.aDeg, format: (v) => `${v}°`, onInput: (v) => { st.aDeg = v; render(); } }),
  slider({ label: "b の向き", min: 0, max: 360, step: 1, value: st.bDeg, format: (v) => `${v}°`, onInput: (v) => { st.bDeg = v; render(); } }));

const scene = createScene(document.getElementById("scene"), { width: 420, height: 360 });
const fmt2 = (v) => (Math.abs(v) < 5e-4 ? "0" : v.toFixed(2));
const explainComp = (prime) => (i, j, v) => h("div", {},
  h("strong", {}, `${prime ? "T'" : "T"}${i + 1}${j + 1} = ${fmt(v)}`),
  h("div", { class: "caption" }, prime
    ? `回した基底のベクトル e'${i + 1}, e'${j + 1} を写像に代入した値。φ を変えると変わる。`
    : `元の基底 e${i + 1}, e${j + 1} を写像に代入した値 T(e${i + 1}, e${j + 1})。`));
const mOrig = createMatrixView(document.getElementById("m-orig"), { labels: ["1", "2"], block: 2, format: fmt2, explain: explainComp(false) });
const mNew = createMatrixView(document.getElementById("m-new"), { labels: ["1′", "2′"], block: 2, format: fmt2, explain: explainComp(true) });
const plot = createPlot(document.getElementById("plot"), { width: 960, height: 300, xlabel: "基底の回転 φ（度）", ylabel: "成分" });

function render() {
  const Q = rotation2(deg(st.phi));
  const T = st.T, Tp = transformTensor(Q, T);
  const a = unit(st.aDeg), b = unit(st.bDeg), ap = transformVector(Q, a), bp = transformVector(Q, b);
  const vmax = Math.max(...T.flat().map(Math.abs), ...Tp.flat().map(Math.abs));
  mOrig.update(T, { vmax });
  mNew.update(Tp, { vmax });

  // scene: both bases, the two vectors (unchanged by the basis rotation)
  const e1p = Q[0], e2p = Q[1];
  scene.draw({
    view: { cx: 0, cz: 0, span: 3.4 },
    axes: { x: "x", up: "y", out: "z⊙" },
    lines: [{ x1: -1.6, z1: 0, x2: 1.6, z2: 0, color: "var(--line)" }, { x1: 0, z1: -1.6, x2: 0, z2: 1.6, color: "var(--line)" }],
    arrows: [
      { x: 0, z: 0, vx: 1.2, vz: 0, color: "var(--muted)", label: "e₁", width: 1.5 },
      { x: 0, z: 0, vx: 0, vz: 1.2, color: "var(--muted)", label: "e₂", width: 1.5 },
      { x: 0, z: 0, vx: 1.2 * e1p[0], vz: 1.2 * e1p[1], color: "var(--c1)", label: "e′₁" },
      { x: 0, z: 0, vx: 1.2 * e2p[0], vz: 1.2 * e2p[1], color: "var(--c1)", label: "e′₂" },
      { x: 0, z: 0, vx: a[0], vz: a[1], color: "var(--c2)", label: "a", width: 3 },
      { x: 0, z: 0, vx: b[0], vz: b[1], color: "var(--c3)", label: "b", width: 3 },
    ],
  });

  const S = symPart(T), ev = eigSym2(S);
  const row = (label, v1, v2, note) => h("tr", {}, h("td", {}, label), h("td", {}, v1), h("td", {}, v2), h("td", { class: "caption" }, note));
  document.getElementById("invariants").replaceChildren(h("table", { class: "data", style: { marginTop: "10px" } },
    h("tr", {}, h("th", {}, ""), h("th", {}, "元の基底"), h("th", {}, "回した基底"), h("th", {}, "")),
    row("a の成分", `(${fmt2(a[0])}, ${fmt2(a[1])})`, `(${fmt2(ap[0])}, ${fmt2(ap[1])})`, "変わる"),
    row("b の成分", `(${fmt2(b[0])}, ${fmt2(b[1])})`, `(${fmt2(bp[0])}, ${fmt2(bp[1])})`, "変わる"),
    row("値 T(a, b) = aᵢTᵢⱼbⱼ", fmt(bilinear(T, a, b)), fmt(bilinear(Tp, ap, bp)), "変わらない"),
    row("トレース Tᵢᵢ", fmt(trace(T)), fmt(trace(Tp)), "変わらない"),
    row("二重縮約 TᵢⱼTᵢⱼ", fmt(doubleContraction(T, T)), fmt(doubleContraction(Tp, Tp)), "変わらない"),
    row("対称部分の固有値", ev.values.map((v) => fmt(v, 3)).join(", "), eigSym2(symPart(Tp)).values.map((v) => fmt(v, 3)).join(", "), "変わらない")));

  const comps = componentsVsAngle(T);
  const colors = ["var(--c1)", "var(--c2)", "var(--c3)", "var(--c4)"];
  plot.update({
    series: [
      ...[[0, 0], [0, 1], [1, 0], [1, 1]].map(([i, j], k) => ({ name: `T′${i + 1}${j + 1}`, color: colors[k], points: phis.map((p, n) => [p, comps[n][i][j]]) })),
      { name: "トレース", color: "var(--muted)", dash: "2 3", points: phis.map((p) => [p, trace(T)]) },
    ],
    vlines: [{ x: st.phi }, ...(ev.degenerate ? [] : [{ x: ((ev.angle * 180) / Math.PI + 180) % 180, color: "var(--c5)", dash: "6 3" }])],
  });
}

// The curves depend only on T; recompute them only when T changes (not on phi, a, b).
const phis = Array.from({ length: 181 }, (_, k) => k);
let compsKey = "", compsCache = null;
function componentsVsAngle(T) {
  const key = JSON.stringify(T);
  if (key !== compsKey) { compsKey = key; compsCache = phis.map((p) => transformTensor(rotation2(deg(p)), T)); }
  return compsCache;
}
render();
