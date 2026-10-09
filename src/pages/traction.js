import { initPage, h, field, segmented, slider, fmt } from "../ui/page.js";
import { createMatrixView } from "../ui/matrixView.js";
import { createPlot } from "../ui/plot.js";
import { createScene } from "../ui/scene2d.js";
import { traction2, eigSym2, stressTorqueZ, mohr, symPart } from "../physics/tensor.js";

initPage("p3");

const PRESETS = {
  p: { name: "圧力だけ", S: [[-1, 0], [0, -1]] },
  newton: { name: "単純せん断（μ=γ̇=1, p=0）", S: [[0, 1], [1, 0]] },
  q21: { name: "問2.1（p=1, τ=0.5）", S: [[-1, 0.5], [0.5, -1]] },
  general: { name: "一般の対称な例", S: [[1, 0.6], [0.6, -0.5]] },
  asym: { name: "非対称（禁止）", S: [[0, 1], [0.3, 0]] },
};
const st = { preset: "q21", S: PRESETS.q21.S.map((r) => r.slice()), thDeg: 30 };
const deg = (d) => (d * Math.PI) / 180;

const sliders = [];
document.getElementById("controls").append(
  field("応力", segmented(Object.entries(PRESETS).map(([k, p]) => [k, p.name]), st.preset, (v) => {
    st.preset = v;
    st.S = PRESETS[v].S.map((r) => r.slice());
    sliders.forEach((s, k) => s.set(st.S[k >> 1][k & 1]));
    render();
  })),
  ...[[0, 0], [0, 1], [1, 0], [1, 1]].map(([i, j]) => {
    const s = slider({ label: `σ${i + 1}${j + 1}`, min: -2, max: 2, step: 0.05, value: st.S[i][j], format: (v) => v.toFixed(2), onInput: (v) => { st.S[i][j] = v; render(); } });
    sliders.push(s);
    return s;
  }),
  slider({ label: "法線の向き θ", min: 0, max: 360, step: 1, value: st.thDeg, format: (v) => `${v}°`, onInput: (v) => { st.thDeg = v; render(); } }));

const scene = createScene(document.getElementById("scene"), { width: 420, height: 400 });
const matrix = createMatrixView(document.getElementById("matrix"), {
  labels: ["x", "y"], block: 2, format: (v) => (Math.abs(v) < 5e-4 ? "0" : v.toFixed(2)),
  explain: (i, j, v) => h("div", {}, h("strong", {}, `σ${i + 1}${j + 1} = ${fmt(v)}`),
    h("div", { class: "caption" }, `${"xy"[j]} に垂直な面に働く、${"xy"[i]} 方向の力（単位面積当たり）。${i === j ? "面に垂直な成分（法線応力）。" : "面に沿った成分（せん断応力）。"}`)),
});
const plot = createPlot(document.getElementById("plot"), { height: 250, xlabel: "法線の向き θ（度）", ylabel: "表面力の成分" });
const mohrPlot = createPlot(document.getElementById("mohr"), { width: 360, height: 300, xlabel: "t_n（法線方向）", ylabel: "t_s（面に沿う方向）" });

function render() {
  const S = st.S;
  matrix.update(S, { vmax: 2 });
  const n = [Math.cos(deg(st.thDeg)), Math.sin(deg(st.thDeg))];
  const r = traction2(S, n);
  const t1 = traction2(S, [1, 0]).t, t2 = traction2(S, [0, 1]).t;
  const sym = Math.abs(stressTorqueZ(S)) < 1e-12;
  const ev = eigSym2(symPart(S));
  const tang = r.tangent;
  const k = 1.6; // arrow length per unit stress
  const axis = (v) => ({ x1: -2.2 * v[0], z1: -2.2 * v[1], x2: 2.2 * v[0], z2: 2.2 * v[1], color: "var(--c4)", dash: "6 4" });
  scene.draw({
    view: { cx: 0, cz: 0, span: 5 },
    axes: { x: "x", up: "y", out: "z⊙" },
    polys: [{ points: [[-0.6, -0.6], [0.6, -0.6], [0.6, 0.6], [-0.6, 0.6]], opacity: 0.35 }],
    lines: [
      ...(sym ? [axis(ev.vectors[0]), axis(ev.vectors[1])] : []),
      { x1: -2 * tang[0], z1: -2 * tang[1], x2: 2 * tang[0], z2: 2 * tang[1], color: "var(--ink)", dash: "0" },
    ],
    arrows: [
      { x: 0, z: 0, vx: 1.4 * n[0], vz: 1.4 * n[1], color: "var(--muted)", label: "n", width: 2 },
      { x: 0, z: 0, vx: k * n[0] * t1[0], vz: k * n[0] * t1[1], color: "var(--c1)", width: 1.2 },
      { x: k * n[0] * t1[0], z: k * n[0] * t1[1], vx: k * n[1] * t2[0], vz: k * n[1] * t2[1], color: "var(--c3)", width: 1.2 },
      { x: 0, z: 0, vx: k * r.t[0], vz: k * r.t[1], color: "var(--c2)", label: "t", width: 3.5 },
    ],
  });

  const torque = stressTorqueZ(S);
  document.getElementById("readout").replaceChildren(h("table", { class: "data", style: { marginTop: "10px" } },
    h("tr", {}, h("td", {}, "法線 n"), h("td", {}, `(${fmt(n[0], 3)}, ${fmt(n[1], 3)})`)),
    h("tr", {}, h("td", {}, "表面力 t = σn"), h("td", {}, `(${fmt(r.t[0], 3)}, ${fmt(r.t[1], 3)})`)),
    h("tr", {}, h("td", {}, "法線成分 t_n"), h("td", {}, fmt(r.normal, 4))),
    h("tr", {}, h("td", {}, "せん断成分 t_s"), h("td", {}, fmt(r.shear, 4))),
    h("tr", {}, h("td", {}, "主応力（対称部分）"), h("td", {}, ev.values.map((v) => fmt(v, 3)).join(", "))),
    h("tr", {}, h("td", {}, "要素に残るモーメント σ₂₁ − σ₁₂"), h("td", { style: { color: sym ? "" : "var(--warn)" } }, fmt(torque, 3)))),
    h("p", { class: "caption" }, sym
      ? "応力が対称なので、小さな要素に働くモーメントは打ち消し合う。"
      : "応力が非対称なので、小さな要素に正味のモーメントが残る。要素を小さくすると角加速度が発散するので、通常の流体ではこの応力は許されない。"));

  const ths = Array.from({ length: 181 }, (_, k) => 2 * k);
  const rs = ths.map((d) => traction2(S, [Math.cos(deg(d)), Math.sin(deg(d))]));
  plot.update({
    series: [
      { name: "t_n", color: "var(--c1)", points: ths.map((d, k) => [d, rs[k].normal]) },
      { name: "t_s", color: "var(--c2)", points: ths.map((d, k) => [d, rs[k].shear]) },
    ],
    vlines: [{ x: st.thDeg }],
    hlines: [{ y: 0 }],
  });
  const m = mohr(symPart(S));
  const span = Math.max(1, m.radius * 1.3 + Math.abs(m.centre));
  mohrPlot.update({
    series: [{ name: sym ? "(t_n, t_s)：円" : "(t_n, t_s)", color: "var(--c1)", points: rs.map((x) => [x.normal, x.shear]) }],
    points: [{ x: r.normal, y: r.shear, color: "var(--c2)" }, ...(sym ? ev.values.map((v) => ({ x: v, y: 0, color: "var(--c4)", r: 3.5 })) : [])],
    xdomain: [m.centre - span, m.centre + span], ydomain: [-span, span],
    hlines: [{ y: 0 }],
  });
}
render();
