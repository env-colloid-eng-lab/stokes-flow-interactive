import { initPage, h, field, segmented, slider, fmt } from "../ui/page.js";
import { createMatrixView } from "../ui/matrixView.js";
import { createScene } from "../ui/scene2d.js";
import { matvec, dot } from "../core/linalg.js";
import { symPart, antiPart, eigSym2, omegaZ, trace, expm2 } from "../physics/tensor.js";

initPage("p2");

const PRESETS = {
  shear: { name: "単純せん断", L: [[0, 1], [0, 0]] },
  ext: { name: "純粋伸長", L: [[0.5, 0], [0, -0.5]] },
  rot: { name: "剛体回転", L: [[0, -0.5], [0.5, 0]] },
  comp: { name: "圧縮性の例", L: [[0.3, 0.6], [-0.2, 0.1]] },
};
const st = { preset: "shear", L: PRESETS.shear.L.map((r) => r.slice()), part: "L", t: 0, lineDeg: 0, anim: null };
const TMAX = 3;

const current = () => (st.part === "E" ? symPart(st.L) : st.part === "W" ? antiPart(st.L) : st.L);

let presetSeg = null;
const entrySliders = [];
const tSlider = slider({ label: "時刻 t", min: 0, max: TMAX, step: 0.01, value: 0, format: (v) => v.toFixed(2), onInput: (v) => { stop(); st.t = v; render(); } });
document.getElementById("controls").append(
  field("流れ", presetSeg = segmented(Object.entries(PRESETS).map(([k, p]) => [k, p.name]), st.preset, (v) => {
    st.preset = v;
    st.L = PRESETS[v].L.map((r) => r.slice());
    entrySliders.forEach((s, k) => s.set(st.L[k >> 1][k & 1]));
    render();
  })),
  ...[[0, 0], [0, 1], [1, 0], [1, 1]].map(([i, j]) => {
    const s = slider({ label: `L${i + 1}${j + 1}`, min: -1, max: 1, step: 0.05, value: st.L[i][j], format: (v) => v.toFixed(2), onInput: (v) => { st.L[i][j] = v; presetSeg.set(null); render(); } });
    entrySliders.push(s);
    return s;
  }),
  field("動かす部分", segmented([["L", "L 全体"], ["E", "E だけ"], ["W", "W だけ"]], st.part, (v) => { st.part = v; render(); })),
  slider({ label: "物質線の向き", min: 0, max: 180, step: 1, value: st.lineDeg, format: (v) => `${v}°`, onInput: (v) => { st.lineDeg = v; render(); } }));

const playBtn = h("button", { class: "action", type: "button", onclick: () => (st.anim ? stop() : play()) }, "再生");
document.getElementById("anim").append(playBtn, tSlider);
function play() {
  if (st.t >= TMAX) st.t = 0;
  playBtn.textContent = "停止";
  let last = performance.now();
  const step = (now) => {
    st.t = Math.min(TMAX, st.t + Math.max(0, now - last) / 1000); // the first frame can predate the click
    last = now;
    tSlider.set(st.t);
    render();
    if (st.t >= TMAX) { stop(); return; }
    st.anim = requestAnimationFrame(step);
  };
  st.anim = requestAnimationFrame(step);
}
function stop() {
  if (st.anim) cancelAnimationFrame(st.anim);
  st.anim = null;
  playBtn.textContent = "再生";
}

const scene = createScene(document.getElementById("scene"), { width: 440, height: 400 });
const fmt2 = (v) => (Math.abs(v) < 5e-4 ? "0" : v.toFixed(2));
const explain = (name) => (i, j, v) => h("div", {}, h("strong", {}, `${name}${i + 1}${j + 1} = ${fmt(v)}`),
  h("div", { class: "caption" }, name === "E" ? `(L${i + 1}${j + 1} + L${j + 1}${i + 1})/2。伸び縮みを表す。` : `(L${i + 1}${j + 1} − L${j + 1}${i + 1})/2。長さを変えない回転を表す。`));
const mE = createMatrixView(document.getElementById("m-e"), { labels: ["x", "y"], block: 2, format: fmt2, explain: explain("E") });
const mW = createMatrixView(document.getElementById("m-w"), { labels: ["x", "y"], block: 2, format: fmt2, explain: explain("W") });

const circle = Array.from({ length: 72 }, (_, k) => [Math.cos((2 * Math.PI * k) / 72), Math.sin((2 * Math.PI * k) / 72)]);
const gridLines = [];
for (const c of [-0.7, -0.35, 0, 0.35, 0.7]) {
  gridLines.push(Array.from({ length: 15 }, (_, k) => [-0.7 + (1.4 * k) / 14, c]));
  gridLines.push(Array.from({ length: 15 }, (_, k) => [c, -0.7 + (1.4 * k) / 14]));
}

function render() {
  const A = current();
  const F = expm2(A.map((r) => r.map((v) => v * st.t))); // flow map x(t) = exp(A t) x0
  const map = (p) => matvec(F, p);
  const E = symPart(st.L), W = antiPart(st.L), ev = eigSym2(E);
  const evA = eigSym2(symPart(A));
  mE.update(E, { vmax: 1 });
  mW.update(W, { vmax: 1 });

  const vel = [];
  for (let x = -2; x <= 2.001; x += 0.5)
    for (let y = -2; y <= 2.001; y += 0.5) {
      const u = matvec(A, [x, y]);
      vel.push({ x, z: y, vx: 0.35 * u[0], vz: 0.35 * u[1], color: "var(--line)", width: 1.2 });
    }
  const d0 = [Math.cos((st.lineDeg * Math.PI) / 180), Math.sin((st.lineDeg * Math.PI) / 180)];
  const d = map(d0);
  const axis = (v, s) => ({ x1: -2.4 * v[0], z1: -2.4 * v[1], x2: 2.4 * v[0], z2: 2.4 * v[1], color: "var(--c4)", dash: s });
  // zoom out when the element has grown beyond the default view
  const extent = Math.max(2.2, ...circle.map(map).map((p) => Math.max(Math.abs(p[0]), Math.abs(p[1]))));
  scene.draw({
    view: { cx: 0, cz: 0, span: 2.3 * extent },
    axes: { x: "x", up: "y", out: "z⊙" },
    lines: evA.degenerate ? [] : [axis(evA.vectors[0], "6 4"), axis(evA.vectors[1], "2 4")],
    arrows: [...vel, { x: -d[0] / 2, z: -d[1] / 2, vx: d[0], vz: d[1], color: "var(--c2)", width: 3 }],
    polys: [{ points: circle.map(map), opacity: 0.6 }],
    paths: gridLines.map((ln) => ({ points: ln.map(map), color: "var(--c1)", width: 1 })),
  });

  // rate of change of the squared length of the material line, now
  const rate = 2 * dot(d, matvec(symPart(A), d));
  const len = Math.hypot(...d);
  const div = trace(A), wz = omegaZ(A);
  const P = st.part === "L" ? "L" : st.part; // name of the part that moves the element
  document.getElementById("readout").replaceChildren(h("table", { class: "data", style: { marginTop: "10px" } },
    h("tr", {}, h("th", {}, `動かしている部分 A = ${P}`), h("th", {}, "")),
    h("tr", {}, h("td", {}, "A の対称部分の固有値（主ひずみ速度）"), h("td", {}, evA.values.map((v) => fmt(v, 3)).join(", "))),
    h("tr", {}, h("td", {}, "局所角速度 ω_z"), h("td", {}, fmt(wz, 3))),
    h("tr", {}, h("td", {}, "発散 tr A（面積の増加率）"), h("td", {}, fmt(div, 3))),
    h("tr", {}, h("td", {}, `物質線の長さ（t = ${st.t.toFixed(2)}）`), h("td", {}, fmt(len, 4))),
    h("tr", {}, h("td", {}, "d|δx|²/dt = 2δxᵀ(A の対称部分)δx"), h("td", {}, fmt(rate, 4))),
    h("tr", {}, h("td", {}, "要素の面積比 det exp(At)"), h("td", {}, fmt(F[0][0] * F[1][1] - F[0][1] * F[1][0], 4)))),
    h("p", { class: "caption" }, st.part === "W"
      ? `W だけで動かすと、要素は角速度 ${fmt(wz, 3)} で回るだけで、円は円のまま、物質線の長さも変わらない。`
      : st.part === "E"
        ? `E だけで動かすと、要素は主軸の方向に伸び縮みするだけで回転しない（主軸は動かない）。${Math.abs(div) > 1e-9 ? "tr E ≠ 0 なので面積も変わる。" : "tr E = 0 なので面積は保たれる。"}`
        : "L 全体では、伸び縮みと回転が同時に起こる。「E だけ」「W だけ」と見比べよう。"));
}
render();
