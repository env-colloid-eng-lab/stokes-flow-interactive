import { initPage, h, field, segmented, slider, fmt, frameThrottle } from "../ui/page.js";
import { createPlot } from "../ui/plot.js";
import { createScene } from "../ui/scene2d.js";
import { matvec, norm } from "../core/linalg.js";
import { translationField, rotationField, strainField, surfaceMoments, sphereRule } from "../physics/sphere.js";

initPage("p6");

const a = 1, mu = 1;
const MOTIONS = {
  trans: { name: "並進 U = e_z", field: (r) => translationField(r, [0, 0, 1], { a, mu }) },
  rot: { name: "回転 Ω = e_y", field: (r) => rotationField(r, [0, 1, 0], { a, mu }) },
  strain: { name: "伸長 E = diag(½, 0, −½)", field: (r) => strainField(r, E0, { a, mu }) },
};
const E0 = [[0.5, 0, 0], [0, 0, 0], [0, 0, -0.5]];
const st = { motion: "trans", nt: 4, np: 8 };

const schedule = frameThrottle(() => render());
document.getElementById("controls").append(
  field("運動", segmented(Object.entries(MOTIONS).map(([k, m]) => [k, m.name]), st.motion, (v) => { st.motion = v; schedule(); })),
  slider({ label: "cos θ 方向の点数 n_θ", min: 2, max: 12, step: 1, value: st.nt, format: (v) => String(v), onInput: (v) => { st.nt = v; schedule(); } }),
  slider({ label: "方位角の点数 n_φ", min: 3, max: 24, step: 1, value: st.np, format: (v) => String(v), onInput: (v) => { st.np = v; schedule(); } }));

const scene = createScene(document.getElementById("scene"), { width: 440, height: 400 });
const plot = createPlot(document.getElementById("plot"), { width: 960, height: 280, xlabel: "極角 θ（度）", ylabel: "表面力の z 成分 t_z" });

// exact values for each motion (a = mu = 1)
function exactAndComputed(m, nt, np) {
  const res = surfaceMoments(MOTIONS[m].field, { a, nt, np });
  if (m === "trans") return { got: res.F[2], exact: -6 * Math.PI, res, label: "F_z（力）", formula: "−6πμa U" };
  if (m === "rot") return { got: res.T[1], exact: -8 * Math.PI, res, label: "T_y（トルク）", formula: "−8πμa³ Ω" };
  return { got: res.S[0][0], exact: (20 * Math.PI) / 3 * 0.5, res, label: "S_xx（ストレスレット）", formula: "(20π/3) μa³ E_xx" };
}

// The flow, the tractions and the traction plot depend only on the motion, not on the quadrature.
let cache = { motion: null };
function motionPicture(m) {
  if (cache.motion === m) return cache;
  const fieldFn = MOTIONS[m].field;
  const arrows = [];
  for (let x = -3; x <= 3.001; x += 0.5)
    for (let z = -3; z <= 3.001; z += 0.5) {
      if (Math.hypot(x, z) < a * 1.05) continue;
      const u = fieldFn([x, 0, z]).u;
      arrows.push({ x, z, vx: 0.35 * u[0], vz: 0.35 * u[2], color: "var(--line)", width: 1.2 });
    }
  for (let k = 0; k < 24; k++) {
    const th = (2 * Math.PI * k) / 24, n = [Math.sin(th), 0, Math.cos(th)];
    const t = matvec(fieldFn(n.map((v) => a * v)).sigma, n);
    arrows.push({ x: a * n[0], z: a * n[2], vx: 0.4 * t[0], vz: 0.4 * t[2], color: "var(--c2)", width: 2 });
  }
  // traction along the great circle phi = 0, split into pressure and viscous parts
  const ths = Array.from({ length: 91 }, (_, k) => 2 * k);
  const parts = ths.map((d) => {
    const th = (d * Math.PI) / 180, n = [Math.sin(th), 0, Math.cos(th)];
    const fl = fieldFn(n.map((v) => a * v));
    const visc = matvec(fl.L.map((row, i) => row.map((v, j) => mu * (v + fl.L[j][i]))), n);
    return { p: -fl.p * n[2], v: visc[2] };
  });
  plot.update({
    series: [
      { name: "圧力の寄与", color: "var(--c1)", points: ths.map((d, k) => [d, parts[k].p]) },
      { name: "粘性応力の寄与", color: "var(--c3)", points: ths.map((d, k) => [d, parts[k].v]) },
      { name: "合計 t_z", color: "var(--c2)", width: 3, points: ths.map((d, k) => [d, parts[k].p + parts[k].v]) },
    ],
    hlines: [{ y: 0 }],
  });
  cache = { motion: m, arrows };
  return cache;
}

function render() {
  const m = st.motion;
  const { arrows } = motionPicture(m);
  // quadrature nodes that lie on this cross-section (azimuth 0, and pi when n_phi is even)
  const nodes = sphereRule(st.nt, st.np).normals.filter((n) => Math.abs(n[1]) < 1e-9)
    .map((n) => ({ x: a * n[0], z: a * n[2], a: 0.06, fill: "var(--ink)", stroke: "none" }));
  scene.draw({ view: { cx: 0, cz: 0, span: 6.6 }, spheres: [{ x: 0, z: 0, a }, ...nodes], arrows });

  const { got, exact, res, label, formula } = exactAndComputed(m, st.nt, st.np);
  const rows = [
    h("tr", {}, h("th", {}, ""), h("th", {}, "数値積分"), h("th", {}, "厳密値"), h("th", {}, "相対誤差")),
    h("tr", {}, h("td", {}, label), h("td", {}, fmt(got, 8)), h("td", {}, `${fmt(exact, 8)}（${formula}）`), h("td", {}, fmt(Math.abs(got - exact) / Math.abs(exact), 2))),
    h("tr", {}, h("td", {}, "力 |F|"), h("td", {}, fmt(norm(res.F), 3)), h("td", {}, m === "trans" ? "6π" : "0"), h("td", {}, "")),
    h("tr", {}, h("td", {}, "トルク |T|"), h("td", {}, fmt(norm(res.T), 3)), h("td", {}, m === "rot" ? "8π" : "0"), h("td", {}, "")),
  ];
  if (m === "trans")
    rows.push(
      h("tr", {}, h("td", {}, "うち圧力 F_z"), h("td", {}, fmt(res.Fpressure[2], 6)), h("td", {}, "−2π（1/3）"), h("td", {}, "")),
      h("tr", {}, h("td", {}, "うち粘性 F_z"), h("td", {}, fmt(res.Fviscous[2], 6)), h("td", {}, "−4π（2/3）"), h("td", {}, "")));
  document.getElementById("result").replaceChildren(
    h("table", { class: "data" }, ...rows),
    h("p", { class: "caption" }, `求積点 ${st.nt} × ${st.np} = ${st.nt * st.np} 点。` +
      (m === "trans" ? "表面力はどの点でも −(3μ/2a)U で同じなので、少ない点でも正確に積分できる。" :
        m === "rot" ? "圧力は 0。表面力は球面に沿った向きだけを持つ。" : "力とトルクは打ち消し合い、ストレスレットだけが残る。")));
}
render();
