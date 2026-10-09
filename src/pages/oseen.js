import { initPage, h, slider, fmt, frameThrottle } from "../ui/page.js";
import { createPlot } from "../ui/plot.js";
import { createScene } from "../ui/scene2d.js";
import { createMatrixView } from "../ui/matrixView.js";
import { matvec, norm, dot } from "../core/linalg.js";
import { oseen, lapOseen, sphereFlow, fourierStokeslet } from "../physics/oseen.js";
import { rotationField, surfaceAverage } from "../physics/sphere.js";

initPage("p5");

const mu = 1, a = 1;
const deg = (d) => (d * Math.PI) / 180;
const inPlane = (th) => [Math.sin(th), 0, Math.cos(th)]; // angle from +z towards +x

// ---------------------------------------------------------------------
// Fourier space: the velocity is the part of f across k
// ---------------------------------------------------------------------
const fo = { k: 60, f: 0 };
const scheduleF = frameThrottle(() => renderFourier());
document.getElementById("fourier-controls").append(
  slider({ label: "波数ベクトル k の向き", min: 0, max: 180, step: 1, value: fo.k, format: (v) => `${v}°`, onInput: (v) => { fo.k = v; scheduleF(); } }),
  slider({ label: "力 f の向き", min: 0, max: 360, step: 1, value: fo.f, format: (v) => `${v}°`, onInput: (v) => { fo.f = v; scheduleF(); } }));
const fScene = createScene(document.getElementById("fourier-scene"), { width: 420, height: 360 });
function renderFourier() {
  const kh = inPlane(deg(fo.k)), f = inPlane(deg(fo.f));
  const { u, pImag } = fourierStokeslet(kh, f, { mu }); // |k| = 1: mu k^2 u = P f
  const par = kh.map((v) => v * dot(kh, f));
  const t = [kh[2], 0, -kh[0]]; // in-plane direction along the wavefronts
  const lines = [];
  for (let s = -3; s <= 3; s++) {
    const c = kh.map((v) => (v * s) / 2.5);
    lines.push({ x1: c[0] - 3 * t[0], z1: c[2] - 3 * t[2], x2: c[0] + 3 * t[0], z2: c[2] + 3 * t[2], color: "var(--line)", dash: "4 4" });
  }
  fScene.draw({
    view: { cx: 0, cz: 0, span: 3.4 },
    lines,
    arrows: [
      { x: 0, z: 0, vx: 1.2 * kh[0], vz: 1.2 * kh[2], color: "var(--muted)", width: 1.5, label: "k" },
      { x: 0, z: 0, vx: f[0], vz: f[2], color: "var(--ink)", width: 3, label: "f" },
      { x: 0, z: 0, vx: par[0], vz: par[2], color: "var(--c2)", width: 2.5, label: "圧力が受け持つ", offset: -5 },
      { x: 0, z: 0, vx: u[0], vz: u[2], color: "var(--c1)", width: 2.5, label: "μk²û", offset: 5 },
    ],
  });
  document.getElementById("fourier-out").replaceChildren(h("table", { class: "data" },
    h("tr", {}, h("td", {}, "k 方向の成分 k·f/|k|（圧力勾配と釣り合う）"), h("td", {}, fmt(dot(kh, f), 3))),
    h("tr", {}, h("td", {}, "波面に沿った成分 |(I − kkᵀ/k²) f|（流れになる）"), h("td", {}, fmt(norm(u), 3))),
    h("tr", {}, h("td", {}, "k·û（非圧縮の条件）"), h("td", {}, fmt(Math.abs(dot(kh, u)) < 1e-14 ? 0 : dot(kh, u), 2))),
    h("tr", {}, h("td", {}, "圧力の振幅 |p̂|·|k|"), h("td", {}, fmt(Math.abs(pImag), 3)))),
    h("p", { class: "caption" }, Math.abs(dot(kh, f)) > 0.999 ? "力が k に平行だと、すべて圧力が受け持ち、流れは生じない。" :
      Math.abs(dot(kh, f)) < 1e-3 ? "力が波面に沿っていると、圧力は生じず、力はすべて流れになる。" :
        "流速は波面（破線）に沿う向きだけを持つ。波面に沿った流れは、どの点でも流入と流出が釣り合う。"));
}

// ---------------------------------------------------------------------
// The stokeslet in real space
// ---------------------------------------------------------------------
const sk = { f: 0, r: 2.5, th: 50 };
const scheduleS = frameThrottle(() => renderStokeslet());
document.getElementById("stokeslet-controls").append(
  slider({ label: "力 f の向き", min: 0, max: 360, step: 1, value: sk.f, format: (v) => `${v}°`, onInput: (v) => { sk.f = v; scheduleS(); } }),
  slider({ label: "観測点の距離 r", min: 0.5, max: 5, step: 0.05, value: sk.r, format: (v) => v.toFixed(2), onInput: (v) => { sk.r = v; scheduleS(); } }),
  slider({ label: "観測点の向き", min: 0, max: 360, step: 1, value: sk.th, format: (v) => `${v}°`, onInput: (v) => { sk.th = v; scheduleS(); } }));
const sScene = createScene(document.getElementById("stokeslet-scene"), { width: 440, height: 440 });
const c8 = 8 * Math.PI * mu;
const gView = createMatrixView(document.getElementById("stokeslet-matrix"), {
  labels: ["x", "y", "z"], block: 3, format: (v) => (Math.abs(v) < 5e-4 ? "0" : v.toFixed(3)),
  explain: (i, j, v) => {
    const r = inPlane(deg(sk.th)).map((x) => sk.r * x);
    return h("div", {}, h("strong", {}, `8πμ G${"xyz"[i]}${"xyz"[j]} = ${fmt(v, 4)}`),
      h("div", { class: "caption" }, `= δ${"xyz"[i]}${"xyz"[j]}/r + r${"xyz"[i]} r${"xyz"[j]}/r³ = ${i === j ? 1 : 0}/${fmt(sk.r, 3)} + (${fmt(r[i], 3)})(${fmt(r[j], 3)})/${fmt(sk.r ** 3, 4)}。` +
        `${"xyz"[j]} 方向の単位の力が、観測点に作る ${"xyz"[i]} 方向の流速。`));
  },
});
function renderStokeslet() {
  const f = inPlane(deg(sk.f));
  const arrows = [];
  const uref = 1 / (c8 * 1.5);
  for (let x = -4.5; x <= 4.501; x += 0.75)
    for (let z = -4.5; z <= 4.501; z += 0.75) {
      if (Math.hypot(x, z) < 0.4) continue;
      const u = matvec(oseen([x, 0, z], { mu }), f), m = norm(u);
      const len = (0.6 * m) / (m + uref); // saturating length so that the far field stays visible
      arrows.push({ x, z, vx: (len * u[0]) / m, vz: (len * u[2]) / m, color: "var(--line)", width: 1.3 });
    }
  const rp = inPlane(deg(sk.th)).map((x) => sk.r * x);
  const G = oseen(rp, { mu });
  const up = matvec(G, f);
  sScene.draw({
    view: { cx: 0, cz: 0, span: 10 },
    arrows: [...arrows, { x: 0, z: 0, vx: 1.2 * f[0], vz: 1.2 * f[2], color: "var(--ink)", width: 3, label: "f" },
      { x: rp[0], z: rp[2], vx: 20 * up[0], vz: 20 * up[2], color: "var(--c2)", width: 3, label: "u" }],
    spheres: [{ x: rp[0], z: rp[2], a: 0.08, fill: "var(--c2)", stroke: "none" }],
  });
  gView.update(G.map((row) => row.map((v) => v * c8)), { vmax: 2 / 0.5 });
  const rh = rp.map((v) => v / sk.r);
  document.getElementById("stokeslet-out").replaceChildren(h("table", { class: "data" },
    h("tr", {}, h("td", {}, "観測点の向き e_r の固有値 × 8πμ"), h("td", {}, `${fmt(2 / sk.r, 4)}（= 2/r）`)),
    h("tr", {}, h("td", {}, "e_r に垂直な 2 方向の固有値 × 8πμ"), h("td", {}, `${fmt(1 / sk.r, 4)}（= 1/r）`)),
    h("tr", {}, h("td", {}, "f と e_r のなす角"), h("td", {}, `${fmt((Math.acos(Math.max(-1, Math.min(1, dot(f, rh)))) * 180) / Math.PI, 3)}°`)),
    h("tr", {}, h("td", {}, "観測点の流速 |u| × 8πμ"), h("td", {}, fmt(norm(up) * c8, 4)))),
    h("p", { class: "caption" }, "G は対称で、観測点の向き e_r が主軸になる。力の向きに沿って並んだ点は、横に並んだ点の 2 倍の速さで動く。"));
}
const decay = createPlot(document.getElementById("decay-plot"), { width: 960, height: 300, xlog: true, ylog: true, xlabel: "力の向きに垂直な方向の距離 r/a", ylabel: "流速 |u|（球の速さ = 1）" });
{
  const rs = Array.from({ length: 41 }, (_, k) => 10 ** ((2 * k) / 40));
  const f = [0, 0, 6 * Math.PI * mu * a]; // the force that moves a sphere of radius a with unit speed
  const at = (r) => [r, 0, 0];
  decay.update({
    series: [
      { name: "点力 G f（1/r）", color: "var(--c1)", points: rs.map((r) => [r, norm(matvec(oseen(at(r), { mu }), f))]) },
      { name: "回転する球（1/r²）", color: "var(--c3)", points: rs.map((r) => [r, norm(rotationField(at(r), [0, 1, 0], { a, mu }).u)]) },
      { name: "有限サイズの項 (a²/6)∇²G f（1/r³）", color: "var(--c2)", points: rs.map((r) => [r, norm(matvec(lapOseen(at(r), { mu }), f)) * (a * a) / 6]) },
    ],
  });
}

// ---------------------------------------------------------------------
// The point force on a sphere surface
// ---------------------------------------------------------------------
const sf = { R: 1 };
const scheduleSurf = frameThrottle(() => renderSurface());
const surfPlot = createPlot(document.getElementById("surface-plot"), { height: 280, xlabel: "極角 θ（度）", ylabel: "半径 R の球面上の流速 ÷ U" });
document.getElementById("surface-controls").append(
  slider({ label: "評価する球面の半径 R/a", min: 1, max: 4, step: 0.05, value: sf.R, format: (v) => v.toFixed(2), onInput: (v) => { sf.R = v; scheduleSurf(); } }));
function renderSurface() {
  const U = [0, 0, 1], f = U.map((v) => 6 * Math.PI * mu * a * v);
  const Gf = (r) => matvec(oseen(r, { mu }), f);
  const full = (r) => sphereFlow(r, U, { a, mu });
  const ths = Array.from({ length: 91 }, (_, k) => 2 * k);
  const on = (fn, d) => fn(inPlane(deg(d)).map((v) => sf.R * v));
  surfPlot.update({
    series: [
      { name: "点力だけ u_z", color: "var(--c1)", points: ths.map((d) => [d, on(Gf, d)[2]]) },
      { name: "点力だけ u_x", color: "var(--c1)", dash: "5 4", points: ths.map((d) => [d, on(Gf, d)[0]]) },
      { name: "点力 + 有限サイズの項 u_z", color: "var(--c2)", width: 3, points: ths.map((d) => [d, on(full, d)[2]]) },
      { name: "点力 + 有限サイズの項 u_x", color: "var(--c2)", dash: "5 4", points: ths.map((d) => [d, on(full, d)[0]]) },
    ],
    hlines: [{ y: 0 }, { y: 1 }], ydomain: [-0.5, 1.6],
  });
  const avgG = surfaceAverage(Gf, [0, 0, 0], { a: sf.R })[2], avgF = surfaceAverage(full, [0, 0, 0], { a: sf.R })[2];
  document.getElementById("surface-out").replaceChildren(h("table", { class: "data" },
    h("tr", {}, h("th", {}, ""), h("th", {}, "球面平均の u_z"), h("th", {}, "u_z の最大 − 最小")),
    h("tr", {}, h("td", {}, "点力だけ"), h("td", {}, fmt(avgG, 6)), h("td", {}, fmt(0.75 / sf.R, 4))),
    h("tr", {}, h("td", {}, "点力 + 有限サイズの項"), h("td", {}, fmt(avgF, 6)), h("td", {}, fmt(Math.abs(on(full, 0)[2] - on(full, 90)[2]), 4)))),
    h("p", { class: "caption" }, sf.R === 1 ? "R = a では、有限サイズの項を加えた流れは球面上のどこでも U に一致し、剛体の並進になる。点力だけでも平均は U に等しいが、場所によって速さが違う。" :
      `R = ${sf.R.toFixed(2)}a では、平均はどちらも a/R = ${fmt(1 / sf.R, 4)}。有限サイズの項は平均に寄与せず、R とともに速く消える。`));
}

renderFourier();
renderStokeslet();
renderSurface();
