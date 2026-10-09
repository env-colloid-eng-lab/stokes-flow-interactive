import { initPage, h, field, segmented, slider, fmt, frameThrottle } from "../ui/page.js";
import { createPlot } from "../ui/plot.js";
import { createScene } from "../ui/scene2d.js";
import { shearReversal, rmsDistance } from "../physics/kinematics.js";
import { dissipationDensity, dissipationIntegral, translationField, rotationField } from "../physics/sphere.js";

initPage("p4");

// ---------------------------------------------------------------------
// Reynolds number of everyday and microscopic motions (SI units)
// ---------------------------------------------------------------------
const FLUIDS = {
  water: { name: "水", mu: 1.0e-3, rho: 998 },
  glycerol: { name: "グリセリン", mu: 1.41, rho: 1261 },
  air: { name: "空気", mu: 1.8e-5, rho: 1.2 },
};
const PRESETS = {
  colloid: { name: "沈降するコロイド", a: 0.5e-6, U: 0.5e-6, fluid: "water" },
  ecoli: { name: "泳ぐ大腸菌", a: 1e-6, U: 30e-6, fluid: "water" },
  fog: { name: "落ちる霧の粒", a: 10e-6, U: 1e-2, fluid: "air" },
  sand: { name: "沈む砂粒", a: 1e-4, U: 1e-2, fluid: "water" },
  swimmer: { name: "泳ぐ人", a: 1, U: 1, fluid: "water" },
};
const RHO_P = 1000; // particle density assumed for the relaxation time
const re = { ...PRESETS.colloid, preset: "colloid" };

const unit = (v, units) => {
  for (const [scale, name] of units) if (Math.abs(v) >= scale * 0.9999) return `${+(v / scale).toPrecision(3)} ${name}`;
  const [scale, name] = units.at(-1);
  return `${+(v / scale).toPrecision(3)} ${name}`;
};
const LEN = [[1, "m"], [1e-3, "mm"], [1e-6, "µm"], [1e-9, "nm"]];
const SPEED = [[1, "m/s"], [1e-3, "mm/s"], [1e-6, "µm/s"], [1e-9, "nm/s"]];
const TIME = [[1, "s"], [1e-3, "ms"], [1e-6, "µs"], [1e-9, "ns"], [1e-12, "ps"]];

let presetSeg, fluidSeg;
const aSlider = slider({ label: "粒子の半径 a", min: -9, max: 0, step: 0.05, value: Math.log10(re.a), format: (v) => unit(10 ** v, LEN), onInput: (v) => { re.a = 10 ** v; presetSeg.set(null); renderRe(); } });
const uSlider = slider({ label: "速さ U", min: -9, max: 1, step: 0.05, value: Math.log10(re.U), format: (v) => unit(10 ** v, SPEED), onInput: (v) => { re.U = 10 ** v; presetSeg.set(null); renderRe(); } });
document.getElementById("re-controls").append(
  field("例", presetSeg = segmented(Object.entries(PRESETS).map(([k, p]) => [k, p.name]), re.preset, (v) => {
    Object.assign(re, PRESETS[v]);
    aSlider.set(Math.log10(re.a)); uSlider.set(Math.log10(re.U)); fluidSeg.set(re.fluid);
    renderRe();
  })),
  field("流体", fluidSeg = segmented(Object.entries(FLUIDS).map(([k, f]) => [k, f.name]), re.fluid, (v) => { re.fluid = v; presetSeg.set(null); renderRe(); })),
  aSlider, uSlider);
const rePlot = createPlot(document.getElementById("re-plot"), { height: 280, xlog: true, ylog: true, xlabel: "粒子の半径 a（m）", ylabel: "Re = ρUa/μ" });

function renderRe() {
  const { mu, rho } = FLUIDS[re.fluid];
  const Re = (rho * re.U * re.a) / mu;
  const tauNu = (rho * re.a ** 2) / mu, tauAdv = re.a / re.U;
  const tauP = (2 * RHO_P * re.a ** 2) / (9 * mu);
  const as = Array.from({ length: 46 }, (_, k) => 10 ** (-9 + (9 * k) / 45));
  rePlot.update({
    series: Object.entries(FLUIDS).map(([k, f], i) => ({ name: `${f.name}（U = ${unit(re.U, SPEED)}）`, color: `var(--c${i + 1})`, width: k === re.fluid ? 3 : 1.5, points: as.map((a) => [a, (f.rho * re.U * a) / f.mu]) })),
    hlines: [{ y: 1 }], points: [{ x: re.a, y: Re, color: "var(--ink)" }],
    ydomain: [1e-16, 1e8],
  });
  const verdict = Re < 0.01 ? "慣性は無視できる。ストークス方程式がよく成り立つ。" : Re < 1 ? "ストークス近似はおおむね使えるが、Re の一次の補正（オゼーン補正）が見え始める。" : "慣性が効く。ストークス方程式は使えない。";
  document.getElementById("re-out").replaceChildren(h("table", { class: "data" },
    h("tr", {}, h("td", {}, "レイノルズ数 Re = ρUa/μ"), h("td", {}, fmt(Re, 3))),
    h("tr", {}, h("td", {}, "流れが a だけ進む時間 a/U"), h("td", {}, unit(tauAdv, TIME))),
    h("tr", {}, h("td", {}, "運動量が a だけ拡散する時間 ρa²/μ"), h("td", {}, unit(tauNu, TIME))),
    h("tr", {}, h("td", {}, "粒子の慣性が消える時間 m/(6πμa)"), h("td", {}, unit(tauP, TIME))),
    h("tr", {}, h("td", {}, "力を切った後に惰性で進む距離 ÷ a"), h("td", {}, fmt((re.U * tauP) / re.a, 3)))),
    h("p", { class: "caption" }, `${verdict} 運動量の拡散時間と a/U の比が Re である。粒子の密度は ${RHO_P} kg/m³ とした。`));
}

// ---------------------------------------------------------------------
// Reversibility: a dyed blob in simple shear, run forward then backward
// ---------------------------------------------------------------------
const blob = [];
for (let i = -12; i <= 12; i++)
  for (let j = -12; j <= 12; j++) {
    const x = (0.35 * i) / 12, z = (0.35 * j) / 12;
    if (x * x + z * z <= 0.35 * 0.35) blob.push([x, z]);
  }
const rev = { strain: 6, D: 0, t: 0, frames: null, anim: null };
const STEPS = 150;
function recompute() { rev.frames = shearReversal(blob, { strain: rev.strain, steps: STEPS, D: rev.D, seed: 7 }); }
const tSlider = slider({ label: "時間（前進 → 逆転）", min: 0, max: 2 * STEPS, step: 1, value: 0, format: (v) => (v <= STEPS ? `前進 ${((rev.strain * v) / STEPS).toFixed(1)}` : `逆転 ${((rev.strain * (2 * STEPS - v)) / STEPS).toFixed(1)}`), onInput: (v) => { stopRev(); rev.t = v; renderRev(); } });
const playRev = h("button", { class: "action", type: "button", onclick: () => (rev.anim ? stopRev() : startRev()) }, "再生");
document.getElementById("rev-controls").append(
  slider({ label: "ひずみの大きさ γ", min: 1, max: 10, step: 0.5, value: rev.strain, format: (v) => v.toFixed(1), onInput: (v) => { rev.strain = v; recompute(); tSlider.set(rev.t); renderRev(); } }),
  field("ブラウン拡散 D/(γ̇L²)", segmented([[0, "なし"], [1e-5, "10⁻⁵"], [1e-4, "10⁻⁴"], [1e-3, "10⁻³"]], rev.D, (v) => { rev.D = v; recompute(); renderRev(); })),
  tSlider, playRev);
const revScene = createScene(document.getElementById("rev-scene"), { width: 520, height: 300 });

function startRev() {
  if (rev.t >= 2 * STEPS) rev.t = 0;
  playRev.textContent = "停止";
  let last = null;
  const step = (now) => {
    if (last !== null) rev.t = Math.min(2 * STEPS, rev.t + Math.max(1, Math.round(((now - last) / 1000) * 60)));
    last = now;
    tSlider.set(rev.t);
    renderRev();
    if (rev.t >= 2 * STEPS) { stopRev(); return; }
    if (rev.anim) rev.anim = requestAnimationFrame(step);
  };
  rev.anim = requestAnimationFrame(step);
}
function stopRev() {
  if (rev.anim) cancelAnimationFrame(rev.anim);
  rev.anim = null;
  playRev.textContent = "再生";
}

function renderRev() {
  const pts = rev.frames[rev.t];
  const forward = rev.t <= STEPS;
  const arrows = [];
  for (let z = -0.9; z <= 0.91; z += 0.3) arrows.push({ x: -3.0, z, vx: (forward ? 1 : -1) * 0.5 * z, vz: 0, color: "var(--line)", width: 1.3 });
  revScene.draw({
    view: { cx: 0, cz: 0, span: 8 },
    lines: [{ x1: -4, z1: 1, x2: 4, z2: 1, color: "var(--ink)", dash: "0" }, { x1: -4, z1: -1, x2: 4, z2: -1, color: "var(--ink)", dash: "0" }],
    spheres: pts.map(([x, z], k) => ({ x, z, a: 0.025, fill: blob[k][0] < 0 ? "var(--c1)" : "var(--c2)", stroke: "none" })),
    arrows: [...arrows, { x: 2.8, z: 1.12, vx: forward ? 0.8 : -0.8, vz: 0, color: "var(--ink)", width: 2 }],
    axes: { x: "x", up: "z", out: null },
  });
  const back = rmsDistance(rev.frames.at(-1), blob);
  document.getElementById("rev-out").replaceChildren(h("p", { class: "caption" },
    `上の板の向き：${forward ? "右へ（前進）" : "左へ（逆転）"}。最後まで戻したときの、最初の位置からのずれ（二乗平均）：${fmt(back, 2)}` +
    (rev.D === 0 ? "（丸め誤差だけ。流れの履歴は完全に巻き戻る）" : "（拡散で失われた分は、流れを逆にしても戻らない）")));
}

// ---------------------------------------------------------------------
// Dissipation: local 2 mu E:E, and the total around a sphere
// ---------------------------------------------------------------------
const DPRESETS = {
  shear: { name: "単純せん断", L: [1, 0, 0] },
  ext: { name: "純粋な伸長", L: [0, 0, 0.5] },
  rot: { name: "剛体回転", L: [1, -1, 0] },
};
const dis = { L: DPRESETS.shear.L.slice() }; // [L12, L21, L11 (= -L22)]
let dSeg;
const dSliders = ["L₁₂ = ∂u₁/∂x₂", "L₂₁ = ∂u₂/∂x₁", "L₁₁ = −L₂₂"].map((label, k) =>
  slider({ label, min: -1.5, max: 1.5, step: 0.05, value: dis.L[k], format: (v) => v.toFixed(2), onInput: (v) => { dis.L[k] = v; dSeg.set(null); renderDis(); } }));
document.getElementById("dis-controls").append(
  field("例", dSeg = segmented(Object.entries(DPRESETS).map(([k, p]) => [k, p.name]), "shear", (v) => { dis.L = DPRESETS[v].L.slice(); dSliders.forEach((s, k) => s.set(dis.L[k])); renderDis(); })),
  ...dSliders);
function renderDis() {
  const [l12, l21, l11] = dis.L;
  const L = [[l11, l12], [l21, -l11]];
  const E = [[l11, (l12 + l21) / 2], [(l12 + l21) / 2, -l11]], w = (l12 - l21) / 2;
  const phi = dissipationDensity(L, 1);
  const mat = (M) => `[[${M.map((r) => r.map((v) => fmt(v, 3)).join(", ")).join("], [")}]]`;
  document.getElementById("dis-out").replaceChildren(h("table", { class: "data" },
    h("tr", {}, h("td", {}, "変形速度 E"), h("td", {}, mat(E))),
    h("tr", {}, h("td", {}, "回転 W の成分 W₁₂"), h("td", {}, fmt(w, 3))),
    h("tr", {}, h("td", {}, "散逸 Φ = 2μ E:E"), h("td", {}, fmt(phi, 4))),
    h("tr", {}, h("td", {}, "回転の部分 2μ W:W（散逸に入らない）"), h("td", {}, fmt(4 * w * w, 4)))),
    h("p", { class: "caption" }, Math.abs(phi) < 1e-12 ? "剛体回転では流体は変形しないので、散逸は 0。" : "散逸は E だけで決まり、W をいくら変えても変わらない。"));
}

const sphereDis = (() => {
  const Rs = Array.from({ length: 31 }, (_, k) => 10 ** ((3 * k) / 30));
  const tr = (r) => translationField(r, [0, 0, 1]), ro = (r) => rotationField(r, [0, 1, 0]);
  const totT = dissipationIntegral(tr), totR = dissipationIntegral(ro);
  return { Rs, totT, totR, fT: Rs.map((R) => dissipationIntegral(tr, { R }) / totT), fR: Rs.map((R) => dissipationIntegral(ro, { R }) / totR) };
})();
const disPlot = createPlot(document.getElementById("dis-plot"), { height: 260, xlog: true, xlabel: "半径 R/a", ylabel: "a < r < R での散逸の割合" });
disPlot.update({
  series: [
    { name: "並進（流れは 1/r で減衰）", color: "var(--c1)", points: sphereDis.Rs.map((R, k) => [R, sphereDis.fT[k]]) },
    { name: "回転（流れは 1/r² で減衰）", color: "var(--c3)", points: sphereDis.Rs.map((R, k) => [R, sphereDis.fR[k]]) },
  ],
  ydomain: [0, 1.02], hlines: [{ y: 1 }],
});
const at10 = dissipationIntegral((r) => translationField(r, [0, 0, 1]), { R: 10 }) / sphereDis.totT;
document.getElementById("dis-sphere").replaceChildren(h("table", { class: "data" },
  h("tr", {}, h("th", {}, ""), h("th", {}, "流体全体の散逸（数値積分）"), h("th", {}, "球がする仕事率")),
  h("tr", {}, h("td", {}, "並進 U = e_z"), h("td", {}, fmt(sphereDis.totT, 8)), h("td", {}, `6πμaU² = ${fmt(6 * Math.PI, 8)}`)),
  h("tr", {}, h("td", {}, "回転 Ω = e_y"), h("td", {}, fmt(sphereDis.totR, 8)), h("td", {}, `8πμa³Ω² = ${fmt(8 * Math.PI, 8)}`))),
  h("p", { class: "caption" }, `並進では、半径 10a の内側に入る散逸は ${(100 * at10).toFixed(1)}% にすぎない。点力の流れが 1/r でしか減衰しないことの表れである（ページ5）。`));

recompute();
renderRe();
renderRev();
renderDis();
