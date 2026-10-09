import { initPage, h, tex, segmented, slider, fmt } from "../ui/page.js";
import { createMatrixView } from "../ui/matrixView.js";
import { createPlot } from "../ui/plot.js";
import { createScene } from "../ui/scene2d.js";
import { matvec } from "../core/linalg.js";
import { joPolynomials, joXA } from "../physics/jo.js";
import { mobilityCoefficients, resistanceCoefficients, pairMatrix, pairFrame, toPairFrame } from "../models/pairA.js";

initPage("p8");

const a = 1, mu = 1;
const unit = 6 * Math.PI * mu * a; // display: M * unit, R / unit
const polys = joPolynomials(40);

const state = {
  r: 4, thetaDeg: 30, kind: "rpy", show: "R", frame: "lab", forcing: "sediment",
  mode: null, anim: null, t: 0, com: [0, 0], trail: [[0, 0]], hist: [], eigenMode: null,
};

const LAB = ["1x", "1y", "1z", "2x", "2y", "2z"];
const PAIR = ["1∥", "1⊥", "1⊥′", "2∥", "2⊥", "2⊥′"];
const DIRS = ["∥", "⊥", "⊥′"];

// ---------- model ----------
function current() {
  const theta = (state.thetaDeg * Math.PI) / 180;
  const frame = pairFrame(theta);
  const mc = mobilityCoefficients(state.r, { a, mu, kind: state.kind });
  const rc = resistanceCoefficients(state.r, { a, mu, kind: state.kind });
  const coeffs = state.show === "M"
    ? { par: { self: mc.par.self * unit, cross: mc.par.cross * unit }, perp: { self: mc.perp.self * unit, cross: mc.perp.cross * unit } }
    : { par: { self: rc.R.par.self / unit, cross: rc.R.par.cross / unit }, perp: { self: rc.R.perp.self / unit, cross: rc.R.perp.cross / unit } };
  const lab = pairMatrix(coeffs, frame.e);
  const shown = state.frame === "lab" ? lab : toPairFrame(lab, frame.axes);
  const Mdim = pairMatrix(mc, frame.e);
  return { theta, frame, mc, rc, coeffs, lab, shown, Mdim };
}

function forces(frame) {
  if (state.eigenMode) {
    const { dir, sign } = state.eigenMode;
    const v = frame.axes[dir];
    return [...v, ...v.map((x) => sign * x)].map((x) => x / Math.SQRT2);
  }
  if (state.forcing === "sediment") return [0, 0, -1, 0, 0, -1];
  const e = frame.e;
  return [...e, ...e.map((x) => -x)]; // push the spheres toward each other
}

// ---------- controls ----------
const controls = document.getElementById("controls");
const rSlider = slider({ label: "中心間距離 r/a", min: 2.02, max: 12, step: 0.01, value: state.r, format: (v) => v.toFixed(2), onInput: (v) => { state.r = v; render(); } });
const thSlider = slider({ label: "対の向き θ", min: 0, max: 180, step: 1, value: state.thetaDeg, format: (v) => `${v}°`, onInput: (v) => { state.thetaDeg = v; render(); } });
controls.append(rSlider, thSlider,
  h("label", {}, "近似", segmented([["rpy", "RPY"], ["oseen", "Oseen"]], state.kind, (v) => { state.kind = v; render(); })),
  h("label", {}, "力", segmented([["sediment", "同じ向き（沈降）"], ["approach", "逆向き（中心線方向）"]], state.forcing, (v) => { state.forcing = v; state.eigenMode = null; render(); })));

const mctl = document.getElementById("matrix-controls");
mctl.append(
  h("label", {}, "表示", segmented([["M", "移動度 M"], ["R", "抵抗 R"]], state.show, (v) => { state.show = v; state.hist = []; render(); })),
  h("label", {}, "座標系", segmented([["lab", "実験室系"], ["pair", "対の局所系"]], state.frame, (v) => { state.frame = v; state.hist = []; render(); })));

const scene = createScene(document.getElementById("scene"));
const matrixView = createMatrixView(document.getElementById("matrix"), {
  labels: LAB,
  explain: (i, j, v) => explainEntry(i, j, v),
  onSelect: () => { state.hist = []; render(); },
});
matrixView.select(0, 3);
const seriesPlot = createPlot(document.getElementById("series"), { height: 240, xlabel: "時間 t（μa²/F₀）", ylabel: "成分（無次元）" });
const distPlot = createPlot(document.getElementById("distance"), { height: 240, xlabel: "中心間距離 r/a", ylabel: "係数（無次元）" });

// ---------- animation ----------
const anim = document.getElementById("anim");
const playButtons = [
  ["sediment", "沈降させる"], ["rotate", "対を回す"], ["approach", "近づける"],
].map(([mode, label]) => h("button", { class: "action", type: "button", onclick: () => toggle(mode) }, label));
const status = h("span", { class: "caption" });
anim.append(...playButtons, h("button", { class: "action", type: "button", onclick: reset }, "リセット"), status);

function toggle(mode) {
  if (state.mode === mode) { stop(); return; }
  stop();
  state.mode = mode;
  state.hist = [];
  state.t = 0;
  if (mode === "sediment") { state.forcing = "sediment"; state.eigenMode = null; }
  if (mode === "approach") { state.forcing = "approach"; state.eigenMode = null; }
  let last = performance.now();
  const step = (now) => {
    const wall = Math.min(0.05, (now - last) / 1000);
    last = now;
    advance(wall);
    render();
    if (state.mode) state.anim = requestAnimationFrame(step);
  };
  state.anim = requestAnimationFrame(step);
  playButtons.forEach((b, i) => (b.textContent = ["sediment", "rotate", "approach"][i] === mode ? "停止" : ["沈降させる", "対を回す", "近づける"][i]));
}
function stop() {
  if (state.anim) cancelAnimationFrame(state.anim);
  state.mode = null; state.anim = null;
  playButtons.forEach((b, i) => (b.textContent = ["沈降させる", "対を回す", "近づける"][i]));
}
function reset() {
  stop();
  Object.assign(state, { r: 4, thetaDeg: 30, t: 0, com: [0, 0], trail: [[0, 0]], hist: [], eigenMode: null });
  rSlider.set(state.r); thSlider.set(state.thetaDeg);
  status.textContent = "";
  render();
}

function advance(wall) {
  const dt = 20 * wall; // simulated time per wall-clock second
  const c = current();
  if (state.mode === "sediment") {
    const U = matvec(c.Mdim, forces(c.frame));
    state.com = [state.com[0] + dt * (U[0] + U[3]) / 2, state.com[1] + dt * (U[2] + U[5]) / 2];
    state.trail.push(state.com.slice());
    if (state.trail.length > 2000) state.trail.shift();
    status.textContent = "配置は変わらないので、R も速度も一定のまま。";
  } else if (state.mode === "rotate") {
    state.thetaDeg = (state.thetaDeg + 30 * wall) % 180;
    thSlider.set(Math.round(state.thetaDeg));
    status.textContent = "距離が一定なので X, Y は一定。実験室系の成分だけが変わる。";
  } else if (state.mode === "approach") {
    const U = matvec(c.Mdim, forces(c.frame));
    const e = c.frame.e;
    const vrel = (U[0] - U[3]) * e[0] + (U[1] - U[4]) * e[1] + (U[2] - U[5]) * e[2];
    state.r = Math.max(2.02, state.r - dt * vrel);
    rSlider.set(state.r);
    status.textContent = state.r <= 2.02 ? "隙間 h = 0.02a で止めた。この先は 9. で扱う。" : "距離が縮まり、X と Y が変わる。";
    if (state.r <= 2.02) { state.t += dt; record(); stop(); return; }
  }
  state.t += dt;
  record();
}

function record() {
  const sel = matrixView.selected;
  if (!sel) return;
  const c = current();
  const [i, j] = sel;
  const key = Math.floor(i / 3) === Math.floor(j / 3) ? "self" : "cross";
  state.hist.push({ t: state.t, v: c.shown[i][j], X: c.coeffs.par[key], Y: c.coeffs.perp[key] });
}

// ---------- explanations ----------
function explainEntry(i, j, v) {
  const c = current();
  const al = Math.floor(i / 3), be = Math.floor(j / 3), k = i % 3, l = j % 3;
  const key = al === be ? "self" : "cross";
  const P = c.coeffs.par[key], Q = c.coeffs.perp[key];
  const isM = state.show === "M";
  const wrap = h("div");
  const sym = isM ? "M" : "R";
  const parName = isM ? (key === "self" ? "m_0" : "c_\\parallel") : (key === "self" ? "X_{11}" : "X_{12}");
  const perpName = isM ? (key === "self" ? "m_0" : "c_\\perp") : (key === "self" ? "Y_{11}" : "Y_{12}");
  const norm = isM ? "6\\pi\\mu a\\,M" : "R/(6\\pi\\mu a)";
  if (state.frame === "lab") {
    const ee = c.frame.e[k] * c.frame.e[l];
    const d = (k === l ? 1 : 0) - ee;
    wrap.append(
      h("div", {}, h("strong", {}, `${LAB[i]}, ${LAB[j]} 成分`), "（", tex(norm), " の値）"),
      tex(`${sym}_{${al + 1}${"xyz"[k]},${be + 1}${"xyz"[l]}} = ${parName}\\,e_${"xyz"[k]}e_${"xyz"[l]} + ${perpName}\\,(\\delta_{${"xyz"[k]}${"xyz"[l]}}-e_${"xyz"[k]}e_${"xyz"[l]})`, { display: true }),
      h("div", { class: "val" }, `= ${fmt(P)} × ${fmt(ee)} + ${fmt(Q)} × ${fmt(d)} = ${fmt(v)}`),
      h("div", { class: "caption" }, sourceText(key)));
  } else {
    const same = k === l;
    wrap.append(h("div", {}, h("strong", {}, `${PAIR[i]}, ${PAIR[j]} 成分`), "（", tex(norm), " の値）"));
    if (!same) wrap.append(h("div", {}, "対の局所系では、方向の違う成分は結合しない（二球を結ぶ軸のまわりの対称性）。値は 0。"));
    else wrap.append(
      tex(`= ${k === 0 ? parName : perpName} = ${fmt(v)}`, { display: true }),
      h("div", { class: "caption" }, `${k === 0 ? "中心線方向" : "中心線に垂直な方向"}の係数。対の向きを回しても変わらず、距離 r だけで決まる。 ${sourceText(key)}`));
  }
  return wrap;
}

function sourceText(key) {
  const approx = state.kind === "rpy" ? "RPY" : "Oseen";
  if (state.show === "M")
    return key === "self" ? "自己ブロック：孤立球の移動度 1/(6πμa)（6. 球の抵抗）。" : `相互ブロック：${approx} 近似の c∥(r), c⊥(r)。`;
  return key === "self"
    ? `自己抵抗：${approx} 移動度の 2×2 逆行列の対角 m₀/(m₀² − c²)。相手がいるだけで 1 からずれる。`
    : `相互抵抗：${approx} 移動度の 2×2 逆行列の非対角 −c/(m₀² − c²)。負の値になる。`;
}

// ---------- rendering ----------
function render() {
  const c = current();
  // scene
  const [cx, cz] = state.com;
  const e = c.frame.e;
  const x1 = [cx - (state.r / 2) * e[0], cz - (state.r / 2) * e[2]];
  const x2 = [cx + (state.r / 2) * e[0], cz + (state.r / 2) * e[2]];
  const f = forces(c.frame);
  const U = matvec(c.Mdim, f);
  const us = 1.4 * unit, fs = 1.4; // arrow scales
  scene.draw({
    view: { cx, cz, span: Math.max(9, state.r + 5) },
    spheres: [{ x: x1[0], z: x1[1], a, label: "1" }, { x: x2[0], z: x2[1], a, label: "2" }],
    arrows: [
      { x: x1[0], z: x1[1], vx: fs * f[0], vz: fs * f[2], color: "var(--c2)", offset: -5 },
      { x: x2[0], z: x2[1], vx: fs * f[3], vz: fs * f[5], color: "var(--c2)", offset: -5 },
      { x: x1[0], z: x1[1], vx: us * U[0], vz: us * U[2], color: "var(--c1)", offset: 5 },
      { x: x2[0], z: x2[1], vx: us * U[3], vz: us * U[5], color: "var(--c1)", offset: 5 },
    ],
    paths: state.trail.length > 1 ? [{ points: state.trail, color: "var(--muted)", dash: "3 3" }] : [],
    lines: [{ x1: 0, z1: -1000, x2: 0, z2: 1000 }],
  });
  document.getElementById("scene-note").textContent =
    `重心速度 (x, z) = (${fmt((U[0] + U[3]) / 2 * unit, 3)}, ${fmt((U[2] + U[5]) / 2 * unit, 3)}) × m₀F。灰色の破線は重心の軌跡と鉛直線。`;

  // matrix
  document.getElementById("matrix-title").textContent =
    state.show === "M" ? "移動度行列 6πμa·M（自己ブロックの対角が 1）" : "抵抗行列 R/(6πμa)（孤立球なら対角が 1）";
  matrixView.update(c.shown, { labels: state.frame === "lab" ? LAB : PAIR });

  // time series
  const hist = state.hist;
  seriesPlot.update({
    series: [
      { name: "選んだ成分", color: "var(--c1)", points: hist.map((p) => [p.t, p.v]) },
      { name: "X（中心線方向）", color: "var(--c2)", dash: "5 4", points: hist.map((p) => [p.t, p.X]) },
      { name: "Y（垂直方向）", color: "var(--c3)", dash: "2 3", points: hist.map((p) => [p.t, p.Y]) },
    ],
  });

  // distance dependence
  const rs = Array.from({ length: 200 }, (_, i) => 2.02 + (i * (12 - 2.02)) / 199);
  const fam = rs.map((r) => {
    if (state.show === "M") {
      const m = mobilityCoefficients(r, { a, mu, kind: state.kind });
      return { r, a1: m.par.cross * unit, a2: m.perp.cross * unit };
    }
    const rr = resistanceCoefficients(r, { a, mu, kind: state.kind }).R;
    const jo = joXA(r / a, 1, polys);
    return { r, X11: rr.par.self / unit, X12: rr.par.cross / unit, Y11: rr.perp.self / unit, Y12: rr.perp.cross / unit, J11: jo.x11, J12: jo.x12 };
  });
  const tag = state.kind === "rpy" ? "RPY" : "Oseen";
  if (state.show === "M") {
    distPlot.update({
      series: [
        { name: `6πμa·c∥（${tag}）`, color: "var(--c2)", points: fam.map((p) => [p.r, p.a1]) },
        { name: `6πμa·c⊥（${tag}）`, color: "var(--c3)", points: fam.map((p) => [p.r, p.a2]) },
      ],
      vlines: [{ x: state.r }],
    });
    document.getElementById("distance-caption").textContent = "相互移動度は 1/r でゆっくり減衰する。遠く離れても相互作用は消えにくい。";
  } else {
    const nearOk = (p) => p.r > 2.3; // JO series (K = 40) is shown only where it has converged well
    distPlot.update({
      series: [
        { name: `X₁₁（${tag}⁻¹）`, color: "var(--c2)", points: fam.map((p) => [p.r, p.X11]) },
        { name: `X₁₂（${tag}⁻¹）`, color: "var(--c2)", dash: "6 4", points: fam.map((p) => [p.r, p.X12]) },
        { name: `Y₁₁（${tag}⁻¹）`, color: "var(--c3)", points: fam.map((p) => [p.r, p.Y11]) },
        { name: `Y₁₂（${tag}⁻¹）`, color: "var(--c3)", dash: "6 4", points: fam.map((p) => [p.r, p.Y12]) },
        { name: "X^A₁₁, X^A₁₂（JO 級数）", color: "var(--ink)", width: 1.2, dash: "1 3", points: fam.map((p) => [p.r, nearOk(p) ? p.J11 : NaN]) },
        { color: "var(--ink)", width: 1.2, dash: "1 3", points: fam.map((p) => [p.r, nearOk(p) ? p.J12 : NaN]) },
      ],
      vlines: [{ x: state.r }],
    });
    document.getElementById("distance-caption").textContent =
      "点線は二球の厳密な軸方向抵抗（JO 級数、r/a > 2.3 で表示）。遠方では RPY の逆行列と重なり、近づくほど離れる。差の大きさは 9. で調べる。";
  }

  renderEigen(c);
}

function renderEigen(c) {
  const box = document.getElementById("eigen");
  const rows = [];
  for (let d = 0; d < 3; d++)
    for (const sign of [1, -1]) {
      const fam = d === 0 ? "par" : "perp";
      const lam = c.coeffs[fam].self + sign * c.coeffs[fam].cross;
      rows.push({ d, sign, lam, label: `${sign > 0 ? "集団" : "相対"}（${DIRS[d]}）` });
    }
  rows.sort((p, q) => p.lam - q.lam);
  const table = h("table", { class: "data" },
    h("tr", {}, h("th", {}, "モード"), h("th", {}, "力の組（球1, 球2）"), h("th", {}, state.show === "M" ? "固有値 ×6πμa" : "固有値 ÷6πμa")));
  for (const row of rows) {
    const active = state.eigenMode && state.eigenMode.dir === row.d && state.eigenMode.sign === row.sign;
    const tr = h("tr", { style: { cursor: "pointer", background: active ? "var(--accent-soft)" : "" } },
      h("td", {}, row.label), h("td", {}, row.sign > 0 ? `(${DIRS[row.d]}, ${DIRS[row.d]})` : `(${DIRS[row.d]}, −${DIRS[row.d]})`), h("td", {}, fmt(row.lam)));
    tr.addEventListener("click", () => { stop(); state.eigenMode = active ? null : { dir: row.d, sign: row.sign }; render(); });
    table.append(tr);
  }
  const min = Math.min(...rows.map((r) => r.lam));
  box.replaceChildren(table, h("p", { class: "caption" }, min > 0 ? "全ての固有値が正（正定値）。" : "負の固有値がある。この近似は散逸の正値性を破っている。"));
}

render();
