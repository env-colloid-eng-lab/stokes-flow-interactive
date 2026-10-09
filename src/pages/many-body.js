import { initPage, h, field, segmented, slider, fmt, frameThrottle } from "../ui/page.js";
import { createMatrixView } from "../ui/matrixView.js";
import { createPlot } from "../ui/plot.js";
import { createScene } from "../ui/scene2d.js";
import { solve, sub, eigvalsSym, vsub, norm, matvec } from "../core/linalg.js";
import { mulberry32 } from "../core/random.js";
import { mobility, mobilityVelocity, flatten, unflatten } from "../physics/rpy.js";
import { heunStep } from "../physics/integrate.js";
import { resistance, pairwiseSumResistance, axialManyBody, invertMobility, axialComponents } from "../models/manyBody.js";

initPage("p10");

const a = 1, mu = 1;
const unit = 6 * Math.PI * mu * a;

// =====================================================================
// Coaxial three spheres: exact (collocation) vs pairwise sum
// =====================================================================
const ax = { g12: 1, g23: 1, L: 24 };
const axialScene = createScene(document.getElementById("axial-scene"), { width: 300, height: 340 });
document.getElementById("axial-controls").append(
  slider({ label: "隙間 h₁₂/a", min: 0.2, max: 4, step: 0.1, value: ax.g12, format: (v) => v.toFixed(1), onInput: (v) => { ax.g12 = v; scheduleAxial(); } }),
  slider({ label: "隙間 h₂₃/a", min: 0.2, max: 4, step: 0.1, value: ax.g23, format: (v) => v.toFixed(1), onInput: (v) => { ax.g23 = v; scheduleAxial(); } }),
  field("次数 L", segmented([[16, "16"], [24, "24"], [32, "32"]], ax.L, (v) => { ax.L = v; scheduleAxial(); })));

const zz = axialComponents;

// Slider input events can arrive faster than the solves; render at most once per frame.
const scheduleAxial = frameThrottle(() => renderAxial());
const pairCache = new Map(); // two-body solutions depend only on the gap and L

function renderAxial() {
  const z = [0, 2 * a + ax.g12, 4 * a + ax.g12 + ax.g23];
  const exact = axialManyBody(z, [a, a, a], { L: ax.L, mu, pairCache });
  const X = z.map((zi) => [0, 0, zi]);
  const rpy = { R: zz(resistance(X)), R2B: zz(pairwiseSumResistance(X)) };
  const f = [-1, -1, -1];
  const vel = (R) => solve(R, f).map((u) => u * unit); // in units of the isolated settling speed m0 F
  const methods = [
    { name: "厳密（境界条件解法）", R: exact.R, R2B: exact.R2B },
    { name: "RPY の逆行列", R: rpy.R, R2B: rpy.R2B },
  ].map((m) => ({ ...m, U: vel(m.R), U2B: vel(m.R2B) }));
  const entries = [["R₁₁", 0, 0], ["R₂₂", 1, 1], ["R₁₂", 0, 1], ["R₁₃", 0, 2]];
  const table = h("table", { class: "data" },
    h("tr", {}, h("th", {}, ""), ...methods.flatMap((m) => [h("th", {}, `${m.name}：R`), h("th", {}, "R²ᴮ"), h("th", {}, "差")])));
  for (const [lab, i, j] of entries)
    table.append(h("tr", {}, h("td", {}, lab), ...methods.flatMap((m) => {
      const v = m.R[i][j] / unit, w = m.R2B[i][j] / unit, d = v - w;
      return [h("td", {}, fmt(v)), h("td", {}, fmt(w)), h("td", { style: { fontWeight: Math.abs(d) > 0.05 * Math.abs(w) ? "700" : "" } }, fmt(d, 3))];
    })));
  const vtable = h("table", { class: "data", style: { marginTop: "10px" } },
    h("tr", {}, h("th", {}, "沈降速度 U/(m₀F)"), ...methods.flatMap((m) => [h("th", {}, `${m.name}：R`), h("th", {}, "R²ᴮ")])));
  for (let k = 0; k < 3; k++)
    vtable.append(h("tr", {}, h("td", {}, `球 ${k + 1}`), ...methods.flatMap((m) => [h("td", {}, fmt(-m.U[k])), h("td", {}, fmt(-m.U2B[k]))])));
  document.getElementById("axial-tables").replaceChildren(h("div", { style: { overflowX: "auto" } }, table, vtable));

  const shield = exact.R[0][2] / exact.R2B[0][2];
  document.getElementById("axial-note").textContent =
    `R₁₃ は二体の和の ${(100 * shield).toFixed(0)}%。間の球 2 が球 1 と球 3 の相互作用を遮っている（遮蔽）。` +
    `境界残差の最大値 ${fmt(exact.boundaryError, 2)}（L = ${ax.L}）。`;

  axialScene.draw({
    view: { cx: 0, cz: z[1], span: Math.max(8, z[2] - z[0] + 5) * (300 / 340) },
    spheres: z.map((zi, k) => ({ x: 0, z: zi, a, label: String(k + 1) })),
    arrows: z.map((zi) => ({ x: 0, z: zi, vx: 0, vz: -1.4, color: "var(--c2)" })),
    lines: [{ x1: 0, z1: -100, x2: 0, z2: 100 }],
  });
}

// =====================================================================
// Sedimentation of N spheres in the x-z plane
// =====================================================================
const SCEN = {
  three: { name: "三球", make: () => [[-3, 0, 0], [3, 0, 0], [0, 0, 4]] },
  row: { name: "横一列（4球）", make: () => [-4.5, -1.5, 1.5, 4.5].map((x) => [x, 0, 0]) },
  cloud: { name: "雲（8球）", make: () => cloud(8, 5, 2.6, 7) },
};
function cloud(n, R, dmin, seed) {
  const rnd = mulberry32(seed), pts = [];
  while (pts.length < n) {
    const p = [(2 * rnd() - 1) * R, 0, (2 * rnd() - 1) * R];
    if (Math.hypot(p[0], p[2]) <= R && pts.every((q) => norm(vsub(p, q)) >= dmin)) pts.push(p);
  }
  return pts;
}

const sed = { scen: "three", model: "rpy", show: "R", X: null, t: 0, trails: [], hist: [], anim: null, stopped: "" };
function resetSed() {
  stopSed();
  sed.X = SCEN[sed.scen].make();
  sed.t = 0;
  sed.trails = sed.X.map((p) => [[p[0], p[2]]]);
  sed.hist = [];
  sed.stopped = "";
  const n = 3 * sed.X.length;
  const lab = sed.X.flatMap((_, p) => ["x", "y", "z"].map((c) => `${p + 1}${c}`));
  sedMatrix.update(Array.from({ length: n }, () => new Array(n).fill(0)), { labels: lab });
  sedMatrix.select(2, 5);
  renderSed();
}

const sc = document.getElementById("sed-controls");
sc.append(
  field("配置", segmented(Object.entries(SCEN).map(([k, s]) => [k, s.name]), sed.scen, (v) => { sed.scen = v; resetSed(); })),
  field("近似", segmented([["rpy", "RPY"], ["oseen", "Oseen"]], sed.model, (v) => { sed.model = v; resetSed(); })));
document.getElementById("sed-matrix-controls").append(
  field("表示", segmented([["R", "R"], ["R2B", "二体の和 R²ᴮ"], ["diff", "差 R − R²ᴮ"]], sed.show, (v) => { sed.show = v; renderSed(); })));

const sedScene = createScene(document.getElementById("sed-scene"), { width: 480, height: 360 });
const sedMatrix = createMatrixView(document.getElementById("sed-matrix"), {
  explain: (i, j) => explainSed(i, j),
  onSelect: () => { sed.hist = []; renderSed(); },
});
const sedSeries = createPlot(document.getElementById("sed-series"), { height: 240, xlabel: "時間 t（μa²/F₀）", ylabel: "成分 ÷ 6πμa" });
const sedEigen = createPlot(document.getElementById("sed-eigen"), { height: 240, xlabel: "時間 t（μa²/F₀）", ylabel: "固有値 × 6πμa" });

const playBtn = h("button", { class: "action", type: "button", onclick: () => (sed.anim ? stopSed() : startSed()) }, "沈降させる");
document.getElementById("sed-anim").append(playBtn, h("button", { class: "action", type: "button", onclick: resetSed }, "リセット"));

const force = () => sed.X.map(() => [0, 0, -1]);
function startSed() {
  if (sed.stopped) {
    document.getElementById("sed-note").textContent = `${sed.stopped} 続けるには「リセット」を押す。`;
    return;
  }
  playBtn.textContent = "停止";
  let last = performance.now();
  const step = (now) => {
    const wall = Math.min(0.05, (now - last) / 1000);
    last = now;
    try {
      advanceSed(30 * wall);
      const c = compute(); // once per frame
      record(c);
      renderSed(c);
    } catch (err) {
      // e.g. a singular Oseen mobility: stop cleanly instead of leaving a dead animation
      sed.stopped = `計算できなくなったので止めた（${err.message}）。`;
      stopSed();
      document.getElementById("sed-note").textContent = sed.stopped;
      return;
    }
    if (sed.anim) sed.anim = requestAnimationFrame(step);
  };
  sed.anim = requestAnimationFrame(step);
}
function stopSed() {
  if (sed.anim) cancelAnimationFrame(sed.anim);
  sed.anim = null;
  playBtn.textContent = "沈降させる";
}

function minGap(X) {
  let g = Infinity;
  for (let p = 0; p < X.length; p++) for (let q = p + 1; q < X.length; q++) g = Math.min(g, norm(vsub(X[p], X[q])) - 2 * a);
  return g;
}

function advanceSed(dt) {
  if (!(dt > 0)) return; // the first animation frame has no elapsed time
  const vfun = (X) => mobilityVelocity(X, force(), { a, mu, model: sed.model });
  const nsub = Math.max(1, Math.ceil(dt / 0.25));
  for (let k = 0; k < nsub; k++) {
    const next = heunStep(sed.X, vfun, dt / nsub);
    const g = minGap(next);
    if (g < 0.05) {
      // keep the last state that was still separated
      sed.stopped = `球どうしの隙間が ${g > 0 ? g.toPrecision(2) + "a" : "なくなる"}ところまで近づいたので、その手前で止めた。ここから先は近接の扱い（9.）が必要になる。`;
      stopSed();
      break;
    }
    sed.X = next;
    sed.t += dt / nsub;
  }
  sed.X.forEach((p, i) => { sed.trails[i].push([p[0], p[2]]); if (sed.trails[i].length > 1500) sed.trails[i].shift(); });
}

let cur = null;
function compute() {
  const M = mobility(sed.X, { a, mu, model: sed.model }); // built once, reused for R, eigenvalues and U
  const R = invertMobility(M);
  const R2B = pairwiseSumResistance(sed.X, { a, mu, model: sed.model });
  const ev = eigvalsSym(M);
  const scale = (A) => A.map((row) => row.map((v) => v / unit));
  cur = {
    R: scale(R), R2B: scale(R2B), diff: scale(sub(R, R2B)), evMin: ev[0] * unit, evMax: ev.at(-1) * unit,
    U: unflatten(matvec(M, flatten(force()))),
  };
  return cur;
}

function record(c) {
  const sel = sedMatrix.selected;
  if (!sel) return;
  const [i, j] = sel;
  sed.hist.push({ t: sed.t, R: c.R[i][j], R2B: c.R2B[i][j], d: c.diff[i][j], evMin: c.evMin, evMax: c.evMax });
  if (sed.hist.length > 2000) sed.hist.shift();
}

function explainSed(i, j) {
  if (!cur) return "";
  const lab = (k) => `${Math.floor(k / 3) + 1}${"xyz"[k % 3]}`;
  const p = Math.floor(i / 3), q = Math.floor(j / 3);
  const wrap = h("div", {},
    h("div", {}, h("strong", {}, `${lab(i)}, ${lab(j)} 成分`), "（÷6πμa）"),
    h("div", { class: "val" }, `R = ${fmt(cur.R[i][j])}　R²ᴮ = ${fmt(cur.R2B[i][j])}　差 = ${fmt(cur.diff[i][j], 3)}`));
  const others = sed.X.map((_, k) => k + 1).filter((k) => k !== p + 1 && k !== q + 1);
  wrap.append(h("div", { class: "caption" }, p === q
    ? `自己ブロック。R²ᴮ は、球 ${p + 1} と他の各球を二球だけにしたときの補正の和。差は ${others.length ? `球 ${others.join("・")} のうち二つ以上を経由する反射` : "なし"}。`
    : `相互ブロック。R²ᴮ は球 ${p + 1}・${q + 1} の二球だけの抵抗。差は球 ${others.join("・") || "（なし）"} を経由する反射の寄与。`));
  return wrap;
}

function renderSed(c = compute()) {
  const X = sed.X;
  const cx = X.reduce((s, p) => s + p[0], 0) / X.length, cz = X.reduce((s, p) => s + p[2], 0) / X.length;
  const span = Math.max(14, ...X.map((p) => 2 * Math.abs(p[0] - cx) + 6), ...X.map((p) => (2 * Math.abs(p[2] - cz) + 6) * (480 / 360)));
  const U = c.U;
  const colors = ["var(--c1)", "var(--c2)", "var(--c3)", "var(--c4)", "var(--c5)", "var(--c6)"];
  sedScene.draw({
    view: { cx, cz, span },
    spheres: X.map((p, k) => ({ x: p[0], z: p[2], a, label: String(k + 1) })),
    arrows: X.map((p, k) => ({ x: p[0], z: p[2], vx: 1.4 * unit * U[k][0], vz: 1.4 * unit * U[k][2], color: "var(--c1)" })),
    paths: sed.trails.map((tr, k) => ({ points: tr, color: colors[k % colors.length], dash: "3 3", width: 1 })),
  });
  const vmean = -U.reduce((s, u) => s + u[2], 0) / U.length * unit;
  document.getElementById("sed-note").textContent = sed.stopped ||
    `青の矢印は速度（孤立球の沈降速度 m₀F を基準）。平均沈降速度は孤立球の ${vmean.toFixed(3)} 倍。t = ${sed.t.toFixed(1)}`;

  const shown = sed.show === "R" ? c.R : sed.show === "R2B" ? c.R2B : c.diff;
  document.getElementById("sed-matrix-title").textContent =
    sed.show === "R" ? "抵抗行列 R/(6πμa)" : sed.show === "R2B" ? "二体の和 R²ᴮ/(6πμa)" : "多体の寄与 (R − R²ᴮ)/(6πμa)";
  sedMatrix.update(shown);

  const hs = sed.hist;
  sedSeries.update({
    series: [
      { name: "R", color: "var(--c1)", points: hs.map((p) => [p.t, p.R]) },
      { name: "R²ᴮ", color: "var(--c2)", dash: "6 4", points: hs.map((p) => [p.t, p.R2B]) },
      { name: "差", color: "var(--c3)", dash: "2 3", points: hs.map((p) => [p.t, p.d]) },
    ],
  });
  sedEigen.update({
    series: [
      { name: "最小固有値", color: "var(--c1)", points: hs.map((p) => [p.t, p.evMin]) },
      { name: "最大固有値", color: "var(--c2)", points: hs.map((p) => [p.t, p.evMax]) },
    ],
    hlines: [{ y: 0 }],
  });
}

renderAxial();
resetSed();
