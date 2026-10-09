import { initPage, h, field, segmented, slider, fmt } from "../ui/page.js";
import { createPlot } from "../ui/plot.js";
import { matmul, transpose } from "../core/linalg.js";
import { mobility } from "../physics/rpy.js";
import { joPolynomials } from "../physics/jo.js";
import { resistanceCoefficients } from "../models/pairA.js";
import { noiseFactor, sampleDisplacements, covariance, hinderedMobility, stepHeights, stationaryDensity } from "../physics/brownian.js";
import { mulberry32, normalSampler } from "../core/random.js";

initPage("p12");

const a = 1, mu = 1, kBT = 1;
const M0 = 1 / (6 * Math.PI * mu * a);

// ---------------------------------------------------------------------
// Correlated displacements of two spheres (RPY, frozen configuration)
// ---------------------------------------------------------------------
const METHODS = {
  cholesky: "Cholesky 分解 L（正しい）",
  independent: "粒子ごとに独立（誤り）",
  elementwise: "成分ごとの平方根 √M_ij（誤り）",
};
const pair = { r: 2.5, method: "cholesky", seed: 1 };
const NS = 3000;
document.getElementById("pair-controls").append(
  slider({ label: "中心間距離 r/a", min: 2.05, max: 10, step: 0.05, value: pair.r, format: (v) => v.toFixed(2), onInput: (v) => { pair.r = v; renderPair(); } }),
  field("雑音の作り方", segmented(Object.entries(METHODS), pair.method, (v) => { pair.method = v; renderPair(); })),
  h("button", { class: "action", type: "button", onclick: () => { pair.seed += 1; renderPair(); } }, "引き直す"));
const scatter = createPlot(document.getElementById("pair-scatter"), { width: 420, height: 420, xlabel: "球 1 の変位 ΔX₁（中心線方向）", ylabel: "球 2 の変位 ΔX₂（中心線方向）" });

function ellipse(C, k = 2) {
  // points x with x^T C^{-1} x = k^2, from the Cholesky factor of the 2x2 covariance
  const l11 = Math.sqrt(C[0][0]), l21 = C[1][0] / l11, l22 = Math.sqrt(Math.max(0, C[1][1] - l21 * l21));
  return Array.from({ length: 121 }, (_, i) => { const t = (2 * Math.PI * i) / 120; return [k * l11 * Math.cos(t), k * (l21 * Math.cos(t) + l22 * Math.sin(t))]; });
}

function renderPair() {
  const M = mobility([[0, 0, 0], [pair.r, 0, 0]], { a, mu });
  const unit = 2 * kBT * M0; // dt = 1: displacements in units of sqrt(2 D0 dt)
  const samples = sampleDisplacements(M, { kBT, dt: 1, n: NS, method: pair.method, seed: pair.seed }).map((v) => v.map((x) => x / Math.sqrt(unit)));
  const { C } = covariance(samples);
  const B = noiseFactor(M, pair.method), gen = matmul(B, transpose(B)).map((row) => row.map((v) => v / M0));
  const target = M.map((row) => row.map((v) => v / M0));
  const sub = (A) => [[A[0][0], A[0][3]], [A[3][0], A[3][3]]]; // x components of spheres 1 and 2
  scatter.update({
    series: [
      { name: `${NS} 回の変位`, color: "var(--c1)", marker: true, r: 1.6, opacity: 0.5, points: samples.map((v) => [v[0], v[3]]) },
      { name: "正しい共分散の 2σ 楕円", color: "var(--c2)", width: 2.5, points: ellipse(sub(target)) },
      ...(pair.method === "cholesky" ? [] : [{ name: "この方法が作る共分散の 2σ 楕円", color: "var(--ink)", dash: "5 4", points: ellipse(sub(gen)) }]),
    ],
    xdomain: [-4, 4], ydomain: [-4, 4],
  });
  const rel = (A) => A[0][0] + A[3][3] - 2 * A[0][3];
  const row = (label, f) => h("tr", {}, h("td", {}, label), h("td", {}, fmt(f(C), 3)), h("td", {}, fmt(f(gen), 3)), h("td", {}, fmt(f(target), 3)));
  document.getElementById("pair-out").replaceChildren(h("table", { class: "data" },
    h("tr", {}, h("th", {}, "（2D₀Δt を単位に）"), h("th", {}, "標本"), h("th", {}, "この方法の B Bᵀ"), h("th", {}, "正しい値 M/M₀")),
    row("球 1 の分散 ⟨ΔX₁²⟩", (A) => A[0][0]),
    row("相関 ⟨ΔX₁ΔX₂⟩", (A) => A[0][3]),
    row("相対変位の分散 ⟨(ΔX₂−ΔX₁)²⟩", rel),
    row("横方向の相関 ⟨ΔY₁ΔY₂⟩", (A) => A[1][4])),
    h("p", { class: "caption" }, pair.method === "cholesky" ? "標本の共分散は、標本数による揺らぎの範囲で M に一致する。二球は同じ向きに動きやすく、相対変位の分散は独立な場合の 2 より小さい。" :
      pair.method === "independent" ? "相関の項が 0 になり、相対変位の分散を過大に見積もる。近づいた二球が、実際より速く離れたり接触したりする。" :
        "√M_ij を並べた行列 B は B Bᵀ = M を満たさない。自己の分散も相関も正しくない。"));
}

// relative diffusion along the line of centres versus the gap
const polys = joPolynomials(60);
const relPlot = createPlot(document.getElementById("rel-plot"), { width: 960, height: 300, xlog: true, xlabel: "すき間 h/a", ylabel: "相対拡散係数 ÷ D₀（中心線方向）" });
{
  const hs = Array.from({ length: 61 }, (_, k) => 10 ** (-2 + (3 * k) / 60));
  const rpy = hs.map((hh) => { const r = 2 + hh; return [hh, 2 * (1 - 3 / (2 * r) + 1 / r ** 3)]; });
  const jo = hs.map((hh) => {
    const c = resistanceCoefficients(2 + hh, { a, mu, axial: "jo", polys });
    // relative mode: R_- = R11 - R12 = 6 pi mu a (X11 - X12), so D_rel / D0 = 2 / (X11 - X12)
    return [hh, c.jo.converged ? 2 / (c.jo.x11 - c.jo.x12) : NaN];
  });
  relPlot.update({
    series: [
      { name: "独立な雑音（相関なし）", color: "var(--muted)", points: hs.map((hh) => [hh, 2]) },
      { name: "RPY", color: "var(--c1)", points: rpy },
      { name: "厳密（JO 級数、収束した範囲）", color: "var(--c2)", width: 3, points: jo },
      { name: "潤滑の最低次 4h/a", color: "var(--c3)", dash: "5 4", points: hs.filter((hh) => hh < 0.3).map((hh) => [hh, 4 * hh]) },
    ],
    ydomain: [0, 2.1],
  });
}

// ---------------------------------------------------------------------
// Sedimentation equilibrium above a wall
// ---------------------------------------------------------------------
const NP = 3000, DT = 2e-3, STEPS_PER_FRAME = 25, ZMAX = 8, BINS = 40;
const sed = { beta: 0.8, drift: true, visc: 1, t: 0, zs: null, anim: null, hist: [], seed: 5, randn: null };
function resetSed() {
  sed.randn = normalSampler(mulberry32(sed.seed));
  const u = mulberry32(sed.seed + 100);
  sed.zs = Array.from({ length: NP }, () => 2 + 2 * u()); // start in a band 2 < z < 4
  sed.t = 0;
  sed.hist = [];
  record();
}
const mob = () => hinderedMobility({ M0: 1 / sed.visc, beta: sed.beta, lambda: 0.5 });

const playBtn = h("button", { class: "action", type: "button", onclick: () => (sed.anim ? stopSed() : startSed()) }, "時間を進める");
document.getElementById("sed-controls").append(
  field("熱ドリフト k_BT ∂M/∂z", segmented([[true, "入れる"], [false, "省く（誤り）"]], sed.drift, (v) => { sed.drift = v; resetSed(); renderSed(); })),
  slider({ label: "壁の近くでの移動度の低下 β", min: 0, max: 0.95, step: 0.05, value: sed.beta, format: (v) => v.toFixed(2), onInput: (v) => { sed.beta = v; resetSed(); renderSed(); } }),
  field("粘度", segmented([[1, "μ"], [10, "10μ"]], sed.visc, (v) => { sed.visc = v; resetSed(); renderSed(); })),
  playBtn,
  h("button", { class: "action", type: "button", onclick: () => { stopSed(); resetSed(); renderSed(); } }, "最初から"));
const histPlot = createPlot(document.getElementById("sed-hist"), { height: 300, xlabel: "壁からの高さ z / ℓ_g", ylabel: "密度 c(z)（全体で 1）" });
const timePlot = createPlot(document.getElementById("sed-time"), { height: 300, xlabel: "時間 t（ℓ_g² / (k_BT M₀)、M₀ は粘度 μ のとき）", ylabel: "平均の高さ ⟨z⟩ / ℓ_g" });

function record() {
  sed.hist.push([sed.t, sed.zs.reduce((s, z) => s + z, 0) / NP]);
}
function startSed() {
  playBtn.textContent = "停止";
  const step = () => {
    stepHeights(sed.zs, { mob: mob(), kBT, mg: 1, dt: DT, steps: STEPS_PER_FRAME, drift: sed.drift, randn: sed.randn });
    sed.t += DT * STEPS_PER_FRAME;
    record();
    renderSed();
    if (sed.anim) sed.anim = requestAnimationFrame(step);
  };
  sed.anim = requestAnimationFrame(step);
}
function stopSed() {
  if (sed.anim) cancelAnimationFrame(sed.anim);
  sed.anim = null;
  playBtn.textContent = "時間を進める";
}

function renderSed() {
  const m = mob();
  const w = ZMAX / BINS, counts = Array(BINS).fill(0);
  for (const z of sed.zs) if (z < ZMAX) counts[Math.floor(z / w)] += 1;
  const bars = [];
  counts.forEach((c, k) => { const y = c / (NP * w); bars.push([k * w, y], [(k + 1) * w, y]); });
  const grid = Array.from({ length: 201 }, (_, k) => (k * ZMAX) / 200);
  const good = stationaryDensity(grid, { mob: m, kBT, mg: 1, drift: true });
  const bad = stationaryDensity(grid, { mob: m, kBT, mg: 1, drift: false });
  histPlot.update({
    series: [
      { name: "粒子のヒストグラム", color: "var(--c1)", width: 2, points: bars },
      { name: "ボルツマン分布 exp(−z/ℓ_g)", color: "var(--c2)", width: 3, points: grid.map((z, k) => [z, good[k]]) },
      { name: "熱ドリフトを省いたときの定常分布 exp(−z/ℓ_g)/M(z)", color: "var(--ink)", dash: "5 4", points: grid.map((z, k) => [z, bad[k]]) },
      { name: "移動度 M(z)/M₀（参考）", color: "var(--c3)", dash: "2 3", points: grid.map((z) => [z, m.M(z) * sed.visc]) },
    ],
    xdomain: [0, ZMAX], ydomain: [0, Math.max(1.2, ...bad.slice(0, 5)) * 1.05],
  });
  // means over the whole normalisation range [0, 12], not just the plotted part
  const wide = Array.from({ length: 1201 }, (_, k) => k / 100);
  const mean = (p) => p.reduce((s, v, k) => s + (k === 0 || k === 1200 ? 0.5 : 1) * v * wide[k], 0) / 100;
  const meanGood = mean(stationaryDensity(wide, { mob: m, kBT, mg: 1, drift: true }));
  const meanBad = mean(stationaryDensity(wide, { mob: m, kBT, mg: 1, drift: false }));
  const tmax = Math.max(5, sed.t);
  timePlot.update({
    series: [
      { name: "粒子の平均の高さ", color: "var(--c1)", points: sed.hist },
      { name: "ボルツマン分布の平均", color: "var(--c2)", dash: "6 4", points: [[0, meanGood], [tmax, meanGood]] },
      { name: "熱ドリフトなしの定常分布の平均", color: "var(--ink)", dash: "2 3", points: [[0, meanBad], [tmax, meanBad]] },
    ],
    xdomain: [0, tmax], ydomain: [0, 3.2],
  });
  document.getElementById("sed-out").replaceChildren(h("p", { class: "caption" },
    `時間 ${fmt(sed.t, 3)}、平均の高さ ${fmt(sed.hist.at(-1)[1], 3)}。ボルツマン分布の平均は ${fmt(meanGood, 3)}（≈ ℓ_g）、熱ドリフトを省いたときの定常分布の平均は ${fmt(meanBad, 3)}。` +
    (sed.visc > 1 ? " 粘度を 10 倍にすると、平衡に近づく速さは 1/10 になるが、行き着く分布は変わらない。" : "")));
}

renderPair();
resetSed();
renderSed();
