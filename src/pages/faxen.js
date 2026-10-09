import { initPage, h, slider, fmt, frameThrottle } from "../ui/page.js";
import { createPlot } from "../ui/plot.js";
import { createScene } from "../ui/scene2d.js";
import { matvec } from "../core/linalg.js";
import { sphereFlow } from "../physics/oseen.js";
import { rpyPair } from "../physics/rpy.js";
import { surfaceAverage, surfaceMoments, strainField } from "../physics/sphere.js";

initPage("p7");

const mu = 1;

// ---------------------------------------------------------------------
// Sphere in a parabolic channel flow u_x = U0 (1 - z^2 / H^2)
// ---------------------------------------------------------------------
const ch = { aH: 0.3, zH: 0.4 };
const U0 = 1, H = 1;
const uChannel = (r) => [U0 * (1 - (r[2] / H) ** 2), 0, 0];
const schedule = frameThrottle(() => renderChannel());
document.getElementById("controls").append(
  slider({ label: "球の半径 a/H", min: 0.05, max: 0.6, step: 0.01, value: ch.aH, format: (v) => v.toFixed(2), onInput: (v) => { ch.aH = v; schedule(); } }),
  slider({ label: "中心の高さ z₀/H", min: -0.6, max: 0.6, step: 0.01, value: ch.zH, format: (v) => v.toFixed(2), onInput: (v) => { ch.zH = v; schedule(); } }));
const scene = createScene(document.getElementById("scene"), { width: 440, height: 300 });

function renderChannel() {
  const a = ch.aH * H, z0 = ch.zH * H;
  const zmax = H - a; // keep the sphere inside the channel picture
  const zc = Math.max(-zmax, Math.min(zmax, z0));
  const uc = uChannel([0, 0, zc])[0];
  const lap = (-2 * U0) / (H * H);
  const Ufax = uc + ((a * a) / 6) * lap; // free sphere: F^H = 0
  const avg = surfaceAverage(uChannel, [0, 0, zc], { a })[0];
  const omegaY = (-2 * U0 * zc) / (H * H) / 2; // (1/2)(curl u)_y = (1/2) du_x/dz
  const arrows = [];
  for (let z = -0.9; z <= 0.91; z += 0.15) arrows.push({ x: -1.6, z, vx: 0.8 * uChannel([0, 0, z])[0], vz: 0, color: "var(--line)", width: 1.3 });
  scene.draw({
    view: { cx: 0, cz: 0, span: 4 },
    lines: [{ x1: -3, z1: H, x2: 3, z2: H, color: "var(--ink)", dash: "0" }, { x1: -3, z1: -H, x2: 3, z2: -H, color: "var(--ink)", dash: "0" }],
    paths: [{ points: Array.from({ length: 41 }, (_, k) => { const z = -H + (2 * H * k) / 40; return [-1.6 + 0.8 * uChannel([0, 0, z])[0], z]; }), color: "var(--muted)" }],
    spheres: [{ x: 0.6, z: zc, a }],
    arrows: [...arrows, { x: 0.6, z: zc, vx: 0.8 * uc, vz: 0, color: "var(--muted)", width: 1.5, offset: 8 }, { x: 0.6, z: zc, vx: 0.8 * Ufax, vz: 0, color: "var(--c1)", width: 3 }],
    axes: { x: "x", up: "z", out: "y⊗" },
  });
  document.getElementById("faxen-out").replaceChildren(h("table", { class: "data" },
    h("tr", {}, h("td", {}, "中心での流速 u∞(z₀)"), h("td", {}, fmt(uc, 5))),
    h("tr", {}, h("td", {}, "Faxén の補正 (a²/6)∇²u∞"), h("td", {}, fmt(((a * a) / 6) * lap, 5))),
    h("tr", {}, h("td", {}, "球の速度 U = u∞ + (a²/6)∇²u∞"), h("td", {}, fmt(Ufax, 5))),
    h("tr", {}, h("td", {}, "球面上の u∞ の数値平均"), h("td", {}, fmt(avg, 5))),
    h("tr", {}, h("td", {}, "球の角速度 Ω_y = ½(∇×u∞)_y"), h("td", {}, fmt(omegaY, 4)))),
    h("p", { class: "caption" }, `球は中心の流体より ${fmt(uc - Ufax, 3)} だけ遅い（相対 ${(100 * (uc - Ufax) / uc).toFixed(1)}%）。` +
      "放物線の流れは双調和なので、Faxén の速度は球面平均と丸め誤差の範囲で一致する。" +
      (zc !== z0 ? " （球が壁にかからないよう、中心の位置を制限した。壁の影響はこの式に含まれない。）" : "")));
}

// ---------------------------------------------------------------------
// Receiver-averaged flow of a sender sphere = RPY
// ---------------------------------------------------------------------
const rp = { r: 3, deg: 30 };
const scheduleRpy = frameThrottle(() => renderRpy());
document.getElementById("rpy-controls").append(
  slider({ label: "中心間距離 r/a", min: 2.05, max: 8, step: 0.05, value: rp.r, format: (v) => v.toFixed(2), onInput: (v) => { rp.r = v; scheduleRpy(); } }),
  slider({ label: "力の向き（中心線から）", min: 0, max: 90, step: 1, value: rp.deg, format: (v) => `${v}°`, onInput: (v) => { rp.deg = v; scheduleRpy(); } }));
function renderRpy() {
  const a = 1, d = (rp.deg * Math.PI) / 180;
  const f = [Math.sin(d), 0, Math.cos(d)]; // sender at origin, receiver at r e_z
  const centre = [0, 0, rp.r];
  const U = (x) => sphereFlow(x, f.map((v) => v / (6 * Math.PI * mu * a)), { a, mu }); // flow of a sphere pushed by f
  const pointValue = U(centre);
  const avg = surfaceAverage(U, centre, { a, nt: 12, np: 24 });
  const rpy = matvec(rpyPair([0, 0, rp.r], { a, mu }), f);
  const row = (label, v, note) => h("tr", {}, h("td", {}, label), h("td", {}, `(${fmt(v[0], 6)}, ${fmt(v[2], 6)})`), h("td", { class: "caption" }, note));
  const err = Math.hypot(...avg.map((v, i) => v - rpy[i])) / Math.hypot(...rpy);
  document.getElementById("rpy-out").replaceChildren(h("table", { class: "data" },
    h("tr", {}, h("th", {}, ""), h("th", {}, "受け手の速度 (x, z) × f⁻¹"), h("th", {}, "")),
    row("送り手の流れを受け手の中心で評価", pointValue, "受け手の有限サイズを無視"),
    row("受け手の球面で平均（数値）", avg, "受け手の Faxén 則"),
    row("RPY の式 M₁₂ f", rpy, "ページ8で使う式")),
    h("p", { class: "caption" }, `球面平均と RPY の相対差 ${fmt(err, 2)}。中心での値との差は、受け手の大きさによる補正 (a²/6)∇² の分である。`));
}

// ---------------------------------------------------------------------
// Einstein viscosity from the stresslet surface integral
// ---------------------------------------------------------------------
const ein = { phi: 0.05 };
const Es = [[0, 0.5, 0], [0.5, 0, 0], [0, 0, 0]]; // simple shear rate 1: E_xy = 1/2
const Sxy = surfaceMoments((r) => strainField(r, Es, { a: 1, mu }), { a: 1, nt: 10, np: 20 }).S[0][1];
const einPlot = createPlot(document.getElementById("einstein-plot"), { height: 260, xlabel: "体積分率 φ", ylabel: "μ_eff / μ" });
const relVisc = (phi) => (mu + (phi / ((4 * Math.PI) / 3)) * Sxy) / mu;
document.getElementById("einstein-controls").append(
  slider({ label: "体積分率 φ", min: 0, max: 0.1, step: 0.002, value: ein.phi, format: (v) => v.toFixed(3), onInput: (v) => { ein.phi = v; renderEinstein(); } }));
function renderEinstein() {
  const phis = Array.from({ length: 51 }, (_, k) => (0.1 * k) / 50);
  einPlot.update({ series: [{ name: "μ_eff/μ = 1 + n S_xy / (2μE_xy)", color: "var(--c1)", points: phis.map((p) => [p, relVisc(p)]) }], points: [{ x: ein.phi, y: relVisc(ein.phi), color: "var(--c2)" }] });
  document.getElementById("einstein-out").replaceChildren(h("table", { class: "data" },
    h("tr", {}, h("td", {}, "ストレスレット S_xy（数値積分）"), h("td", {}, fmt(Sxy, 8))),
    h("tr", {}, h("td", {}, "厳密値 (20π/3) μ a³ E_xy"), h("td", {}, fmt((20 * Math.PI) / 3 * 0.5, 8))),
    h("tr", {}, h("td", {}, "μ_eff/μ（φ の値で）"), h("td", {}, fmt(relVisc(ein.phi), 6))),
    h("tr", {}, h("td", {}, "係数 (μ_eff/μ − 1)/φ"), h("td", {}, ein.phi > 0 ? fmt((relVisc(ein.phi) - 1) / ein.phi, 8) : "—"))),
    h("p", { class: "caption" }, "係数は φ によらず 5/2。粒子どうしの相互作用を入れると φ² 以上の項が加わる。"));
}

renderChannel();
renderRpy();
renderEinstein();
