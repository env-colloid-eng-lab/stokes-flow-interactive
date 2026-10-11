import { initPage, h, field, segmented, fmt, renderMath } from "../ui/page.js";
import { createScene } from "../ui/scene2d.js";
import { createPlot } from "../ui/plot.js";
import { mulberry32 } from "../core/random.js";
import { FLOWS, sphereMotion, velocity, moments, streamline } from "../models/linearFlow.js";

initPage("p15");

const a = 1, mu = 1, BOX = 9;
const W = 640, H = 460;
const st = { flow: "shear", mode: "free", part: "total", sign: 1, tracers: true, anim: null, phi: 0, pts: [] };

// ---- table of the four flows ---------------------------------------------------------------
const tex = (s) => `$${s}$`;
const FLOW_TEX = {
  uniform: ["Ge_x", "0", "0"],
  shear: ["0", "\\tfrac{G}{2}e_y", "\\tfrac{G}{2}(e_xe_z+e_ze_x)"],
  extension: ["0", "0", "G(e_xe_x-e_ze_z)"],
  rotation: ["0", "Ge_y", "0"],
};
const UINF_TEX = { uniform: "G\\,e_x", shear: "Gz\\,e_x", extension: "G(x\\,e_x-z\\,e_z)", rotation: "G(z\\,e_x-x\\,e_z)" };
document.getElementById("flow-table").append(h("div", { style: { overflowX: "auto" } }, h("table", { class: "data" },
  h("tr", {}, ...["流れ", "u^∞(x)", "U", "Ω^∞", "E"].map((t) => h("th", {}, t))),
  ...Object.entries(FLOWS).map(([k, f]) => h("tr", {}, h("td", {}, f.name), h("td", {}, tex(UINF_TEX[k])), ...FLOW_TEX[k].map((t) => h("td", {}, tex(t))))))));

// ---- controls ------------------------------------------------------------------------------
document.getElementById("lf-controls").append(
  field("流れ", segmented(Object.entries(FLOWS).map(([k, f]) => [k, f.name]), st.flow, (v) => { st.flow = v; update(); })),
  field("球", segmented([["free", "自由に動く"], ["fixed", "固定する（回転も止める）"]], st.mode, (v) => { st.mode = v; update(); })),
  field("表示", segmented([["total", "流れ全体"], ["disturbance", "乱れだけ"]], st.part, (v) => { st.part = v; update(); })));
const playBtn = h("button", { class: "action", type: "button", onclick: () => (st.anim ? stop() : start()) }, "止める");
document.getElementById("lf-anim").append(
  playBtn,
  h("button", { class: "action", type: "button", onclick: () => { st.sign = -st.sign; update(false); } }, "流れを逆にする"),
  field("トレーサー", segmented([[true, "表示"], [false, "隠す"]], st.tracers, (v) => { st.tracers = v; draw(); })));

const scene = createScene(document.getElementById("lf-scene"), { width: W, height: H });
const view = { cx: 0, cz: 0, span: (2 * BOX * W) / H }; // the box -BOX..BOX fills the height
const caption = document.getElementById("lf-caption");
const table = document.getElementById("lf-table");
const decayNote = h("p", { class: "caption" });
document.getElementById("lf-decay").after(decayNote);
const decay = createPlot(document.getElementById("lf-decay"), { width: 960, height: 280, xlog: true, ylog: true, xlabel: "距離 r/a", ylabel: "乱れの速さ（平均）/ Ga" });

// ---- state derived from the controls ---------------------------------------------------------
let flow, motion, lines = [];
const frame = () => (st.part === "total" ? motion.Up : [0, 0, 0]); // seen from the sphere's centre
const vel = (x, z) => {
  const u = velocity([x, 0, z], flow, motion, { a, part: st.part });
  const f = frame();
  return [u[0] - f[0], u[2] - f[2]];
};

function update(reseed = true) {
  flow = FLOWS[st.flow].make(st.sign);
  motion = sphereMotion(flow, st.mode);
  lines = streamlines();
  if (reseed || !st.pts.length) seedTracers();
  renderTable();
  renderDecay();
  caption.textContent = captionText();
  draw();
}

function captionText() {
  const still = Math.hypot(...vel(3, 2)) < 1e-12 && Math.hypot(...vel(-2, 5)) < 1e-12;
  const frameNote = st.part === "total" && Math.hypot(...motion.Up) > 0 ? "球と一緒に動く人から見た流れ（球の速度を引いた）。" : "";
  if (still) return `${frameNote}流れは止まって見える：この条件では球は流れを乱さない。`;
  return `${frameNote}細線は流線、点はトレーサー。球の上の線は球の回転を表す。`;
}

// streamlines: candidate seeds on the axes and the diagonals, nearest first, integrated both
// ways; a seed that lies close to a line already drawn is skipped, so the lines stay apart
function streamlines() {
  const radii = [1.04, 1.12, 1.3, 1.6, 2.1, 2.8, 3.8, 5, 6.5, 8.5];
  const angles = [90, 270, 0, 180, 45, 135, 225, 315, 22.5, 67.5, 112.5, 157.5, 202.5, 247.5, 292.5, 337.5];
  const seeds = radii.flatMap((r) => angles.map((d) => [r * Math.cos((d * Math.PI) / 180), r * Math.sin((d * Math.PI) / 180)]));
  const out = [], drawn = [];
  const opts = { a, box: BOX * 1.4, h: 0.06, frame: [frame()[0], 0, frame()[2]], part: st.part };
  const near = (p, d) => drawn.some((q) => Math.abs(q[0] - p[0]) < d && Math.abs(q[1] - p[1]) < d && Math.hypot(q[0] - p[0], q[1] - p[1]) < d);
  for (const p of seeds) {
    if (Math.hypot(...vel(p[0], p[1])) < 1e-9) continue;
    if (near(p, 0.18 + 0.12 * Math.hypot(...p))) continue; // looser spacing further out
    const fwd = streamline(p, flow, motion, { ...opts, dir: 1 });
    const closed = fwd.length > 2 && fwd[fwd.length - 1] === p;
    const pts = closed ? fwd : [...streamline(p, flow, motion, { ...opts, dir: -1 }).reverse(), ...fwd.slice(1)];
    out.push(pts);
    for (let i = 0; i < pts.length; i += 3) drawn.push(pts[i]);
    if (out.length >= 28) break;
  }
  return out;
}

// ---- tracers ---------------------------------------------------------------------------------
const rnd = mulberry32(5);
function randomPoint() {
  for (;;) {
    const p = [(2 * rnd() - 1) * view.span / 2, (2 * rnd() - 1) * BOX];
    if (Math.hypot(...p) > a * 1.02) return p;
  }
}
function seedTracers() {
  st.pts = Array.from({ length: 260 }, () => ({ p: randomPoint(), age: Math.floor(rnd() * 400) }));
}
function moveTracers(dt) {
  const half = view.span / 2;
  for (const t of st.pts) {
    const [x, z] = t.p, k1 = vel(x, z), k2 = vel(x + dt * k1[0], z + dt * k1[1]);
    t.p = [x + (dt / 2) * (k1[0] + k2[0]), z + (dt / 2) * (k1[1] + k2[1])];
    t.age++;
    // tracers that leave, touch the sphere or have lived long are put back at random
    if (Math.abs(t.p[0]) > half || Math.abs(t.p[1]) > BOX || Math.hypot(...t.p) < a || t.age > 900) { t.p = randomPoint(); t.age = 0; }
  }
}

function start() {
  playBtn.textContent = "止める";
  let last = performance.now();
  const step = (now) => {
    const dt = Math.min(0.05, (now - last) / 1000) * 1.2;
    last = now;
    moveTracers(dt);
    st.phi -= motion.Wp[1] * dt; // rotation about +y (into the screen) turns clockwise on the screen
    draw();
    st.anim = requestAnimationFrame(step);
  };
  st.anim = requestAnimationFrame(step);
}
function stop() {
  if (st.anim) cancelAnimationFrame(st.anim);
  st.anim = null;
  playBtn.textContent = "動かす";
}

// ---- drawing ---------------------------------------------------------------------------------
// the sphere is drawn as a polygon so that the rotation mark (a path) lies on top of it
const disk = Array.from({ length: 72 }, (_, k) => [a * Math.cos((k * Math.PI) / 36), a * Math.sin((k * Math.PI) / 36)]);
function draw() {
  const mark = [Math.cos(st.phi) * a, Math.sin(st.phi) * a];
  const paths = [
    ...lines.map((pts) => ({ points: pts, color: "var(--muted)", width: 1 })),
    { points: [[-mark[0], -mark[1]], mark], color: "var(--accent)", width: 2.5 },
  ];
  const spheres = st.tracers ? st.pts.map((t) => ({ x: t.p[0], z: t.p[1], a: 0.07, fill: "var(--c2)", stroke: "none" })) : [];
  scene.draw({ view, polys: [{ points: disk }], spheres, paths });
}

function renderTable() {
  const mo = moments(flow, motion, { a, mu });
  const six = 6 * Math.PI * mu * a, eight = 8 * Math.PI * mu * a ** 3, twenty = (20 / 3) * Math.PI * mu * a ** 3;
  const v2 = (v) => `(${fmt(v[0], 3)}, ${fmt(v[2], 3)})`;
  const rows = [
    ["球の速度 (U_x, U_z)", v2(motion.Up)],
    ["球の角速度 Ω_y", fmt(motion.Wp[1], 3)],
    ["力 (F_x, F_z) ÷ 6πμa", v2(mo.F.map((v) => v / six))],
    ["トルク L_y ÷ 8πμa³", fmt(mo.L[1] / eight, 3)],
    ["ストレスレット (S_xx, S_xz, S_zz) ÷ (20/3)πμa³", `(${fmt(mo.S[0][0] / twenty, 3)}, ${fmt(mo.S[0][2] / twenty, 3)}, ${fmt(mo.S[2][2] / twenty, 3)})`],
  ];
  table.replaceChildren(h("table", { class: "data" }, ...rows.map(([k, v]) => h("tr", {}, h("td", {}, k), h("td", {}, v)))));
}

function renderDecay() {
  const rs = Array.from({ length: 40 }, (_, i) => 1.05 * 10 ** ((i / 39) * Math.log10(60 / 1.05)));
  const n = 48;
  const mean = (r) => {
    let s = 0;
    for (let k = 0; k < n; k++) {
      const th = (2 * Math.PI * (k + 0.5)) / n;
      s += Math.hypot(...velocity([r * Math.cos(th), 0, r * Math.sin(th)], flow, motion, { a, part: "disturbance" }));
    }
    return s / n;
  };
  const data = rs.map((r) => [r, mean(r)]);
  // no disturbance at all (only rounding noise): nothing to draw
  if (Math.max(...data.map(([, v]) => v)) < 1e-9) { decay.update({ series: [] }); decayNote.textContent = "この条件では乱れがない。"; return; }
  decayNote.textContent = "";
  // reference slopes through the farthest point, drawn only where they stay near the data
  const [r0, v0] = data[data.length - 1], vmax = Math.max(...data.map(([, v]) => v));
  const ref = (p) => {
    const r1 = Math.max(1.05, r0 * (v0 / (2 * vmax)) ** (1 / p));
    return [[r1, v0 * (r0 / r1) ** p], [r0, v0]];
  };
  decay.update({
    series: [
      { name: "乱れの速さ", color: "var(--c1)", points: data },
      { name: "∝ 1/r", color: "var(--muted)", dash: "6 4", points: ref(1) },
      { name: "∝ 1/r²", color: "var(--c3)", dash: "2 3", points: ref(2) },
    ],
  });
}

update();
start();
renderMath();
