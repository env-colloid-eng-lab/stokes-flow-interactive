import { initPage, h, field, segmented, slider, fmt } from "../ui/page.js";
import { createScene } from "../ui/scene2d.js";
import { norm, vsub } from "../core/linalg.js";
import { mulberry32 } from "../core/random.js";
import { velocities, advance, minGap, repulsion, flowAt } from "../models/live.js";

initPage("p14");

const a = 1, mu = 1;
const F0 = 6 * Math.PI * mu * a; // unit force: an isolated sphere moves at speed 1 (a per time unit)
const SPRING = 1, FMAX = 8; // pull = F0 * SPRING * (pointer - centre) / a, at most FMAX * F0
const MAXN = { none: 16, rpy: 16, sd: 10 }; // SD solves a 6N system with N(N-1)/2 exact pairs

const PRESETS = {
  pair: { name: "二球", make: () => [[-2, 0, 0], [2, 0, 0]] },
  tri: { name: "三角", make: () => [[-2.5, 0, -1.5], [2.5, 0, -1.5], [0, 0, 2.8]] },
  row: { name: "一列", make: () => [-8, -4, 0, 4, 8].map((x) => [x, 0, 0]) },
  cloud: { name: "雲", make: () => cloud(8, 4.2, 2.3, 11) },
  empty: { name: "空", make: () => [] },
};
function cloud(n, R, dmin, seed) {
  const rnd = mulberry32(seed), pts = [];
  while (pts.length < n) {
    const p = [(2 * rnd() - 1) * R, 0, (2 * rnd() - 1) * R];
    if (Math.hypot(p[0], p[2]) <= R && pts.every((q) => norm(vsub(p, q)) >= dmin)) pts.push(p);
  }
  return pts;
}

const st = {
  preset: "pair", model: "rpy", click: "add", gravity: 0, speed: 1, field: true, trails: true, follow: false,
  X: [], trails0: [], drag: null, anim: null, U: null, F: null, note: "",
  view: { cx: 0, cz: 0, span: 26 },
};
// a narrow screen shows less of the plane, so that the spheres stay large enough to grab
const baseSpan = () => (scene.svg.clientWidth && scene.svg.clientWidth < 520 ? 20 : 26);
const homeView = () => ({ cx: 0, cz: 0, span: baseSpan() });

const scene = createScene(document.getElementById("live-scene"), { width: 640, height: 460 });
scene.svg.style.touchAction = "none"; // let the pointer drag spheres instead of scrolling the page
const note = document.getElementById("live-note");
const table = document.getElementById("live-table");

const modelSeg = segmented([["none", "なし"], ["rpy", "RPY（案A）"], ["sd", "SD（案B）"]], st.model, (v) => setModel(v));
document.getElementById("live-controls").append(
  field("配置", segmented(Object.entries(PRESETS).map(([k, p]) => [k, p.name]), st.preset, (v) => { st.preset = v; reset(); })),
  field("相互作用", modelSeg),
  field("クリック", segmented([["add", "球を追加"], ["remove", "削除"]], st.click, (v) => { st.click = v; })),
  slider({ label: "重力 F/F₀", min: 0, max: 1, step: 0.05, value: st.gravity, format: (v) => v.toFixed(2), onInput: (v) => { st.gravity = v; wake(); } }),
  slider({ label: "再生速度", min: 0.25, max: 4, step: 0.25, value: st.speed, format: (v) => `×${v}`, onInput: (v) => { st.speed = v; } }));
document.getElementById("live-view").append(
  field("流れ", segmented([[true, "表示"], [false, "隠す"]], st.field, (v) => { st.field = v; draw(); })),
  field("軌跡", segmented([[true, "表示"], [false, "隠す"]], st.trails, (v) => { st.trails = v; draw(); })),
  field("視点", segmented([[false, "固定"], [true, "重心を追う"]], st.follow, (v) => { st.follow = v; if (!v) st.view = homeView(); draw(); })),
  h("button", { class: "action", type: "button", onclick: reset }, "配置を戻す"));

function setModel(m) {
  if (st.X.length > MAXN[m]) {
    modelSeg.set(st.model);
    note.textContent = `SD（案B）で扱える球は ${MAXN.sd} 個まで。球を減らしてから切り替える。`;
    return;
  }
  st.model = m;
  wake();
}

function reset() {
  st.X = PRESETS[st.preset].make();
  st.trails0 = st.X.map((p) => [[p[0], p[2]]]);
  st.drag = null;
  st.view = homeView();
  st.note = "";
  if (st.X.length > MAXN[st.model]) { st.model = "rpy"; modelSeg.set("rpy"); }
  wake();
}

// ---- forces --------------------------------------------------------------------------------
function externalForces(X) {
  return X.map((p, i) => {
    const f = [0, 0, -st.gravity * F0];
    if (st.drag && st.drag.index === i) {
      let dx = (st.drag.target[0] - p[0]) * SPRING * F0 / a, dz = (st.drag.target[1] - p[2]) * SPRING * F0 / a;
      const m = Math.hypot(dx, dz), cap = FMAX * F0;
      if (m > cap) { dx *= cap / m; dz *= cap / m; }
      f[0] += dx; f[2] += dz;
    }
    return f;
  });
}
const totalForces = (X) => {
  const Fe = externalForces(X), Fr = repulsion(X, { a, F0 });
  return Fe.map((f, p) => f.map((v, k) => v + Fr[p][k]));
};

// ---- pointer -------------------------------------------------------------------------------
const hit = ([x, z]) => st.X.findIndex((p) => Math.hypot(p[0] - x, p[2] - z) <= a);
scene.svg.addEventListener("pointerdown", (e) => {
  const w = scene.toWorld(e, st.view), i = hit(w);
  if (i >= 0 && st.click === "remove") {
    st.X.splice(i, 1); st.trails0.splice(i, 1);
    wake();
    return;
  }
  if (i >= 0) {
    st.drag = { index: i, target: w, id: e.pointerId };
    scene.svg.setPointerCapture(e.pointerId);
    wake();
    return;
  }
  if (st.click !== "add") return;
  const p = [w[0], 0, w[1]];
  if (st.X.length >= MAXN[st.model]) {
    note.textContent = `この相互作用で置ける球は ${MAXN[st.model]} 個まで。`;
    return;
  }
  if (st.X.some((q) => norm(vsub(p, q)) < 2 * a + 0.2)) {
    note.textContent = "ほかの球と重なる位置には置けない。";
    return;
  }
  st.X.push(p); st.trails0.push([[p[0], p[2]]]);
  wake();
});
scene.svg.addEventListener("pointermove", (e) => {
  if (!st.drag || e.pointerId !== st.drag.id) return;
  st.drag.target = scene.toWorld(e, st.view);
});
const release = (e) => { if (st.drag && e.pointerId === st.drag.id) { st.drag = null; wake(); } };
scene.svg.addEventListener("pointerup", release);
scene.svg.addEventListener("pointercancel", release);

// ---- time loop -----------------------------------------------------------------------------
// Runs while something can move (a drag, gravity, or spheres still pushed apart by repulsion).
function wake() {
  if (!st.anim) {
    let last = performance.now();
    const step = (now) => {
      const wall = Math.min(0.05, (now - last) / 1000);
      last = now;
      const moving = tick(wall * st.speed);
      st.anim = moving ? requestAnimationFrame(step) : null;
      draw();
    };
    st.anim = requestAnimationFrame(step);
  }
}

function tick(dt) {
  if (!st.X.length) { st.U = null; return false; }
  try {
    if (dt > 0) {
      const s = advance(st.X, externalForces, st.model, dt, { a, mu, F0, maxSub: st.model === "sd" ? 3 : 12 });
      st.X = s.X;
      st.note = s.t < 0.5 * dt ? "球が近いので、時間をゆっくり進めている。" : "";
      st.X.forEach((p, i) => { st.trails0[i].push([p[0], p[2]]); if (st.trails0[i].length > 600) st.trails0[i].shift(); });
    }
    st.F = totalForces(st.X);
    st.U = velocities(st.X, st.F, st.model, { a, mu });
  } catch (err) {
    st.note = `計算できなかった（${err.message}）。`;
    st.drag = null;
    return false;
  }
  if (st.follow) {
    const n = st.X.length;
    st.view = { ...st.view, cx: st.X.reduce((s, p) => s + p[0], 0) / n, cz: st.X.reduce((s, p) => s + p[2], 0) / n };
  }
  const vmax = Math.max(...st.U.map((u) => norm(u)));
  return !!st.drag || st.gravity > 0 || vmax > 1e-3;
}

// ---- drawing -------------------------------------------------------------------------------
function draw() {
  const { X, U, F, view } = st;
  const spheres = X.map((p, i) => ({
    x: p[0], z: p[2], a, label: `${i + 1}`,
    fill: st.drag?.index === i ? "var(--warn-soft)" : undefined,
    stroke: st.drag?.index === i ? "var(--warn)" : undefined,
  }));
  const arrows = [], lines = [];
  const vmax = U ? Math.max(1e-9, ...U.map((u) => norm(u))) : 1;
  const vscale = 2.5 / Math.max(vmax, 0.5); // the fastest sphere gets an arrow of 2.5a (or less)
  if (st.field && U && X.length) {
    const sp = 2, half = view.span / 2, halfZ = (view.span * 460) / 640 / 2;
    for (let x = Math.ceil((view.cx - half) / sp) * sp; x <= view.cx + half; x += sp)
      for (let z = Math.ceil((view.cz - halfZ) / sp) * sp; z <= view.cz + halfZ; z += sp) {
        const u = flowAt([x, 0, z], X, F, { a, mu });
        if (!u) continue;
        let vx = u[0] * vscale, vz = u[2] * vscale;
        const m = Math.hypot(vx, vz);
        if (m > 0.9 * sp) { vx *= (0.9 * sp) / m; vz *= (0.9 * sp) / m; }
        arrows.push({ x, z, vx, vz, color: "var(--muted)", width: 1.2, opacity: 0.7 });
      }
  }
  if (U) X.forEach((p, i) => arrows.push({ x: p[0], z: p[2], vx: U[i][0] * vscale, vz: U[i][2] * vscale, color: "var(--c2)", width: 2.5 }));
  if (st.drag) {
    const p = X[st.drag.index];
    if (p) lines.push({ x1: p[0], z1: p[2], x2: st.drag.target[0], z2: st.drag.target[1], color: "var(--warn)", dash: "4 3" });
  }
  const paths = st.trails ? st.trails0.map((t) => ({ points: t, color: "var(--c1)", width: 1, dash: "2 3" })) : [];
  scene.draw({ view, spheres, arrows, paths, lines });
  note.textContent = st.note || (X.length ? "" : "空いた所をクリックして球を置く。");
  renderTable();
}

let tableFrame = 0;
function renderTable() {
  if (tableFrame++ % 6 && st.anim) return; // while running, the numbers need not change every frame
  const { X, U } = st;
  if (!U || !X.length) { table.replaceChildren(); return; }
  const ref = st.drag ? st.drag.index : null;
  const uref = ref != null ? norm(U[ref]) : null;
  const rows = X.map((_, i) => h("tr", {},
    h("td", {}, `${i + 1}${i === ref ? "（引っ張り中）" : ""}`),
    h("td", {}, fmt(norm(U[i]), 3)),
    h("td", {}, uref > 1e-9 ? fmt(norm(U[i]) / uref, 3) : "—")));
  const gap = minGap(X, a);
  table.replaceChildren(
    h("table", { class: "data" }, h("tr", {}, h("th", {}, "球"), h("th", {}, "速さ"), h("th", {}, "比")), ...rows),
    h("p", { class: "caption" }, `最小の隙間 h/a = ${Number.isFinite(gap) ? fmt(gap, 3) : "—"}`));
}

addEventListener("resize", () => { st.view = { ...st.view, span: baseSpan() }; draw(); });
reset();
