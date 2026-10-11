import { initPage, h, field, segmented, slider, fmt } from "../ui/page.js";
import { createScene } from "../ui/scene2d.js";
import { createPlot } from "../ui/plot.js";
import { FLOWS } from "../models/linearFlow.js";
import { bgTable, batchelorGreen, relativeTrajectory } from "../models/pairLinear.js";

initPage("p16");

const W = 640, H = 460, BOX = 8;
const st = { flow: "shear", model: "exact", ratio: 1, y0: 0, pick: null, anim: null, k: 0 };
const tables = new Map(); // A, B per model and ratio (the exact one takes a moment to build)
const tableFor = () => {
  const key = `${st.model}:${st.ratio}`;
  if (!tables.has(key)) tables.set(key, bgTable(st.model, { a1: 1, a2: st.ratio }));
  return tables.get(key);
};

const MODEL_NAMES = [["none", "相互作用なし"], ["far", "遠方近似"], ["exact", "厳密（案B）"]];
document.getElementById("pl-controls").append(
  field("流れ", segmented([["shear", "単純せん断"], ["extension", "平面伸長"]], st.flow, (v) => { st.flow = v; st.pick = null; update(); })),
  field("相互作用", segmented(MODEL_NAMES, st.model, (v) => { st.model = v; update(); })),
  field("半径比 a₂/a₁", segmented([[1, "1"], [0.5, "1/2"], [0.25, "1/4"]], st.ratio, (v) => { st.ratio = v; update(); })),
  slider({ label: "渦度方向のずれ y₀", min: 0, max: 2, step: 0.1, value: st.y0, format: (v) => v.toFixed(1), onInput: (v) => { st.y0 = v; schedule(); } }));
const playBtn = h("button", { class: "action", type: "button", onclick: () => (st.anim ? stop() : start()) }, "動かす");
// a closed orbit starts within ~1e-5 a of contact in the plane of shear: too close to click
const closedBtn = h("button", { class: "action", type: "button", onclick: () => pickClosed() }, "閉じた軌道を選ぶ");
document.getElementById("pl-anim").append(playBtn, closedBtn);
function pickClosed() {
  const c = family.filter((t) => t.closed).at(-1);
  if (!c) { caption.textContent = "この条件では閉じた軌道が見つからない（せん断・厳密で試す）。"; return; }
  st.pick = { p0: c.pts[0], traj: c };
  st.k = 0;
  renderInfo();
  draw();
}

const scene = createScene(document.getElementById("pl-scene"), { width: W, height: H });
const view = { cx: 0, cz: 0, span: (2 * BOX * W) / H };
const caption = document.getElementById("pl-caption"), info = document.getElementById("pl-info");
const abPlot = createPlot(document.getElementById("pl-ab"), { width: 960, height: 300, xlog: true, ylog: true, xlabel: "隙間 r/a − 2（等しい球）", ylabel: "A、B" });

let family = [], pending = false;
const schedule = () => { if (!pending) { pending = true; requestAnimationFrame(() => { pending = false; update(); }); } };

// seeds: upstream of the reference sphere (and next to it, where closed orbits may be)
function seeds() {
  const zs = [0.05, 0.25, 0.5, 0.8, 1.2, 1.7, 2.3, 3, 4, 5.5];
  if (st.flow === "shear") return [...zs.map((z) => [-BOX, st.y0, z]), ...zs.map((z) => [BOX, st.y0, -z]), ...[1e-6, 1e-5, 3e-5, 1e-3, 1e-2, 3e-2].map(topSeed).filter(Boolean)];
  return [...zs.map((x) => [x, st.y0, BOX]), ...zs.map((x) => [-x, st.y0, -BOX])];
}
// above sphere 1 (x = 0, y = y0) at the centre distance (a1 + a2)(1 + g)
function topSeed(g) {
  const R = (1 + st.ratio) * (1 + g);
  return R > st.y0 ? [0, st.y0, Math.sqrt(R * R - st.y0 * st.y0)] : null;
}

function trajectory(p0, tab, flow) {
  const fwd = relativeTrajectory(p0, flow, tab, { box: BOX * 1.05 });
  if (fwd.closed) return fwd;
  const bwd = relativeTrajectory(p0, flow, tab, { box: BOX * 1.05, dir: -1 });
  return { pts: [...bwd.pts.reverse(), ...fwd.pts.slice(1)], rmin: Math.min(fwd.rmin, bwd.rmin), closed: false };
}

function update() {
  const flow = FLOWS[st.flow].make(1), tab = tableFor();
  const t0 = performance.now();
  family = seeds().map((p) => trajectory(p, tab, flow));
  if (st.pick) st.pick.traj = trajectory(st.pick.p0, tab, flow);
  st.k = 0;
  caption.textContent = `${MODEL_NAMES.find(([k]) => k === st.model)[1]}：${family.length} 本の軌道（${Math.round(performance.now() - t0)} ms）。` +
    (st.model === "exact" ? "厳密な関数では、球どうしは重ならない（y₀ = 0 では、どの軌道も破線の円の外に留まる）。" : "この近似では、球どうしが重なる軌道がある（y₀ = 0 では破線の円の中に入る）。");
  renderInfo();
  draw();
  renderAB();
}

scene.svg.addEventListener("pointerdown", (e) => {
  const [x, z] = scene.toWorld(e, view), a12 = 1 + st.ratio;
  if (Math.hypot(x, st.y0, z) <= a12 && st.model === "exact") return; // inside the excluded sphere
  st.pick = { p0: [x, st.y0, z] };
  st.pick.traj = trajectory(st.pick.p0, tableFor(), FLOWS[st.flow].make(1));
  st.k = 0;
  renderInfo();
  draw();
});

function renderInfo() {
  const t = st.pick?.traj;
  if (!t) { info.replaceChildren(h("p", { class: "caption" }, "図をクリックして軌道を選ぶ。")); return; }
  const gap = t.rmin - (1 + st.ratio);
  info.replaceChildren(h("table", { class: "data" },
    h("tr", {}, h("td", {}, "出発点 (x, y, z)"), h("td", {}, st.pick.p0.map((v) => v.toFixed(2)).join(", "))),
    h("tr", {}, h("td", {}, "種類"), h("td", {}, t.closed ? "閉じた軌道（回り続ける）" : "開いた軌道（通り過ぎる）")),
    h("tr", {}, h("td", {}, "最も近いときの隙間 |r| − (a₁+a₂)"), h("td", {}, gap > 0 ? fmt(gap, 3) : `${fmt(gap, 3)}（重なり）`))));
}

function draw() {
  const paths = family.map((t) => ({ points: t.pts.map((p) => [p[0], p[2]]), color: t.closed ? "var(--c3)" : "var(--muted)", width: 1 }));
  const a2 = st.ratio, spheres = [{ x: 0, z: 0, a: 1, label: "1" }];
  if (st.pick?.traj) {
    paths.push({ points: st.pick.traj.pts.map((p) => [p[0], p[2]]), color: "var(--c2)", width: 2.5 });
    const p = st.pick.traj.pts[Math.min(st.k, st.pick.traj.pts.length - 1)];
    spheres.push({ x: p[0], z: p[2], a: a2, label: "2", fill: "var(--warn-soft)", stroke: "var(--warn)" });
  }
  const ring = Array.from({ length: 73 }, (_, i) => [(1 + a2) * Math.cos((i * Math.PI) / 36), (1 + a2) * Math.sin((i * Math.PI) / 36)]);
  paths.push({ points: ring, color: "var(--accent)", dash: "4 4", width: 1 });
  scene.draw({ view, spheres, paths });
}

function start() {
  if (!st.pick?.traj) return;
  playBtn.textContent = "止める";
  const step = () => {
    st.k = (st.k + 3) % st.pick.traj.pts.length;
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

let abDone = false;
function renderAB() {
  if (abDone) return; // equal spheres only, the same for every flow
  abDone = true;
  const xs = Array.from({ length: 48 }, (_, i) => 10 ** (-5 + (i * 7) / 47));
  const ex = xs.map((x) => batchelorGreen(2 + x)), fa = xs.map((x) => batchelorGreen(2 + x, "far"));
  const r = (x) => 2 + x;
  abPlot.update({
    series: [
      { name: "A（厳密）", color: "var(--c1)", points: xs.map((x, i) => [x, ex[i].A]) },
      { name: "B（厳密）", color: "var(--c2)", points: xs.map((x, i) => [x, ex[i].B]) },
      { name: "A（遠方近似）", color: "var(--c1)", dash: "6 4", points: xs.map((x, i) => [x, fa[i].A]) },
      { name: "B（遠方近似）", color: "var(--c2)", dash: "6 4", points: xs.map((x, i) => [x, fa[i].B]) },
      { name: "遠方の展開", color: "var(--muted)", dash: "2 3", points: xs.filter((x) => x > 1).map((x) => [x, 5 / r(x) ** 3 - 8 / r(x) ** 5]) },
      { name: "", color: "var(--muted)", dash: "2 3", points: xs.filter((x) => x > 1).map((x) => [x, 16 / (3 * r(x) ** 5)]) },
    ],
  });
}

update();
