// Minimal SVG line plot with linear or logarithmic axes.
import { h } from "./page.js";

const NS = "http://www.w3.org/2000/svg";
const s = (tag, attrs = {}, text) => {
  const el = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) if (v != null) el.setAttribute(k, v);
  if (text != null) el.textContent = text;
  return el;
};

function niceTicks(lo, hi, n = 5) {
  const span = hi - lo || 1;
  const step0 = span / n;
  const mag = 10 ** Math.floor(Math.log10(step0));
  const step = [1, 2, 2.5, 5, 10].map((m) => m * mag).find((st) => span / st <= n + 0.5) ?? 10 * mag;
  const ticks = [];
  for (let t = Math.ceil(lo / step) * step; t <= hi + 1e-9 * span; t += step) ticks.push(+t.toPrecision(12));
  return ticks;
}
function logTicks(lo, hi) {
  const e0 = Math.ceil(Math.log10(lo) - 1e-9), e1 = Math.floor(Math.log10(hi) + 1e-9);
  const step = Math.max(1, Math.ceil((e1 - e0 + 1) / 10)); // at most about ten labelled decades
  const out = [];
  for (let e = Math.ceil(e0 / step) * step; e <= e1; e += step) out.push(10 ** e);
  return out;
}
const tickLabel = (t, log) => {
  if (log) { const e = Math.round(Math.log10(t)); return e === 0 ? "1" : e === 1 ? "10" : `1e${e}`; }
  return Math.abs(t) >= 1e4 || (Math.abs(t) < 1e-3 && t !== 0) ? t.toExponential(0) : String(+t.toPrecision(4));
};

export function createPlot(container, { width = 520, height = 300, xlabel = "", ylabel = "", xlog = false, ylog = false, title } = {}) {
  const svg = s("svg", { class: "plot", viewBox: `0 0 ${width} ${height}`, role: "img" });
  const legend = h("div", { class: "legend" });
  if (title) container.append(h("div", { class: "caption" }, title));
  container.append(svg, legend);
  const m = { l: 56, r: 12, t: 10, b: 40 };
  const W = width - m.l - m.r, H = height - m.t - m.b;
  let state = {};

  function update({ series = [], vlines = [], hlines = [], xdomain, ydomain, points = [] } = {}) {
    state = { series, vlines, hlines, xdomain, ydomain, points };
    svg.replaceChildren();
    const ok = (x, y) => Number.isFinite(x) && Number.isFinite(y) && (!xlog || x > 0) && (!ylog || y > 0);
    const all = series.flatMap((se) => se.points.filter(([x, y]) => ok(x, y))).concat(points.filter((p) => ok(p.x, p.y)).map((p) => [p.x, p.y]));
    if (!all.length) return;
    let [x0, x1] = xdomain ?? [Math.min(...all.map((p) => p[0])), Math.max(...all.map((p) => p[0]))];
    let [y0, y1] = ydomain ?? [Math.min(...all.map((p) => p[1])), Math.max(...all.map((p) => p[1]))];
    if (!ydomain && !ylog) { const pad = 0.06 * (y1 - y0 || Math.abs(y1) || 1); y0 -= pad; y1 += pad; }
    if (y0 === y1) { if (ylog) { y0 /= 2; y1 *= 2; } else { y0 -= 1; y1 += 1; } }
    if (x0 === x1) { if (xlog) { x0 /= 2; x1 *= 2; } else { const d = Math.abs(x0) || 1; x0 -= d; x1 += d; } }
    const fx = xlog ? (x) => (Math.log10(x) - Math.log10(x0)) / (Math.log10(x1) - Math.log10(x0)) : (x) => (x - x0) / (x1 - x0);
    const fy = ylog ? (y) => (Math.log10(y) - Math.log10(y0)) / (Math.log10(y1) - Math.log10(y0)) : (y) => (y - y0) / (y1 - y0);
    const X = (x) => m.l + fx(x) * W, Y = (y) => m.t + (1 - fy(y)) * H;

    const g = s("g", { class: "axis" });
    const xt = xlog ? logTicks(x0, x1) : niceTicks(x0, x1);
    const yt = ylog ? logTicks(y0, y1) : niceTicks(y0, y1);
    for (const t of xt) {
      g.append(s("line", { class: "grid-line", x1: X(t), x2: X(t), y1: m.t, y2: m.t + H }));
      g.append(s("text", { x: X(t), y: m.t + H + 15, "text-anchor": "middle" }, tickLabel(t, xlog)));
    }
    for (const t of yt) {
      g.append(s("line", { class: "grid-line", x1: m.l, x2: m.l + W, y1: Y(t), y2: Y(t) }));
      g.append(s("text", { x: m.l - 6, y: Y(t) + 4, "text-anchor": "end" }, tickLabel(t, ylog)));
    }
    g.append(s("rect", { x: m.l, y: m.t, width: W, height: H, fill: "none", stroke: "var(--line)" }));
    g.append(s("text", { x: m.l + W / 2, y: height - 6, "text-anchor": "middle" }, xlabel));
    g.append(s("text", { x: 12, y: m.t + H / 2, "text-anchor": "middle", transform: `rotate(-90 12 ${m.t + H / 2})` }, ylabel));
    svg.append(g);

    const clip = s("clipPath", { id: `clip${Math.random().toString(36).slice(2)}` });
    clip.append(s("rect", { x: m.l, y: m.t, width: W, height: H }));
    svg.append(clip);
    const body = s("g", { "clip-path": `url(#${clip.id})` });
    for (const v of vlines) body.append(s("line", { x1: X(v.x), x2: X(v.x), y1: m.t, y2: m.t + H, stroke: v.color ?? "var(--muted)", "stroke-dasharray": v.dash ?? "4 3" }));
    for (const v of hlines) body.append(s("line", { x1: m.l, x2: m.l + W, y1: Y(v.y), y2: Y(v.y), stroke: v.color ?? "var(--muted)", "stroke-dasharray": v.dash ?? "4 3" }));
    for (const se of series) {
      const pts = se.points.filter(([x, y]) => ok(x, y));
      if (!pts.length) continue;
      if (se.marker) {
        for (const [x, y] of pts) body.append(s("circle", { cx: X(x), cy: Y(y), r: se.r ?? 3.2, fill: se.color, stroke: "var(--panel)" }));
      } else {
        // break the polyline at NaN gaps
        let d = "", pen = false;
        for (const [x, y] of se.points) {
          if (!ok(x, y)) { pen = false; continue; }
          d += `${pen ? "L" : "M"}${X(x).toFixed(2)},${Y(y).toFixed(2)}`;
          pen = true;
        }
        body.append(s("path", { d, fill: "none", stroke: se.color, "stroke-width": se.width ?? 2, "stroke-dasharray": se.dash ?? null }));
      }
    }
    for (const p of points.filter((q) => ok(q.x, q.y))) {
      body.append(s("circle", { cx: X(p.x), cy: Y(p.y), r: p.r ?? 4.5, fill: p.color ?? "var(--ink)", stroke: "var(--panel)", "stroke-width": 1.5 }));
    }
    svg.append(body);
    legend.replaceChildren(...series.filter((se) => se.name).map((se) =>
      h("span", { class: se.dash ? "dash" : "", style: { "--sw": se.color } }, se.name)));
  }
  return { update, svg };
}
