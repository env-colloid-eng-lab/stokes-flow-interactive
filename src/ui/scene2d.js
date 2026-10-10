// 2D (x, z) view of spheres with force / velocity arrows. z points up on screen.
const NS = "http://www.w3.org/2000/svg";
const s = (tag, attrs = {}, text) => {
  const el = document.createElementNS(NS, tag);
  for (const [k, v] of Object.entries(attrs)) if (v != null) el.setAttribute(k, v);
  if (text != null) el.textContent = text;
  return el;
};

export function createScene(container, { width = 480, height = 340 } = {}) {
  const svg = s("svg", { class: "scene", viewBox: `0 0 ${width} ${height}`, role: "img" });
  container.append(svg);
  const defs = s("defs");
  svg.append(defs);
  const markers = new Map();
  function marker(color) {
    if (markers.has(color)) return markers.get(color);
    const id = `ah${markers.size}${Math.random().toString(36).slice(2, 6)}`;
    const m = s("marker", { id, viewBox: "0 0 10 10", refX: 8, refY: 5, markerWidth: 7, markerHeight: 7, orient: "auto-start-reverse" });
    m.append(s("path", { d: "M0,0 L10,5 L0,10 z", fill: color }));
    defs.append(m);
    markers.set(color, id);
    return id;
  }
  const layer = s("g");
  svg.append(layer);

  /**
   * view: { cx, cz, span }  world box centred at (cx, cz), `span` world units across the width
   * spheres: [{ x, z, a, label, fill }]
   * arrows:  [{ x, z, vx, vz, color, label, width }]   (vectors already in world units)
   * paths:   [{ points: [[x, z], ...], color, dash }]
   * axes:    true to draw the x / z axis triad, or { x, up, out } labels for another plane
   * polys:   [{ points: [[x, z], ...], fill, stroke }] closed polygons
   */
  function draw({ view, spheres = [], arrows = [], paths = [], polys = [], axes = true, lines = [], texts = [] }) {
    layer.replaceChildren();
    const k = width / view.span;
    const X = (x) => width / 2 + (x - view.cx) * k;
    const Z = (z) => height / 2 - (z - view.cz) * k;
    for (const l of lines)
      layer.append(s("line", { x1: X(l.x1), y1: Z(l.z1), x2: X(l.x2), y2: Z(l.z2), stroke: l.color ?? "var(--line)", "stroke-dasharray": l.dash ?? "5 4" }));
    for (const pg of polys) {
      const d = pg.points.map(([x, z], i) => `${i ? "L" : "M"}${X(x).toFixed(1)},${Z(z).toFixed(1)}`).join("") + "Z";
      layer.append(s("path", { d, fill: pg.fill ?? "var(--accent-soft)", stroke: pg.stroke ?? "var(--accent)", "stroke-width": pg.width ?? 1.5, "fill-opacity": pg.opacity ?? 1 }));
    }
    for (const p of paths) {
      const d = p.points.map(([x, z], i) => `${i ? "L" : "M"}${X(x).toFixed(1)},${Z(z).toFixed(1)}`).join("");
      layer.append(s("path", { d, fill: "none", stroke: p.color ?? "var(--muted)", "stroke-width": p.width ?? 1.5, "stroke-dasharray": p.dash ?? null }));
    }
    for (const sp of spheres) {
      layer.append(s("circle", { cx: X(sp.x), cy: Z(sp.z), r: sp.a * k, fill: sp.fill ?? "var(--accent-soft)", stroke: sp.stroke ?? "var(--accent)", "stroke-width": sp.strokeWidth ?? 1.5 }));
      if (sp.label) layer.append(s("text", { x: X(sp.x), y: Z(sp.z) + 4, "text-anchor": "middle", style: "font-size:13px;font-weight:700;fill:var(--ink)" }, sp.label));
    }
    for (const a of arrows) {
      const len = Math.hypot(a.vx, a.vz) * k;
      if (len < 2) continue;
      const color = a.color ?? "var(--ink)";
      // optional sideways offset (screen px) so that parallel arrows stay visible
      const off = a.offset ?? 0;
      const ox = (-a.vz * k / len) * off, oy = (-a.vx * k / len) * off;
      layer.append(s("line", {
        x1: X(a.x) + ox, y1: Z(a.z) + oy, x2: X(a.x + a.vx) + ox, y2: Z(a.z + a.vz) + oy,
        stroke: color, "stroke-width": a.width ?? 2.5, "marker-end": `url(#${marker(color)})`, "stroke-opacity": a.opacity ?? null,
      }));
      if (a.label) layer.append(s("text", { x: X(a.x + a.vx) + 6, y: Z(a.z + a.vz) - 4, style: `fill:${color};font-size:12px` }, a.label));
    }
    for (const t of texts) layer.append(s("text", { x: X(t.x), y: Z(t.z), "text-anchor": t.anchor ?? "middle", style: t.style ?? "" }, t.text));
    if (axes) {
      const ox = 26, oz = height - 22;
      layer.append(s("line", { x1: ox, y1: oz, x2: ox + 30, y2: oz, stroke: "var(--muted)", "marker-end": `url(#${marker("var(--muted)")})` }));
      layer.append(s("line", { x1: ox, y1: oz, x2: ox, y2: oz - 30, stroke: "var(--muted)", "marker-end": `url(#${marker("var(--muted)")})` }));
      // default: x right, z up, y into the screen; pass { x, up, out } for other planes
      const lab = typeof axes === "object" ? axes : { x: "x", up: "z", out: "y⊗" };
      layer.append(s("text", { x: ox + 36, y: oz + 4 }, lab.x));
      layer.append(s("text", { x: ox - 4, y: oz - 34 }, lab.up));
      if (lab.out) layer.append(s("text", { x: ox - 14, y: oz + 14 }, lab.out));
    }
  }
  // world coordinates (x, z) of a pointer event, for the view last passed to draw()
  function toWorld(evt, view) {
    const pt = svg.createSVGPoint();
    pt.x = evt.clientX; pt.y = evt.clientY;
    const p = pt.matrixTransform(svg.getScreenCTM().inverse());
    const k = width / view.span;
    return [view.cx + (p.x - width / 2) / k, view.cz - (p.y - height / 2) / k];
  }
  return { draw, svg, toWorld };
}
