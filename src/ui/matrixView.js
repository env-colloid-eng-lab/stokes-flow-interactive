// Heat-map view of a (block) matrix with per-entry explanations.
// Each entry can be hovered (explanation) and clicked (selection, e.g. for a time series).
import { h, cssVar } from "./page.js";

function hexToRgb(hex) {
  const m = hex.replace("#", "");
  const n = parseInt(m.length === 3 ? m.split("").map((c) => c + c).join("") : m, 16);
  return [(n >> 16) & 255, (n >> 8) & 255, n & 255];
}
const mix = (a, b, t) => a.map((x, i) => Math.round(x + (b[i] - x) * t));

export function colorFor(v, vmax, gamma = 0.5) {
  const zero = hexToRgb(cssVar("--zero") || "#f3f1ec");
  const end = hexToRgb(cssVar(v >= 0 ? "--pos" : "--neg") || (v >= 0 ? "#c2410c" : "#1d4ed8"));
  const t = vmax > 0 ? Math.min(1, Math.abs(v) / vmax) ** gamma : 0;
  return { rgb: `rgb(${mix(zero, end, t).join(",")})`, t };
}

/**
 * opts.labels      row/column labels (e.g. ["1x","1y",...])
 * opts.block       block size for separators (3 for translation)
 * opts.explain     (i, j, value) => Node | string, shown under the matrix
 * opts.onSelect    (i, j) => void, called on click
 * opts.format      value formatter for the cell text
 */
export function createMatrixView(container, opts = {}) {
  const { block = 3, format = (v) => (Math.abs(v) < 5e-4 ? "0" : v.toFixed(Math.abs(v) >= 10 ? 1 : 3)) } = opts;
  const grid = h("div", { class: "matrix", role: "grid" });
  const legend = h("div", { class: "matrix-legend" });
  const explain = h("div", { class: "explain", "aria-live": "polite" }, "成分にカーソルを合わせると、その値を与える式と代入した数値を表示します。クリックで選択します。");
  container.append(h("div", { class: "matrix-wrap" }, grid), legend, explain);

  let M = null, cells = [], selected = opts.selected ?? null, hover = null, vmax = 1, labels = opts.labels ?? [];

  function showExplain(i, j) {
    explain.replaceChildren();
    if (!M || i == null) return;
    const content = opts.explain ? opts.explain(i, j, M[i][j]) : `${labels[i]}, ${labels[j]}: ${M[i][j]}`;
    explain.append(content instanceof Node ? content : document.createTextNode(content));
  }

  function build(n) {
    grid.replaceChildren();
    const small = n > 12;
    grid.style.gridTemplateColumns = `auto repeat(${n}, auto)`;
    grid.append(h("div", { class: "head" }));
    for (let j = 0; j < n; j++) {
      const lab = h("div", { class: "head" + ((j + 1) % block === 0 && j < n - 1 ? " cell block-r" : "") }, small ? "" : labels[j] ?? "");
      if (!small) lab.style.width = "46px";
      grid.append(lab);
    }
    cells = [];
    for (let i = 0; i < n; i++) {
      grid.append(h("div", { class: "head", style: { paddingRight: "4px" } }, small && i % block ? "" : labels[i] ?? ""));
      const row = [];
      for (let j = 0; j < n; j++) {
        let cls = "cell" + (small ? " small" : "");
        if ((j + 1) % block === 0 && j < n - 1) cls += " block-r";
        if ((i + 1) % block === 0 && i < n - 1) cls += " block-b";
        const c = h("div", { class: cls, role: "gridcell", tabindex: small ? -1 : 0 });
        c.addEventListener("mouseenter", () => { hover = [i, j]; showExplain(i, j); });
        c.addEventListener("mouseleave", () => { hover = null; if (selected) showExplain(...selected); });
        c.addEventListener("focus", () => showExplain(i, j));
        c.addEventListener("click", () => {
          selected = [i, j];
          paintSelection();
          showExplain(i, j);
          opts.onSelect?.(i, j);
        });
        row.push(c);
        grid.append(c);
      }
      cells.push(row);
    }
  }

  function paintSelection() {
    cells.forEach((row, i) => row.forEach((c, j) => c.classList.toggle("hl", !!selected && selected[0] === i && selected[1] === j)));
  }

  function update(matrix, { labels: newLabels, vmax: fixedMax } = {}) {
    if (newLabels) labels = newLabels;
    const n = matrix.length;
    if (!M || M.length !== n || newLabels) build(n);
    M = matrix;
    vmax = fixedMax ?? Math.max(...matrix.flat().map(Math.abs));
    for (let i = 0; i < n; i++)
      for (let j = 0; j < n; j++) {
        const v = matrix[i][j];
        const { rgb, t } = colorFor(v, vmax);
        const c = cells[i][j];
        c.style.background = rgb;
        c.style.color = t > 0.6 ? "#fff" : "";
        if (!c.classList.contains("small")) c.textContent = format(v);
        c.title = `${labels[i] ?? i}, ${labels[j] ?? j}: ${v.toPrecision(6)}`;
      }
    paintSelection();
    legend.replaceChildren(
      h("span", {}, `−${format(vmax)}`),
      h("span", { class: "bar", style: { background: `linear-gradient(90deg, ${colorFor(-vmax, vmax).rgb}, ${colorFor(0, vmax).rgb}, ${colorFor(vmax, vmax).rgb})` } }),
      h("span", {}, `+${format(vmax)}`),
      h("span", {}, "（色は |値|^0.5 で強調。小さな非対角成分も見えるようにしている）"));
    const target = hover ?? selected;
    if (target) showExplain(...target);
  }

  return {
    update,
    select(i, j) { selected = [i, j]; paintSelection(); showExplain(i, j); },
    get selected() { return selected; },
  };
}
