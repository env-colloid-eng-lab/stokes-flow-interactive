// Common page chrome: header navigation, topic-link boxes, math rendering.
import { pages, links, pageById } from "../topics.js";

export function h(tag, attrs = {}, ...children) {
  const el = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === "class") el.className = v;
    else if (k === "style" && typeof v === "object")
      for (const [sk, sv] of Object.entries(v)) sk.startsWith("--") ? el.style.setProperty(sk, sv) : (el.style[sk] = sv);
    else if (k.startsWith("on")) el.addEventListener(k.slice(2), v);
    else if (v !== false && v != null) el.setAttribute(k, v === true ? "" : v);
  }
  for (const c of children.flat()) if (c != null) el.append(c instanceof Node ? c : document.createTextNode(c));
  return el;
}

const katex = () => window.katex;

export function tex(src, { display = false } = {}) {
  const span = document.createElement(display ? "div" : "span");
  if (katex()) katex().render(src, span, { displayMode: display, throwOnError: false });
  else span.textContent = src;
  return span;
}

export function renderMath(root = document.body) {
  if (window.renderMathInElement)
    window.renderMathInElement(root, {
      delimiters: [
        { left: "$$", right: "$$", display: true },
        { left: "\\(", right: "\\)", display: false },
        { left: "$", right: "$", display: false },
      ],
      throwOnError: false,
    });
}

const base = () => (location.pathname.includes("/pages/") ? "../" : "./");

export function header(currentId) {
  const ready = pages.filter((p) => p.status === "ready");
  return h("header", { class: "site-header" },
    h("a", { class: "brand", href: `${base()}index.html` }, "テンソルからストークス動力学へ"),
    h("nav", {}, ready.map((p) =>
      h("a", { href: `${base()}pages/${p.file}`, "aria-current": p.id === currentId ? "page" : false }, `${p.num}. ${p.title}`))));
}

function pageLink(p) {
  return p.status === "ready"
    ? h("a", { href: `${base()}pages/${p.file}` }, `${p.num}. ${p.title}`)
    : h("span", {}, `${p.num}. ${p.title}`, h("span", { class: "src-tag" }, "準備中"));
}

export function topicLinks(id) {
  const uses = links.filter((l) => l.to === id);
  const usedBy = links.filter((l) => l.from === id);
  const list = (items, key) => h("ul", {}, items.map((l) =>
    h("li", {}, pageLink(pageById(l[key])), "：", key === "from" ? l.result : l.use)));
  return h("div", { class: "links-box" },
    uses.length ? h("section", {}, h("h4", {}, "このページで使う結果"), list(uses, "from")) : null,
    usedBy.length ? h("section", {}, h("h4", {}, "この結果を使う先"), list(usedBy, "to")) : null);
}

// Call from each page: inserts the header and the link boxes (into #links), then renders math.
export function initPage(id) {
  document.body.prepend(header(id));
  const slot = document.getElementById("links");
  if (slot) slot.replaceWith(topicLinks(id));
  const p = pageById(id);
  const secEl = document.getElementById("sec");
  if (secEl && p) secEl.textContent = `教科書の対応箇所：${p.sec}`;
  renderMath();
}

export function cssVar(name) {
  return getComputedStyle(document.documentElement).getPropertyValue(name).trim();
}

export function segmented(options, value, onChange) {
  const wrap = h("span", { class: "seg", role: "group" });
  const buttons = options.map(([v, label]) => {
    const b = h("button", { type: "button", "aria-pressed": String(v === value) }, label);
    b.addEventListener("click", () => {
      buttons.forEach((x) => x.setAttribute("aria-pressed", "false"));
      b.setAttribute("aria-pressed", "true");
      onChange(v);
    });
    return b;
  });
  wrap.append(...buttons);
  // reflect a value changed elsewhere; a value not among the options presses none
  wrap.set = (v) => buttons.forEach((b, i) => b.setAttribute("aria-pressed", String(options[i][0] === v)));
  return wrap;
}

export function slider({ label, min, max, step, value, format = (v) => v, onInput }) {
  const out = h("output", {}, format(value));
  const input = h("input", { type: "range", min, max, step, value });
  input.addEventListener("input", () => { out.textContent = format(+input.value); onInput(+input.value); });
  const el = h("label", {}, label, input, out);
  el.set = (v) => { input.value = v; out.textContent = format(v); };
  return el;
}

export const fmt = (v, d = 4) => {
  if (!Number.isFinite(v)) return String(v);
  if (v === 0) return "0";
  const a = Math.abs(v);
  return a >= 1e4 || a < 1e-3 ? v.toExponential(d - 1) : v.toPrecision(d);
};
