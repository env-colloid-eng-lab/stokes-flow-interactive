import { h, header, renderMath } from "../ui/page.js";
import { parts, pages, links, pageById } from "../topics.js";

document.body.prepend(header(null));
const map = document.getElementById("map");
for (const part of parts) {
  map.append(h("div", { class: "part-title" }, part.title));
  map.append(h("div", { class: "map" }, pages.filter((p) => p.part === part.id).map((p) =>
    h("a", { class: "topic", href: p.status === "ready" ? `pages/${p.file}` : null, "data-status": p.status },
      h("span", { class: "num" }, `${p.num}`, p.status === "ready" ? "" : h("span", { class: "src-tag" }, "準備中")),
      h("span", { class: "t" }, p.title),
      h("span", { class: "d" }, p.desc)))));
}

// Connections as a table grouped by the page that produces the result.
const graph = document.getElementById("graph");
const table = h("table", { class: "data" },
  h("tr", {}, h("th", {}, "結果を作るページ"), h("th", {}, "結果"), h("th", {}, "使うページ"), h("th", {}, "使い方")));
for (const l of links) {
  const from = pageById(l.from), to = pageById(l.to);
  table.append(h("tr", {}, h("td", {}, `${from.num}. ${from.title}`), h("td", {}, l.result), h("td", {}, `${to.num}. ${to.title}`), h("td", {}, l.use)));
}
graph.append(h("div", { style: { overflowX: "auto" } }, table));
renderMath();
