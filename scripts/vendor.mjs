// Copies pinned third-party files from node_modules into vendor/ so that the site works
// on GitHub Pages and offline without any CDN. Run after `npm install`: npm run vendor
import { cpSync, mkdirSync, rmSync, readFileSync, readdirSync } from "node:fs";

const files = {
  "three/build/three.module.min.js": "three/three.module.min.js",
  "three/build/three.core.min.js": "three/three.core.min.js",
  "three/examples/jsm/controls/OrbitControls.js": "three/addons/controls/OrbitControls.js",
  "three/LICENSE": "three/LICENSE",
  "katex/dist/katex.min.js": "katex/katex.min.js",
  "katex/dist/katex.min.css": "katex/katex.min.css",
  "katex/dist/contrib/auto-render.min.js": "katex/contrib/auto-render.min.js",
  "katex/dist/fonts": "katex/fonts",
  "katex/LICENSE": "katex/LICENSE",
};
rmSync("vendor", { recursive: true, force: true });
for (const [from, to] of Object.entries(files)) {
  mkdirSync(`vendor/${to}`.replace(/\/[^/]+$/, ""), { recursive: true });
  cpSync(`node_modules/${from}`, `vendor/${to}`, { recursive: true });
}
// KaTeX ships woff2, woff and ttf; every browser we target reads woff2.
for (const f of readdirSync("vendor/katex/fonts")) if (!f.endsWith(".woff2")) rmSync(`vendor/katex/fonts/${f}`);
const v = (p) => JSON.parse(readFileSync(`node_modules/${p}/package.json`, "utf8")).version;
console.log(`vendored three ${v("three")}, katex ${v("katex")}`);
