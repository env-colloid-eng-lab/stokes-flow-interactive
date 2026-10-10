// Check that every local file the published site refers to is inside it.
// Usage: node scripts/check-site.mjs _site
// Scans src/href attributes of the HTML files and the static/dynamic imports of the JS files.
import { readdirSync, readFileSync, statSync, existsSync } from "node:fs";
import { join, dirname, resolve, relative, sep } from "node:path";

const root = resolve(process.argv[2] ?? "_site");
const walk = (d) => readdirSync(d).flatMap((f) => (statSync(join(d, f)).isDirectory() ? walk(join(d, f)) : [join(d, f)]));
const isLocal = (u) => !/^([a-z]+:|\/\/|#)/i.test(u);
const patterns = {
  ".html": /\b(?:src|href)="([^"]+)"/g,
  // our modules import each other by relative paths only
  ".js": /(?:\bfrom\s*|\bimport\s*\(\s*|\bimport\s+)["'](\.{1,2}\/[^"']+)["']/g,
};

const missing = [];
let checked = 0;
for (const file of walk(root)) {
  const re = patterns[file.slice(file.lastIndexOf("."))];
  // vendored libraries are copied as published; only our own files are checked
  if (!re || relative(root, file).startsWith("vendor")) continue;
  for (const [, ref] of readFileSync(file, "utf8").matchAll(re)) {
    if (!isLocal(ref)) continue;
    checked++;
    const target = resolve(dirname(file), ref.split(/[?#]/)[0]);
    if (!target.startsWith(root + sep) || !existsSync(target)) missing.push(`${relative(root, file)} -> ${ref}`);
  }
}
if (missing.length) {
  console.error(`missing from the site (${missing.length}):\n  ${missing.join("\n  ")}`);
  process.exit(1);
}
console.log(`site check: ${checked} local references, all present`);
