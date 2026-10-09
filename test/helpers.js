import { readFileSync } from "node:fs";
import assert from "node:assert/strict";

const root = new URL("../", import.meta.url);

// Reads reference/output/<name>.csv into { header, rows } with numeric rows.
export function readCsv(name) {
  const text = readFileSync(new URL(`reference/output/${name}.csv`, root), "utf8").trim();
  const [head, ...lines] = text.split(/\r?\n/);
  return { header: head.split(","), rows: lines.map((l) => l.split(",").map(Number)) };
}

export function close(actual, expected, { rtol = 1e-12, atol = 0 } = {}, msg = "") {
  const ok = Math.abs(actual - expected) <= atol + rtol * Math.abs(expected);
  assert.ok(ok, `${msg} expected ${expected}, got ${actual} (rtol ${rtol}, atol ${atol})`);
}

export function closeArray(actual, expected, tol = {}, msg = "") {
  const a = actual.flat(Infinity), e = expected.flat(Infinity);
  assert.equal(a.length, e.length, `${msg} length`);
  a.forEach((v, i) => close(v, e[i], tol, `${msg}[${i}]`));
}
