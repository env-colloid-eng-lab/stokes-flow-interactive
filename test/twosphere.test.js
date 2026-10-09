// Mirrors the JO / collocation parts of reference/julia/regression_foundation.jl
// and the values quoted in chapters 12-13.
import { test } from "node:test";
import assert from "node:assert/strict";
import { fnorm, sub, transpose, eigvalsSym } from "../src/core/linalg.js";
import { Q } from "../src/core/rational.js";
import { joPolynomials, joXA, joResistance, pairMotion, evalPoly } from "../src/physics/jo.js";
import { axialCollocation } from "../src/physics/collocation.js";
import { reflectionCoefficients } from "../src/physics/bridge.js";
import { readCsv, close, closeArray } from "./helpers.js";

const exact = joPolynomials(10, { exact: true });
const polys = joPolynomials(60);

test("JO coefficients f0-f5 and f10 are exact rationals", () => {
  const targets = [[1], [0, 3], [0, 9, 0], [0, -4, 27, -4], [0, -24, 81, 36, 0], [0, 0, 72, 243, 72, 0]];
  targets.forEach((t, k) => assert.deepEqual(exact[k].map(String), t.map(String), `f${k}`));
  assert.deepEqual(exact[10].map(String), [0, 0, 2304, 20736, 42804, 115849, 76176, 39264, 20736, 2304, 0].map(String));
  // float generation agrees with the exact one
  exact.forEach((row, k) => closeArray(polys[k], row.map((c) => c.toNumber()), { rtol: 1e-13, atol: 1e-9 }));
});

test("low-order reflections reproduce f0-f4; dropping the stresslet changes only f4", () => {
  for (const lam of [Q(1), Q(1, 2), Q(3)]) {
    const f = reflectionCoefficients(lam);
    f.forEach((fk, k) => {
      const expected = exact[k].reduce((acc, c, q) => acc.add(c.mul(new (c.constructor)(lam.num ** BigInt(q), lam.den ** BigInt(q)))), Q(0));
      assert.ok(fk.equals(expected), `lambda=${lam} f${k}: ${fk} vs ${expected}`);
    });
    const g = reflectionCoefficients(lam, { stresslet: false });
    for (let k = 0; k < 4; k++) assert.ok(g[k].equals(f[k]));
    assert.ok(!g[4].equals(f[4]));
  }
});

test("X^A at s = 3, lambda = 1 (section 13.6) and jo_convergence.csv", () => {
  const { x11, x12 } = joXA(3, 1, polys);
  close(x11, 1.3684796, { atol: 5e-8 });
  close(x12, -0.670175, { atol: 5e-8 });
  for (const [K, X11, X12] of readCsv("jo_convergence").rows) {
    const r = joXA(3, 1, joPolynomials(K));
    close(r.x11, X11, { rtol: 1e-13 }, `K=${K}`);
    close(r.x12, X12, { rtol: 1e-13 }, `K=${K}`);
  }
});

test("JO series agrees with independent boundary collocation", () => {
  for (const [a1, a2, r] of [[1, 1, 4], [1, 1, 3], [1, 0.5, 3]]) {
    const jo = joResistance(a1, a2, r, polys);
    const c = axialCollocation([0, r], [a1, a2], { L: 24, ncheck: 151, withCondition: false });
    assert.ok(fnorm(sub(jo, c.R)) / fnorm(c.R) < 2e-9, `${a1},${a2},${r}`);
    assert.ok(c.boundaryError < 2e-9);
    assert.ok(fnorm(sub(jo, transpose(jo))) / fnorm(jo) < 2e-13);
    assert.ok(eigvalsSym(jo)[0] > 0);
  }
  const one = axialCollocation([0], [1.3], { L: 4 });
  close(one.R[0][0], 6 * Math.PI * 1.3, { rtol: 1e-12 });
});

test("collocation convergence table matches axial_convergence.csv", () => {
  for (const [gap, L, X11, X12, resid, condition] of readCsv("axial_convergence").rows) {
    const c = axialCollocation([0, 2 + gap], [1, 1], { L, ncheck: 151 });
    // Least-squares solutions agree closely; the residual and condition are diagnostics.
    close(c.R[0][0] / (6 * Math.PI), X11, { rtol: 1e-9 }, `gap=${gap} L=${L} X11`);
    close(c.R[0][1] / (6 * Math.PI), X12, { rtol: 1e-9 }, `gap=${gap} L=${L} X12`);
    close(c.boundaryError, resid, { rtol: 1e-3, atol: 1e-13 }, `gap=${gap} L=${L} residual`);
    close(c.condition, condition, { rtol: 1e-6 }, `gap=${gap} L=${L} cond`);
  }
});

test("pair approach under opposite forces matches jo_pair_motion.csv", () => {
  const rows = pairMotion([0, 4], [1, -1], [1, 1], polys, { dt: 0.1, steps: 100, minGap: 1 });
  const ref = readCsv("jo_pair_motion").rows;
  assert.equal(rows.length, ref.length);
  rows.forEach((row, i) => closeArray(row, ref[i], { rtol: 1e-11, atol: 1e-14 }, `row ${i}`));
});

test("evalPoly handles exact coefficients", () => {
  close(evalPoly(exact[3], 1), 19, { rtol: 0 });
});
