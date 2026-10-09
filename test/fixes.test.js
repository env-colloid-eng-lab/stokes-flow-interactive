// Regression tests for the issues found in the first review of PR #1.
import { test } from "node:test";
import assert from "node:assert/strict";
import { cond, singularValues, choleskySolve, cholesky, matmul, transpose, matvec } from "../src/core/linalg.js";
import { Rational } from "../src/core/rational.js";
import { joPolynomials } from "../src/physics/jo.js";
import { resistanceCoefficients } from "../src/models/pairA.js";
import { approachTrajectory } from "../src/models/approach.js";
import { close, closeArray } from "./helpers.js";

const polys60 = joPolynomials(60);
const unit = 6 * Math.PI;

test("JO convergence check works for every offered order, including K = 10", () => {
  for (const K of [10, 20, 40, 60]) {
    const far = resistanceCoefficients(12, { axial: "jo", polys: polys60.slice(0, K + 1) }).jo;
    assert.ok(far.converged, `K=${K} far field: change ${far.change}`);
    assert.ok(far.Kc < far.K && far.Kc >= far.K - Math.max(2, Math.floor(K / 4)));
    assert.equal(resistanceCoefficients(2.01, { axial: "jo", polys: polys60.slice(0, K + 1) }).jo.converged, false);
  }
});

function rel(model, K = 40) {
  return (h) => {
    const c = resistanceCoefficients(2 + h, model === "lub" ? { lubrication: true, hc: 0.2 }
      : model === "jo" ? { axial: "jo", polys: polys60.slice(0, K + 1) } : {});
    return { value: (c.R.par.self - c.R.par.cross) / unit, ok: model !== "jo" || c.jo.converged };
  };
}

test("lubricated approach never reports contact and decays at the rate K = 2F/(3 pi mu a^2)", () => {
  const tr = approachTrajectory(rel("lub"), { h0: 2, Tmax: 120 });
  assert.equal(tr.stop, null);
  close(tr.ts.at(-1), 120, { rtol: 1e-12 });
  const hEnd = tr.hs.at(-1);
  assert.ok(hEnd > 0 && hEnd < 1e-4, `h(120) = ${hEnd}`);
  const n = tr.ts.length, i = n - 200;
  const rate = -Math.log(tr.hs[n - 1] / tr.hs[i]) / (tr.ts[n - 1] - tr.ts[i]);
  close(rate, 2 / (3 * Math.PI), { rtol: 1e-3 });
});

test("RPY alone reaches contact; the JO series stops where it stops converging", () => {
  const rpy = approachTrajectory(rel("rpy"), { h0: 2, Tmax: 120, contactBelow: 1e-4 });
  assert.equal(rpy.stop, "contact");
  assert.ok(rpy.ts.at(-1) > 20 && rpy.ts.at(-1) < 60);
  const jo = approachTrajectory(rel("jo"), { h0: 2, Tmax: 120 });
  assert.equal(jo.stop, "series");
  const hLast = jo.hs.at(-1);
  assert.ok(hLast > 0.1 && hLast < 2);
  // every accepted state is inside the converged range
  for (const hh of jo.hs) assert.ok(rel("jo")(hh).ok);
});

test("Rational.toNumber keeps precision when one part is small", () => {
  close(new Rational(2n ** 1010n, 3n).toNumber(), 2 ** 1010 / 3, { rtol: 1e-15 });
  close(new Rational(2n ** 1005n + 1n, 2049n).toNumber(), 2 ** 1005 / 2049, { rtol: 1e-15 });
  close(new Rational(7n, 2n ** 1100n).toNumber(), 7 * 2 ** -1100, { rtol: 1e-15 });
  close(new Rational(-(3n ** 700n), 5n ** 400n).toNumber(), -Math.exp(700 * Math.log(3) - 400 * Math.log(5)), { rtol: 1e-12 });
});

test("cond uses singular values and survives an ill-conditioned matrix", () => {
  const c = Math.cos(0.4), s = Math.sin(0.4);
  const Q = [[c, -s, 0], [s, c, 0], [0, 0, 1]];
  const D = [[1, 0, 0], [0, 1e-10, 0], [0, 0, 0.5]];
  const A = matmul(matmul(Q, D), transpose(Q));
  close(cond(A), 1e10, { rtol: 1e-5 });
  closeArray(singularValues([[3, 0], [0, 4], [0, 0]]), [4, 3], { rtol: 1e-15 });
});

test("choleskySolve solves with the Cholesky factor", () => {
  const A = [[4, 1, 0.5], [1, 3, 0.2], [0.5, 0.2, 2]];
  const x = [1, -2, 0.5];
  closeArray(choleskySolve(cholesky(A), matvec(A, x)), x, { rtol: 1e-13 });
});
