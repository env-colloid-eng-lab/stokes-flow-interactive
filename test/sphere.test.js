// Mirrors the surface-integral parts of reference/julia/regression_foundation.jl.
import { test } from "node:test";
import assert from "node:assert/strict";
import { matvec } from "../src/core/linalg.js";
import {
  sphereRule, translationField, rotationField, strainField, surfaceMoments, stf,
  faxenStresslet, surfaceAverage, faxenForce,
} from "../src/physics/sphere.js";
import { close, closeArray } from "./helpers.js";

const a = 1.3, mu = 0.7;
const U = [0.4, -0.2, 0.7], Om = [0.2, 0.8, -0.3];
const E = [[0.2, 0.4, -0.1], [0.4, -0.3, 0.2], [-0.1, 0.2, 0.1]];

test("quadrature weights integrate the solid angle", () => {
  close(sphereRule().weights.reduce((s, w) => s + w, 0), 4 * Math.PI, { rtol: 1e-14 });
});

test("surface traction gives 6 pi mu a U, 8 pi mu a^3 Omega and 20 pi/3 mu a^3 E", () => {
  const mt = surfaceMoments((r) => translationField(r, U, { a, mu }), { a, mu });
  const mr = surfaceMoments((r) => rotationField(r, Om, { a, mu }), { a, mu });
  const ms = surfaceMoments((r) => strainField(r, E, { a, mu }), { a, mu });
  closeArray(mt.F, U.map((x) => -6 * Math.PI * mu * a * x), { rtol: 1e-12 });
  closeArray(mr.T, Om.map((x) => -8 * Math.PI * mu * a ** 3 * x), { rtol: 1e-12 });
  closeArray(ms.S, E.map((row) => row.map((x) => ((20 * Math.PI) / 3) * mu * a ** 3 * x)), { rtol: 1e-12, atol: 1e-14 });
  assert.ok(Math.hypot(...ms.F) < 1e-12 && Math.hypot(...ms.T) < 1e-12);
  // pressure drag is one third, viscous drag two thirds
  closeArray(mt.Fpressure, U.map((x) => -2 * Math.PI * mu * a * x), { rtol: 1e-12 });
  closeArray(mt.Fviscous, U.map((x) => -4 * Math.PI * mu * a * x), { rtol: 1e-12 });
});

test("boundary conditions on the surface; the translating sphere has uniform traction", () => {
  const { normals } = sphereRule();
  for (const n of normals.slice(0, 40)) {
    const r = n.map((x) => a * x);
    closeArray(translationField(r, U, { a, mu }).u, U, { rtol: 1e-12, atol: 1e-14 });
    assert.ok(Math.hypot(...strainField(r, E, { a, mu }).u) < 1e-12);
    closeArray(matvec(strainField(r, E, { a, mu }).sigma, n), matvec(E, n).map((x) => 5 * mu * x), { rtol: 1e-10, atol: 1e-12 });
    closeArray(matvec(translationField(r, U, { a, mu }).sigma, n), U.map((x) => ((-3 * mu) / (2 * a)) * x), { rtol: 1e-10, atol: 1e-12 });
  }
});

test("Faxen curvature term and the Einstein coefficient 5/2 from a surface integral", () => {
  const A = 1.2, M = 0.8;
  // surface mean of (a y)^2 equals a^2/3
  closeArray(surfaceAverage((r) => [r[1] ** 2, 0, 0], [0, 0, 0], { a: A }), [A * A / 3, 0, 0], { rtol: 1e-13, atol: 1e-15 });
  // a cubic flow: stresslet from the traction 5 mu ... equals the Faxen stresslet with lap E
  const { normals, weights } = sphereRule(10, 20);
  const cubic = (r) => {
    const r2 = r[0] ** 2 + r[1] ** 2 + r[2] ** 2;
    return [5 * r2 * 2 * r[0], 5 * r2 * -2 * r[1], 0].map((v, i) => (v - 4 * r[i] * (r[0] ** 2 - r[1] ** 2)) / (42 * M));
  };
  let S = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  normals.forEach((n, k) => {
    const v = cubic(n.map((x) => A * x));
    const st = stf(v.map((vi) => n.map((nj) => vi * nj)));
    S = S.map((row, i) => row.map((s, j) => s + 5 * M * A * A * weights[k] * st[i][j]));
  });
  const lapE = [[2 / M, 0, 0], [0, -2 / M, 0], [0, 0, 0]];
  closeArray(S, faxenStresslet([[0, 0, 0], [0, 0, 0], [0, 0, 0]], lapE, { a: A, mu: M }), { rtol: 1e-12, atol: 1e-14 });
  // Einstein: relative viscosity 1 + 5/2 phi
  const Es = [[0, 0.5, 0], [0.5, 0, 0], [0, 0, 0]];
  const ms = surfaceMoments((r) => strainField(r, Es, { a: A, mu: M }), { a: A, mu: M, nt: 10, np: 20 });
  const phi = 0.02, num = phi / ((4 * Math.PI * A ** 3) / 3);
  close(((M + num * ms.S[0][1]) / M - 1) / phi, 2.5, { rtol: 1e-12 });
  // a free sphere (F = 0) in u = u0 + quadratic moves with u0 + a^2/6 lap u
  closeArray(faxenForce([1 + A * A / 6 * -2, 0, 0], [1, 0, 0], [-2, 0, 0], { a: A, mu: M }), [0, 0, 0], { atol: 1e-14 });
});
