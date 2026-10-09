import { test } from "node:test";
import assert from "node:assert/strict";
import { matmul, transpose, fnorm, sub } from "../src/core/linalg.js";
import { mulberry32, normalSampler } from "../src/core/random.js";
import { stf, stf5, unstf5, ddot, naive5, stfBasis } from "../src/physics/stf.js";
import { schurReduce, fullSolve, sdResistance, isolatedResistance11 } from "../src/models/sd.js";
import { axialManyBody, resistance, pairwiseSumResistance, axialComponents } from "../src/models/manyBody.js";
import { close, closeArray } from "./helpers.js";

const randn = normalSampler(mulberry32(14));

test("STF basis is orthonormal; S:E = s.e, but not the naive 5 numbers (regression 01)", () => {
  const B = stfBasis();
  B.forEach((Bi, i) => B.forEach((Bj, j) => close(ddot(Bi, Bj), i === j ? 1 : 0, { atol: 1e-15 })));
  const A = [[2, 1, 3], [1, -1, 4], [3, 4, -1]];
  const v = stf5(A);
  closeArray(unstf5(v), A, { atol: 1e-14 });
  close(v.reduce((s, x) => s + x * x, 0), A.flat().reduce((s, x) => s + x * x, 0), { rtol: 1e-14 });
  const S = stf([[0.3, -1, 2], [0.4, 1.2, 0.5], [-0.7, 0.1, -0.2]]), E = stf([[1, 0.2, 0], [0.6, -0.5, 0.3], [0.1, 0.9, 0.8]]);
  const s5 = stf5(S), e5 = stf5(E);
  close(s5.reduce((acc, x, k) => acc + x * e5[k], 0), ddot(S, E), { rtol: 1e-13 });
  const n5 = naive5(S).reduce((acc, x, k) => acc + x * naive5(E)[k], 0);
  assert.ok(Math.abs(n5 - ddot(S, E)) > 1e-3);
});

test("Schur complement reproduces the full solve; induced stresslet vanishes without coupling", () => {
  const n = 11, nv = 6;
  const Bm = Array.from({ length: n }, () => Array.from({ length: n }, randn));
  const M = matmul(Bm, transpose(Bm)).map((r, i) => r.map((v, j) => v + (i === j ? 0.5 : 0)));
  const gv = Array.from({ length: nv }, randn), einf = Array.from({ length: n - nv }, randn);
  const a = schurReduce(M, nv, { gv, einf }), b = fullSolve(M, nv, { gv, einf });
  closeArray(a.qv, b.qv, { rtol: 1e-10, atol: 1e-12 });
  closeArray(a.gs, b.gs, { rtol: 1e-10, atol: 1e-12 });
  assert.ok(Math.hypot(...a.gs) > 1e-3);
  const decoupled = M.map((r, i) => r.map((v, j) => ((i < nv) !== (j < nv) ? 0 : v)));
  const c = schurReduce(decoupled, nv, { gv, einf: new Array(n - nv).fill(0) });
  assert.ok(Math.hypot(...c.gs) < 1e-14);
  // background strain alone (no coupling): g_s = -Mss^-1 e_inf, from both routes
  const d1 = schurReduce(decoupled, nv, { gv, einf }), d2 = fullSolve(decoupled, nv, { gv, einf });
  closeArray(d1.gs, d2.gs, { rtol: 1e-10, atol: 1e-12 });
  assert.ok(Math.hypot(...d1.gs) > 1e-3);
  closeArray(isolatedResistance11().map((r, i) => r[i]), [...Array(3).fill(6 * Math.PI), ...Array(3).fill(8 * Math.PI), ...Array(5).fill(20 * Math.PI / 3)], { rtol: 1e-15 });
});

const zz = axialComponents;
function axialCase(z) {
  const ex = axialManyBody(z, z.map(() => 1), { L: 24 });
  const X = z.map((v) => [0, 0, v]);
  const Rr = zz(resistance(X)), R2r = zz(pairwiseSumResistance(X));
  return { ex, Rr, sd: sdResistance(Rr, ex.R2B, R2r) };
}

test("SD-type construction is exact for two spheres and beats RPY and the pairwise sum for three", () => {
  const two = axialCase([0, 3]);
  closeArray(two.sd, two.ex.R, { rtol: 1e-8 });
  const { ex, Rr, sd } = axialCase([0, 3, 6]);
  const err = (A) => fnorm(sub(A, ex.R)) / fnorm(ex.R);
  assert.ok(err(sd) < err(Rr) && err(sd) < err(ex.R2B), `${err(sd)} ${err(Rr)} ${err(ex.R2B)}`);
});
