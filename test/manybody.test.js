import { test } from "node:test";
import assert from "node:assert/strict";
import { fnorm, sub, transpose, inv, asymmetry } from "../src/core/linalg.js";
import { mobility } from "../src/physics/rpy.js";
import { joPolynomials, joResistance } from "../src/physics/jo.js";
import { resistance, pairwiseSumResistance, axialManyBody } from "../src/models/manyBody.js";
import { closeArray } from "./helpers.js";

test("for two particles the pairwise sum is the exact resistance", () => {
  const X = [[0, 0, 0], [2.7, 0.4, -1.1]];
  closeArray(pairwiseSumResistance(X), resistance(X), { rtol: 1e-12, atol: 1e-12 });
});

test("RPY resistance of three spheres is not pairwise additive, but stays symmetric", () => {
  const X = [[-3, 0, 0], [3, 0, 0], [0, 0, 4]];
  const R = resistance(X), R2B = pairwiseSumResistance(X);
  closeArray(R, inv(mobility(X)), { rtol: 1e-12, atol: 1e-12 });
  assert.ok(asymmetry(R) < 1e-13 && asymmetry(R2B) < 1e-13);
  const diff = fnorm(sub(R, R2B)) / fnorm(R);
  assert.ok(diff > 1e-4, `many-body part ${diff}`);
  // the many-body part decays as the third sphere is taken far away
  const far = [[-3, 0, 0], [3, 0, 0], [0, 0, 400]];
  const R12 = (M) => M.slice(0, 6).map((r) => r.slice(0, 6));
  assert.ok(fnorm(sub(R12(resistance(far)), R12(pairwiseSumResistance(far)))) < 1e-3 * fnorm(sub(R12(R), R12(R2B))));
});

test("coaxial three spheres: exact resistance is symmetric and differs from the exact pairwise sum", () => {
  const { R, R2B, boundaryError } = axialManyBody([0, 3, 6], [1, 1, 1], { L: 16 });
  assert.ok(boundaryError < 1e-6, `residual ${boundaryError}`);
  assert.ok(fnorm(sub(R, transpose(R))) / fnorm(R) < 1e-8);
  // pairwise sum built from two-body collocation agrees with JO for each isolated pair
  const polys = joPolynomials(60);
  const jo13 = joResistance(1, 1, 6, polys);
  assert.ok(Math.abs(R2B[0][2] - jo13[0][1]) / Math.abs(jo13[0][1]) < 1e-8);
  // shielding: the middle sphere changes the 1-3 coupling
  assert.ok(Math.abs(R[0][2] - R2B[0][2]) / Math.abs(R2B[0][2]) > 1e-2);
});
