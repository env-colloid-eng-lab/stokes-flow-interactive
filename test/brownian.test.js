// Mirrors reference/julia/regression_multibody.jl testsets 07 and 08, plus the wrong noise models.
import { test } from "node:test";
import assert from "node:assert/strict";
import { matmul, transpose, norm } from "../src/core/linalg.js";
import { mobility } from "../src/physics/rpy.js";
import { noiseFactor, sampleDisplacements, covariance, hinderedMobility, stepHeights, stationaryDensity, flux } from "../src/physics/brownian.js";
import { mulberry32, normalSampler } from "../src/core/random.js";
import { close } from "./helpers.js";

const fro = (A) => Math.sqrt(A.flat().reduce((s, v) => s + v * v, 0));
const diff = (A, B) => A.map((row, i) => row.map((v, j) => v - B[i][j]));
const M = mobility([[0, 0, 0], [4, 0, 0]]);

test("Cholesky noise has covariance 2 kBT M dt; the other two do not", () => {
  const kBT = 1, dt = 0.01;
  const { mean, C } = covariance(sampleDisplacements(M, { kBT, dt, n: 100000, seed: 3 }));
  const target = M.map((row) => row.map((v) => 2 * kBT * dt * v));
  assert.ok(fro(diff(C, target)) / fro(target) < 0.025);
  assert.ok(norm(mean) < 1e-3);
  const BBt = (m) => { const B = noiseFactor(M, m); return matmul(B, transpose(B)); };
  assert.ok(fro(diff(BBt("cholesky"), M)) / fro(M) < 1e-14);
  for (const m of ["independent", "elementwise"]) assert.ok(fro(diff(BBt(m), M)) / fro(M) > 0.05, m);
  // independent noise overestimates the relative diffusion along the line of centres
  const rel = (A) => A[0][0] + A[3][3] - 2 * A[0][3];
  assert.ok(rel(BBt("independent")) > 1.3 * rel(M));
});

test("zero flux: Boltzmann with the thermal drift, Boltzmann / M without it", () => {
  // the textbook check: M(z) = 1 + z, kBT = 1, m_b g = 2
  const mob = { M: (z) => 1 + z, dM: () => 1 };
  const z = 0.7, hh = 1e-6;
  const c = (x) => Math.exp(-2 * x), wrong = (x) => Math.exp(-2 * x) / mob.M(x);
  const d = (f) => (f(z + hh) - f(z - hh)) / (2 * hh);
  assert.ok(Math.abs(flux(c(z), d(c), z, { mob, mg: 2 })) < 1e-9);
  assert.ok(Math.abs(flux(wrong(z), d(wrong), z, { mob, mg: 2, drift: false })) < 1e-9);
  assert.ok(Math.abs(flux(c(z), d(c), z, { mob, mg: 2, drift: false })) > 0.1);
  // normalisation
  const mobH = hinderedMobility({ beta: 0.8 });
  const zs = Array.from({ length: 1201 }, (_, k) => k / 100);
  for (const drift of [true, false]) {
    const p = stationaryDensity(zs, { mob: mobH, drift });
    close(p.reduce((s, v, k) => s + (k === 0 || k === 1200 ? 0.5 : 1) * v, 0) / 100, 1, { rtol: 1e-3 }); // trapezoid, h = 0.01
  }
});

test("an ensemble relaxes to Boltzmann only with the thermal drift", () => {
  const mob = hinderedMobility({ beta: 0.8, lambda: 0.5 });
  const meanOf = (drift) => {
    const randn = normalSampler(mulberry32(11));
    const zs = Array.from({ length: 3000 }, () => 1);
    stepHeights(zs, { mob, drift, dt: 2e-3, steps: 3000, randn });
    return zs.reduce((s, z) => s + z, 0) / zs.length;
  };
  const grid = Array.from({ length: 4001 }, (_, k) => (k * 12) / 4000);
  const theory = (drift) => { const p = stationaryDensity(grid, { mob, drift }); return p.reduce((s, v, k) => s + v * grid[k], 0) * (12 / 4000); };
  close(theory(true), 1, { rtol: 1e-4 }); // mean height of e^{-z} is l_g = 1
  assert.ok(Math.abs(meanOf(true) - theory(true)) < 0.05);
  assert.ok(Math.abs(meanOf(false) - theory(false)) < 0.05);
  assert.ok(theory(true) - theory(false) > 0.15); // the two really differ
});
