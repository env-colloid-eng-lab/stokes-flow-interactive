// Pages 4, 5 and 7: dissipation, the point force, and the Faxen route to RPY.
import { test } from "node:test";
import assert from "node:assert/strict";
import { dot, matvec } from "../src/core/linalg.js";
import { oseen, lapOseen, sphereFlow, fourierStokeslet, projector } from "../src/physics/oseen.js";
import { rpyPair } from "../src/physics/rpy.js";
import { translationField, rotationField, surfaceMoments, surfaceAverage, dissipationDensity, dissipationIntegral, gaussLegendre } from "../src/physics/sphere.js";
import { shearReversal, rmsDistance } from "../src/physics/kinematics.js";
import { close, closeArray } from "./helpers.js";

test("Gauss-Legendre integrates polynomials of degree 2n-1", () => {
  const { nodes, weights } = gaussLegendre(5);
  close(nodes.reduce((s, x, k) => s + weights[k] * x ** 8, 0), 2 / 9, { rtol: 1e-13 });
});

test("dissipation: the spin does not dissipate, and the fluid dissipates exactly F . U", () => {
  close(dissipationDensity([[0, 1, 0], [-1, 0, 0], [0, 0, 0]], 2), 0, { atol: 1e-15 });
  close(dissipationDensity([[0, 1, 0], [0, 0, 0], [0, 0, 0]], 2), 2, { rtol: 1e-15 }); // 2 mu (2 * (1/2)^2)
  const a = 1.3, mu = 0.7, U = [0.3, -0.5, 0.8], Om = [0.4, 0.1, -0.6];
  const U2 = dot(U, U), O2 = dot(Om, Om);
  close(dissipationIntegral((r) => translationField(r, U, { a, mu }), { a }), 6 * Math.PI * mu * a * U2, { rtol: 1e-11 });
  close(dissipationIntegral((r) => rotationField(r, Om, { a, mu }), { a }), 8 * Math.PI * mu * a ** 3 * O2, { rtol: 1e-11 });
  // a finite shell a < r < R holds the fraction 1 - 3s/2 + s^3 - s^5/2 (translation) and
  // 1 - s^3 (rotation) of it, s = a/R: the translating sphere's dissipation reaches far out
  for (const R of [1.5 * a, 4 * a, 10 * a, 300 * a]) {
    const s = a / R;
    close(dissipationIntegral((r) => translationField(r, U, { a, mu }), { a, R, nr: 3, nt: 3, np: 6 }) / (6 * Math.PI * mu * a * U2), 1 - 1.5 * s + s ** 3 - s ** 5 / 2, { rtol: 1e-12 });
    close(dissipationIntegral((r) => rotationField(r, Om, { a, mu }), { a, R }) / (8 * Math.PI * mu * a ** 3 * O2), 1 - s ** 3, { rtol: 1e-12 });
  }
});

test("Oseen tensor: twice as fast along the force as across it", () => {
  const mu = 0.6, r = 2.5, f = [0, 0, 1];
  closeArray(matvec(oseen([0, 0, r], { mu }), f), [0, 0, 2 / (8 * Math.PI * mu * r)], { rtol: 1e-14, atol: 1e-16 });
  closeArray(matvec(oseen([r, 0, 0], { mu }), f), [0, 0, 1 / (8 * Math.PI * mu * r)], { rtol: 1e-14, atol: 1e-16 });
});

test("Fourier space: the velocity is the transverse part of the force", () => {
  const k = [0.3, -1.2, 0.5], f = [1, 0.4, -0.7], mu = 1.7;
  const { u, pImag } = fourierStokeslet(k, f, { mu });
  close(dot(k, u), 0, { atol: 1e-15 });
  const kk = dot(k, k);
  closeArray(u, matvec(projector(k), f).map((v) => v / (mu * kk)), { rtol: 1e-14 });
  // momentum balance: mu k^2 u + i k p = f, with i k p = -k pImag carrying the part along k
  closeArray(u.map((v, i) => mu * kk * v - k[i] * pImag), f, { rtol: 1e-14 });
});

test("sphere surface: the point force alone averages to U, the dipole makes it uniform", () => {
  const a = 1.4, mu = 0.8, U = [0.2, -0.1, 0.9];
  const f = U.map((v) => 6 * Math.PI * mu * a * v);
  const Gf = (r) => matvec(oseen(r, { mu }), f);
  closeArray(surfaceAverage(Gf, [0, 0, 0], { a }), U, { rtol: 1e-13 });
  for (const n of [[0, 0, 1], [1, 0, 0], [0.6, 0, 0.8]]) {
    const r = n.map((v) => a * v);
    closeArray(sphereFlow(r, U, { a, mu }), U, { rtol: 1e-13 });
    const lap = matvec(lapOseen(r, { mu }), f).map((v) => ((a * a) / 6) * v);
    closeArray(Gf(r).map((v, i) => v + lap[i]), U, { rtol: 1e-13 });
  }
});

test("the singularity construction and the closed-form field agree; force split uses the field's mu", () => {
  const a = 1.1, mu = 0.6, U = [0.3, 0.2, -0.9];
  for (const r of [[1.1, 0, 0], [0.5, -1.4, 2], [0, 0, 7]]) closeArray(sphereFlow(r, U, { a, mu }), translationField(r, U, { a, mu }).u, { rtol: 1e-13, atol: 1e-16 });
  const m = surfaceMoments((r) => translationField(r, U, { a, mu }), { a }); // no mu passed here
  closeArray(m.Fpressure.map((v, i) => v + m.Fviscous[i]), m.F, { rtol: 1e-13 });
  closeArray(m.Fviscous, U.map((x) => -4 * Math.PI * mu * a * x), { rtol: 1e-12 });
});

test("receiver's surface average of the sender's flow is the RPY mobility", () => {
  const a = 1, mu = 1, f = [0.6, -0.3, 0.75];
  for (const R of [[0, 0, 2.2], [1.5, 2, -1], [0, 4, 0]]) {
    const avg = surfaceAverage((x) => sphereFlow(x, f.map((v) => v / (6 * Math.PI * mu * a)), { a, mu }), R, { a, nt: 12, np: 24 });
    // exact in principle (mean-value property); the quadrature error is largest near contact
    closeArray(avg, matvec(rpyPair(R, { a, mu }), f), { rtol: 1e-7 });
  }
});

test("shear flow is undone exactly by reversing it; diffusion is not", () => {
  const pts = Array.from({ length: 50 }, (_, k) => [0.3 * Math.cos(k), 0.3 * Math.sin(2 * k)]);
  const det = shearReversal(pts, { strain: 8, steps: 100 });
  assert.ok(rmsDistance(det[100], pts) > 0.5);
  assert.ok(rmsDistance(det.at(-1), pts) < 1e-12);
  const noisy = shearReversal(pts, { strain: 8, steps: 100, D: 1e-3 });
  assert.ok(rmsDistance(noisy.at(-1), pts) > 0.1);
  assert.deepEqual(noisy, shearReversal(pts, { strain: 8, steps: 100, D: 1e-3 })); // seeded
});
