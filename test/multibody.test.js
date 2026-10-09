// Mirrors reference/julia/regression_multibody.jl (sections 01-08) and problem 9.2.
import { test } from "node:test";
import assert from "node:assert/strict";
import {
  matmul, transpose, eye, sub, fnorm, norm, cross, vsub, eigvalsSym, isPosDef, matvec, dot, zeros,
} from "../src/core/linalg.js";
import { mulberry32, normalSampler } from "../src/core/random.js";
import { projector, oseen, gradOseen, lapOseen, pressure, dipole, sphereFlow } from "../src/physics/oseen.js";
import { rpyPair, mobility, flatten, rpyPairUnequal } from "../src/physics/rpy.js";
import { velocity, normalResistance } from "../src/physics/lubrication.js";
import { heunStep, brownianIncrement } from "../src/physics/integrate.js";
import { readCsv, close, closeArray } from "./helpers.js";

const randn = normalSampler(mulberry32(20260913));

test("01 projection", () => {
  const q = [1, 2, 3], P = projector(q);
  assert.ok(fnorm(sub(matmul(P, P), P)) < 1e-14);
  assert.ok(norm(matvec(P, q)) < 1e-14);
  closeArray(eigvalsSym(P), [0, 1, 1], { atol: 1e-14 });
});

test("02 Oseen derivatives, pressure and rotation covariance", () => {
  const r = [1.3, -0.7, 2.1], f = [0.4, 0.8, -0.3];
  const H = gradOseen(r), h = 1e-5;
  for (let k = 0; k < 3; k++) {
    const e = [0, 0, 0]; e[k] = h;
    const Gp = oseen(r.map((x, i) => x + e[i])), Gm = oseen(r.map((x, i) => x - e[i]));
    for (let i = 0; i < 3; i++)
      for (let j = 0; j < 3; j++) close((Gp[i][j] - Gm[i][j]) / (2 * h), H[i][j][k], { rtol: 1e-8, atol: 1e-11 });
  }
  // incompressibility: partial_i G_ij = 0
  for (let j = 0; j < 3; j++) assert.ok(Math.abs(H[0][j][0] + H[1][j][1] + H[2][j][2]) < 1e-14);
  // Stokes equation: mu lap(G f) = grad p
  const gradp = [0, 1, 2].map((k) => {
    const e = [0, 0, 0]; e[k] = h;
    return (pressure(r.map((x, i) => x + e[i]), f) - pressure(r.map((x, i) => x - e[i]), f)) / (2 * h);
  });
  assert.ok(norm(vsub(matvec(lapOseen(r), f), gradp)) < 1e-10);
  // G(Q r) = Q G Q^T for a rotation Q
  const c = Math.cos(0.7), s = Math.sin(0.7);
  const Q = [[c, -s, 0], [s, c, 0], [0, 0, 1]];
  closeArray(oseen(matvec(Q, r)), matmul(matmul(Q, oseen(r)), transpose(Q)), { atol: 1e-15 });
});

test("03 dipole decomposition and moving sphere", () => {
  const r = [1.3, -0.7, 2.1], f = [0.4, 0.8, -0.3], d = [0.2, 0.4, 0.1];
  const D = f.map((fi) => d.map((dj) => fi * dj));
  const e = 1e-3;
  const Gm = oseen(r.map((x, i) => x - (e * d[i]) / 2)), Gp = oseen(r.map((x, i) => x + (e * d[i]) / 2));
  const finite = matvec(sub(Gm, Gp), f).map((v) => v / e);
  closeArray(finite, dipole(r, D), { rtol: 1e-6 });
  const trD = D[0][0] + D[1][1] + D[2][2];
  const S = D.map((row, i) => row.map((v, j) => (v + D[j][i]) / 2 - (i === j ? trD / 3 : 0)));
  const C = D.map((row, i) => row.map((v, j) => (v - D[j][i]) / 2));
  closeArray(dipole(r, D), dipole(r, S).map((v, i) => v + dipole(r, C)[i]), { atol: 1e-15 });
  const torque = cross(d, f), R3 = norm(r) ** 3;
  closeArray(dipole(r, C), cross(torque, r).map((v) => v / (8 * Math.PI * R3)), { rtol: 1e-12 });
  assert.ok(norm(dipole(r, eye(3))) < 1e-14);
  const U = [0.2, 0.4, 0.6];
  closeArray(sphereFlow(r.map((x) => x / norm(r)), U), U, { rtol: 1e-12 });
});

test("04 two spheres: RPY eigenvalues, continuity, positivity", () => {
  const X = [[0, 0, 0], [4, 0, 0]];
  const M = mobility(X);
  const m0 = 1 / (6 * Math.PI);
  const pair = rpyPair([4, 0, 0]);
  const cp = (1 - 2 / (3 * 16)) / (4 * Math.PI * 4), ct = (1 + 2 / (3 * 16)) / (8 * Math.PI * 4);
  closeArray([pair[0][0], pair[1][1], pair[2][2]], [cp, ct, ct], { rtol: 1e-13 });
  closeArray(eigvalsSym(M), [m0 + cp, m0 - cp, m0 + ct, m0 - ct, m0 + ct, m0 - ct].sort((a, b) => a - b), { rtol: 1e-12 });
  assert.ok(isPosDef(M));
  closeArray(rpyPair([2 - 1e-8, 0, 0]), rpyPair([2 + 1e-8, 0, 0]), { rtol: 1e-7 });
  const F = [[0, 0, -1], [0, 0, -1]];
  const U = velocity(X, F);
  closeArray(U[0], U[1], { rtol: 1e-14 });
  assert.ok(Math.abs(U[0][2]) > m0);
  assert.ok(dot(flatten(F), matvec(M, flatten(F))) > 0);
  // unequal-radius formula reduces to the equal-radius one
  closeArray(rpyPairUnequal([3, 1, -2], 1, 1), rpyPair([3, 1, -2]), { rtol: 1e-13 });
});

test("problem 9.2: Oseen loses positive definiteness, RPY does not", () => {
  const X = Array.from({ length: 10 }, (_, i) => [2.01 * i, 0, 0]);
  assert.ok(eigvalsSym(mobility(X, { model: "oseen" }))[0] < 0);
  assert.ok(eigvalsSym(mobility(X, { model: "rpy" }))[0] > 0);
});

test("05 three-body sedimentation: Heun convergence and three_spheres.csv", () => {
  const X0 = [[-3, 0, 0], [3, 0, 0], [0, 0, 4]];
  const F = [[0, 0, -1], [0, 0, -1], [0, 0, -1]];
  const vfun = (X) => velocity(X, F);
  const integrate = (n) => { let X = X0; for (let k = 0; k < n; k++) X = heunStep(X, vfun, 1 / n); return X.flat(); };
  const x1 = integrate(10), x2 = integrate(20), x3 = integrate(40);
  const ratio = norm(vsub(x1, x2)) / norm(vsub(x2, x3));
  assert.ok(ratio > 3.5 && ratio < 4.5, `ratio ${ratio}`);
  const { rows } = readCsv("three_spheres");
  let X = X0;
  rows.forEach((row, k) => {
    close(row[0], k * 0.05, { rtol: 1e-12, atol: 1e-15 });
    closeArray(X.flat(), row.slice(1), { rtol: 1e-12, atol: 1e-14 }, `t=${row[0]}`);
    if (k < rows.length - 1) X = heunStep(X, vfun, 0.05);
  });
});

test("06 normal lubrication model reproduces lubrication.csv", () => {
  const F = [[1, 0, 0], [-1, 0, 0]];
  const { rows } = readCsv("lubrication");
  for (const [h, vrpy, vlub, vasym] of rows) {
    const X = [[0, 0, 0], [2 + h, 0, 0]];
    const Ur = velocity(X, F), Ul = velocity(X, F, { lubrication: true });
    close(Ur[0][0] - Ur[1][0], vrpy, { rtol: 1e-9 }, `rpy h=${h}`);
    close(Ul[0][0] - Ul[1][0], vlub, { rtol: 1e-7 }, `lub h=${h}`);
    close((2 * h) / (3 * Math.PI), vasym, { rtol: 1e-13 });
    if (h <= 1e-4) close(Ul[0][0] - Ul[1][0], vasym, { rtol: 0.003 });
  }
  const { pairs } = normalResistance([[0, 0, 0], [2.05, 0, 0], [10, 0, 0]]);
  assert.deepEqual(pairs.map(({ p, q }) => [p, q]), [[0, 1]]);
});

test("07 correlated Brownian increments have covariance 2 kBT M dt", () => {
  const M = mobility([[0, 0, 0], [4, 0, 0]]);
  const kBT = 1, dt = 0.01, n = 100000;
  const C = zeros(6), mean = new Array(6).fill(0);
  const samples = [];
  for (let s = 0; s < n; s++) {
    const y = brownianIncrement(M, kBT, dt, randn);
    samples.push(y);
    y.forEach((v, i) => (mean[i] += v / n));
  }
  for (const y of samples)
    for (let i = 0; i < 6; i++) for (let j = 0; j < 6; j++) C[i][j] += ((y[i] - mean[i]) * (y[j] - mean[j])) / (n - 1);
  const target = M.map((row) => row.map((v) => 2 * kBT * dt * v));
  assert.ok(fnorm(sub(C, target)) / fnorm(target) < 0.025);
  assert.ok(norm(mean) < 0.001);
  // the reference CSV holds the same theoretical covariance
  for (const [i, j, , theory] of readCsv("brownian_covariance").rows) close(target[i - 1][j - 1], theory, { rtol: 1e-12, atol: 1e-18 });
});

test("08 thermal drift: zero flux with drift is Boltzmann, without drift is exp(-V)/M", () => {
  const z = 0.7, h = 1e-6, c = (x) => Math.exp(-2 * x), mob = (x) => 1 + x;
  const cp = (c(z + h) - c(z - h)) / (2 * h);
  assert.ok(Math.abs(mob(z) * (-2 * c(z) - cp)) < 1e-9);
  const wrong = (x) => Math.exp(-2 * x) / mob(x);
  const dc = (mob(z + h) * wrong(z + h) - mob(z - h) * wrong(z - h)) / (2 * h);
  assert.ok(Math.abs(-2 * mob(z) * wrong(z) - dc) < 1e-9);
});
