// Basics pages: values from problems 1.6, 1.7, 1.8, 2.1 and the simple-shear example.
import { test } from "node:test";
import assert from "node:assert/strict";
import { matmul, transpose, eye } from "../src/core/linalg.js";
import {
  rotation2, transformVector, transformTensor, bilinear, trace, doubleContraction,
  symPart, antiPart, eigSym2, omegaZ, traction2, stressTorqueZ, mohr, expm2,
} from "../src/physics/tensor.js";
import { close, closeArray } from "./helpers.js";

test("problem 1.6: bilinear in each argument, the quadratic form is not linear", () => {
  const B = [[2, 1], [-3, 4]];
  const a = [0.7, -1.2], b = [1.5, 0.4];
  close(bilinear(B, a.map((x) => 2 * x), b.map((x) => 3 * x)), 6 * bilinear(B, a, b), { rtol: 1e-14 });
  closeArray([[1, 0], [0, 1]].map((ei) => [[1, 0], [0, 1]].map((ej) => bilinear(B, ei, ej))), B, { rtol: 0 });
  const q = (x) => bilinear(B, x, x);
  close(q(a.map((x) => 2 * x)), 4 * q(a), { rtol: 1e-14 });
});

test("problem 1.7: components change, the value of the form does not", () => {
  const T = [[2, 0, 0], [0, 1, 0], [0, 0, 3]];
  const Q = [[0, 1, 0], [-1, 0, 0], [0, 0, 1]];
  closeArray(matmul(Q, transpose(Q)), eye(3), { rtol: 0 });
  const Tp = transformTensor(Q, T);
  closeArray(Tp, [[1, 0, 0], [0, 2, 0], [0, 0, 3]], { atol: 0 });
  const a = [1, 0, 0];
  const ap = transformVector(Q, a);
  close(bilinear(Tp, ap, ap), 2, { rtol: 0 });
  close(bilinear(T, a, a), 2, { rtol: 0 });
});

test("rotation leaves trace, double contraction, eigenvalues and form values unchanged", () => {
  const T = [[2, 1], [-3, 4]], S = [[0.5, -1], [2, 0.3]];
  const a = [0.3, 1.1], b = [-0.8, 0.6];
  for (const phi of [0.1, 0.9, 2.5]) {
    const Q = rotation2(phi);
    closeArray(matmul(Q, transpose(Q)), eye(2), { atol: 1e-15 });
    const Tp = transformTensor(Q, T), Sp = transformTensor(Q, S);
    close(trace(Tp), trace(T), { rtol: 1e-14 });
    close(doubleContraction(Tp, Sp), doubleContraction(T, S), { rtol: 1e-13 });
    close(bilinear(Tp, transformVector(Q, a), transformVector(Q, b)), bilinear(T, a, b), { rtol: 1e-13 });
    const sym = symPart(T), symp = transformTensor(Q, sym);
    closeArray(eigSym2(symp).values, eigSym2(sym).values, { rtol: 1e-13 });
  }
});

test("problem 1.8 and simple shear: W does not stretch, E has eigenvalues +-gamma/2 at 45 degrees", () => {
  const g = 1.6;
  const L = [[0, g], [0, 0]]; // u_x = g y
  const E = symPart(L), W = antiPart(L);
  closeArray(E, [[0, g / 2], [g / 2, 0]], { rtol: 0 });
  closeArray(W, [[0, g / 2], [-g / 2, 0]], { rtol: 0 });
  close(omegaZ(L), -g / 2, { rtol: 0 });
  for (const x of [[1, 0], [0.3, -2]]) close(bilinear(W, x, x), 0, { atol: 1e-15 });
  assert.equal(eigSym2([[0.3, 0], [0, 0.3]]).degenerate, true); // isotropic: every direction is principal
  assert.equal(eigSym2([[0, 0], [0, -0.5]]).degenerate, false);
  const { values, angle } = eigSym2(E);
  closeArray(values, [g / 2, -g / 2], { rtol: 1e-14 });
  close(angle, Math.PI / 4, { rtol: 1e-14 });
  // the flow map of simple shear is x -> x + g t y
  closeArray(expm2(L.map((r) => r.map((v) => v * 0.7))), [[1, g * 0.7], [0, 1]], { atol: 1e-13 });
  // a rigid rotation keeps lengths
  const R = expm2([[0, -1.3], [1.3, 0]]);
  closeArray(R, [[Math.cos(1.3), -Math.sin(1.3)], [Math.sin(1.3), Math.cos(1.3)]], { atol: 1e-13 });
});

test("problem 2.1: traction is linear in the normal; symmetric stress exerts no torque", () => {
  const p = 1.2, tau = 0.5;
  const sigma = [[-p, tau], [tau, -p]];
  const tx = traction2(sigma, [1, 0]).t, ty = traction2(sigma, [0, 1]).t;
  closeArray(tx, [-p, tau], { rtol: 0 });
  const n = [Math.SQRT1_2, Math.SQRT1_2];
  closeArray(traction2(sigma, n).t, tx.map((v, i) => (v + ty[i]) * Math.SQRT1_2), { rtol: 1e-15 });
  close(traction2(sigma, n).normal, -p + tau, { rtol: 1e-14 });
  close(traction2(sigma, n).shear, 0, { atol: 1e-15 });
  close(stressTorqueZ(sigma), 0, { atol: 0 });
  close(stressTorqueZ([[0, 1], [0.4, 0]]), -0.6, { rtol: 1e-15 });
  // Mohr circle: every (normal, shear) pair lies on it
  const sg = [[1, 0.7], [0.7, -0.4]], { centre, radius } = mohr(sg);
  for (const th of [0.2, 1.1, 2.9]) {
    const r = traction2(sg, [Math.cos(th), Math.sin(th)]);
    close(Math.hypot(r.normal - centre, r.shear), radius, { rtol: 1e-13 });
  }
});
