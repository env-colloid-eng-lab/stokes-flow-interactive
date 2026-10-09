import { test } from "node:test";
import assert from "node:assert/strict";
import {
  matmul, transpose, solve, inv, eye, cholesky, lstsq, symEig, cond, fnorm, sub, asymmetry, isPosDef, matvec,
} from "../src/core/linalg.js";
import { mulberry32, normalSampler } from "../src/core/random.js";
import { closeArray, close } from "./helpers.js";

const randn = normalSampler(mulberry32(1));
const randMat = (m, n) => Array.from({ length: m }, () => Array.from({ length: n }, randn));

test("LU solve and inverse", () => {
  const A = randMat(7, 7);
  const x = Array.from({ length: 7 }, randn);
  closeArray(solve(A, matvec(A, x)), x, { atol: 1e-11 });
  closeArray(matmul(A, inv(A)), eye(7), { atol: 1e-11 });
});

test("Cholesky reproduces an SPD matrix and rejects an indefinite one", () => {
  const B = randMat(6, 6);
  const A = matmul(B, transpose(B)).map((r, i) => r.map((v, j) => v + (i === j ? 0.1 : 0)));
  const L = cholesky(A);
  closeArray(matmul(L, transpose(L)), A, { atol: 1e-12 });
  assert.equal(isPosDef([[1, 2], [2, 1]]), false);
});

test("least squares by QR matches the normal equations", () => {
  const A = randMat(30, 6), b = randMat(30, 2);
  const X = lstsq(A, b);
  const At = transpose(A);
  closeArray(X, solve(matmul(At, A), matmul(At, b)), { atol: 1e-10 });
});

test("symmetric eigen-decomposition", () => {
  const B = randMat(8, 8);
  const A = matmul(B, transpose(B));
  const { values, vectors } = symEig(A);
  for (let i = 1; i < values.length; i++) assert.ok(values[i] >= values[i - 1]);
  const D = values.map((v, i) => values.map((_, j) => (i === j ? v : 0)));
  const rec = matmul(matmul(vectors, D), transpose(vectors));
  assert.ok(fnorm(sub(rec, A)) / fnorm(A) < 1e-12);
  closeArray(matmul(transpose(vectors), vectors), eye(8), { atol: 1e-12 });
  close(cond([[2, 0], [0, 0.5]]), 4, { rtol: 1e-12 });
  assert.equal(asymmetry(A), 0);
});
