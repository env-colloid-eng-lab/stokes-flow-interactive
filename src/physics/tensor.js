// Small tensor utilities for the basics pages (chapter 1-2): basis rotation, invariants,
// symmetric / antisymmetric split, traction on a plane, linear flows.
// Conventions follow the textbook: L_ij = d_j u_i, E = (L + L^T)/2, W = (L - L^T)/2,
// W r = omega x r, t_i = sigma_ij n_j (first index: force direction, second: plane normal).
import { matmul, transpose, matvec, dot } from "../core/linalg.js";

// Rows of Q are the new basis vectors: Q_ij = e'_i . e_j. In 2D the basis is turned by phi.
export const rotation2 = (phi) => [[Math.cos(phi), Math.sin(phi)], [-Math.sin(phi), Math.cos(phi)]];

// Components in the new basis: v' = Q v, T' = Q T Q^T (eq. 1.2).
export const transformVector = (Q, v) => matvec(Q, v);
export const transformTensor = (Q, T) => matmul(matmul(Q, T), transpose(Q));

// Value of the bilinear form T(a, b) = a_i T_ij b_j.
export const bilinear = (T, a, b) => dot(a, matvec(T, b));

export const trace = (T) => T.reduce((s, row, i) => s + row[i], 0);
export const doubleContraction = (A, B) => A.reduce((s, row, i) => s + row.reduce((t, v, j) => t + v * B[i][j], 0), 0);

export const symPart = (L) => L.map((row, i) => row.map((v, j) => (v + L[j][i]) / 2));
export const antiPart = (L) => L.map((row, i) => row.map((v, j) => (v - L[j][i]) / 2));

// Eigen-decomposition of a symmetric 2x2 matrix: values descending, unit vectors.
// degenerate: the two eigenvalues coincide (isotropic), so every direction is principal and
// the returned vectors are arbitrary.
export function eigSym2(S) {
  const [[a, b], [, d]] = S;
  const m = (a + d) / 2, r = Math.hypot((a - d) / 2, b);
  const theta = 0.5 * Math.atan2(2 * b, a - d); // direction of the larger eigenvalue
  const v1 = [Math.cos(theta), Math.sin(theta)], v2 = [-Math.sin(theta), Math.cos(theta)];
  const degenerate = r <= 1e-12 * (Math.abs(m) + 1);
  return { values: [m + r, m - r], vectors: [v1, v2], angle: theta, degenerate };
}

// Planar rotation rate omega_z = W_21 (so that W r = omega x r).
export const omegaZ = (L) => antiPart(L)[1][0];

// Traction t = sigma n and its normal / tangential parts. The tangent is n rotated by +90 deg.
export function traction2(sigma, n) {
  const t = matvec(sigma, n);
  const s = [-n[1], n[0]];
  return { t, normal: dot(t, n), shear: dot(t, s), tangent: s };
}

// Net torque per unit volume on a small element from the stress (z component): sigma_21 - sigma_12.
export const stressTorqueZ = (sigma) => sigma[1][0] - sigma[0][1];

// Mohr circle of a symmetric 2D stress: centre and radius.
export function mohr(sigma) {
  const [s1, s2] = eigSym2(sigma).values;
  return { centre: (s1 + s2) / 2, radius: (s1 - s2) / 2 };
}

// Exact flow map of the linear field u = L x: x(t) = exp(L t) x0 (2x2, by scaling and squaring).
export function expm2(A) {
  const nrm = Math.max(...A.flat().map(Math.abs));
  const s = Math.max(0, Math.ceil(Math.log2(nrm + 1e-300)) + 4);
  const B = A.map((row) => row.map((v) => v / 2 ** s));
  // Taylor series of exp(B) to high order (B is small)
  let E = [[1, 0], [0, 1]], term = [[1, 0], [0, 1]];
  for (let k = 1; k <= 12; k++) {
    term = matmul(term, B).map((row) => row.map((v) => v / k));
    E = E.map((row, i) => row.map((v, j) => v + term[i][j]));
  }
  for (let k = 0; k < s; k++) E = matmul(E, E);
  return E;
}
