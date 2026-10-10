// Stokesian-dynamics-type resistance with forces and torques (model B in the many-body setting):
//   R = (M_far)^-1 + sum_pairs (R2B_exact - R2B_far)          (eq. 14.2, 6 unknowns per sphere)
// M_far is the RPY mobility with rotation, R2B_exact the Jeffrey-Onishi pair resistance
// (src/physics/joFull.js) and R2B_far the inverse of the two-sphere RPY mobility. For two spheres
// the construction is exact. Freely rotating (torque-free) spheres: the rotations are eliminated
// with the Schur complement, leaving a 3N x 3N translational resistance.
// Layout of the 6N system: per sphere [U_p (3), Omega_p (3)].
import { zeros, eye, inv, getBlock, setBlock, solve, sub, matmul } from "../core/linalg.js";
import { joScalarsAtGap } from "../physics/joFull.js";
import { resistance12, rpyMobility12 } from "./pairB.js";
import { pairwiseSum } from "./manyBody.js";
import { sdResistance } from "./sd.js";

// [U1, U2, W1, W2] (pairB layout) -> [U1, W1, U2, W2] (per-sphere layout)
const PERM = [0, 1, 2, 6, 7, 8, 3, 4, 5, 9, 10, 11];
const toPerSphere = (A) => PERM.map((i) => PERM.map((j) => A[i][j]));

function pairVector(X, p, q) {
  return [0, 1, 2].map((k) => X[q][k] - X[p][k]);
}

// two-sphere matrices in the per-sphere layout
export function pairExactFT(X, p, q, { a = 1, mu = 1 } = {}) {
  const rv = pairVector(X, p, q), r = Math.hypot(...rv);
  if (!(r > 2 * a)) throw new RangeError(`spheres ${p + 1} and ${q + 1} touch or overlap`);
  return toPerSphere(resistance12(joScalarsAtGap((r - 2 * a) / a), rv.map((x) => x / r), { a, mu }));
}
export function pairFarFT(X, p, q, { a = 1, mu = 1 } = {}) {
  return toPerSphere(inv(rpyMobility12(pairVector(X, p, q), { a, mu })));
}

// isolated sphere: diag(6 pi mu a I, 8 pi mu a^3 I)
function selfFT({ a = 1, mu = 1 } = {}) {
  const S = zeros(6);
  for (let i = 0; i < 3; i++) { S[i][i] = 6 * Math.PI * mu * a; S[i + 3][i + 3] = 8 * Math.PI * mu * a ** 3; }
  return S;
}

// RPY mobility with rotation for N spheres (6N x 6N, per-sphere layout)
export function rpyMobilityFT(X, { a = 1, mu = 1 } = {}) {
  const N = X.length, M = zeros(6 * N);
  const s = inv(selfFT({ a, mu }));
  for (let p = 0; p < N; p++) setBlock(M, 6 * p, 6 * p, s);
  for (let p = 0; p < N; p++)
    for (let q = p + 1; q < N; q++) {
      const M2 = toPerSphere(rpyMobility12(pairVector(X, p, q), { a, mu }));
      setBlock(M, 6 * p, 6 * q, getBlock(M2, 0, 6, 6));
      setBlock(M, 6 * q, 6 * p, getBlock(M2, 6, 0, 6));
    }
  return M;
}

/** Force-torque SD resistance (6N x 6N, per-sphere layout). */
export function sdResistanceFT(X, { a = 1, mu = 1 } = {}) {
  const N = X.length, self = () => selfFT({ a, mu });
  const Rfar = inv(rpyMobilityFT(X, { a, mu }));
  const exact = pairwiseSum(N, 6, (p, q) => pairExactFT(X, p, q, { a, mu }), self);
  const far = pairwiseSum(N, 6, (p, q) => pairFarFT(X, p, q, { a, mu }), self);
  return sdResistance(Rfar, exact, far);
}

/** Eliminate the rotations of a per-sphere 6N matrix: translational resistance of torque-free spheres. */
export function torqueFreeN(R6) {
  const N = R6.length / 6;
  const ui = [], wi = [];
  for (let p = 0; p < N; p++) for (let k = 0; k < 3; k++) { ui.push(6 * p + k); wi.push(6 * p + 3 + k); }
  const pick = (rows, cols) => rows.map((i) => cols.map((j) => R6[i][j]));
  const A = pick(ui, ui), Bt = pick(ui, wi), B = pick(wi, ui), C = pick(wi, wi);
  return sub(A, matmul(Bt, solve(C, B)));
}

/** Translational resistance (3N x 3N) of freely rotating spheres, SD with forces and torques. */
export const sdTranslationalFT = (X, opts) => torqueFreeN(sdResistanceFT(X, opts));

/**
 * Everything page 10 needs for one configuration, with each exact pair built once:
 * R (3N, torque free, SD with forces and torques) and R2B (pairwise sum of torque-free exact pairs).
 */
export function sdBundle(X, { a = 1, mu = 1 } = {}) {
  const N = X.length, opts = { a, mu };
  const exactPairs = new Map();
  const exactPair = (p, q) => {
    const k = p * N + q;
    if (!exactPairs.has(k)) exactPairs.set(k, pairExactFT(X, p, q, opts));
    return exactPairs.get(k);
  };
  const self6 = () => selfFT(opts);
  const R6 = sdResistance(inv(rpyMobilityFT(X, opts)),
    pairwiseSum(N, 6, exactPair, self6), pairwiseSum(N, 6, (p, q) => pairFarFT(X, p, q, opts), self6));
  const iso = eye(3).map((row) => row.map((d) => d * 6 * Math.PI * mu * a));
  const R2B = pairwiseSum(N, 3, (p, q) => torqueFreeN(exactPair(p, q)), () => iso);
  return { R: torqueFreeN(R6), R2B };
}

/** Pairwise sum of torque-free exact two-sphere resistances (model B), 3N x 3N. */
export const pairwiseSumB = (X, opts) => sdBundle(X, opts).R2B;
