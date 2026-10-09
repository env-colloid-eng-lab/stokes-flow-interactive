// Brownian displacements and sedimentation equilibrium (overdamped, isothermal).
// Units: kBT, the buoyant weight m_b g and the mobility scale are parameters; nothing is fixed to SI.
import { cholesky, matvec, zeros } from "../core/linalg.js";
import { mulberry32, normalSampler } from "../core/random.js";

/**
 * Noise factor B used as dX = sqrt(2 kBT dt) B xi.
 *   "cholesky":    B = L with L L^T = M (correct: B B^T = M)
 *   "independent": B = diag(sqrt(M_ii)) (each particle on its own; drops the hydrodynamic correlation)
 *   "elementwise": B_ij = sign(M_ij) sqrt|M_ij| (the element-wise square root; B B^T != M)
 */
export function noiseFactor(M, method = "cholesky") {
  if (method === "cholesky") return cholesky(M);
  if (method === "independent") return M.map((row, i) => row.map((_, j) => (i === j ? Math.sqrt(M[i][i]) : 0)));
  if (method === "elementwise") return M.map((row) => row.map((v) => Math.sign(v) * Math.sqrt(Math.abs(v))));
  throw new RangeError(`unknown method ${method}`);
}

// One displacement sqrt(2 kBT dt) B xi for a given noise factor B.
export const applyNoise = (B, kBT, dt, randn) => matvec(B, B.map(() => randn())).map((v) => Math.sqrt(2 * kBT * dt) * v);

// n samples of dX at a frozen configuration (no force, no drift); B = noiseFactor(M, method).
export function sampleDisplacements(B, { kBT = 1, dt = 1, n = 1000, seed = 1 } = {}) {
  const randn = normalSampler(mulberry32(seed));
  return Array.from({ length: n }, () => applyNoise(B, kBT, dt, randn));
}

// Sample covariance (about the sample mean).
export function covariance(samples) {
  const n = samples.length, d = samples[0].length;
  const mean = Array(d).fill(0);
  for (const s of samples) s.forEach((v, i) => (mean[i] += v / n));
  const C = zeros(d);
  for (const s of samples)
    for (let i = 0; i < d; i++) for (let j = 0; j < d; j++) C[i][j] += ((s[i] - mean[i]) * (s[j] - mean[j])) / (n - 1);
  return { mean, C };
}

// ---------------------------------------------------------------------
// One dimension: height z > 0 above a wall, potential V = m_b g z, mobility M(z).
// ---------------------------------------------------------------------

// Teaching model of wall hindrance: M(z) = M0 (1 - beta e^{-z/lambda}), 0 <= beta < 1.
export function hinderedMobility({ M0 = 1, beta = 0.8, lambda = 0.5 } = {}) {
  if (!(beta >= 0 && beta < 1 && lambda > 0 && M0 > 0)) throw new RangeError("0 <= beta < 1, lambda > 0, M0 > 0 required");
  return {
    M: (z) => M0 * (1 - beta * Math.exp(-z / lambda)),
    dM: (z) => ((M0 * beta) / lambda) * Math.exp(-z / lambda),
  };
}

/**
 * Euler-Maruyama (Ito) steps for an ensemble of heights, reflecting at z = 0:
 *   dz = [-M m_b g + (drift ? kBT M' : 0)] dt + sqrt(2 kBT M dt) xi
 * Mutates and returns zs.
 */
export function stepHeights(zs, { mob, kBT = 1, mg = 1, dt = 1e-3, steps = 1, drift = true, randn }) {
  for (let s = 0; s < steps; s++)
    for (let k = 0; k < zs.length; k++) {
      const z = zs[k], M = mob.M(z);
      const zn = z + (-M * mg + (drift ? kBT * mob.dM(z) : 0)) * dt + Math.sqrt(2 * kBT * M * dt) * randn();
      zs[k] = Math.abs(zn);
    }
  return zs;
}

// Zero-flux densities on a grid, normalised on [0, zmax]: with the thermal drift the Boltzmann
// distribution e^{-m_b g z / kBT}; without it e^{-m_b g z / kBT} / M(z).
export function stationaryDensity(zs, { mob, kBT = 1, mg = 1, drift = true, zmax = 12 } = {}) {
  const shape = (z) => Math.exp((-mg * z) / kBT) / (drift ? 1 : mob.M(z));
  // normalise with the composite Simpson rule on [0, zmax]
  const N = 2000, h = zmax / N;
  let I = shape(0) + shape(zmax);
  for (let k = 1; k < N; k++) I += (k % 2 ? 4 : 2) * shape(k * h);
  I *= h / 3;
  return zs.map((z) => shape(z) / I);
}

// Particle flux J = A c - d(D c)/dz with D = kBT M and A = -M m_b g + (drift ? kBT M' : 0).
export function flux(c, dc, z, { mob, kBT = 1, mg = 1, drift = true }) {
  const M = mob.M(z), dM = mob.dM(z);
  const A = -M * mg + (drift ? kBT * dM : 0);
  return A * c - kBT * (dM * c + M * dc);
}
