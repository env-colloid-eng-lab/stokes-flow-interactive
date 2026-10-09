// Pieces of the Stokesian Dynamics construction (chapter 14).

import { solve, inv, matmul, sub, getBlock, zeros, matvec } from "../core/linalg.js";

/**
 * Eliminate the rigid-particle constraint with a Schur complement (eq. 14.1).
 * Generalised mobility in blocks (v: motion, s: strain/stresslet):
 *   [ q_v      ]   [ Mvv Mvs ] [ g_v ]
 *   [ -e_inf   ] = [ Msv Mss ] [ g_s ]
 * Returns the reduced operator Mvv - Mvs Mss^-1 Msv, the ambient term -Mvs Mss^-1 e_inf,
 * and, for a given g_v, the induced g_s and q_v.
 */
export function schurReduce(M, nv, { gv, einf }) {
  const n = M.length, ns = n - nv;
  const Mvv = getBlock(M, 0, 0, nv), Mvs = getBlock(M, 0, nv, nv, ns);
  const Msv = getBlock(M, nv, 0, ns, nv), Mss = getBlock(M, nv, nv, ns);
  // one factorisation of Mss for both right-hand sides
  const rhs = Msv.map((row, i) => [...row, einf[i]]);
  const sol = solve(Mss, rhs);
  const MssInvMsv = sol.map((row) => row.slice(0, nv)), MssInvE = sol.map((row) => row[nv]);
  const reduced = sub(Mvv, matmul(Mvs, MssInvMsv));
  const ambient = matvec(Mvs, MssInvE).map((x) => -x);
  // lower row: Msv g_v + Mss g_s = -e_inf  =>  g_s = -(Mss^-1 Msv g_v + Mss^-1 e_inf)
  const gs = matvec(MssInvMsv, gv).map((x, i) => -x - MssInvE[i]);
  const qv = matvec(reduced, gv).map((x, i) => x + ambient[i]);
  return { reduced, ambient, gs, qv };
}

// The same problem solved independently through the resistance R = M^-1 (the form SD uses):
// g = R q with q = (q_v, -e_inf) and g_v given, so
//   q_v = R_vv^-1 (g_v + R_vs e_inf),   g_s = R_sv q_v - R_ss e_inf.
export function fullSolve(M, nv, { gv, einf }) {
  const n = M.length, ns = n - nv;
  const R = inv(M);
  const Rvv = getBlock(R, 0, 0, nv), Rvs = getBlock(R, 0, nv, nv, ns);
  const Rsv = getBlock(R, nv, 0, ns, nv), Rss = getBlock(R, nv, nv, ns);
  const RvsE = matvec(Rvs, einf), RssE = matvec(Rss, einf);
  const qv = solve(Rvv, gv.map((g, i) => g + RvsE[i]));
  const gs = matvec(Rsv, qv).map((x, i) => x - RssE[i]);
  return { gs, qv };
}

/**
 * The resistance correction of eq. (14.2), for any matrices of matching size:
 *   R_SD = R_far + sum_pairs (R2B_exact - R2B_far)
 * Here R2B_* are already the pairwise sums (see manyBody.pairwiseSum), so the correction is
 * their difference; it vanishes for well separated pairs and is exact for two particles.
 */
export function sdResistance(Rfar, R2Bexact, R2Bfar) {
  return Rfar.map((row, i) => row.map((v, j) => v + R2Bexact[i][j] - R2Bfar[i][j]));
}

// Isolated-sphere 11x11 generalised resistance diag(6 pi mu a I3, 8 pi mu a^3 I3, 20 pi/3 mu a^3 I5).
export function isolatedResistance11({ a = 1, mu = 1 } = {}) {
  const R = zeros(11);
  const d = [...Array(3).fill(6 * Math.PI * mu * a), ...Array(3).fill(8 * Math.PI * mu * a ** 3), ...Array(5).fill((20 * Math.PI / 3) * mu * a ** 3)];
  d.forEach((v, i) => (R[i][i] = v));
  return R;
}
