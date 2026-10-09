// Many-body resistance and its pairwise-additive approximation.
//
// The mobility of RPY / Oseen is pairwise additive: M_ab depends only on particles a and b.
// The resistance R = M^-1 is not: each block of R depends on every particle. To make this
// visible we compare R with the "pairwise sum" built from isolated pairs,
//
//   R2B_aa = R0 + sum_{c != a} ( R^(a,c)_aa - R0 ),     R2B_ab = R^(a,b)_ab ,
//
// where R^(a,b) is the resistance of particles a and b alone and R0 is the isolated-sphere
// resistance. R - R2B is the part that needs three or more bodies (multi-body reflections).
import { zeros, solve, eye, getBlock, setBlock, cholesky, choleskySolve } from "../core/linalg.js";
import { mobility } from "../physics/rpy.js";
import { normalResistance } from "../physics/lubrication.js";
import { axialCollocation } from "../physics/collocation.js";

// Translational resistance R = M^-1 of N equal spheres (3N x 3N).
export function resistance(X, { a = 1, mu = 1, model = "rpy", lubrication = false, hc = 0.2 * a } = {}) {
  if (lubrication) return normalResistance(X, { a, mu, hc }).R;
  const M = mobility(X, { a, mu, model });
  // Oseen can be indefinite, so fall back to LU when Cholesky fails.
  try { return choleskySolve(cholesky(M), eye(M.length)); } catch { return solve(M, eye(M.length)); }
}

/**
 * Pairwise-additive resistance from isolated pairs, for any block size b.
 * pairResistance(p, q) -> (2b x 2b) resistance of particles p and q alone, ordered (p, q)
 * self -> (b x b) isolated-particle resistance
 */
export function pairwiseSum(N, b, pairResistance, self) {
  const R = zeros(b * N);
  for (let p = 0; p < N; p++) setBlock(R, b * p, b * p, self);
  for (let p = 0; p < N; p++)
    for (let q = p + 1; q < N; q++) {
      const R2 = pairResistance(p, q);
      for (const [i, I] of [[0, p], [1, q]]) {
        const blk = getBlock(R2, b * i, b * i, b);
        const cur = getBlock(R, b * I, b * I, b);
        setBlock(R, b * I, b * I, cur.map((row, k) => row.map((v, l) => v + blk[k][l] - self[k][l])));
      }
      setBlock(R, b * p, b * q, getBlock(R2, 0, b, b));
      setBlock(R, b * q, b * p, getBlock(R2, b, 0, b));
    }
  return R;
}

// 3D translational version for RPY / Oseen (optionally with the teaching lubrication model).
export function pairwiseSumResistance(X, opts = {}) {
  const { a = 1, mu = 1 } = opts;
  const self = eye(3).map((row) => row.map((d) => d * 6 * Math.PI * mu * a));
  return pairwiseSum(X.length, 3, (p, q) => resistance([X[p], X[q]], opts), self);
}

// Coaxial spheres moving along the axis: exact (collocation) N-body resistance and the
// pairwise sum built from exact two-body solutions with the same order L.
export function axialManyBody(centers, radii, { L = 16, mu = 1 } = {}) {
  const full = axialCollocation(centers, radii, { L, mu, withCondition: false });
  const N = radii.length;
  let worst = full.boundaryError;
  const pairR = (p, q) => {
    const res = axialCollocation([centers[p], centers[q]], [radii[p], radii[q]], { L, mu, withCondition: false });
    worst = Math.max(worst, res.boundaryError);
    return res.R;
  };
  const R2B = pairwiseSumAxial(N, radii, pairR, mu);
  return { R: full.R, R2B, boundaryError: worst };
}

function pairwiseSumAxial(N, radii, pairR, mu) {
  // block size 1, but the isolated value differs per particle (radii may differ)
  const R = zeros(N);
  for (let p = 0; p < N; p++) R[p][p] = 6 * Math.PI * mu * radii[p];
  for (let p = 0; p < N; p++)
    for (let q = p + 1; q < N; q++) {
      const R2 = pairR(p, q);
      R[p][p] += R2[0][0] - 6 * Math.PI * mu * radii[p];
      R[q][q] += R2[1][1] - 6 * Math.PI * mu * radii[q];
      R[p][q] = R2[0][1];
      R[q][p] = R2[1][0];
    }
  return R;
}
