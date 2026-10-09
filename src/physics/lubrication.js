// Teaching model: ONLY a normal leading-order lubrication correction (StokesLab.jl).
// This is not the matched two-body correction of full Stokesian Dynamics.
import { cholesky, solve, eye, norm, vsub, matvec } from "../core/linalg.js";
import { mobility, flatten, unflatten } from "./rpy.js";

// zeta_n = (3 pi mu a^2 / 2)(1/h - 1/hc) for h < hc, added as zeta_n b b^T.
export const normalLubricationZeta = (h, { a = 1, mu = 1, hc = 0.2 * a } = {}) =>
  h < hc ? ((3 * Math.PI * mu * a * a) / 2) * (1 / h - 1 / hc) : 0;

// Returns the resistance matrix and a list of the pairs that received a correction,
// so that the UI can say which components came from which term.
export function normalResistance(X, { a = 1, mu = 1, hc = 0.2 * a } = {}) {
  if (!(hc > 0)) throw new RangeError("hc > 0 required");
  const M = mobility(X, { a, mu });
  cholesky(M); // throws if M is not positive definite
  const n = M.length;
  const R = solve(M, eye(n));
  const pairs = [];
  for (let p = 0; p < X.length; p++) {
    for (let q = p + 1; q < X.length; q++) {
      const r = vsub(X[p], X[q]);
      const gap = norm(r) - 2 * a;
      if (!(gap > 0)) throw new RangeError("lubrication model requires gap > 0");
      if (gap < hc) {
        const e = r.map((x) => x / norm(r));
        const b = new Array(n).fill(0);
        for (let k = 0; k < 3; k++) { b[3 * p + k] = e[k]; b[3 * q + k] = -e[k]; }
        const zeta = normalLubricationZeta(gap, { a, mu, hc });
        for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) R[i][j] += zeta * b[i] * b[j];
        pairs.push({ p, q, gap, zeta });
      }
    }
  }
  return { R, pairs };
}

export function velocity(X, F, { a = 1, mu = 1, lubrication = false, hc = 0.2 * a } = {}) {
  const f = flatten(F);
  const u = lubrication
    ? solve(normalResistance(X, { a, mu, hc }).R, f)
    : matvec(mobility(X, { a, mu }), f);
  return unflatten(u);
}
