// Rotne-Prager-Yamakawa translational mobility in free space. Port of StokesLab.jl.
// Positions X are an array of points [[x,y,z], ...]; force/velocity vectors are
// flattened particle by particle: [x1,y1,z1, x2,y2,z2, ...] (same as Julia's vec).
import { eye, zeros, norm, vsub, setBlock, transpose, matvec } from "../core/linalg.js";
import { oseen, lapOseen } from "./oseen.js";

export const selfMobility = (a = 1, mu = 1) => 1 / (6 * Math.PI * mu * a);

// Equal radii. Overlapping pairs (r < 2a) use the standard regularised form, which is a
// way to keep M positive definite, NOT contact physics.
export function rpyPair(r, { a = 1, mu = 1 } = {}) {
  if (!(a > 0 && mu > 0)) throw new RangeError("a, mu > 0 required");
  const R = norm(r);
  if (R === 0) return eye(3).map((row) => row.map((d) => d * selfMobility(a, mu)));
  const e = r.map((x) => x / R);
  let cI, cE;
  if (R >= 2 * a) {
    const c = 1 / (8 * Math.PI * mu * R);
    cI = c * (1 + (2 * a * a) / (3 * R * R));
    cE = c * (1 - (2 * a * a) / (R * R));
  } else {
    const c = selfMobility(a, mu);
    cI = c * (1 - (9 * R) / (32 * a));
    cE = c * ((3 * R) / (32 * a));
  }
  return eye(3).map((row, i) => row.map((d, j) => cI * d + cE * e[i] * e[j]));
}

// Unequal radii, non-overlapping: (1 + a1^2/6 lap)(1 + a2^2/6 lap) G = G + (a1^2+a2^2)/6 lap G.
export function rpyPairUnequal(r, a1, a2, { mu = 1 } = {}) {
  if (!(norm(r) > a1 + a2)) throw new RangeError("spheres must not overlap");
  const G = oseen(r, { mu }), L = lapOseen(r, { mu });
  const c = (a1 * a1 + a2 * a2) / 6;
  return G.map((row, i) => row.map((g, j) => g + c * L[i][j]));
}

// Pair coefficients along / across the line of centres for equal radii (r >= 2a).
export function rpyParallelPerp(r, { a = 1, mu = 1 } = {}) {
  return {
    par: (1 / (4 * Math.PI * mu * r)) * (1 - (2 * a * a) / (3 * r * r)),
    perp: (1 / (8 * Math.PI * mu * r)) * (1 + (2 * a * a) / (3 * r * r)),
  };
}

// model: "rpy" | "oseen" | "self"
export function mobility(X, { a = 1, mu = 1, model = "rpy" } = {}) {
  if (!(a > 0 && mu > 0)) throw new RangeError("a, mu > 0 required");
  if (!["rpy", "oseen", "self"].includes(model)) throw new RangeError("unknown model");
  const N = X.length;
  const M = zeros(3 * N);
  const m0 = selfMobility(a, mu);
  for (let p = 0; p < N; p++) {
    for (let k = 0; k < 3; k++) M[3 * p + k][3 * p + k] = m0;
    for (let q = p + 1; q < N; q++) {
      const r = vsub(X[p], X[q]);
      const block =
        model === "rpy" ? rpyPair(r, { a, mu }) : model === "oseen" ? oseen(r, { mu }) : zeros(3);
      setBlock(M, 3 * p, 3 * q, block);
      setBlock(M, 3 * q, 3 * p, transpose(block));
    }
  }
  return M;
}

export const flatten = (P) => P.flat();
export const unflatten = (v) => Array.from({ length: v.length / 3 }, (_, p) => v.slice(3 * p, 3 * p + 3));

export function mobilityVelocity(X, F, opts = {}) {
  return unflatten(matvec(mobility(X, opts), flatten(F)));
}
