// Jeffrey-Onishi (1984) axisymmetric resistance X^A. Port of TwoSphereLab.jl.
// The coefficients are generated from the recurrence, eqs. (3.6)-(3.9), (3.15);
// no table is read. Exact mode uses BigInt rationals.
import { Rational, Q } from "../core/rational.js";
import { solve } from "../core/linalg.js";

function binomial(n, k) {
  let r = 1n;
  for (let i = 1n; i <= BigInt(k); i++) r = (r * (BigInt(n) - BigInt(k) + i)) / i;
  return r;
}

// Arithmetic back-ends: exact rationals or Float64.
const exactOps = {
  zero: Q(0), one: Q(1),
  frac: (n, d) => Q(n, d), int: (n) => new Rational(n),
  add: (a, b) => a.add(b), sub: (a, b) => a.sub(b), mul: (a, b) => a.mul(b),
};
const floatOps = {
  zero: 0, one: 1,
  frac: (n, d) => n / d, int: (n) => Number(n),
  add: (a, b) => a + b, sub: (a, b) => a - b, mul: (a, b) => a * b,
};

// Returns polys[k] = coefficients of f_k(lambda) in increasing powers of lambda, k = 0..K.
export function joPolynomials(K, { exact = false } = {}) {
  if (!(Number.isInteger(K) && K >= 0)) throw new RangeError("K must be a nonnegative integer");
  const T = exact ? exactOps : floatOps;
  const pc = new Map(), vc = new Map();
  const choose = (n, s) => T.int(binomial(n + s, n));
  function P(n, p, q) {
    if (n < 1 || p < 0 || q < 0) return T.zero;
    if (p === 0 && q === 0) return n === 1 ? T.one : T.zero;
    const key = `${n},${p},${q}`;
    if (pc.has(key)) return pc.get(key);
    let v = T.zero;
    for (let s = 1; s <= q; s++) {
      const c1 = T.frac(n * (2 * n + 1) * (2 * n * s - n - s + 2), 2 * (n + 1) * (2 * s - 1) * (n + s));
      const c2 = T.frac(n * (2 * n - 1), 2 * (n + 1));
      const c3 = T.frac(n * (4 * n * n - 1), 2 * (n + 1) * (2 * s + 1));
      const inner = T.sub(T.sub(T.mul(c1, P(s, q - s, p - n + 1)), T.mul(c2, P(s, q - s, p - n - 1))),
        T.mul(c3, V(s, q - s - 2, p - n + 1)));
      v = T.add(v, T.mul(choose(n, s), inner));
    }
    pc.set(key, v);
    return v;
  }
  function V(n, p, q) {
    if (n < 1 || p < 0 || q < 0) return T.zero;
    if (p === 0 && q === 0) return n === 1 ? T.one : T.zero;
    const key = `${n},${p},${q}`;
    if (vc.has(key)) return vc.get(key);
    let v = P(n, p, q);
    const fac = T.frac(2 * n, (n + 1) * (2 * n + 3));
    for (let s = 1; s <= q; s++) v = T.sub(v, T.mul(T.mul(fac, choose(n, s)), P(s, q - s, p - n - 1)));
    vc.set(key, v);
    return v;
  }
  const polys = [];
  for (let k = 0; k <= K; k++) {
    const row = [];
    for (let q = 0; q <= k; q++) row.push(T.mul(T.int(2n ** BigInt(k)), P(1, k - q, q)));
    polys.push(row);
  }
  return polys;
}

const toFloat = (c) => (c instanceof Rational ? c.toNumber() : c);
export const evalPoly = (coeffs, x) => coeffs.reduceRight((acc, c) => acc * x + toFloat(c), 0);

// s = 2r/(a1+a2) > 2, lambda = a2/a1. Truncated FAR-FIELD series; near contact the
// convergence must be checked. Returns also the individual terms for display.
export function joXA(s, lambda, polys) {
  if (!(s > 2 && lambda > 0)) throw new RangeError("s > 2, lambda > 0 required");
  let x11 = 0, x12 = 0;
  const terms = [];
  for (let k = 0; k < polys.length; k++) {
    const term = evalPoly(polys[k], lambda) / ((1 + lambda) * s) ** k;
    terms.push(term);
    if (k % 2 === 0) x11 += term;
    else x12 -= (2 / (1 + lambda)) * term;
  }
  return { x11, x12, terms };
}

// Dimensional 2x2 axial resistance, eq. (13.1): f = R U (force on fluid, along the axis).
export function joResistance(a1, a2, r, polys, { mu = 1 } = {}) {
  if (!(a1 > 0 && a2 > 0 && mu > 0)) throw new RangeError("positive radii and viscosity required");
  const s = (2 * r) / (a1 + a2), lambda = a2 / a1;
  const { x11, x12 } = joXA(s, lambda, polys);
  const { x11: x22, x12: x21 } = joXA(s, 1 / lambda, polys);
  return [
    [6 * Math.PI * mu * a1 * x11, 3 * Math.PI * mu * (a1 + a2) * x12],
    [3 * Math.PI * mu * (a1 + a2) * x21, 6 * Math.PI * mu * a2 * x22],
  ];
}

// Deterministic midpoint stepping along the axis, stopping before gap < minGap
// (outside the range where the series was checked). Rows: [t, z1, z2, U1, U2, gap].
export function pairMotion(z0, F, radii, polys, { dt = 0.1, steps = 100, minGap = 1, mu = 1 } = {}) {
  let z = z0.slice();
  const [a1, a2] = radii;
  if (!(z[1] > z[0])) throw new RangeError("centers must be ordered");
  const speed = (x) => solve(joResistance(a1, a2, x[1] - x[0], polys, { mu }), F);
  const rows = [];
  for (let j = 0; j <= steps; j++) {
    const gap = z[1] - z[0] - a1 - a2;
    if (gap < minGap) break;
    const u = speed(z);
    rows.push([j * dt, z[0], z[1], u[0], u[1], gap]);
    if (j === steps) break;
    const mid = z.map((zi, i) => zi + (dt / 2) * u[i]);
    if (mid[1] - mid[0] - a1 - a2 < minGap) break;
    const um = speed(mid);
    z = z.map((zi, i) => zi + dt * um[i]);
  }
  return rows;
}
