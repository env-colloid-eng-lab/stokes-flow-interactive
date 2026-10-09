// Jeffrey-Onishi (1984) axisymmetric resistance X^A. Port of TwoSphereLab.jl.
// The coefficients are generated from the recurrence, eqs. (3.6)-(3.9), (3.15) (in joRecurrence.js);
// no table is read. Exact mode uses BigInt rationals.
import { Rational } from "../core/rational.js";
import { solve } from "../core/linalg.js";
import { joFamilyPolynomials } from "./joRecurrence.js";

// Returns polys[k] = coefficients of f_k(lambda) in increasing powers of lambda, k = 0..K.
export const joPolynomials = (K, opts) => joFamilyPolynomials("XA", K, opts);

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
