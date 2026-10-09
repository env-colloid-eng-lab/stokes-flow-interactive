// Model A (textbook scope) for two equal spheres.
// Every pair matrix is written as  A_par e e^T + A_perp (I - e e^T)  per 3x3 block,
// so each lab-frame entry can be explained as  A_par e_k e_l + A_perp (delta_kl - e_k e_l).
import { joXA } from "../physics/jo.js";
import { normalLubricationZeta } from "../physics/lubrication.js";
import { zeros } from "../core/linalg.js";
import { selfMobility, rpyParallelPerp } from "../physics/rpy.js";

const PI = Math.PI;

// Mobility coefficients along (par) and across (perp) the line of centres.
export function mobilityCoefficients(r, { a = 1, mu = 1, kind = "rpy" } = {}) {
  const m0 = selfMobility(a, mu);
  const c = kind === "oseen"
    ? { par: 1 / (4 * PI * mu * r), perp: 1 / (8 * PI * mu * r) }
    : rpyParallelPerp(r, { a, mu });
  return { par: { self: m0, cross: c.par }, perp: { self: m0, cross: c.perp } };
}

// Inverse of [[m0, c], [c, m0]]: self = m0/(m0^2 - c^2), cross = -c/(m0^2 - c^2).
const invert2 = ({ self: m0, cross: c }) => {
  const det = m0 * m0 - c * c;
  return { self: m0 / det, cross: -c / det };
};

/**
 * Resistance coefficients for model A.
 * axial: "rpy"  -> inverse of the RPY mobility
 *        "jo"   -> JO X^A series (needs polys); perpendicular entries still from RPY
 * lubrication: add the teaching normal model zeta_n to the axial relative mode.
 * Returns also the provenance of each family for display.
 */
export function resistanceCoefficients(r, { a = 1, mu = 1, kind = "rpy", axial = "rpy", polys = null, lubrication = false, hc = 0.2 * a } = {}) {
  const M = mobilityCoefficients(r, { a, mu, kind });
  const perp = invert2(M.perp);
  let par = invert2(M.par);
  const source = { par: kind === "oseen" ? "Oseen 移動度の逆行列" : "RPY 移動度の逆行列", perp: kind === "oseen" ? "Oseen 移動度の逆行列" : "RPY 移動度の逆行列（トルクゼロ）" };
  let jo = null;
  if (axial === "jo") {
    if (!polys) throw new Error("JO polynomials required");
    const s = r / a; // 2r/(a1+a2) for equal radii
    const res = joXA(s, 1, polys);
    // Convergence estimate: compare with the series truncated about a quarter of the
    // orders earlier (at least one even and one odd term fewer). Equal radii: lambda = 1.
    const K = res.terms.length - 1;
    const Kc = Math.max(0, K - Math.max(2, Math.floor(K / 4)));
    let x11c = 0, x12c = 0;
    for (let k = 0; k <= Kc; k++) (k % 2 === 0 ? (x11c += res.terms[k]) : (x12c -= res.terms[k]));
    const rel = Math.max(Math.abs(res.x11 - x11c) / Math.abs(res.x11), Math.abs(res.x12 - x12c) / Math.abs(res.x12));
    jo = { x11: res.x11, x12: res.x12, K, Kc, change: rel, converged: rel < 1e-4 };
    par = { self: 6 * PI * mu * a * res.x11, cross: 6 * PI * mu * a * res.x12 };
    source.par = `JO 級数 X^A（K = ${jo.K}）`;
  }
  let zeta = 0;
  if (lubrication) {
    const h = r - 2 * a;
    if (!(h > 0)) throw new RangeError("lubrication model requires gap > 0");
    zeta = normalLubricationZeta(h, { a, mu, hc });
    if (zeta > 0) {
      par = { self: par.self + zeta, cross: par.cross - zeta };
      source.par += " ＋ 教材潤滑 ζₙ";
    }
  }
  return { M, R: { par, perp }, source, jo, zeta };
}

// Build the 6x6 lab-frame matrix from block coefficients and the unit vector e (1 -> 2).
export function pairMatrix(coeffs, e) {
  const A = zeros(6);
  for (let al = 0; al < 2; al++)
    for (let be = 0; be < 2; be++) {
      const key = al === be ? "self" : "cross";
      for (let k = 0; k < 3; k++)
        for (let l = 0; l < 3; l++) {
          const ee = e[k] * e[l];
          A[3 * al + k][3 * be + l] = coeffs.par[key] * ee + coeffs.perp[key] * ((k === l ? 1 : 0) - ee);
        }
    }
  return A;
}

// Orthonormal pair frame (e, n1, n2) for e in the x-z plane: n1 in-plane, n2 = y.
export function pairFrame(theta) {
  const e = [Math.cos(theta), 0, Math.sin(theta)];
  const n1 = [-Math.sin(theta), 0, Math.cos(theta)];
  const n2 = [0, 1, 0];
  return { e, n1, n2, axes: [e, n1, n2] };
}

// Q^T A Q with Q = blockdiag(F, F), F having the frame vectors as columns.
export function toPairFrame(A, axes) {
  const n = A.length, B = zeros(n);
  const Q = (i, a) => (Math.floor(i / 3) === Math.floor(a / 3) ? axes[a % 3][i % 3] : 0);
  for (let a = 0; a < n; a++)
    for (let b = 0; b < n; b++) {
      let s = 0;
      for (let i = 0; i < n; i++) {
        const qa = Q(i, a);
        if (!qa) continue;
        for (let j = 0; j < n; j++) {
          const qb = Q(j, b);
          if (qb) s += qa * A[i][j] * qb;
        }
      }
      B[a][b] = s;
    }
  return B;
}
