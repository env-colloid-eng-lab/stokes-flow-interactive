// Single sphere in Stokes flow: exterior fields, surface quadrature and force moments.
// Port of FoundationLab.jl. Conventions: L_ij = d_j u_i, n points from the particle into the
// fluid, t = sigma n is the traction the fluid exerts on the particle, and S is the stresslet
// that enters the suspension stress with a plus sign.
import { eye, zeros, norm, dot, cross, symEig, matvec } from "../core/linalg.js";

const I3 = () => eye(3);
const outer = (u, v) => u.map((ui) => v.map((vj) => ui * vj));
const addM = (...Ms) => Ms[0].map((row, i) => row.map((_, j) => Ms.reduce((s, M) => s + M[i][j], 0)));
const scaleM = (M, c) => M.map((row) => row.map((v) => v * c));
const transposeM = (M) => M[0].map((_, j) => M.map((row) => row[j]));
const stressFrom = (p, L, mu) => addM(scaleM(I3(), -p), scaleM(addM(L, transposeM(L)), mu));

// Symmetric traceless part (A + A^T)/2 - tr(A)/3 I.
export function stf(A) {
  const tr = (A[0][0] + A[1][1] + A[2][2]) / 3;
  return A.map((row, i) => row.map((v, j) => (v + A[j][i]) / 2 - (i === j ? tr : 0)));
}

// Gauss-Legendre nodes and weights on [-1, 1] (Golub-Welsch).
export function gaussLegendre(n) {
  if (!(n >= 1)) throw new RangeError("n >= 1 required");
  const J = zeros(n);
  for (let k = 1; k < n; k++) J[k - 1][k] = J[k][k - 1] = k / Math.sqrt(4 * k * k - 1);
  const { values, vectors } = symEig(J);
  return { nodes: values, weights: values.map((_, j) => 2 * vectors[0][j] ** 2) };
}

// Gauss-Legendre nodes in t = cos(theta) times the trapezoid rule in azimuth.
// The weights integrate dOmega (they sum to 4 pi).
export function sphereRule(nt = 8, np = 16) {
  if (!(nt >= 2 && np >= 3)) throw new RangeError("increase quadrature order");
  const gl = gaussLegendre(nt);
  const normals = [], weights = [];
  gl.nodes.forEach((t, j) => {
    const w = gl.weights[j];
    const s = Math.sqrt(Math.max(0, 1 - t * t));
    for (let k = 0; k < np; k++) {
      const phi = (2 * Math.PI * k) / np;
      normals.push([s * Math.cos(phi), s * Math.sin(phi), t]);
      weights.push((w * 2 * Math.PI) / np);
    }
  });
  return { normals, weights };
}

function exterior(r, a) {
  const R = norm(r);
  if (!(R >= a * (1 - 1e-12) && a > 0)) throw new RangeError("exterior points required");
  return R;
}

// Sphere translating with U in fluid at rest.
export function translationField(r, U, { a = 1, mu = 1 } = {}) {
  const R = exterior(r, a), n = r.map((x) => x / R), q = dot(r, U);
  const f = (3 * a) / (4 * R) + a ** 3 / (4 * R ** 3), g = (3 * a) / (4 * R ** 3) - (3 * a ** 3) / (4 * R ** 5);
  const fp = (-3 * a) / (4 * R ** 2) - (3 * a ** 3) / (4 * R ** 4), gp = (-9 * a) / (4 * R ** 4) + (15 * a ** 3) / (4 * R ** 6);
  const u = U.map((Ui, i) => f * Ui + g * q * r[i]);
  const L = addM(scaleM(outer(U, n), fp), scaleM(outer(r, n), gp * q), scaleM(I3(), g * q), scaleM(outer(r, U), g));
  const p = (3 * mu * a * q) / (2 * R ** 3);
  return { u, p, L, sigma: stressFrom(p, L, mu) };
}

// Force-free, torque-free sphere held fixed in a pure straining flow E r (E symmetric traceless).
export function strainField(r, E, { a = 1, mu = 1 } = {}) {
  const R = exterior(r, a), n = r.map((x) => x / R), v = matvec(E, r), q = dot(r, v);
  const f = 1 - a ** 5 / R ** 5, g = (-5 * a ** 3) / (2 * R ** 5) + (5 * a ** 5) / (2 * R ** 7);
  const fp = (5 * a ** 5) / R ** 6, gp = (25 * a ** 3) / (2 * R ** 6) - (35 * a ** 5) / (2 * R ** 8);
  const u = v.map((vi, i) => f * vi + g * q * r[i]);
  const L = addM(scaleM(outer(v, n), fp), scaleM(E, f), scaleM(outer(r, n), gp * q), scaleM(addM(scaleM(outer(r, v), 2), scaleM(I3(), q)), g));
  const p = (-5 * mu * a ** 3 * q) / R ** 5;
  return { u, p, L, sigma: stressFrom(p, L, mu) };
}

// Sphere rotating with Omega in fluid at rest.
export function rotationField(r, Omega, { a = 1, mu = 1 } = {}) {
  const R = exterior(r, a), v = cross(Omega, r);
  const C = [[0, -Omega[2], Omega[1]], [Omega[2], 0, -Omega[0]], [-Omega[1], Omega[0], 0]];
  const u = v.map((x) => (a ** 3 * x) / R ** 3);
  const L = addM(scaleM(C, a ** 3 / R ** 3), scaleM(outer(v, r), (-3 * a ** 3) / R ** 5));
  return { u, p: 0, L, sigma: stressFrom(0, L, mu) };
}

// Force, torque and stresslet from the traction on the sphere surface, with the pressure and
// viscous parts of the force kept separately.
export function surfaceMoments(field, { a = 1, nt = 8, np = 16, mu = 1 } = {}) {
  const { normals, weights } = sphereRule(nt, np);
  const F = [0, 0, 0], T = [0, 0, 0], Fp = [0, 0, 0], Fv = [0, 0, 0];
  let S = zeros(3);
  normals.forEach((n, k) => {
    const r = n.map((x) => a * x), fl = field(r), dA = a * a * weights[k];
    const t = matvec(fl.sigma, n), torque = cross(r, t);
    const tv = matvec(addM(fl.L, transposeM(fl.L)), n).map((x) => mu * x);
    for (let i = 0; i < 3; i++) {
      F[i] += dA * t[i]; T[i] += dA * torque[i];
      Fp[i] += dA * -fl.p * n[i]; Fv[i] += dA * tv[i];
    }
    S = addM(S, scaleM(stf(outer(r, t)), dA));
  });
  return { F, T, S, Fpressure: Fp, Fviscous: Fv };
}

// Faxen relations (force and torque on the fluid side: F^H, T^H; stresslet S).
export const faxenForce = (U, uamb, lapu, { a = 1, mu = 1 } = {}) => U.map((Ui, i) => -6 * Math.PI * mu * a * (Ui - uamb[i] - ((a * a) / 6) * lapu[i]));
export const faxenTorque = (Omega, curlu, { a = 1, mu = 1 } = {}) => Omega.map((Oi, i) => -8 * Math.PI * mu * a ** 3 * (Oi - curlu[i] / 2));
export const faxenStresslet = (E, lapE, { a = 1, mu = 1 } = {}) => E.map((row, i) => row.map((v, j) => ((20 * Math.PI) / 3) * mu * a ** 3 * (v + ((a * a) / 10) * lapE[i][j])));

// Average of a vector field over the sphere surface |r - centre| = a.
export function surfaceAverage(fn, centre, { a = 1, nt = 10, np = 20 } = {}) {
  const { normals, weights } = sphereRule(nt, np);
  let acc = null; // the field need not be defined at the centre (e.g. a point force there)
  normals.forEach((n, k) => {
    const v = fn(centre.map((c, i) => c + a * n[i]));
    acc ??= v.map(() => 0);
    v.forEach((x, i) => (acc[i] += (weights[k] * x) / (4 * Math.PI)));
  });
  return acc;
}

// Viscous dissipation per unit volume, 2 mu E:E with E = (L + L^T)/2. The spin W does not dissipate.
export function dissipationDensity(L, mu = 1) {
  let s = 0;
  for (let i = 0; i < L.length; i++) for (let j = 0; j < L.length; j++) s += ((L[i][j] + L[j][i]) / 2) ** 2;
  return 2 * mu * s;
}

// Dissipation in the fluid between the sphere surface r = a and r = R (R = Infinity for all of it).
// With s = a/r the radial integral runs over s in [a/R, 1] and dV = a^3 s^-4 ds dOmega.
export function dissipationIntegral(field, { a = 1, mu = 1, R = Infinity, nr = 12, nt = 8, np = 16 } = {}) {
  const s0 = Math.min(1, a / R);
  const gl = gaussLegendre(nr), { normals, weights } = sphereRule(nt, np);
  let total = 0;
  gl.nodes.forEach((x, m) => {
    const sv = s0 + ((1 - s0) * (x + 1)) / 2, ws = ((1 - s0) / 2) * gl.weights[m];
    const r = a / sv, jac = (a ** 3 / sv ** 4) * ws;
    normals.forEach((n, k) => (total += jac * weights[k] * dissipationDensity(field(n.map((v) => r * v)).L, mu)));
  });
  return total;
}
