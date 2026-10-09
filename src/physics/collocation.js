// N coaxial spheres translating along the axis: boundary collocation with
// axisymmetric Stokes streamfunction modes. Port of TwoSphereLab.jl (chapter 12).
// Non-axisymmetric motion is excluded.
import { zeros, lstsq, cond } from "../core/linalg.js";

export function legendreAndDerivative(n, t) {
  if (n === 0) return [1, 0];
  let p0 = 1, p1 = t, d0 = 0, d1 = 1;
  for (let k = 1; k < n; k++) {
    const p2 = ((2 * k + 1) * t * p1 - k * p0) / (k + 1);
    const d2 = ((2 * k + 1) * (p1 + t * d1) - k * d0) / (k + 1);
    p0 = p1; p1 = p2; d0 = d1; d1 = d2;
  }
  return [p1, d1];
}

// Exterior mode around a sphere of radius a; k = 0 pressure mode, k = 1 potential mode.
// psi = c a^2 (r/a)^q (1 - t^2) P_n'(t). Returns [u_rho, u_z].
export function axialMode(rho, z, a, n, k) {
  const R = Math.hypot(rho, z), t = z / R, s = rho / R;
  const [P, D] = legendreAndDerivative(n, t);
  const q = k === 0 ? 2 - n : -n;
  const f = a * a * (R / a) ** q;
  const ur = (n * (n + 1) * f) / (R * R) * P;
  const ut = (-q * f) / (R * R) * s * D;
  return [ur * s + ut * t, ur * t - ut * s];
}

function axialRows(rho, z, centers, radii, L) {
  const rowRho = [], rowZ = [];
  for (let b = 0; b < radii.length; b++)
    for (let n = 1; n <= L; n++)
      for (let k = 0; k <= 1; k++) {
        const [ur, uz] = axialMode(rho, z - centers[b], radii[b], n, k);
        rowRho.push(ur); rowZ.push(uz);
      }
  return [rowRho, rowZ];
}

// Solves all unit-velocity problems at once without imposing symmetry on R.
// Returns R (N x N, f = R U along the axis), the coefficients, an independent
// boundary residual and the condition number of the scaled system.
export function axialCollocation(centers, radii, { L = 16, nc = 3 * L + 4, mu = 1, ncheck = 121, withCondition = true } = {}) {
  const N = radii.length;
  if (!(centers.length === N && N >= 1 && L >= 1 && nc >= L + 1 && mu > 0)) throw new RangeError("invalid dimensions/orders");
  if (!radii.every((a) => a > 0)) throw new RangeError("positive radii required");
  for (let i = 0; i < N; i++)
    for (let j = i + 1; j < N; j++)
      if (!(Math.abs(centers[i] - centers[j]) > radii[i] + radii[j])) throw new RangeError("spheres overlap or touch");

  const A = [], rhs = [];
  for (let al = 0; al < N; al++) {
    for (let j = 1; j <= nc; j++) {
      const t = Math.cos((Math.PI * (j - 0.5)) / nc);
      const rho = radii[al] * Math.sqrt(1 - t * t), z = centers[al] + radii[al] * t;
      const [rr, rz] = axialRows(rho, z, centers, radii, L);
      A.push(rr, rz);
      const bRho = new Array(N).fill(0), bZ = new Array(N).fill(0);
      bZ[al] = 1;
      rhs.push(bRho, bZ);
    }
  }
  const ncol = A[0].length;
  const sc = Array.from({ length: ncol }, (_, c) => Math.sqrt(A.reduce((s, row) => s + row[c] * row[c], 0)));
  const As = A.map((row) => row.map((v, c) => v / sc[c]));
  const C = lstsq(As, rhs).map((row, c) => row.map((v) => v / sc[c]));

  const R = zeros(N);
  for (let al = 0; al < N; al++)
    for (let b = 0; b < N; b++) R[al][b] = 8 * Math.PI * mu * radii[al] * C[2 * L * al][b];

  let err = 0;
  for (let al = 0; al < N; al++) {
    for (let i = 0; i < ncheck; i++) {
      const t = -0.999999 + (1.999998 * i) / (ncheck - 1);
      const rho = radii[al] * Math.sqrt(1 - t * t), z = centers[al] + radii[al] * t;
      const [rr, rz] = axialRows(rho, z, centers, radii, L);
      for (let b = 0; b < N; b++) {
        let ur = 0, uz = 0;
        for (let c = 0; c < ncol; c++) { ur += rr[c] * C[c][b]; uz += rz[c] * C[c][b]; }
        err = Math.max(err, Math.abs(ur), Math.abs(uz - (b === al ? 1 : 0)));
      }
    }
  }
  return { R, C, boundaryError: err, condition: withCondition ? cond(As) : NaN, L, nc, centers, radii };
}
