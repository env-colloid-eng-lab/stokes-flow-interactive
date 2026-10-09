// Point-force (Oseen) tensor and its derivatives. Port of StokesLab.jl.
// Convention: f is the force the particle exerts on the fluid; u = G f.
import { eye, dot, norm, zeros } from "../core/linalg.js";

export const eye3 = () => eye(3);

// I - q q^T / |q|^2 : removes the component along q (Fourier-space pressure projection).
export function projector(q) {
  const qq = dot(q, q);
  if (!(qq > 0)) throw new RangeError("q must be nonzero");
  return eye(3).map((row, i) => row.map((d, j) => d - (q[i] * q[j]) / qq));
}

function check(r, mu) {
  const R = norm(r);
  if (!(R > 0 && mu > 0)) throw new RangeError("r != 0, mu > 0 required");
  return R;
}

// G_ij = (delta_ij / r + r_i r_j / r^3) / (8 pi mu)
export function oseen(r, { mu = 1 } = {}) {
  const R = check(r, mu);
  const c = 1 / (8 * Math.PI * mu);
  return eye(3).map((row, i) => row.map((d, j) => c * (d / R + (r[i] * r[j]) / R ** 3)));
}

// H[i][j][k] = partial_k G_ij (not partial_i G_jk).
export function gradOseen(r, { mu = 1 } = {}) {
  const R = check(r, mu);
  const c = 1 / (8 * Math.PI * mu);
  const d = (a, b) => (a === b ? 1 : 0);
  const H = [0, 1, 2].map(() => zeros(3));
  for (let i = 0; i < 3; i++)
    for (let j = 0; j < 3; j++)
      for (let k = 0; k < 3; k++)
        H[i][j][k] =
          c *
          ((-d(i, j) * r[k] + d(i, k) * r[j] + d(j, k) * r[i]) / R ** 3 -
            (3 * r[i] * r[j] * r[k]) / R ** 5);
  return H;
}

// nabla^2 G = (2 I / r^3 - 6 r r^T / r^5) / (8 pi mu)
export function lapOseen(r, { mu = 1 } = {}) {
  const R = check(r, mu);
  const c = 1 / (8 * Math.PI * mu);
  return eye(3).map((row, i) => row.map((dd, j) => c * ((2 * dd) / R ** 3 - (6 * r[i] * r[j]) / R ** 5)));
}

// Pressure of a point force: p = f . r / (4 pi r^3). r != 0.
export const pressure = (r, f) => dot(f, r) / (4 * Math.PI * norm(r) ** 3);

// Force dipole: +f at +d/2 and -f at -d/2, D = f d^T, u_i = -H_ijk D_jk.
export function dipole(r, D, { mu = 1 } = {}) {
  const H = gradOseen(r, { mu });
  return [0, 1, 2].map((i) => {
    let s = 0;
    for (let j = 0; j < 3; j++) for (let k = 0; k < 3; k++) s -= H[i][j][k] * D[j][k];
    return s;
  });
}

// Flow outside a sphere translating with U: (G + a^2/6 nabla^2 G)(6 pi mu a U).
export function sphereFlow(r, U, { a = 1, mu = 1 } = {}) {
  if (!(a > 0 && norm(r) >= a)) throw new RangeError("r must be outside sphere");
  const G = oseen(r, { mu }), L = lapOseen(r, { mu });
  const f = U.map((u) => 6 * Math.PI * mu * a * u);
  return [0, 1, 2].map((i) => [0, 1, 2].reduce((s, j) => s + (G[i][j] + (a * a / 6) * L[i][j]) * f[j], 0));
}
