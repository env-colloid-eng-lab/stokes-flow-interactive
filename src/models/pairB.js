// Two equal spheres with translation and rotation: the 12 x 12 resistance matrix
//   (F1, F2, L1, L2) = R (U1, U2, Omega1, Omega2)
// (force and torque that the spheres exert on the fluid), in the block layout of Jeffrey & Onishi
// (1984, eq. 1.3). Each 3 x 3 block is fixed by the line of centres e = (x2 - x1)/r and scalars:
//   A_ab = 6 pi mu a [X^A_ab e e + Y^A_ab (I - e e)]
//   B_ab = 4 pi mu a^2 Y^B_ab (eps . e)          (torque on a from the translation of b)
//   C_ab = 8 pi mu a^3 [X^C_ab e e + Y^C_ab (I - e e)]
// with the equal-sphere relations A22 = A11, B22 = -B11, B21 = -B12, C22 = C11 (JO 1.9).
// Model A inverts the RPY mobility with rotation (outside the textbook; tested separately);
// model B uses the Jeffrey-Onishi functions of src/physics/joFull.js.
import { zeros, inv, getBlock, setBlock, sub, matmul, transpose, solve } from "../core/linalg.js";
import { rpyPair, selfMobility } from "../physics/rpy.js";
import { joScalarsAtGap } from "../physics/joFull.js";

export const NAMES = ["XA11", "XA12", "YA11", "YA12", "YB11", "YB12", "XC11", "XC12", "YC11", "YC12"];

// eps_ijk e_k: (eps . e) v = v x e ... as a matrix M_ij = eps_ijk e_k
const epsE = (e) => [[0, e[2], -e[1]], [-e[2], 0, e[0]], [e[1], -e[0], 0]];
const axisym = (X, Y, e) => [0, 1, 2].map((i) => [0, 1, 2].map((j) => X * e[i] * e[j] + Y * ((i === j ? 1 : 0) - e[i] * e[j])));
const scaled = (M, c) => M.map((row) => row.map((v) => v * c));

/** 12 x 12 resistance matrix from dimensionless scalars (JO normalisation). */
export function resistance12(sc, e, { a = 1, mu = 1 } = {}) {
  const cA = 6 * Math.PI * mu * a, cB = 4 * Math.PI * mu * a * a, cC = 8 * Math.PI * mu * a ** 3;
  const A11 = scaled(axisym(sc.XA11, sc.YA11, e), cA), A12 = scaled(axisym(sc.XA12, sc.YA12, e), cA);
  const B11 = scaled(epsE(e), cB * sc.YB11), B12 = scaled(epsE(e), cB * sc.YB12);
  const B21 = scaled(B12, -1), B22 = scaled(B11, -1);
  const C11 = scaled(axisym(sc.XC11, sc.YC11, e), cC), C12 = scaled(axisym(sc.XC12, sc.YC12, e), cC);
  const R = zeros(12);
  setBlock(R, 0, 0, A11); setBlock(R, 0, 3, A12); setBlock(R, 3, 0, A12); setBlock(R, 3, 3, A11);
  setBlock(R, 6, 6, C11); setBlock(R, 6, 9, C12); setBlock(R, 9, 6, C12); setBlock(R, 9, 9, C11);
  // lower-left: torque rows, velocity columns; upper-right is its transpose (reciprocal theorem)
  [[6, 0, B11], [6, 3, B12], [9, 0, B21], [9, 3, B22]].forEach(([r0, c0, B]) => {
    setBlock(R, r0, c0, B);
    setBlock(R, c0, r0, transpose(B));
  });
  return R;
}

/** Scalars back from a 12 x 12 matrix (works for any model), with e along z. */
export function scalarsOf(R, { a = 1, mu = 1 } = {}) {
  const cA = 6 * Math.PI * mu * a, cB = 4 * Math.PI * mu * a * a, cC = 8 * Math.PI * mu * a ** 3;
  return {
    XA11: R[2][2] / cA, XA12: R[2][5] / cA, YA11: R[0][0] / cA, YA12: R[0][3] / cA,
    YB11: R[6][1] / cB, YB12: R[6][4] / cB, // B_xy = Y^B eps_xyz e_z with e = z
    XC11: R[8][8] / cC, XC12: R[8][11] / cC, YC11: R[6][6] / cC, YC12: R[6][9] / cC,
  };
}

/**
 * RPY mobility with rotation for two equal spheres, r >= 2a, same layout (U1, U2, Omega1, Omega2).
 * Translation-translation: rpyPair. Rotation-rotation: I/(8 pi mu a^3) and (3ee - I)/(16 pi mu r^3).
 * Rotation from a force on the other sphere: Omega_b = F_a x (x_b - x_a) / (8 pi mu r^3).
 */
export function rpyMobility12(rvec, { a = 1, mu = 1 } = {}) {
  const r = Math.hypot(...rvec);
  if (!(r >= 2 * a)) throw new RangeError("non-overlapping spheres required (r >= 2a)");
  const e = rvec.map((x) => x / r);
  const M = zeros(12);
  const self = selfMobility(a, mu), cross = rpyPair(rvec, { a, mu });
  const eye3 = (c) => [0, 1, 2].map((i) => [0, 1, 2].map((j) => (i === j ? c : 0)));
  setBlock(M, 0, 0, eye3(self)); setBlock(M, 3, 3, eye3(self)); setBlock(M, 0, 3, cross); setBlock(M, 3, 0, cross);
  const rr = scaled(axisym(2, -1, e), 1 / (16 * Math.PI * mu * r ** 3)); // (3ee - I) = 2 ee - (I - ee)
  setBlock(M, 6, 6, eye3(1 / (8 * Math.PI * mu * a ** 3))); setBlock(M, 9, 9, eye3(1 / (8 * Math.PI * mu * a ** 3)));
  setBlock(M, 6, 9, rr); setBlock(M, 9, 6, rr);
  // Omega_2 = F_1 x e / (8 pi mu r^2): (F x e)_i = eps_ijk F_j e_k  ->  matrix eps_ijk e_k
  const w21 = scaled(epsE(e), 1 / (8 * Math.PI * mu * r * r)), w12 = scaled(w21, -1);
  setBlock(M, 9, 0, w21); setBlock(M, 6, 3, w12);
  setBlock(M, 0, 9, transpose(w21)); setBlock(M, 3, 6, transpose(w12));
  return M;
}

/** Model A and model B scalars at the gap h = r - 2a. */
export function pairScalarsAtGap(h, model, { a = 1, mu = 1 } = {}) {
  if (model !== "A" && model !== "B") throw new RangeError(`model must be "A" or "B", got ${model}`);
  if (model === "B") return joScalarsAtGap(h / a);
  return scalarsOf(inv(rpyMobility12([0, 0, 2 * a + h], { a, mu })), { a, mu });
}

/** Torque-free translational resistance: eliminate the rotations (Schur complement A - B~ C^-1 B). */
export function torqueFree(R) {
  const A = getBlock(R, 0, 0, 6), Bt = getBlock(R, 0, 6, 6, 6), B = getBlock(R, 6, 0, 6, 6), C = getBlock(R, 6, 6, 6);
  return sub(A, matmul(Bt, solve(C, B)));
}

/**
 * Model B in the form used by pages 8 and 9: torque-free pair coefficients along (par) and across
 * (perp) the line of centres, self and cross, for the resistance R and its inverse M (dimensional).
 */
export function coefficientsBAtGap(h, { a = 1, mu = 1 } = {}) {
  if (!(h > 0)) throw new RangeError("the spheres must not touch (gap > 0)");
  const T = torqueFree(resistance12(joScalarsAtGap(h / a), [0, 0, 1], { a, mu }));
  const R = { par: { self: T[2][2], cross: T[2][5] }, perp: { self: T[0][0], cross: T[0][3] } };
  // self ~ -cross near contact: factor the determinant instead of forming self^2 - cross^2
  const invert2 = ({ self, cross }) => { const d = (self - cross) * (self + cross); return { self: self / d, cross: -cross / d }; };
  return { R, M: { par: invert2(R.par), perp: invert2(R.perp) } };
}
// The same at the centre distance r (pass the gap when it is tiny, to keep its digits).
export const coefficientsB = (r, opts = {}) => coefficientsBAtGap(r - 2 * (opts.a ?? 1), opts);
