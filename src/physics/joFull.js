// Model B: the ten resistance scalars of two equal spheres (Jeffrey & Onishi 1984, lambda = 1),
// valid for every separation s = r/a > 2.
//
// Each function is written as  F(s) = S(s) + sum_m (c_m - sigma_m) s^-m,
// where c_m are the far-field series coefficients (from the recurrences, precomputed in
// src/data/joLambda1.js) and S(s) is a closed-form function that carries the near-contact
// singularities (1/xi, ln 1/xi, xi ln 1/xi; xi = s - 2). sigma_m are the series coefficients of S,
// computed here, so the identity holds exactly for s > 2 whatever S is; S only speeds up the
// convergence near contact. The leading singular coefficients are those of JO (1984), eqs. (3.19),
// (4.15)-(4.16), (5.5)-(5.6), (6.9)-(6.10), (7.9)-(7.10), at lambda = 1.
//
// Normalisation (JO 1.7): A11 = 6 pi mu a XA11, A12 = 6 pi mu a XA12 (equal spheres), B = 4 pi mu a^2 Y^B,
// C = 8 pi mu a^3 X^C, Y^C. Forces and torques are those the spheres exert on the fluid.
import { F, K as KMAX } from "../data/joLambda1.js";

const lnOm = (xi) => Math.log((xi * (4 + xi)) / (2 + xi) ** 2); // ln(1 - 4/s^2)
const lnRatio = (xi) => Math.log((4 + xi) / xi);                  // ln((s+2)/(s-2))
const om = (xi) => (xi * (4 + xi)) / (2 + xi) ** 2;               // 1 - 4/s^2

// series coefficients (powers of 1/s) of the basic singular functions
const cE2 = (m) => (m >= 2 && m % 2 === 0 ? 2 ** (m + 1) / m : 0); // -ln(1 - 4/s^2)
const cO2 = (m) => (m >= 1 && m % 2 === 1 ? 2 ** (m + 1) / m : 0); // ln((s+2)/(s-2))
const BASIS = {
  E1: { value: (xi) => 1 / om(xi), coeff: (m) => (m % 2 === 0 ? 2 ** m : 0) },                      // (1 - 4/s^2)^-1
  E2: { value: (xi) => -lnOm(xi), coeff: cE2 },
  E3: { value: (xi) => -om(xi) * lnOm(xi), coeff: (m) => cE2(m) - 4 * cE2(m - 2) },                 // -(1 - 4/s^2) ln(1 - 4/s^2)
  O1: { value: (xi) => 2 / ((2 + xi) * om(xi)), coeff: (m) => (m % 2 === 1 ? 2 ** m : 0) },         // 2 s^-1 (1 - 4/s^2)^-1
  O2: { value: lnRatio, coeff: cO2 },
  O3: { value: (xi) => om(xi) * lnRatio(xi), coeff: (m) => cO2(m) - 4 * cO2(m - 2) },               // (1 - 4/s^2) ln((s+2)/(s-2))
  C1: { value: (xi) => 0.5 * lnOm(xi) + lnRatio(xi) / (2 + xi), coeff: (m) => -0.5 * cE2(m) + cO2(m - 1) }, // JO (6.12)
  C2: { value: (xi) => lnRatio(xi) + (2 / (2 + xi)) * lnOm(xi), coeff: (m) => cO2(m) - 2 * cE2(m - 1) },     // JO (6.13)
};

// [family, parity of the powers, sign of the series, singular part [basis, coefficient]]
const even = 0, odd = 1;
export const FUNCTIONS = {
  XA11: ["XA", even, 1, [["E1", 1 / 4], ["E2", 9 / 40], ["E3", 3 / 112]]],
  XA12: ["XA", odd, -1, [["O1", -1 / 4], ["O2", -9 / 40], ["O3", -3 / 112]]],
  YA11: ["YA", even, 1, [["E2", 1 / 6]]],
  YA12: ["YA", odd, -1, [["O2", -1 / 6]]],
  YB11: ["YB", odd, 1, [["O2", -1 / 4], ["O3", -1 / 8]]],
  YB12: ["YB", even, -1, [["E2", 1 / 4], ["E3", 1 / 8]]],
  XC11: ["XC", even, 1, [["C1", 1 / 2]]],
  XC12: ["XC", odd, -1, [["C2", 1 / 4]]],
  YC11: ["YC", even, 1, [["E2", 1 / 5], ["E3", 47 / 250]]],
  YC12: ["YC", odd, 1, [["O2", 1 / 20], ["O3", 31 / 250]]],
};

// remainder coefficients (c_m - sigma_m), cached per function and order
const cache = new Map();
function remainder(name, K) {
  const key = `${name}:${K}`;
  if (cache.has(key)) return cache.get(key);
  const [fam, parity, sign, sing] = FUNCTIONS[name];
  const r = [];
  for (let m = 0; m <= K; m++) {
    if (m % 2 !== parity) { r.push(0); continue; }
    const c = (sign * F[fam][m]) / 2 ** m; // f_m (1 + lambda)^-m with lambda = 1
    r.push(c - sing.reduce((acc, [b, g]) => acc + g * BASIS[b].coeff(m), 0));
  }
  cache.set(key, r);
  return r;
}

/**
 * One scalar at the gap xi = s - 2 = (r - 2a)/a > 0. Passing the gap keeps 1/xi accurate near
 * contact (2 + xi - 2 would lose digits). { plain: true } sums the far-field series alone.
 */
export function joScalarAtGap(name, xi, { K = KMAX, plain = false } = {}) {
  if (!FUNCTIONS[name]) throw new RangeError(`unknown function ${name}`);
  if (!(xi > 0)) throw new RangeError("a positive gap is required (the spheres must not touch)");
  if (!(K <= KMAX)) throw new RangeError(`at most ${KMAX} terms are available`);
  const s = 2 + xi;
  const [fam, parity, sign, sing] = FUNCTIONS[name];
  let v = 0, p = 1;
  if (plain) {
    for (let m = 0; m <= K; m++, p /= s) if (m % 2 === parity) v += ((sign * F[fam][m]) / 2 ** m) * p;
    return v;
  }
  v = sing.reduce((acc, [b, g]) => acc + g * BASIS[b].value(xi), 0);
  const r = remainder(name, K);
  for (let m = 0; m <= K; m++, p /= s) v += r[m] * p;
  return v;
}

// The same at s = r/a.
export const joScalar = (name, s, opts) => joScalarAtGap(name, s - 2, opts);

export function joScalarsAtGap(xi, opts) {
  return Object.fromEntries(Object.keys(FUNCTIONS).map((n) => [n, joScalarAtGap(n, xi, opts)]));
}
export const joScalars = (s, opts) => joScalarsAtGap(s - 2, opts);

// Leading near-contact behaviour (lambda = 1), for display next to the full functions.
export const nearContact = {
  XA11: (xi) => 1 / (4 * xi) + (9 / 40) * Math.log(1 / xi),
  XA12: (xi) => -1 / (4 * xi) - (9 / 40) * Math.log(1 / xi),
  YA11: (xi) => Math.log(1 / xi) / 6,
  YA12: (xi) => -Math.log(1 / xi) / 6,
  YB11: (xi) => -Math.log(1 / xi) / 4,
  YB12: (xi) => Math.log(1 / xi) / 4,
  YC11: (xi) => Math.log(1 / xi) / 5,
  YC12: (xi) => Math.log(1 / xi) / 20,
};
