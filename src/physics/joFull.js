// Model B: the resistance scalars of two spheres (Jeffrey & Onishi 1984) for the size ratios
// lambda = a2/a1 in LAMBDAS (1, 1/2, 2, 1/4, 4), valid for every separation s = 2r/(a1 + a2) > 2.
//
// Each function is written as  F(s) = S(s) + sum_m (c_m - sigma_m) s^-m,
// where c_m are the far-field series coefficients (from the recurrences, precomputed in
// src/data/joLambda1.js) and S(s) is a closed-form function that carries the near-contact
// singularities (1/xi, ln 1/xi, xi ln 1/xi; xi = s - 2). sigma_m are the series coefficients of S,
// computed here, so the identity holds exactly for s > 2 whatever S is; S only speeds up the
// convergence near contact. The singular coefficients are those of JO (1984), eqs. (3.19),
// (4.15)-(4.16), (5.5)-(5.6), (6.9)-(6.10), (7.9)-(7.10), with two corrections for Y^C_12:
//  - g5 = (2/125) lambda (43 - 24 lambda + 43 lambda^2)(1 + lambda)^-4 (Townsend, arXiv:1802.08226,
//    eq. 30; Phys. Fluids 35, 127126, 2023); JO print twice this. His corrected recurrence for V
//    (eq. 18) is the one used in joRecurrence.js.
//  - the singular part multiplies Y^C_12 itself, Y^C_12 = g4 ln(1/xi) + C + g5 xi ln(1/xi), as in JO (7.10).
//    The factor (1 + lambda)^3/8 of JO (7.15), and 8/(1 + lambda)^3 of Townsend (34), would leave a
//    ln(1/xi) singularity in the remainder for lambda != 1 (both equal 1 at lambda = 1).
// Every choice is checked by the size of the series tail (test/joFull.test.js) and by the contact
// constants tabulated by JO for lambda = 1, 1/2, 1/4.
//
// Normalisation (JO 1.7): A_ab = 3 pi mu (a_a + a_b) X^A_ab, B_ab = pi mu (a_a + a_b)^2 Y^B_ab,
// C_ab = pi mu (a_a + a_b)^3 X^C_ab, Y^C_ab. Forces and torques are those the spheres exert on the fluid.
// Functions "22" and "21" follow from lambda -> 1/lambda (JO 1.9), see joPairScalarsAtGap.
import { F as DATA, K as KMAX, LAMBDAS } from "../data/joValues.js";

export { LAMBDAS };

// 1 - 4/s^2 and its logarithm: the product form keeps digits near contact, log1p far away
const om = (xi) => (xi < 1 ? (xi * (4 + xi)) / (2 + xi) ** 2 : 1 - 4 / (2 + xi) ** 2);
const lnOm = (xi) => (xi < 1 ? Math.log(xi) + Math.log1p(xi / 4) - 2 * Math.log1p(xi / 2) : Math.log1p(-4 / (2 + xi) ** 2));
const lnRatio = (xi) => Math.log1p(4 / xi);                        // ln((s+2)/(s-2))

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

// For each scalar: [family, parity of the powers, prefactor of the series, singular part [[basis, g]]].
// The series is prefactor * sum_m f_m(lambda) (1 + lambda)^-m s^-m over the given parity.
const even = 0, odd = 1;
export const NAMES10 = ["XA11", "XA12", "YA11", "YA12", "YB11", "YB12", "XC11", "XC12", "YC11", "YC12"];
const defsCache = new Map();
export function functionsFor(lambda = 1) {
  if (defsCache.has(lambda)) return defsCache.get(lambda);
  const l = lambda, L = 1 + l;
  const A = { g1: (2 * l * l) / L ** 3, g2: (l * (1 + 7 * l + l * l)) / (5 * L ** 3), g3: (1 + 18 * l - 29 * l * l + 18 * l ** 3 + l ** 4) / (42 * L ** 3) };
  const YA = { g2: (4 * l * (2 + l + 2 * l * l)) / (15 * L ** 3), g3: (2 * (16 - 45 * l + 58 * l * l - 45 * l ** 3 + 16 * l ** 4)) / (375 * L ** 3) };
  const YB = { g2: -(l * (4 + l)) / (5 * L * L), g3: -(32 - 33 * l + 83 * l * l + 43 * l ** 3) / (250 * L * L) };
  const YC = { g2: (2 * l) / (5 * L), g3: (8 + 6 * l + 33 * l * l) / (125 * L), g4: (4 * l * l) / (5 * L ** 4), g5: (2 * l * (43 - 24 * l + 43 * l * l)) / (125 * L ** 4) };
  const k12 = -2 / L, kB = -4 / L ** 2;
  const defs = {
    XA11: ["XA", even, 1, [["E1", A.g1], ["E2", A.g2], ["E3", A.g3]]],
    XA12: ["XA", odd, k12, [["O1", k12 * A.g1], ["O2", k12 * A.g2], ["O3", k12 * A.g3]]],
    YA11: ["YA", even, 1, [["E2", YA.g2], ["E3", YA.g3]]],
    YA12: ["YA", odd, k12, [["O2", k12 * YA.g2], ["O3", k12 * YA.g3]]],
    YB11: ["YB", odd, 1, [["O2", YB.g2], ["O3", YB.g3]]],
    YB12: ["YB", even, kB, [["E2", kB * YB.g2], ["E3", kB * YB.g3]]],
    XC11: ["XC", even, 1, [["C1", (l * l) / L]]],
    XC12: ["XC", odd, -8 / L ** 3, [["C2", (4 * l * l) / L ** 4]]],
    YC11: ["YC", even, 1, [["E2", YC.g2], ["E3", YC.g3]]],
    YC12: ["YC", odd, 8 / L ** 3, [["O2", YC.g4], ["O3", YC.g5]]],
  };
  // leading behaviour near contact (the 1/xi and ln(1/xi) parts), for display
  const lead = Object.fromEntries(Object.entries(defs).map(([n, [, , , sing]]) => {
    const c1 = sing.find(([b]) => b === "E1" || b === "O1")?.[1] ?? 0, c2 = sing.find(([b]) => b === "E2" || b === "O2")?.[1] ?? 0;
    // E1 ~ 1/xi and O1 ~ 1/xi near contact; E2 ~ ln(1/xi) and O2 ~ ln(1/xi) + ln 4 (the constant is not part of the leading term)
    return [n, (xi) => c1 / xi + c2 * Math.log(1 / xi)];
  }));
  const out = { defs, lead };
  defsCache.set(lambda, out);
  return out;
}
export const FUNCTIONS = functionsFor(1).defs;

function dataFor(lambda) {
  const d = DATA[String(lambda)];
  if (!d) throw new RangeError(`size ratio ${lambda} is not available (use one of ${LAMBDAS.join(", ")})`);
  return d;
}

// remainder coefficients (c_m - sigma_m), cached per function, ratio and order
const cache = new Map();
export function remainder(name, K = KMAX, lambda = 1) {
  const key = `${name}:${K}:${lambda}`;
  if (cache.has(key)) return cache.get(key);
  const [fam, parity, pre, sing] = functionsFor(lambda).defs[name];
  const f = dataFor(lambda)[fam], L = 1 + lambda;
  const r = [];
  for (let m = 0; m <= K; m++) {
    if (m % 2 !== parity) { r.push(0); continue; }
    const c = pre * f[m] * L ** -m;
    r.push(c - sing.reduce((acc, [b, g]) => acc + g * BASIS[b].coeff(m), 0));
  }
  cache.set(key, r);
  return r;
}

/**
 * One scalar at the gap xi = s - 2 = (r - 2a)/a > 0. Passing the gap keeps 1/xi accurate near
 * contact (2 + xi - 2 would lose digits). { plain: true } sums the far-field series alone.
 */
export function joScalarAtGap(name, xi, { K = KMAX, plain = false, lambda = 1 } = {}) {
  if (!NAMES10.includes(name)) throw new RangeError(`unknown function ${name}`);
  if (!(xi > 0)) throw new RangeError("a positive gap is required (the spheres must not touch)");
  if (!(Number.isInteger(K) && K >= 0 && K <= KMAX)) throw new RangeError(`K must be an integer in [0, ${KMAX}]`);
  const s = 2 + xi;
  const [fam, parity, pre, sing] = functionsFor(lambda).defs[name];
  let v = 0, p = 1;
  if (plain) {
    const f = dataFor(lambda)[fam], L = 1 + lambda;
    for (let m = 0; m <= K; m++, p /= s) if (m % 2 === parity) v += pre * f[m] * L ** -m * p;
    return v;
  }
  v = sing.reduce((acc, [b, g]) => acc + g * BASIS[b].value(xi), 0);
  const r = remainder(name, K, lambda);
  for (let m = 0; m <= K; m++, p /= s) v += r[m] * p;
  return v;
}

// The same at s = r/a.
export const joScalar = (name, s, opts) => joScalarAtGap(name, s - 2, opts);

export function joScalarsAtGap(xi, opts) {
  return Object.fromEntries(NAMES10.map((n) => [n, joScalarAtGap(n, xi, opts)]));
}

/**
 * All scalars of an unequal pair at lambda = a2/a1: the ten above plus, from JO (1.9),
 * X_22(lambda) = X_11(1/lambda), Y^B_22(lambda) = -Y^B_11(1/lambda), Y^B_21(lambda) = -Y^B_12(1/lambda).
 * (X^A_21 = X^A_12 etc. for the symmetric families.)
 */
export function joPairScalarsAtGap(xi, { lambda = 1, K } = {}) {
  const one = joScalarsAtGap(xi, { lambda, K });
  const inv = lambda === 1 ? one : joScalarsAtGap(xi, { lambda: 1 / lambda, K });
  return {
    ...one,
    XA22: inv.XA11, YA22: inv.YA11, XC22: inv.XC11, YC22: inv.YC11,
    YB22: -inv.YB11, YB21: -inv.YB12,
  };
}
export const joScalars = (s, opts) => joScalarsAtGap(s - 2, opts);

// Leading near-contact behaviour (1/xi and ln(1/xi) parts), for display next to the full functions.
export const nearContactAt = (name, xi, lambda = 1) => functionsFor(lambda).lead[name](xi);
export const nearContact = Object.fromEntries(NAMES10.filter((n) => n[0] !== "X" || n[1] !== "C")
  .map((n) => [n, (xi) => nearContactAt(n, xi, 1)]));
