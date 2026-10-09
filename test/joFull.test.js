// Model B: Jeffrey & Onishi (1984) resistance functions for two equal spheres, with rotation.
// Reference values: the polynomials and the contact constants printed in JO (1984) at lambda = 1
// (Tables 2, 4, 5, 6), exact zeta(3) values for X^C at contact, and the collocation solver.
import { test } from "node:test";
import assert from "node:assert/strict";
import { Q } from "../src/core/rational.js";
import { inv, eigvalsSym, asymmetry, matvec } from "../src/core/linalg.js";
import { joAllPolynomials, joAllValues } from "../src/physics/joRecurrence.js";
import { joPolynomials } from "../src/physics/jo.js";
import { F as DATA } from "../src/data/joLambda1.js";
import { joScalar, joScalarAtGap, joScalars, nearContact, remainder, FUNCTIONS } from "../src/physics/joFull.js";
import { resistance12, scalarsOf, rpyMobility12, pairScalarsAtGap, torqueFree, NAMES } from "../src/models/pairB.js";
import { mobility } from "../src/physics/rpy.js";
import { axialCollocation } from "../src/physics/collocation.js";
import { readCsv, close } from "./helpers.js";

// coefficients printed in JO (1984) for k = 0..8, increasing powers of lambda
const P = (...xs) => xs.map((x) => (Array.isArray(x) ? Q(x[0], x[1]) : Q(x)));
const PRINTED = {
  YA: [P(1), P(0, [3, 2]), P(0, [9, 4]), P(0, 2, [27, 8], 2), P(0, 6, [81, 16], 18), P(0, 0, [63, 2], [243, 32], [63, 2]),
    P(0, 4, 54, [1241, 64], 81, 72), P(0, 0, 144, [1053, 8], [19083, 128], [1053, 8], 144),
    P(0, 0, 279, [4261, 8], [126369, 256], [-117, 8], 648, 288)],
  YB: [P(0), P(0, 0), P(0, -6, 0), P(0, -9, 0, 0), P(0, 0, [-27, 2], 0, 0), P(0, -12, [-81, 4], -36, 0, 0),
    P(0, 0, -108, [-243, 8], -72, 0, 0), P(0, 0, -189, [-8409, 16], -243, -144, 0, 0),
    P(0, 0, -432, -486, [-77451, 32], -405, -288, 0, 0)],
  XC: [P(1), P(0, 0), P(0, 0, 0), P(0, 0, 0, 8), P(0, 0, 0, 0, 0), P(0, 0, 0, 0, 0, 0), P(0, 0, 0, 64, 0, 0, 0),
    P(0, 0, 0, 0, 0, 0, 0, 0), P(0, 0, 0, 0, 0, 768, 0, 0, 0)],
  YC: [P(1), P(0, 0), P(0, 0, 0), P(0, 0, 0, 4), P(0, 12, 0, 0, 0), P(0, 0, 0, 0, 18, 0), P(0, 0, 27, 256, 0, 0, 0),
    P(0, 0, 0, 0, 72, [81, 2], 72, 0), P(0, 0, 216, [243, 4], 216, 2496, 0, 0, 0)],
};

test("the recurrences reproduce the polynomials printed in JO (1984)", () => {
  const ex = joAllPolynomials(8, { exact: true });
  for (const [fam, rows] of Object.entries(PRINTED))
    rows.forEach((row, k) => {
      const got = ex[fam][k];
      row.forEach((c, q) => assert.ok(got[q].equals(c), `${fam} f_${k} lambda^${q}: ${got[q]} vs ${c}`));
      for (let q = row.length; q < got.length; q++) assert.ok(got[q].isZero(), `${fam} f_${k} lambda^${q}`);
    });
  // X^A agrees with the existing implementation in jo.js (used by model A on pages 8-11)
  joPolynomials(8, { exact: true }).forEach((row, k) => row.forEach((c, q) => assert.ok(c.equals(ex.XA[k][q]), `XA f_${k}`)));
});

test("the generated data equal a fresh Float64 run", () => {
  const v = joAllValues(40, 1);
  for (const fam of Object.keys(v)) v[fam].forEach((x, k) => close(DATA[fam][k], x, { rtol: 1e-14, atol: 1e-300 }, `${fam} f_${k}`));
});

test("contact constants agree with JO (1984) tables at lambda = 1", () => {
  const xi = 1e-7;
  const rest = (name) => joScalarAtGap(name, xi) - nearContact[name](xi);
  close(rest("XA11"), 0.99536, { atol: 5e-5 }); // Cooley & O'Neill 0.99536
  close(rest("XA12"), -0.35022, { atol: 5e-5 }); // -0.35022
  close(rest("YA11"), 0.9983, { atol: 1e-4 });
  close(rest("YA12"), -0.2737, { atol: 1e-4 });
  close(rest("YB11"), 0.2390, { atol: 1e-4 });
  close(rest("YB12"), -0.0017, { atol: 1e-4 });
  close(rest("YC11"), 0.7028, { atol: 1e-4 });
  close(rest("YC12"), -0.0274, { atol: 1e-4 });
  // X^C stays finite: (7/8) zeta(3) and -(1/8) zeta(3)
  const z3 = 1.2020569031595942;
  close(joScalarAtGap("XC11", xi), (7 / 8) * z3, { atol: 2e-4 });
  close(joScalarAtGap("XC12", xi), -z3 / 8, { atol: 2e-4 });
});

test("the singular coefficients make every remainder series converge (tail ~ m^-2 with a small constant)", () => {
  // A wrong 1/xi or ln coefficient would make r_m 2^-m grow or stay O(1/m); a wrong xi ln xi coefficient
  // leaves r_m 2^-m m^2 of the order of that error times 4 (YC12 with the printed g5 gives about 0.25).
  for (const name of Object.keys(FUNCTIONS)) {
    const r = remainder(name);
    for (let m = 190; m <= 200; m++) if (r[m] !== 0) assert.ok(Math.abs((r[m] / 2 ** m) * m * m) < 0.1, `${name} m=${m}: ${(r[m] / 2 ** m) * m * m}`);
  }
});

test("invalid arguments are rejected", () => {
  assert.throws(() => joScalarAtGap("XA11", 0.5, { K: -1 }), RangeError);
  assert.throws(() => joScalarAtGap("XA11", 0.5, { K: 2.5 }), RangeError);
  assert.throws(() => joScalarAtGap("XA11", 0), RangeError);
  assert.throws(() => pairScalarsAtGap(0.1, "b"), RangeError);
  assert.throws(() => joAllPolynomials(1000), RangeError);
  // far away the functions tend to the isolated-sphere values without NaN or round-off growth
  close(joScalarAtGap("XA11", 1e160), 1, { rtol: 1e-15 });
  assert.ok(Math.abs(joScalarAtGap("YC12", 1e8)) < 1e-20);
});

test("the near-contact form equals the plain series, and X^A matches collocation", () => {
  for (const s of [2.5, 3, 4, 8])
    for (const n of NAMES) close(joScalar(n, s), joScalar(n, s, { plain: true }), { rtol: 1e-12, atol: 1e-15 }, `${n} s=${s}`);
  // collocation reference (independent method), L = 32; at gap 0.2 it is itself converged to ~1e-7
  for (const [gap, L, X11, X12] of readCsv("axial_convergence").rows.filter((r) => r[1] === 32)) {
    const rtol = gap < 0.5 ? 1e-7 : 1e-10;
    close(joScalarAtGap("XA11", gap), X11, { rtol }, `XA11 gap=${gap}`);
    close(joScalarAtGap("XA12", gap), X12, { rtol }, `XA12 gap=${gap}`);
  }
  const c = axialCollocation([0, 2.05], [1, 1], { L: 64, withCondition: false });
  close(joScalar("XA11", 2.05), c.R[0][0] / (6 * Math.PI), { rtol: 1e-5 });
});

test("12 x 12 matrices: symmetric, positive definite, and scalars round-trip", () => {
  for (const h of [0.001, 0.01, 0.1, 0.5, 2, 8]) {
    for (const model of ["A", "B"]) {
      const sc = pairScalarsAtGap(h, model);
      const R = resistance12(sc, [0, 0, 1]);
      assert.ok(asymmetry(R) < 1e-12, `${model} h=${h} symmetric`);
      assert.ok(Math.min(...eigvalsSym(R)) > 0, `${model} h=${h} positive definite`);
      const back = scalarsOf(R);
      for (const n of NAMES) close(back[n], sc[n], { rtol: 1e-12, atol: 1e-14 }, `${model} ${n}`);
    }
  }
  // model A: the RPY mobility with rotation inverted gives the same matrix as resistance12
  const RA = inv(rpyMobility12([0, 0, 3]));
  const RA2 = resistance12(pairScalarsAtGap(1, "A"), [0, 0, 1]);
  RA.forEach((row, i) => row.forEach((v, j) => close(RA2[i][j], v, { rtol: 1e-10, atol: 1e-12 })));
});

test("model A eliminating rotation gives the translational RPY inverse used on pages 8-11", () => {
  const rv = [0.6, -1.1, 2.7]; // |r| = 3 a
  const RA = inv(rpyMobility12(rv));
  const tf = torqueFree(RA), ref = inv(mobility([[0, 0, 0], rv]));
  tf.forEach((row, i) => row.forEach((v, j) => close(v, ref[i][j], { rtol: 1e-10, atol: 1e-12 })));
});

test("far field: the exact self mobility differs from RPY at (a/r)^4, rotation couples with the right sign", () => {
  const a = 1, mu = 1, m0 = 1 / (6 * Math.PI);
  for (const r of [8, 12]) {
    const MB = inv(resistance12(joScalars(r), [0, 0, 1], { a, mu }));
    // JO mobility x^a_11 = 1 - 60 lambda^3 / ((1+lambda)^4 s^4) + ... = 1 - (15/4)(a/r)^4 at lambda = 1
    close(MB[2][2] / m0 - 1, -(15 / 4) / r ** 4, { rtol: 0.05 });
  }
  // translation-rotation coupling: Omega_2 = F_1 x e / (8 pi mu r^2)
  const r = 15, MB = inv(resistance12(joScalars(r), [0, 0, 1])), MA = rpyMobility12([0, 0, r]);
  const F1 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
  const OmB = matvec(MB, F1).slice(9, 12), OmA = matvec(MA, F1).slice(9, 12);
  const expected = [0, -1 / (8 * Math.PI * r * r), 0]; // e_x x e_z = -e_y
  OmA.forEach((v, i) => close(v, expected[i], { rtol: 1e-12, atol: 1e-15 }));
  OmB.forEach((v, i) => close(v, expected[i], { rtol: 0.02, atol: 1e-7 }));
});
