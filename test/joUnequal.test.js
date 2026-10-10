// Model B for unequal spheres (lambda = a2/a1 in 1/4, 1/2, 2, 4). References: the contact constants
// tabulated by Jeffrey & Onishi (1984, Tables 2, 4, 5, 6), exact X^C at contact, the collocation
// solver with unequal radii, and the far-field mobility of JO (section 8).
import { test } from "node:test";
import assert from "node:assert/strict";
import { Rational } from "../src/core/rational.js";
import { inv, eigvalsSym, asymmetry } from "../src/core/linalg.js";
import { joAllPolynomials } from "../src/physics/joRecurrence.js";
import { F as DATA, LAMBDAS } from "../src/data/joValues.js";
import { joScalarAtGap, joPairScalarsAtGap, remainder, NAMES10, functionsFor } from "../src/physics/joFull.js";
import { resistance12, scalarsOf, rpyMobility12, pairScalarsUnequal, torqueFree, NAMES16 } from "../src/models/pairB.js";
import { rpyPairUnequal, selfMobility } from "../src/physics/rpy.js";
import { axialCollocation } from "../src/physics/collocation.js";
import { close } from "./helpers.js";

test("the stored f_k(lambda) equal the exact polynomials evaluated in rational arithmetic", () => {
  const ex = joAllPolynomials(30, { exact: true });
  for (const [num, den] of [[1, 2], [2, 1], [1, 4], [4, 1]]) {
    const l = new Rational(BigInt(num), BigInt(den));
    for (const fam of Object.keys(ex))
      ex[fam].forEach((coeffs, k) => {
        const exact = coeffs.reduceRight((acc, c) => acc.mul(l).add(c), new Rational(0n)).toNumber();
        close(DATA[String(num / den)][fam][k], exact, { rtol: 1e-13, atol: 1e-300 }, `${fam} f_${k} lambda=${num}/${den}`);
      });
  }
});

test("every remainder series converges for every stored ratio", () => {
  for (const l of LAMBDAS)
    for (const name of NAMES10) {
      const r = remainder(name, undefined, l);
      for (let m = 190; m <= 200; m++)
        if (r[m] !== 0) assert.ok(Math.abs((r[m] / 2 ** m) * m * m) < 0.15, `lambda=${l} ${name} m=${m}: ${(r[m] / 2 ** m) * m * m}`);
    }
});

test("contact constants agree with JO (1984) tables 2, 4, 5, 6 for lambda = 1/2 and 1/4", () => {
  const xi = 1e-7;
  // tabulated: A^X 11, 12, 22; A^Y 11, 12, 22; B^Y 11, 12, 21, 22; C^Y 11, 12, 22
  const T = {
    0.5: { AX: [1.0881, -0.2957, 0.8083], AY: [1.0193, -0.2246, 0.9009], BY: [0.1201, 0.0817, 0.0686, -0.4016], CY: [0.8489, -0.0349, 0.4839] },
    0.25: { AX: [1.0836, -0.1880, 0.6253], AY: [1.0073, -0.1181, 0.6871], BY: [0.0620, 0.0592, 0.0594, -0.5664], CY: [0.9280, -0.0349, 0.2097] },
  };
  for (const [ls, t] of Object.entries(T)) {
    const l = +ls, s = joPairScalarsAtGap(xi, { lambda: l });
    const lead = (n, lam = l) => functionsFor(lam).lead[n](xi);
    // the "22" and "21" constants are those of 11 and 12 at 1/lambda (JO 1.9)
    const AX = [s.XA11 - lead("XA11"), s.XA12 - lead("XA12"), s.XA22 - lead("XA11", 1 / l)];
    const AY = [s.YA11 - lead("YA11"), s.YA12 - lead("YA12"), s.YA22 - lead("YA11", 1 / l)];
    const BY = [s.YB11 - lead("YB11"), s.YB12 - lead("YB12"), s.YB21 + lead("YB12", 1 / l), s.YB22 + lead("YB11", 1 / l)];
    const CY = [s.YC11 - lead("YC11"), s.YC12 - lead("YC12"), s.YC22 - lead("YC11", 1 / l)];
    const got = { AX, AY, BY, CY };
    for (const k of Object.keys(t)) t[k].forEach((v, i) => close(got[k][i], v, { atol: 4e-4 }, `lambda=${l} ${k}[${i}]`));
  }
});

test("X^C at contact: lambda^3/(1+lambda)^3 zeta(3, lambda/(1+lambda)) and -8 lambda^3/(1+lambda)^6 zeta(3)", () => {
  const hurwitz = (z, q) => { let s = 0; for (let k = 0; k < 20000; k++) s += (k + q) ** -z; return s + (20000 + q) ** (1 - z) / (z - 1); };
  const z3 = 1.2020569031595942;
  for (const l of [0.5, 2, 0.25]) {
    const L = 1 + l;
    close(joScalarAtGap("XC11", 1e-7, { lambda: l }), (l ** 3 / L ** 3) * hurwitz(3, l / L), { rtol: 3e-4 });
    close(joScalarAtGap("XC12", 1e-7, { lambda: l }), (-8 * l ** 3 / L ** 6) * z3, { rtol: 3e-4 });
  }
});

test("lambda and 1/lambda describe the same pair seen from the other sphere", () => {
  for (const xi of [0.01, 0.3, 2]) {
    const a = joPairScalarsAtGap(xi, { lambda: 0.5 }), b = joPairScalarsAtGap(xi, { lambda: 2 });
    for (const [x, y, sgn] of [["XA12", "XA12", 1], ["YA12", "YA12", 1], ["XC12", "XC12", 1], ["YC12", "YC12", 1], ["XA22", "XA11", 1], ["YB22", "YB11", -1], ["YB21", "YB12", -1]])
      close(a[x], sgn * b[y], { rtol: 1e-11, atol: 1e-14 }, `xi=${xi} ${x}`);
  }
});

test("X^A for unequal spheres matches the collocation solver", () => {
  for (const [a2, gap] of [[0.5, 0.3], [0.5, 1], [0.25, 0.5]]) {
    const c = axialCollocation([0, 1 + a2 + gap], [1, a2], { L: 48, withCondition: false });
    const xi = (2 * gap) / (1 + a2), s = joPairScalarsAtGap(xi, { lambda: a2 });
    close(c.R[0][0] / (6 * Math.PI), s.XA11, { rtol: 1e-6 }, `XA11 a2=${a2} gap=${gap}`);
    close(c.R[1][1] / (6 * Math.PI * a2), s.XA22, { rtol: 1e-6 }, `XA22 a2=${a2} gap=${gap}`);
    close(c.R[0][1] / (3 * Math.PI * (1 + a2)), s.XA12, { rtol: 1e-6 }, `XA12 a2=${a2} gap=${gap}`);
  }
});

test("12 x 12 for unequal spheres: symmetric, positive definite, round trip, far-field mobility", () => {
  for (const a2 of [0.5, 2, 0.25]) {
    for (const xi of [0.001, 0.05, 0.5, 3]) {
      for (const model of ["A", "B"]) {
        const sc = pairScalarsUnequal(xi, model, { a1: 1, a2 });
        const R = resistance12(sc, [0, 0, 1], { a1: 1, a2 });
        assert.ok(asymmetry(R) < 1e-12 * Math.max(...R.flat().map(Math.abs)), `${model} a2=${a2} xi=${xi}`);
        assert.ok(Math.min(...eigvalsSym(R)) > 0, `${model} a2=${a2} xi=${xi} positive definite`);
        const back = scalarsOf(R, { a1: 1, a2 });
        for (const n of NAMES16) close(back[n], sc[n], { rtol: 1e-11, atol: 1e-13 }, `${model} ${n}`);
      }
    }
    // JO mobility x^a_11 = 1 - 60 lambda^3 / ((1+lambda)^4 s^4) + O(s^-6)
    const s = 14, l = a2, r = (s * (1 + a2)) / 2;
    const MB = inv(resistance12(joPairScalarsAtGap(s - 2, { lambda: l }), [0, 0, 1], { a1: 1, a2 }));
    close(MB[2][2] / selfMobility(1) - 1, (-60 * l ** 3) / ((1 + l) ** 4 * s ** 4), { rtol: 0.05 });
    // model A eliminating rotation gives the translational RPY inverse for unequal radii
    const tf = torqueFree(inv(rpyMobility12([0.3, -0.4, r], { a1: 1, a2 })));
    const m11 = selfMobility(1), m22 = selfMobility(a2), c = rpyPairUnequal([0.3, -0.4, r], 1, a2);
    const Mt = Array.from({ length: 6 }, (_, i) => Array.from({ length: 6 }, (_, j) => {
      if ((i < 3) === (j < 3)) return i === j ? (i < 3 ? m11 : m22) : 0;
      return c[i % 3][j % 3];
    }));
    const ref = inv(Mt);
    tf.forEach((row, i) => row.forEach((v, j) => close(v, ref[i][j], { rtol: 1e-10, atol: 1e-12 })));
  }
});

test("ratios that were not generated are rejected", () => {
  assert.throws(() => joScalarAtGap("XA11", 0.5, { lambda: 3 }), RangeError);
  assert.throws(() => pairScalarsUnequal(0.5, "B", { a1: 1, a2: 3 }), RangeError);
});

test("unequal radii without the scalars of sphere 2, or a ratio given as a string, are rejected", () => {
  assert.throws(() => resistance12({ XA11: 1, XA12: 0, YA11: 1, YA12: 0, YB11: 0, YB12: 0, XC11: 1, XC12: 0, YC11: 1, YC12: 0 }, [0, 0, 1], { a1: 1, a2: 0.5 }), RangeError);
  assert.throws(() => joScalarAtGap("XA11", 0.1, { lambda: "0.5" }), RangeError);
  assert.throws(() => rpyMobility12([0, 0, 1.5], { a1: 1, a2: 0.5 }), RangeError); // touching unequal spheres
  // plain series and leading terms for the whole pair
  const p = joPairScalarsAtGap(3, { lambda: 0.5, plain: true }), f = joPairScalarsAtGap(3, { lambda: 0.5 });
  for (const n of NAMES16) close(p[n], f[n], { rtol: 1e-11, atol: 1e-14 }, n);
});
