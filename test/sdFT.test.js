// Many-body SD with forces and torques (model B) and its pair coefficients for pages 8-10.
import { test } from "node:test";
import assert from "node:assert/strict";
import { eigvalsSym, asymmetry, inv } from "../src/core/linalg.js";
import { sdResistanceFT, sdTranslationalFT, torqueFreeN, rpyMobilityFT, pairwiseSumB, sdBundle } from "../src/models/sdFT.js";
import { resistance12, torqueFree, coefficientsB, coefficientsBAtGap, pairScalarsAtGap } from "../src/models/pairB.js";
import { joScalarsAtGap } from "../src/physics/joFull.js";
import { mobility } from "../src/physics/rpy.js";
import { sdResistance } from "../src/models/sd.js";
import { axialManyBody, resistance, pairwiseSumResistance, axialComponents } from "../src/models/manyBody.js";
import { close } from "./helpers.js";

test("two spheres: the SD construction reproduces model B exactly", () => {
  const rv = [0.4, -0.7, 2.05];
  const r = Math.hypot(...rv), e = rv.map((x) => x / r);
  const exact = torqueFree(resistance12(joScalarsAtGap(r - 2), e));
  const sd = sdTranslationalFT([[0, 0, 0], rv]);
  sd.forEach((row, i) => row.forEach((v, j) => close(v, exact[i][j], { rtol: 1e-9, atol: 1e-9 })));
  // for two spheres the pairwise sum of exact pairs is the same matrix
  const two = pairwiseSumB([[0, 0, 0], rv]);
  two.forEach((row, i) => row.forEach((v, j) => close(v, exact[i][j], { rtol: 1e-9, atol: 1e-9 })));
});

test("the RPY mobility with rotation, inverted and reduced, gives the translational RPY inverse", () => {
  const X = [[0, 0, 0], [2.6, 0, 0.8], [-0.5, 1.9, 2.4]];
  const tf = torqueFreeN(inv(rpyMobilityFT(X)));
  const ref = inv(mobility(X));
  tf.forEach((row, i) => row.forEach((v, j) => close(v, ref[i][j], { rtol: 1e-10, atol: 1e-12 })));
});

test("coaxial spheres: axial components agree with the SD-type matrix of page 11", () => {
  const z = [0, 2.3, 4.4], X = z.map((zz) => [0, 0, zz]);
  const R = sdTranslationalFT(X);
  // page 11: RPY far field and exact (collocation) pairs, axial only
  const ex = axialManyBody(z, z.map(() => 1), { L: 56 }); // gap 0.1 needs a high order
  const Rr = axialComponents(resistance(X)), R2r = axialComponents(pairwiseSumResistance(X));
  const page11 = sdResistance(Rr, ex.R2B, R2r);
  const axial = axialComponents(R);
  axial.forEach((row, i) => row.forEach((v, j) => close(v, page11[i][j], { rtol: 1e-6, atol: 1e-8 })));
});

test("SD with forces and torques is symmetric and positive definite, also near contact", () => {
  const X = [[0, 0, 0], [2.01, 0, 0], [1.0, 1.74, 0.1], [3.1, 1.8, -0.4]];
  const R6 = sdResistanceFT(X), R3 = torqueFreeN(R6);
  assert.ok(asymmetry(R6) < 1e-9 * Math.max(...R6.flat().map(Math.abs)));
  assert.ok(Math.min(...eigvalsSym(R6)) > 0);
  assert.ok(Math.min(...eigvalsSym(R3)) > 0);
  const b = sdBundle(X);
  b.R.forEach((row, i) => row.forEach((v, j) => close(v, R3[i][j], { rtol: 1e-12, atol: 1e-12 })));
  assert.equal(b.R2B.length, 12);
});

test("model B coefficients for pages 8-9: M is the inverse of R, the axial part is X^A", () => {
  for (const r of [2.05, 3, 6]) {
    const c = coefficientsB(r);
    const s = joScalarsAtGap(r - 2);
    close(c.R.par.self / (6 * Math.PI), s.XA11, { rtol: 1e-12 });
    close(c.R.par.cross / (6 * Math.PI), s.XA12, { rtol: 1e-12 });
    for (const d of ["par", "perp"]) {
      const { self, cross } = c.R[d], m = c.M[d];
      close(self * m.self + cross * m.cross, 1, { rtol: 1e-12 });
      close(self * m.cross + cross * m.self, 0, { atol: 1e-12 });
    }
    // torque-free transverse resistance lies between model A's and the rotation-fixed Y^A
    const tfA = torqueFree(resistance12(pairScalarsAtGap(r - 2, "A"), [0, 0, 1]))[0][0];
    assert.ok(c.R.perp.self / (6 * Math.PI) < s.YA11 + 1e-12);
    assert.ok(c.R.perp.self > tfA * (1 - 1e-12));
  }
  assert.throws(() => coefficientsB(2), RangeError);
  // near contact the gap is passed directly and the 2x2 inverse is factored
  const h = 1e-8, cg = coefficientsBAtGap(h), par = cg.R.par;
  close((par.self - par.cross) / (6 * Math.PI), joScalarsAtGap(h).XA11 - joScalarsAtGap(h).XA12, { rtol: 1e-12 });
  // mobility from the two eigenmodes, m_self = (1/(self - cross) + 1/(self + cross)) / 2, to full precision
  close(cg.M.par.self, (1 / (par.self - par.cross) + 1 / (par.self + par.cross)) / 2, { rtol: 1e-12 });
});
