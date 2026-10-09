// The block decomposition used for display must reproduce the direct matrices.
import { test } from "node:test";
import assert from "node:assert/strict";
import { inv } from "../src/core/linalg.js";
import { mobility } from "../src/physics/rpy.js";
import { normalResistance } from "../src/physics/lubrication.js";
import { joPolynomials, joResistance } from "../src/physics/jo.js";
import { mobilityCoefficients, resistanceCoefficients, pairMatrix, pairFrame, toPairFrame } from "../src/models/pairA.js";
import { closeArray, close } from "./helpers.js";

const theta = 0.6;
const { e, axes } = pairFrame(theta);
const place = (r) => [[0.3, -0.2, 0.5], [0.3 + r * e[0], -0.2 + r * e[1], 0.5 + r * e[2]]];

test("pair decomposition reproduces the RPY and Oseen mobility matrices", () => {
  for (const r of [2.4, 3, 7]) {
    closeArray(pairMatrix(mobilityCoefficients(r), e), mobility(place(r)), { rtol: 1e-13, atol: 1e-16 });
    closeArray(pairMatrix(mobilityCoefficients(r, { kind: "oseen" }), e), mobility(place(r), { model: "oseen" }), { rtol: 1e-13, atol: 1e-16 });
  }
  closeArray(pairMatrix(mobilityCoefficients(1.5), e), mobility(place(1.5)), { rtol: 1e-13, atol: 1e-16 });
});

test("resistance coefficients equal the inverse matrix, with and without lubrication", () => {
  for (const r of [2.05, 2.15, 3]) {
    const c = resistanceCoefficients(r);
    closeArray(pairMatrix(c.R, e), inv(mobility(place(r))), { rtol: 1e-10, atol: 1e-12 });
    const cl = resistanceCoefficients(r, { lubrication: true });
    closeArray(pairMatrix(cl.R, e), normalResistance(place(r)).R, { rtol: 1e-10, atol: 1e-12 });
  }
});

test("JO axial entries equal joResistance; pair frame is block diagonal", () => {
  const polys = joPolynomials(60);
  const c = resistanceCoefficients(3, { axial: "jo", polys });
  const R2 = joResistance(1, 1, 3, polys);
  close(c.R.par.self, R2[0][0], { rtol: 1e-14 });
  close(c.R.par.cross, R2[0][1], { rtol: 1e-14 });
  assert.ok(c.jo.converged);
  const B = toPairFrame(pairMatrix(c.R, e), axes);
  const expected = [
    [c.R.par.self, 0, 0, c.R.par.cross, 0, 0],
    [0, c.R.perp.self, 0, 0, c.R.perp.cross, 0],
    [0, 0, c.R.perp.self, 0, 0, c.R.perp.cross],
  ];
  const full = expected.concat(expected.map((row) => row.slice(3).concat(row.slice(0, 3))));
  closeArray(B, full, { rtol: 1e-12, atol: 1e-12 });
  // near contact the series is flagged as not converged
  assert.equal(resistanceCoefficients(2.02, { axial: "jo", polys }).jo.converged, false);
});
