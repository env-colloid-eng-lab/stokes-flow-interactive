// Stresslet couplings X^G, Y^G, Y^H (joRecurrence families XG, YG, YH; joFull) and the Cartesian
// two-sphere problem in a linear flow (models/pairLinear.js).
import { test } from "node:test";
import assert from "node:assert/strict";
import { joScalarAtGap, joStressletScalarsAtGap, remainder, NAMES_GH, LAMBDAS } from "../src/physics/joFull.js";
import { axialSolve } from "../src/physics/collocation.js";
import { pairInFlow, batchelorGreen, bgTable, relativeTrajectory } from "../src/models/pairLinear.js";
import { FLOWS } from "../src/models/linearFlow.js";

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg}: ${a} vs ${b}`);

// X^G exactly: surface velocity E (x - x_a) with E = diag(-1/2, -1/2, 1) on one sphere; the axial forces
// on the fluid are G~ : E = pi mu (a_a + a_b)^2 X^G_ab e (e from sphere 1 to sphere 2)
function xgCollocation(s, lam, L = 36) {
  const a1 = 1, a2 = lam, r = (s * (a1 + a2)) / 2;
  const strain = (rho, z) => [-rho / 2, z], zero = () => [0, 0];
  const { F } = axialSolve([0, r], [a1, a2], [(al) => (al === 0 ? strain : zero), (al) => (al === 1 ? strain : zero)], { L });
  return { XG11: F[0][0] / (4 * Math.PI), XG12: F[1][0] / (Math.PI * (a1 + a2) ** 2), XG22: F[1][1] / (4 * Math.PI * a2 * a2), XG21: F[0][1] / (Math.PI * (a1 + a2) ** 2) };
}

test("X^G (all four blocks) equals boundary collocation for lambda = 1, 1/2, 2", () => {
  for (const lam of [1, 0.5, 2])
    for (const s of [2.5, 3, 4]) {
      const ref = xgCollocation(s, lam), got = joStressletScalarsAtGap(s - 2, { lambda: lam });
      for (const k of ["XG11", "XG12", "XG22", "XG21"]) close(got[k], ref[k], 2e-8, `${k} lambda=${lam} s=${s}`);
    }
});

test("far field: Faxen laws for the coupling to the other sphere (equal spheres)", () => {
  // sphere 1 stresslet from sphere 2: strain of a Stokeslet (X: -15/4 s^-2, with the first reflection of
  // the forces -8.4375 s^-4), a^2 corrections (Y: -2 s^-4), strain of a rotlet (Y^H: 5/4 s^-3)
  for (const s of [20, 40]) {
    const X12 = joScalarAtGap("XG12", s - 2), Y12 = joScalarAtGap("YG12", s - 2), H12 = joScalarAtGap("YH12", s - 2);
    close(X12 * s ** 2, -3.75 - 2.4375 / s ** 2, 60 / s ** 4, `X^G_12 s=${s}`);
    close(Y12 * s ** 4, -2, 30 / s ** 2, `Y^G_12 s=${s}`);
    close(H12 * s ** 3, 1.25, 30 / s ** 4, `Y^H_12 s=${s}`);
    // the self couplings start with the first reflection: Y^G_11 ~ 27/8 s^-5, Y^H_11 ~ -3/2 s^-6
    close(joScalarAtGap("YG11", s - 2) * s ** 5, 3.375, 30 / s ** 2, `Y^G_11 s=${s}`);
    close(joScalarAtGap("YH11", s - 2) * s ** 6, -1.5, 30 / s ** 2, `Y^H_11 s=${s}`);
  }
});

test("near contact: 150 and 200 terms agree (the singular parts carry the divergence)", () => {
  // the singular coefficients (Jeffrey 1992) only set the speed of convergence
  for (const lam of LAMBDAS)
    for (const n of NAMES_GH)
      for (const xi of [1e-6, 1e-4, 1e-2, 0.3]) {
        const a = joScalarAtGap(n, xi, { lambda: lam }), b = joScalarAtGap(n, xi, { lambda: lam, K: 150 });
        assert.ok(Math.abs(a - b) <= 1e-4 * Math.max(1, Math.abs(a)), `${n} lambda=${lam} xi=${xi}: ${a} vs ${b}`);
      }
  // and the remainders shrink like m^-2 (a wrong ln coefficient would leave m^-1)
  for (const lam of LAMBDAS)
    for (const n of NAMES_GH) {
      const r = remainder(n, undefined, lam), k = r[r.length - 1] !== 0 ? r.length - 1 : r.length - 2; // last term of the parity
      assert.ok(Math.abs((r[k] / 2 ** k) * k) < 0.02, `${n} lambda=${lam}: r_m 2^-m m = ${(r[k] / 2 ** k) * k}`);
    }
});

test("Y^G, Y^H equal the output of Jeffrey's programs ryg.f, ryh.f (4 decimals)", () => {
  const S = [2.01, 2.02, 2.05, 2.1, 2.5], lams = [1, 2, 4, 0.5];
  const T = {
    YG11: [[0.4390, 0.7058, 0.9518, 0.2881], [0.3560, 0.5611, 0.7294, 0.2358], [0.2506, 0.3808, 0.4625, 0.1688], [0.1768, 0.2583, 0.2922, 0.1210], [0.0457, 0.0574, 0.0500, 0.0327]],
    YG12: [[-0.4773, -0.3356, -0.1615, -0.5717], [-0.3938, -0.2711, -0.1258, -0.4782], [-0.2872, -0.1902, -0.0828, -0.3572], [-0.2114, -0.1347, -0.0551, -0.2694], [-0.0678, -0.0381, -0.0133, -0.0928]],
    YH11: [[0.0446, -0.0513, -0.1302, 0.0997], [0.0299, -0.0477, -0.1051, 0.0780], [0.0132, -0.0400, -0.0712, 0.0509], [0.0038, -0.0315, -0.0464, 0.0326], [-0.0033, -0.0084, -0.0067, 0.0045]],
    YH12: [[0.4341, 0.3743, 0.1813, 0.2647], [0.3670, 0.3095, 0.1444, 0.2253], [0.2806, 0.2279, 0.0999, 0.1748], [0.2180, 0.1711, 0.0710, 0.1383], [0.0901, 0.0651, 0.0245, 0.0610]],
  };
  for (const [name, rows] of Object.entries(T))
    rows.forEach((row, i) => lams.forEach((lam, j) => close(joStressletScalarsAtGap(S[i] - 2, { lambda: lam })[name], row[j], 5.01e-5, `${name} s=${S[i]} lambda=${lam}`)));
  // the 22 and 21 functions of the same program at lambda = 2 (= 1/2 of the 11, 12 ones, with the signs of G)
  const v = joStressletScalarsAtGap(0.01, { lambda: 2 });
  close(v.YG22, -0.2881, 5.01e-5, "YG22"); close(v.YG21, 0.5717, 5.01e-5, "YG21");
  close(v.YH22, 0.0997, 5.01e-5, "YH22"); close(v.YH21, 0.2647, 5.01e-5, "YH21");
});

test("X^G is 3/2 times X^A at the 1/xi order; Y^G, Y^H stay finite times ln(1/xi)", () => {
  const xi = 1e-6;
  close(joScalarAtGap("XG11", xi) / joScalarAtGap("XA11", xi), 1.5, 1e-4, "X^G/X^A");
  close(joScalarAtGap("YG11", xi) / Math.log(1 / xi), 0.125, 0.02, "Y^G_11 / ln");
  close(joScalarAtGap("YH12", xi) / Math.log(1 / xi), 0.1, 0.02, "Y^H_12 / ln");
});

// A(r) exactly: two force-free spheres in the axisymmetric strain E = diag(-1/2, -1/2, 1)
function aCollocation(r, L = 36) {
  const strain = (rho, z) => [rho / 2, -z]; // disturbance data -E (x - x_a)
  const unit = (b) => (al) => (al === b ? () => [0, 1] : () => [0, 0]);
  const { F } = axialSolve([0, r], [1, 1], [unit(0), unit(1), () => strain], { L });
  const R = [[F[0][0], F[0][1]], [F[1][0], F[1][1]]], f = [F[0][2], F[1][2]];
  const det = R[0][0] * R[1][1] - R[0][1] * R[1][0];
  const V1 = (-f[0] * R[1][1] + f[1] * R[0][1]) / det, V2 = (-f[1] * R[0][0] + f[0] * R[1][0]) / det;
  return -(V2 - V1) / r; // relative velocity (1 - A) r along e  ->  V2 - V1 = -A r
}

test("Batchelor-Green A(r) from the Cartesian solve equals the axisymmetric collocation", () => {
  for (const r of [2.2, 2.5, 3, 4]) close(batchelorGreen(r).A, aCollocation(r), 1e-7, `A(${r})`);
});

test("Batchelor-Green functions: far field and contact", () => {
  // A ~ 5/r^3 - 8/r^5 + 25/r^6, B ~ 16/(3 r^5) far away; at contact the radial relative motion stops
  for (const r of [15, 25]) {
    const { A, B } = batchelorGreen(r);
    close(A * r ** 3, 5 - 8 / r ** 2, 60 / r ** 3, `A r^3 at ${r}`);
    close(B * r ** 5, 16 / 3, 40 / r ** 2, `B r^5 at ${r}`);
  }
  // near contact (Batchelor & Green 1972): 1 - A ~ 4.077 xi and B ~ 0.406 - 0.78/ln(1/xi), xi = r/a - 2
  for (const xi of [1e-6, 1e-5]) {
    const { A, B } = batchelorGreen(2 + xi);
    close((1 - A) / xi, 4.077, 0.005, `(1 - A)/xi at ${xi}`);
    close(B, 0.406 - 0.78 / Math.log(1 / xi), 0.01, `B at ${xi}`);
  }
  // the far-field model has the right leading term but not the contact behaviour
  const far = batchelorGreen(20, "far");
  close(far.A * 20 ** 3, 5, 0.1, "far model A");
  assert.ok(batchelorGreen(2.001, "far").A < 0.9, "far model lets the spheres meet");
});

test("pair in a linear flow: free spheres rotate, the pair centre follows the flow (equal spheres)", () => {
  const flow = { U: [0, 0, 0], W: [0, 0.5, 0], E: [[0, 0, 0.5], [0, 0, 0], [0.5, 0, 0]] }; // simple shear u = z e_x
  const x1 = [-1.3, 0.2, -0.6], x2 = [1.1, -0.3, 0.9];
  const v = pairInFlow(x1, x2, flow, "exact");
  const mid = x1.map((c, i) => (c + x2[i]) / 2), um = [mid[2], 0, 0];
  [0, 1, 2].forEach((i) => close((v.U1[i] + v.U2[i]) / 2, um[i], 1e-10, `centre velocity ${i}`));
  // symmetry of equal spheres: the angular velocities are equal
  [0, 1, 2].forEach((i) => close(v.W1[i], v.W2[i], 1e-10, `rotation ${i}`));
});

test("relative trajectories in shear: closed next to contact, open further out; overlap without lubrication", () => {
  const flow = FLOWS.shear.make(1), exact = bgTable("exact", { n: 120 });
  const top = (g) => [0, 0, 2 * (1 + g)];
  assert.ok(relativeTrajectory(top(1e-5), flow, exact, { box: 8 }).closed, "closed orbit 2e-5 above contact");
  const open = relativeTrajectory(top(1e-3), flow, exact, { box: 8 });
  assert.ok(!open.closed && open.rmin > 2, "open orbit, no contact");
  const up = relativeTrajectory([-8, 0, 0.05], flow, exact, { box: 8.4 });
  assert.ok(up.rmin > 2 && Math.abs(up.pts.at(-1)[0]) > 8, "an upstream sphere passes round without touching");
  // extension, near and on the compression axis: the gap shrinks below the table but stays positive
  const ext = FLOWS.extension.make(1), finite = (t) => t.pts.every((p) => p.every(Number.isFinite));
  const near = relativeTrajectory([0.01, 0, 7], ext, exact, { box: 8.4 });
  assert.ok(finite(near) && near.rmin > 2 && near.pts.at(-1)[0] > 8, "leaves along the extension axis");
  const axis = relativeTrajectory([0, 0, 7], ext, exact, { box: 8.4 });
  assert.ok(finite(axis) && axis.pts.length < 8000, "on the axis it stops at contact");
  assert.throws(() => relativeTrajectory([0, 0, 1.5], flow, exact), RangeError);
  // the table keeps 1 - A ~ 4.08 xi below its grid
  close(exact.at(1e-12).A1 / 1e-12, 4.0775, 1e-3, "1 - A at xi = 1e-12");
  // without lubrication the path crosses the excluded sphere and stays finite
  for (const model of ["none", "far"]) {
    const t = relativeTrajectory([-8, 0, 0.3], flow, bgTable(model, { n: 60 }), { box: 8.4 });
    assert.ok(t.rmin < 2 && t.pts.every((p) => p.every(Number.isFinite)) && t.pts.at(-1)[0] > 8, `${model}: overlaps and passes`);
  }
});
