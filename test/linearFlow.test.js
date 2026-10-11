// Page 15: one sphere in a linear flow (exact solution assembled from sphere.js).
import { test } from "node:test";
import assert from "node:assert/strict";
import { FLOWS, ambient, sphereMotion, velocity, moments, streamline } from "../src/models/linearFlow.js";
import { translationField, rotationField, strainField, surfaceMoments } from "../src/physics/sphere.js";
import { cross } from "../src/core/linalg.js";

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol, `${msg}: ${a} vs ${b}`);
// a general linear flow (all three parts at once)
const general = { U: [0.3, -0.2, 0.5], W: [0.1, 0.4, -0.3], E: [[0.2, 0.1, -0.3], [0.1, -0.5, 0.05], [-0.3, 0.05, 0.3]] };
const flows = [general, ...Object.values(FLOWS).map((f) => f.make(1.3))];
const dirs = [[1, 0, 0], [0, 0, 1], [0.6, 0, -0.8], [0.48, 0.6, 0.64], [0, -1, 0]];

test("closed form = ambient flow + the sphere.js fields", () => {
  for (const flow of flows)
    for (const mode of ["free", "fixed"]) {
      const m = sphereMotion(flow, mode);
      for (const n of dirs)
        for (const r of [1, 1.3, 4, 25]) {
          const x = n.map((v) => r * v), ex = flow.E.map((row) => row[0] * x[0] + row[1] * x[1] + row[2] * x[2]);
          const ref = [strainField(x, flow.E), translationField(x, m.Up.map((u, i) => u - flow.U[i])), rotationField(x, m.Wp.map((w, i) => w - flow.W[i]))];
          const ua = ambient(flow, x);
          const want = [0, 1, 2].map((i) => ua[i] - ex[i] + ref.reduce((s, f) => s + f.u[i], 0));
          velocity(x, flow, m).forEach((v, i) => close(v, want[i], 1e-12, `${mode} r=${r} ${i}`));
        }
    }
});

test("no slip on the sphere: u = U_p + Omega_p x x", () => {
  for (const flow of flows)
    for (const mode of ["free", "fixed"]) {
      const m = sphereMotion(flow, mode);
      for (const n of dirs) {
        const u = velocity(n, flow, m), w = cross(m.Wp, n);
        u.forEach((v, i) => close(v, m.Up[i] + w[i], 1e-12, `${mode} ${n} ${i}`));
      }
    }
});

test("Stokes equations: div u = 0 and curl(lap u) = 0 (finite differences)", () => {
  const h = 1e-3, m = sphereMotion(general, "fixed"), x0 = [1.7, -0.9, 1.2];
  const u = (x) => velocity(x, general, m);
  const at = (dx) => u(x0.map((v, i) => v + dx[i]));
  const e = (k, s) => [0, 1, 2].map((i) => (i === k ? s : 0));
  const d = (k) => at(e(k, h)).map((v, i) => (v - at(e(k, -h))[i]) / (2 * h));
  close(d(0)[0] + d(1)[1] + d(2)[2], 0, 1e-7, "div");
  const lap = (x) => {
    const c = u(x);
    return [0, 1, 2].map((i) => [0, 1, 2].reduce((s, k) => s + (u(x.map((v, j) => v + (j === k ? h : 0)))[i] - 2 * c[i] + u(x.map((v, j) => v - (j === k ? h : 0)))[i]) / (h * h), 0));
  };
  const H = 2e-2, L = (dx) => lap(x0.map((v, i) => v + dx[i]));
  const dl = (k) => L(e(k, H)).map((v, i) => (v - L(e(k, -H))[i]) / (2 * H));
  const [d0, d1, d2] = [dl(0), dl(1), dl(2)];
  const curl = [d1[2] - d2[1], d2[0] - d0[2], d0[1] - d1[0]];
  curl.forEach((c, i) => close(c, 0, 2e-3, `curl lap u ${i}`));
});

test("force, torque and stresslet agree with the surface traction", () => {
  for (const flow of flows)
    for (const mode of ["free", "fixed"]) {
      const m = sphereMotion(flow, mode), mo = moments(flow, m);
      const field = (r) => {
        const fs = [strainField(r, flow.E), translationField(r, m.Up.map((u, i) => u - flow.U[i])), rotationField(r, m.Wp.map((w, i) => w - flow.W[i]))];
        return { p: fs.reduce((s, f) => s + f.p, 0), sigma: [0, 1, 2].map((i) => [0, 1, 2].map((j) => fs.reduce((s, f) => s + f.sigma[i][j], 0))) };
      };
      const sm = surfaceMoments(field, { nt: 10, np: 20 });
      // traction integrals give what the fluid exerts on the sphere: minus what the sphere exerts;
      // the stresslet is defined from the same traction (sphere.js), with no sign change
      [0, 1, 2].forEach((i) => {
        close(sm.F[i], -mo.F[i], 1e-9, `${mode} F${i}`);
        close(sm.T[i], -mo.L[i], 1e-9, `${mode} L${i}`);
        [0, 1, 2].forEach((j) => close(sm.S[i][j], mo.S[i][j], 1e-9, `${mode} S${i}${j}`));
      });
    }
});

test("a free sphere does not disturb a uniform flow or a rigid rotation", () => {
  for (const k of ["uniform", "rotation"]) {
    const flow = FLOWS[k].make(2), m = sphereMotion(flow, "free");
    for (const n of dirs) velocity(n.map((v) => 1.7 * v), flow, m, { part: "disturbance" }).forEach((v) => close(v, 0, 1e-14, k));
  }
});

test("far field: stresslet for a free sphere in strain, Stokeslet for a fixed one in uniform flow", () => {
  const ext = FLOWS.extension.make(1), x = [30, 0, 40], r = 50;
  const q = x[0] * x[0] - x[2] * x[2]; // x.E.x
  const d = velocity(x, ext, sphereMotion(ext, "free"), { part: "disturbance" });
  d.forEach((v, i) => close(v, (-2.5 * q * x[i]) / r ** 5, 1e-6, `stresslet ${i}`));
  const uni = FLOWS.uniform.make(1), R = 1e4;
  const s = velocity([0, 0, R], uni, sphereMotion(uni, "fixed"), { part: "disturbance" });
  close(s[0], -0.75 / R, 1e-8, "Stokeslet across the flow");
});

test("shear: closed streamlines next to a free sphere, open ones further out", () => {
  const flow = FLOWS.shear.make(1), m = sphereMotion(flow, "free");
  for (const z0 of [1.02, 1.05, 1.1]) {
    const near = streamline([0, z0], flow, m), last = near[near.length - 1];
    assert.ok(last[0] === 0 && last[1] === z0 && near.length > 20 && near.length < 600, `near ${z0}: closed after one turn (${near.length} points)`);
  }
  const far = streamline([0, 4], flow, m);
  assert.ok(Math.abs(far[far.length - 1][0]) > 12, "far: leaves the box");
});

test("velocity: an unknown part is refused", () => {
  const flow = FLOWS.shear.make(1);
  assert.throws(() => velocity([2, 0, 0], flow, sphereMotion(flow, "free"), { part: "disturb" }), RangeError);
});

test("ambient flows: shear u = G z e_x, extension G (x, 0, -z), rotation G (z, 0, -x)", () => {
  const x = [0.7, 0.2, -1.1], G = 1.3;
  assert.deepEqual(ambient(FLOWS.shear.make(G), x).map((v) => +v.toFixed(12)), [+(G * x[2]).toFixed(12), 0, 0]);
  assert.deepEqual(ambient(FLOWS.extension.make(G), x).map((v) => +v.toFixed(12)), [+(G * x[0]).toFixed(12), 0, +(-G * x[2]).toFixed(12)]);
  assert.deepEqual(ambient(FLOWS.rotation.make(G), x).map((v) => +v.toFixed(12)), [+(G * x[2]).toFixed(12), 0, +(-G * x[0]).toFixed(12)]);
});
