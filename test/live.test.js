// Page 14 (interactive sandbox): velocities, flow field, time stepping.
import { test } from "node:test";
import assert from "node:assert/strict";
import { velocities, flowAt, advance, minGap, repulsion } from "../src/models/live.js";
import { rpyParallelPerp } from "../src/physics/rpy.js";
import { coefficientsB } from "../src/models/pairB.js";

const close = (a, b, tol, msg) => assert.ok(Math.abs(a - b) <= tol * Math.max(1, Math.abs(b)), `${msg}: ${a} vs ${b}`);
const six = 6 * Math.PI;

test("dragging one sphere: the other follows with the pair mobility", () => {
  for (const r of [2.5, 4, 8]) {
    const X = [[0, 0, 0], [r, 0, 0]];
    for (const [dir, key] of [[[1, 0, 0], "par"], [[0, 0, 1], "perp"]]) {
      const F = [dir, [0, 0, 0]];
      const k = dir[0] ? 0 : 2;
      assert.equal(velocities(X, F, "none")[1][k], 0);
      close(velocities(X, F, "rpy")[1][k], rpyParallelPerp(r)[key], 1e-12, `RPY ${key} r=${r}`);
      // model B: torque-free exact pair, U = M F with M the inverse of the 2x2 (self, cross) block
      close(velocities(X, F, "sd")[1][k], coefficientsB(r).M[key].cross, 1e-9, `SD ${key} r=${r}`);
      close(velocities(X, F, "sd")[0][k], coefficientsB(r).M[key].self, 1e-9, `SD self ${key} r=${r}`);
    }
  }
});

test("flow field of one sphere: uniform F/(6 pi mu a) on the surface, Stokeslet far away", () => {
  const F = [[0.3, -0.2, 1]];
  for (const n of [[1, 0, 0], [0, 0, 1], [0.6, 0, 0.8], [0, 1, 0]]) {
    const u = flowAt(n.map((v) => v * (1 + 1e-12)), [[0, 0, 0]], F);
    u.forEach((v, i) => close(v, F[0][i] / six, 1e-9, `surface ${n}`));
  }
  assert.equal(flowAt([0.5, 0, 0], [[0, 0, 0]], F), null);
  const R = 1e4, u = flowAt([0, 0, R], [[0, 0, 0]], [[0, 0, 1]]);
  close(u[2], 1 / (4 * Math.PI * R), 1e-6, "far field along the force");
});

test("time stepping: strong pushes never make spheres overlap", () => {
  for (const model of ["none", "rpy", "sd"]) {
    let X = [[0, 0, 0], [2.6, 0, 0], [1.3, 0, 2.3]];
    const push = (Y) => Y.map((p) => [-8 * six * Math.sign(p[0] - 1.3) || 0, 0, -8 * six * Math.sign(p[2] - 0.8)]);
    let t = 0;
    for (let i = 0; i < 40; i++) {
      const s = advance(X, push, model, 0.05, { F0: six, maxSub: 10 });
      X = s.X; t += s.t;
    }
    assert.ok(minGap(X) > 0, `${model}: gap ${minGap(X)}`);
    assert.ok(t > 0, `${model}: time advanced`);
  }
});

test("repulsion: equal and opposite, short ranged", () => {
  const F = repulsion([[0, 0, 0], [2.02, 0, 0], [10, 0, 0]], { F0: 1 });
  close(F[0][0], -F[1][0], 1e-15, "pair");
  assert.ok(F[0][0] < 0 && F[2][0] === 0);
});

test("one sphere under a force moves at F / (6 pi mu a) in every model", () => {
  for (const model of ["none", "rpy", "sd"]) close(velocities([[1, 2, 3]], [[0, 0, -six]], model)[0][2], -1, 1e-12, model);
});

test("time stepping from a very small gap: no overlap and no failed solve", () => {
  for (const model of ["rpy", "sd"])
    for (const g of [0.003, 0.006, 0.009]) {
      const X0 = [[0, 0, 0], [2 + g, 0, 0]];
      const s = advance(X0, () => [[30 * six, 0, 0], [0, 0, 0]], model, 0.05, { F0: six, maxSub: 10 });
      assert.ok(minGap(s.X) > 0, `${model} gap ${g}: ${minGap(s.X)}`);
      assert.ok(s.t > 0, `${model} gap ${g}: time advanced`);
    }
});
