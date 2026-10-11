// One sphere (centre at the origin) in a linear ambient flow u_inf(x) = U + W x x + E x
// (W: angular velocity of the fluid = half the vorticity, E: rate of strain). Page 15.
// The exact Stokes solution is the ambient flow plus three disturbances, those of
// src/physics/sphere.js (strainField minus E x, translationField(U_p - U), rotationField(Omega_p - W);
// the tests check the closed forms below against them):
//   the sphere resists the strain E (it cannot deform); it moves or turns unlike the fluid
// A free sphere (no force, no torque) has U_p = U, Omega_p = W; a fixed one U_p = Omega_p = 0.
import { matvec, cross, segmentDistance } from "../core/linalg.js";
import { faxenForce, faxenTorque, faxenStresslet } from "../physics/sphere.js";

// The flows of the page, in the x-z plane (x right, z up, y into the screen), rate G.
export const FLOWS = {
  uniform: { name: "一様流", make: (G) => ({ U: [G, 0, 0], W: [0, 0, 0], E: zero3() }) },
  shear: { name: "単純せん断", make: (G) => ({ U: [0, 0, 0], W: [0, G / 2, 0], E: [[0, 0, G / 2], [0, 0, 0], [G / 2, 0, 0]] }) },
  extension: { name: "平面伸長", make: (G) => ({ U: [0, 0, 0], W: [0, 0, 0], E: [[G, 0, 0], [0, 0, 0], [0, 0, -G]] }) },
  rotation: { name: "純回転", make: (G) => ({ U: [0, 0, 0], W: [0, G, 0], E: zero3() }) },
};
function zero3() { return [[0, 0, 0], [0, 0, 0], [0, 0, 0]]; }

export const ambient = (flow, x) => {
  const ex = matvec(flow.E, x), wx = cross(flow.W, x);
  return flow.U.map((u, i) => u + wx[i] + ex[i]);
};

/** Sphere velocity and angular velocity: "free" (force- and torque-free) or "fixed". */
export function sphereMotion(flow, mode) {
  if (mode === "free") return { Up: [...flow.U], Wp: [...flow.W] };
  if (mode === "fixed") return { Up: [0, 0, 0], Wp: [0, 0, 0] };
  throw new RangeError(`mode must be "free" or "fixed", got ${mode}`);
}

/**
 * Velocity at x (|x| >= a) in the lab frame: total, or only the disturbance the sphere adds.
 * Inside the sphere, the rigid-body motion.
 */
export function velocity(x, flow, motion, { a = 1, part = "total" } = {}) {
  if (part !== "total" && part !== "disturbance") throw new RangeError(`part must be "total" or "disturbance", got ${part}`);
  const r = Math.hypot(...x);
  if (r < a) {
    const w = cross(motion.Wp, x), ua = part === "total" ? null : ambient(flow, x);
    return motion.Up.map((u, i) => u + w[i] - (ua ? ua[i] : 0));
  }
  const V = motion.Up.map((u, i) => u - flow.U[i]), dW = motion.Wp.map((w, i) => w - flow.W[i]);
  const ex = matvec(flow.E, x), q = x[0] * ex[0] + x[1] * ex[1] + x[2] * ex[2], xv = x[0] * V[0] + x[1] * V[1] + x[2] * V[2];
  const a3 = a ** 3, a5 = a ** 5, r2 = r * r, r3 = r2 * r, r5 = r3 * r2;
  const ft = (3 * a) / (4 * r) + a3 / (4 * r3), gt = (3 * a) / (4 * r3) - (3 * a3) / (4 * r5); // translation
  const gs = (-5 * a3) / (2 * r5) + (5 * a5) / (2 * r5 * r2); // strain
  const wx = cross(dW, x);
  const d = [0, 1, 2].map((i) => ft * V[i] + gt * xv * x[i] + (a3 / r3) * wx[i] - (a5 / r5) * ex[i] + gs * q * x[i]);
  if (part === "disturbance") return d;
  const w = cross(flow.W, x);
  return d.map((v, i) => flow.U[i] + w[i] + ex[i] + v);
}

/**
 * Force and torque the sphere exerts on the fluid (= the external force and torque that hold it,
 * in Stokes flow): F = 6 pi mu a (U_p - U), L = 8 pi mu a^3 (Omega_p - W); minus the Faxen force
 * and torque on the sphere (the Laplacian terms vanish in a linear flow). S is the stresslet in
 * the sense of sphere.js (it enters the suspension stress with a plus sign):
 * S = (20/3) pi mu a^3 E, the same for a free and a fixed sphere.
 */
export function moments(flow, motion, { a = 1, mu = 1 } = {}) {
  const zero = [0, 0, 0], zeroM = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  return {
    F: faxenForce(motion.Up, flow.U, zero, { a, mu }).map((v) => -v),
    L: faxenTorque(motion.Wp, flow.W.map((w) => 2 * w), { a, mu }).map((v) => -v),
    S: faxenStresslet(flow.E, zeroM, { a, mu }),
  };
}

/**
 * Streamline in the x-z plane through p0 (classical RK4 on the in-plane velocity, which is the
 * whole velocity there: every flow of the page is symmetric about y = 0). frame: velocity
 * subtracted (e.g. the sphere's), so that the lines are those seen by an observer moving with it.
 * Stops on leaving the box, entering the sphere, or closing on itself.
 */
export function streamline(p0, flow, motion, { a = 1, box = 12, h = 0.08, nmax = 1500, frame = [0, 0, 0], part = "total", dir = 1 } = {}) {
  const v = ([x, z]) => {
    const u = velocity([x, 0, z], flow, motion, { a, part });
    const ux = u[0] - frame[0], uz = u[2] - frame[2], m = Math.hypot(ux, uz);
    return m < 1e-12 ? [0, 0] : [(dir * ux) / m, (dir * uz) / m]; // unit speed: equal arc length steps
  };
  const pts = [p0];
  let p = p0;
  for (let n = 0; n < nmax; n++) {
    const k1 = v(p), k2 = v([p[0] + (h / 2) * k1[0], p[1] + (h / 2) * k1[1]]);
    const k3 = v([p[0] + (h / 2) * k2[0], p[1] + (h / 2) * k2[1]]), k4 = v([p[0] + h * k3[0], p[1] + h * k3[1]]);
    if (!k1[0] && !k1[1]) break; // stagnation point
    const prev = p;
    p = [p[0] + (h / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]), p[1] + (h / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])];
    // closed: the path comes back past p0 (within a fraction of a step, allowing for drift)
    if (n > 20 && segmentDistance(p0, prev, p) < 0.25 * h) { pts.push(p0); break; }
    pts.push(p);
    if (Math.abs(p[0]) > box || Math.abs(p[1]) > box || Math.hypot(p[0], p[1]) < a) break;
  }
  return pts;
}
