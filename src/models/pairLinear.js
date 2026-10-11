// Two force- and torque-free spheres in a linear flow u_inf(x) = U + W x x + E x, solved as a
// resistance problem in Cartesian components (page 16).
//
// On sphere a the disturbance velocity is V_a + O_a x (x - x_a) - E (x - x_a), with
// V_a = U_a - u_inf(x_a) and O_a = Omega_a - W. The forces and torques the spheres exert on the fluid are
//   (F, L) = R_FU (V, O) - (G~ : E, H~ : E),
// where R_FU is the 12 x 12 resistance matrix of src/models/pairB.js and G~, H~ the forces and torques
// caused by surface velocities E (x - x_a) (reciprocal theorem: G~^(ba)_ijk = G^(ab)_jki, H~ likewise).
// For one block (a -> b) along e:
//   G~ : E = pi mu (a_a + a_b)^2 [X^G (e.E.e) e + 2 Y^G (E.e - (e.E.e) e)]
//   H~ : E = pi mu (a_a + a_b)^3 2 Y^H e x (E.e)
// Force- and torque-free spheres: R_FU (V, O) = (G~ : E, H~ : E). No polar coordinates are needed;
// Batchelor & Green's A(r), B(r) follow from the relative velocity, for comparison.
import { solve, matvec, cross, norm } from "../core/linalg.js";
import { resistance12 } from "./pairB.js";
import { joPairScalarsAtGap, joStressletScalarsAtGap } from "../physics/joFull.js";
import { ambient } from "./linearFlow.js";
import { strainField } from "../physics/sphere.js";

export const MODELS = ["none", "far", "exact"];

const dot = (u, v) => u[0] * v[0] + u[1] * v[1] + u[2] * v[2];

/** G~:E and H~:E for every sphere (forces F1, F2 and torques L1, L2), stacked as [F1, F2, L1, L2]. */
export function strainForcing(e, E, sc, { a1 = 1, a2 = 1, mu = 1 } = {}) {
  const Ee = matvec(E, e), eEe = dot(e, Ee), t = Ee.map((v, i) => v - eEe * e[i]), exEe = cross(e, Ee);
  const rad = { 1: a1, 2: a2 };
  const F = (al, be) => {
    const X = sc[`XG${al}${be}`], Y = sc[`YG${al}${be}`], c = Math.PI * mu * (rad[al] + rad[be]) ** 2;
    return e.map((v, i) => c * (X * eEe * v + 2 * Y * t[i]));
  };
  const L = (al, be) => {
    const c = Math.PI * mu * (rad[al] + rad[be]) ** 3 * 2 * sc[`YH${al}${be}`];
    return exEe.map((v) => c * v);
  };
  const add = (p, q) => p.map((v, i) => v + q[i]);
  // the force on b collects the strain data of both spheres a = 1, 2 (block G^(ab))
  return [...add(F(1, 1), F(2, 1)), ...add(F(1, 2), F(2, 2)), ...add(L(1, 1), L(2, 1)), ...add(L(1, 2), L(2, 2))];
}

/**
 * Velocities and angular velocities of two free spheres at x1, x2 in the flow
 * { U, W, E } (src/models/linearFlow.js).
 *   "none":  no interaction, each sphere follows the flow at its centre
 *   "far":   first reflection: each sphere moves with the Faxen average of the other's stresslet flow
 *   "exact": the two-sphere resistance problem above (model B functions)
 */
export function pairInFlow(x1, x2, flow, model, { a1 = 1, a2 = 1, mu = 1 } = {}) {
  if (!MODELS.includes(model)) throw new RangeError(`unknown model ${model}`);
  const rv = x2.map((v, i) => v - x1[i]), r = norm(rv), e = rv.map((v) => v / r);
  const base = { U1: ambient(flow, x1), U2: ambient(flow, x2), W1: [...flow.W], W2: [...flow.W] };
  if (model === "none") return base;
  if (model === "far") {
    // disturbance of a free sphere b in strain E, and its Laplacian (= grad p / mu)
    const kick = (xa, xb, ab, aa) => {
      // inside the other sphere (only reached when the model lets the spheres overlap): use its surface
      const y0 = xa.map((v, i) => v - xb[i]), R0 = norm(y0), y = R0 > ab ? y0 : y0.map((v) => (v * ab * (1 + 1e-9)) / Math.max(R0, 1e-12));
      const R = norm(y), Ey = matvec(flow.E, y), q = dot(y, Ey);
      const u = strainField(y, flow.E, { a: ab, mu }).u.map((v, i) => v - Ey[i]);
      const lap = y.map((v, i) => -5 * ab ** 3 * ((2 * Ey[i]) / R ** 5 - (5 * q * v) / R ** 7));
      return u.map((v, i) => v + ((aa * aa) / 6) * lap[i]);
    };
    const v1 = kick(x1, x2, a2, a1), v2 = kick(x2, x1, a1, a2);
    return { ...base, U1: base.U1.map((v, i) => v + v1[i]), U2: base.U2.map((v, i) => v + v2[i]) };
  }
  if (!(r > a1 + a2)) throw new RangeError("the spheres must not touch");
  const xi = (2 * r) / (a1 + a2) - 2, lambda = a2 / a1, opts = { a1, a2, mu };
  const R = resistance12(joPairScalarsAtGap(xi, { lambda }), e, opts);
  const rhs = strainForcing(e, flow.E, joStressletScalarsAtGap(xi, { lambda }), opts);
  const s = solve(R, rhs); // [V1, V2, O1, O2]
  const pick = (k) => s.slice(3 * k, 3 * k + 3);
  return {
    U1: base.U1.map((v, i) => v + pick(0)[i]), U2: base.U2.map((v, i) => v + pick(1)[i]),
    W1: base.W1.map((v, i) => v + pick(2)[i]), W2: base.W2.map((v, i) => v + pick(3)[i]),
  };
}

/**
 * Batchelor & Green's functions from the Cartesian solution: for a pure strain the relative velocity
 * is U2 - U1 = E r - [A e e + B (I - e e)] E r. A from an axisymmetric strain along e, B from a
 * strain whose E e is perpendicular to e (equal radii a = 1 unless given).
 */
export function batchelorGreen(r, model = "exact", opts = {}) {
  const x1 = [0, 0, 0], x2 = [0, 0, r];
  const flowA = { U: [0, 0, 0], W: [0, 0, 0], E: [[-0.5, 0, 0], [0, -0.5, 0], [0, 0, 1]] };
  const flowB = { U: [0, 0, 0], W: [0, 0, 0], E: [[0, 0, 0.5], [0, 0, 0], [0.5, 0, 0]] };
  const va = pairInFlow(x1, x2, flowA, model, opts), vb = pairInFlow(x1, x2, flowB, model, opts);
  const relA = va.U2[2] - va.U1[2], relB = vb.U2[0] - vb.U1[0];
  // E r along e is r for flowA; E r = (r/2) e_x for flowB
  return { A: 1 - relA / r, B: 1 - relB / (r / 2) };
}

/**
 * A(r), B(r) tabulated from the Cartesian solve on a logarithmic grid of the gap
 * xi = 2r/(a1 + a2) - 2 (1e-6 .. 2e3) and interpolated in ln xi; beyond the grid they decay as r^-3
 * and r^-5. Model "none": A = B = 0. Model "far" lets the spheres overlap; inside it keeps the contact
 * values. `lubricated` (model "exact") tells relativeTrajectory that the spheres never touch.
 */
export function bgTable(model = "exact", { a1 = 1, a2 = 1, n = 260 } = {}) {
  if (!MODELS.includes(model)) throw new RangeError(`unknown model ${model}`);
  const h = (a1 + a2) / 2;
  if (model === "none") return { A: () => 0, B: () => 0, h, lubricated: false };
  const u0 = Math.log(1e-6), u1 = Math.log(2e3), du = (u1 - u0) / (n - 1);
  const As = [], Bs = [];
  for (let k = 0; k < n; k++) {
    const r = h * (2 + Math.exp(u0 + k * du)), v = batchelorGreen(r, model, { a1, a2 });
    As.push(v.A); Bs.push(v.B);
  }
  const rmax = h * (2 + Math.exp(u1));
  const interp = (tab, p) => (r) => {
    if (r >= rmax) return tab[n - 1] * (rmax / r) ** p;
    const xi = r / h - 2;
    if (!(xi > 0)) return tab[0]; // overlap (far model only): keep the contact values
    const u = Math.min(Math.max((Math.log(xi) - u0) / du, 0), n - 1.000001), k = Math.floor(u), f = u - k;
    return tab[k] * (1 - f) + tab[k + 1] * f;
  };
  return { A: interp(As, 3), B: interp(Bs, 5), h, lubricated: model === "exact" };
}

/** Relative velocity of sphere 2 with respect to sphere 1 at separation r: W x r + E r - [A ee + B(I - ee)] E r. */
export function relativeVelocity(rv, flow, tab) {
  const r = norm(rv), e = rv.map((v) => v / r), Er = matvec(flow.E, rv), w = cross(flow.W, rv);
  const A = tab.A(r), B = tab.B(r), eEr = dot(e, Er);
  return rv.map((_, i) => w[i] + Er[i] - (A * eEr * e[i] + B * (Er[i] - eEr * e[i])));
}

/**
 * Relative trajectory from r0 (RK4). Stops on leaving the box |x|, |z| < box, on closing on itself,
 * or after nmax steps. Returns the points, the smallest centre distance and whether it closed.
 * With lubrication (tab.lubricated, model "exact") the state is the direction n and ln(gap): near contact
 * the gap changes on its own scale, which a Cartesian step along the curved path would swamp.
 * The velocity itself is always the Cartesian one above.
 */
export function relativeTrajectory(r0, flow, tab, { box = 10, nmax = 8000, dir = 1, lubricated = tab.lubricated } = {}) {
  const rc = 2 * tab.h, v = (x) => relativeVelocity(x, flow, tab).map((c) => dir * c);
  const gap0 = norm(r0) - rc;
  const useLog = lubricated && gap0 > 0;
  // state -> position, and the time derivative of the state
  // the direction is normalised wherever it is used, so that the RK stages stay at the right gap
  const unit = (n) => { const m = norm(n); return n.map((c) => c / m); };
  const toPos = (S) => (useLog ? unit(S.n).map((c) => c * (rc + Math.exp(S.u))) : S.x);
  const deriv = (S) => {
    const x = toPos(S), vel = v(x);
    if (!useLog) return { x: vel };
    const n = unit(S.n), R = rc + Math.exp(S.u), vr = dot(vel, n);
    return { n: vel.map((c, i) => (c - vr * n[i]) / R), u: vr / Math.exp(S.u), speed: norm(vel) };
  };
  const axpy = (S, dt, D) => (useLog ? { n: S.n.map((c, i) => c + dt * D.n[i]), u: S.u + dt * D.u } : { x: S.x.map((c, i) => c + dt * D.x[i]) });
  let S = useLog ? { n: r0.map((c) => c / norm(r0)), u: Math.log(gap0) } : { x: [...r0] };
  const pts = [r0];
  let rmin = norm(r0), t = 0, far = false, prev = null;
  for (let n = 0; n < nmax; n++) {
    const k1 = deriv(S);
    // at most 0.08 along the path and a change of 0.25 in ln(gap) per step
    const dt = useLog
      ? Math.min(0.08 / Math.max(k1.speed, 1e-12), 0.25 / Math.max(Math.abs(k1.u), 1e-12), 50)
      : Math.min(0.08 / Math.max(norm(k1.x), 1e-12), 50);
    const k2 = deriv(axpy(S, dt / 2, k1)), k3 = deriv(axpy(S, dt / 2, k2)), k4 = deriv(axpy(S, dt, k3));
    const comb = (key) => (useLog && key === "u"
      ? S.u + (dt / 6) * (k1.u + 2 * k2.u + 2 * k3.u + k4.u)
      : S[key].map((c, i) => c + (dt / 6) * (k1[key][i] + 2 * k2[key][i] + 2 * k3[key][i] + k4[key][i])));
    S = useLog ? { n: comb("n"), u: comb("u") } : { x: comb("x") };
    if (useLog) S.n = unit(S.n);
    t += dt;
    const p = toPos(S);
    pts.push(p);
    rmin = Math.min(rmin, norm(p));
    if (Math.abs(p[0]) > box || Math.abs(p[2]) > box) break;
    // closed: back past the start (distance from r0 to the last segment), after having left it
    const away = Math.hypot(p[0] - r0[0], p[1] - r0[1], p[2] - r0[2]) > 0.5;
    if (away && !far) { far = true; prev = p; continue; }
    if (far && prev && segDist(r0, prev, p) < 0.02) { pts.push(r0); return { pts, rmin, closed: true, t }; }
    prev = p;
  }
  return { pts, rmin, closed: false, t };
}

function segDist(p0, b, c) {
  const d = c.map((v, i) => v - b[i]), L2 = dot(d, d) || 1;
  const s = Math.max(0, Math.min(1, dot(p0.map((v, i) => v - b[i]), d) / L2));
  return norm(b.map((v, i) => v + s * d[i] - p0[i]));
}
