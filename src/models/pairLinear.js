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
import { solve, matvec, cross, norm, dot, vadd, vsub, vscale, fnorm, segmentDistance } from "../core/linalg.js";
import { resistance12 } from "./pairB.js";
import { joPairScalarsAtGap, joStressletScalarsAtGap } from "../physics/joFull.js";
import { ambient } from "./linearFlow.js";
import { strainField } from "../physics/sphere.js";

export const MODELS = ["none", "far", "exact"];

/** G~:E and H~:E for every sphere (forces F1, F2 and torques L1, L2), stacked as [F1, F2, L1, L2]. */
export function strainForcing(e, E, sc, { a1 = 1, a2 = 1, mu = 1 } = {}) {
  const Ee = matvec(E, e), eEe = dot(e, Ee), t = vsub(Ee, vscale(e, eEe)), exEe = cross(e, Ee);
  const rad = { 1: a1, 2: a2 };
  const F = (al, be) => {
    const X = sc[`XG${al}${be}`], Y = sc[`YG${al}${be}`], c = Math.PI * mu * (rad[al] + rad[be]) ** 2;
    return vadd(vscale(e, c * X * eEe), vscale(t, 2 * c * Y));
  };
  const L = (al, be) => vscale(exEe, Math.PI * mu * (rad[al] + rad[be]) ** 3 * 2 * sc[`YH${al}${be}`]);
  // the force on b collects the strain data of both spheres a = 1, 2 (block G^(ab))
  return [...vadd(F(1, 1), F(2, 1)), ...vadd(F(1, 2), F(2, 2)), ...vadd(L(1, 1), L(2, 1)), ...vadd(L(1, 2), L(2, 2))];
}

// [V1, V2, O1, O2] of the exact problem for each strain in Es (one factorisation of R_FU)
function exactSolve(e, r, Es, opts) {
  const { a1, a2 } = opts;
  if (!(r > a1 + a2)) throw new RangeError("the spheres must not touch");
  const xi = (2 * r) / (a1 + a2) - 2, lambda = a2 / a1;
  const R = resistance12(joPairScalarsAtGap(xi, { lambda }), e, opts), sc = joStressletScalarsAtGap(xi, { lambda });
  const cols = Es.map((E) => strainForcing(e, E, sc, opts));
  const S = solve(R, cols[0].map((_, i) => cols.map((c) => c[i])));
  return Es.map((_, j) => S.map((row) => row[j]));
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
  const rv = vsub(x2, x1), r = norm(rv), e = vscale(rv, 1 / r);
  const base = { U1: ambient(flow, x1), U2: ambient(flow, x2), W1: [...flow.W], W2: [...flow.W] };
  if (model === "none") return base;
  if (model === "far") {
    // disturbance of a free sphere b in strain E, and its Laplacian (= grad p / mu)
    const kick = (xa, xb, ab, aa) => {
      // inside the other sphere (only reached when the model lets the spheres overlap): use its surface
      const y0 = vsub(xa, xb), R0 = norm(y0), y = R0 > ab ? y0 : vscale(y0, (ab * (1 + 1e-9)) / Math.max(R0, 1e-12));
      const R = norm(y), Ey = matvec(flow.E, y), q = dot(y, Ey);
      const u = vsub(strainField(y, flow.E, { a: ab, mu }).u, Ey);
      const lap = y.map((v, i) => -5 * ab ** 3 * ((2 * Ey[i]) / R ** 5 - (5 * q * v) / R ** 7));
      return vadd(u, vscale(lap, (aa * aa) / 6));
    };
    return { ...base, U1: vadd(base.U1, kick(x1, x2, a2, a1)), U2: vadd(base.U2, kick(x2, x1, a1, a2)) };
  }
  const [s] = exactSolve(e, r, [flow.E], { a1, a2, mu });
  return {
    U1: vadd(base.U1, s.slice(0, 3)), U2: vadd(base.U2, s.slice(3, 6)),
    W1: vadd(base.W1, s.slice(6, 9)), W2: vadd(base.W2, s.slice(9, 12)),
  };
}

const STRAIN_A = [[-0.5, 0, 0], [0, -0.5, 0], [0, 0, 1]]; // axisymmetric about e = e_z: E r = r e_z
const STRAIN_B = [[0, 0, 0.5], [0, 0, 0], [0.5, 0, 0]]; // E r = (r/2) e_x, perpendicular to e

/**
 * Batchelor & Green's functions from the Cartesian solution: for a pure strain the relative velocity
 * is U2 - U1 = E r - [A e e + B (I - e e)] E r. A from an axisymmetric strain along e, B from a
 * strain whose E e is perpendicular to e (equal radii a = 1 unless given).
 */
export function batchelorGreen(r, model = "exact", { a1 = 1, a2 = 1, mu = 1 } = {}) {
  const x1 = [0, 0, 0], x2 = [0, 0, r];
  if (model === "exact") {
    // only the disturbance V2 - V1 is needed: -A r along e and -B r/2 across
    const [sA, sB] = exactSolve([0, 0, 1], r, [STRAIN_A, STRAIN_B], { a1, a2, mu });
    return { A: -(sA[5] - sA[2]) / r, B: -(sB[3] - sB[0]) / (r / 2) };
  }
  const flow = (E) => ({ U: [0, 0, 0], W: [0, 0, 0], E });
  const va = pairInFlow(x1, x2, flow(STRAIN_A), model, { a1, a2, mu }), vb = pairInFlow(x1, x2, flow(STRAIN_B), model, { a1, a2, mu });
  return { A: 1 - (va.U2[2] - va.U1[2]) / r, B: 1 - (vb.U2[0] - vb.U1[0]) / (r / 2) };
}

/**
 * A, B tabulated from the Cartesian solve on a logarithmic grid of the gap xi = 2r/(a1 + a2) - 2
 * (1e-6 .. 2e3), interpolated in ln xi; beyond the grid A, B decay as r^-3, r^-5. tab.at(xi) gives
 * { A, B, A1 } with A1 = 1 - A kept separately (it is ~ 4 xi near contact, below the rounding of A). Closer than the grid the exact functions follow the lubrication forms 1 - A ~ xi and
 * B ~ c0 - c1/ln(1/xi) (matched to the grid), so the radial motion keeps slowing down with the gap.
 * Model "none": A = B = 0. Model "far" lets the spheres overlap and keeps the contact values there.
 * `lubricated` (model "exact") tells relativeTrajectory that the spheres never touch.
 */
export function bgTable(model = "exact", { a1 = 1, a2 = 1, n = 260 } = {}) {
  if (!MODELS.includes(model)) throw new RangeError(`unknown model ${model}`);
  const h = (a1 + a2) / 2;
  if (model === "none") return { at: () => ({ A: 0, B: 0, A1: 1 }), h, lubricated: false };
  const u0 = Math.log(1e-6), u1 = Math.log(2e3), du = (u1 - u0) / (n - 1);
  const As = [], Bs = [];
  for (let k = 0; k < n; k++) {
    const v = batchelorGreen(h * (2 + Math.exp(u0 + k * du)), model, { a1, a2 });
    As.push(v.A); Bs.push(v.B);
  }
  const lubricated = model === "exact", xi0 = Math.exp(u0), k1 = Math.min(10, n - 1), xi1 = Math.exp(u0 + k1 * du);
  const c1 = (Bs[0] - Bs[k1]) / (1 / Math.log(1 / xi1) - 1 / Math.log(1 / xi0)), c0 = Bs[0] + c1 / Math.log(1 / xi0);
  const out = (A, B, A1 = 1 - A) => ({ A, B, A1 });
  const at = (xi) => {
    if (xi >= Math.exp(u1)) {
      const f = (2 + Math.exp(u1)) / (2 + xi);
      return out(As[n - 1] * f ** 3, Bs[n - 1] * f ** 5);
    }
    if (!(xi > xi0)) {
      if (!lubricated) return out(As[0], Bs[0]);
      const x = Math.max(xi, 0), A1 = ((1 - As[0]) * x) / xi0;
      return out(1 - A1, x > 0 ? c0 - c1 / Math.log(1 / x) : c0, A1);
    }
    const u = Math.min((Math.log(xi) - u0) / du, n - 1.000001), k = Math.floor(u), f = u - k;
    return out(As[k] * (1 - f) + As[k + 1] * f, Bs[k] * (1 - f) + Bs[k + 1] * f);
  };
  return { at, h, lubricated };
}

/**
 * Relative velocity of sphere 2 with respect to sphere 1 at separation r: W x r + E r - [A ee + B(I - ee)] E r.
 * The gap xi = |r|/h - 2 may be given separately, when |r| would round it away near contact.
 */
export function relativeVelocity(rv, flow, tab, xi = norm(rv) / tab.h - 2) {
  const r = norm(rv);
  if (r === 0) return [0, 0, 0]; // the stagnation point (models that let the spheres overlap)
  const e = vscale(rv, 1 / r), Er = matvec(flow.E, rv), eEr = dot(e, Er);
  const { B, A1 } = tab.at(xi);
  // = W x r + (1 - A)(e.E.r) e + (1 - B)(E r - (e.E.r) e)
  return vadd(cross(flow.W, rv), vadd(vscale(e, A1 * eEr), vscale(vsub(Er, vscale(e, eEr)), 1 - B)));
}

/**
 * Relative trajectory from r0 (RK4). Stops on leaving the box |x|, |z| < box, on closing on itself,
 * or after nmax steps. Returns the points, the smallest centre distance and whether it closed.
 * With lubrication (tab.lubricated, model "exact") the state is the direction n and ln(gap): near
 * contact the gap changes on its own scale, which a Cartesian step along the curved path would swamp.
 * The velocity itself is always the Cartesian one above.
 */
export function relativeTrajectory(r0, flow, tab, { box = 10, nmax = 8000, dir = 1, lubricated = tab.lubricated } = {}) {
  const rc = 2 * tab.h, v = (x, xi) => vscale(relativeVelocity(x, flow, tab, xi), dir);
  const gap0 = norm(r0) - rc;
  if (lubricated && !(gap0 > 0)) throw new RangeError("the spheres must not touch at the start");
  // the step resolves the flow's time scale (RK4 is unstable near stagnation points otherwise)
  const tmax = 1 / Math.max(fnorm(flow.E) + norm(flow.W), 1e-12);
  // state -> position, and the time derivative of the state
  // the direction is normalised wherever it is used, so that the RK stages stay at the right gap
  const unit = (n) => vscale(n, 1 / norm(n));
  const toPos = (S) => (lubricated ? vscale(unit(S.n), rc + Math.exp(S.u)) : S.x);
  const deriv = (S) => {
    if (!lubricated) { const vel = v(S.x); return { x: vel, speed: norm(vel) }; }
    const n = unit(S.n), g = Math.exp(S.u), vel = v(vscale(n, rc + g), g / tab.h), vr = dot(vel, n);
    return { n: vscale(vsub(vel, vscale(n, vr)), 1 / (rc + g)), u: vr / g, speed: norm(vel) };
  };
  const axpy = (S, dt, D) => (lubricated ? { n: vadd(S.n, vscale(D.n, dt)), u: S.u + dt * D.u } : { x: vadd(S.x, vscale(D.x, dt)) });
  let S = lubricated ? { n: unit(r0), u: Math.log(gap0) } : { x: [...r0] };
  const pts = [r0];
  let rmin = norm(r0), t = 0, far = false, prev = null;
  for (let n = 0; n < nmax; n++) {
    const k1 = deriv(S);
    if (!(k1.speed > 0)) break; // at a stagnation point
    // at most 0.08 along the path and a change of 0.25 in ln(gap) per step
    const dt = Math.min(0.08 / Math.max(k1.speed, 1e-12), lubricated ? 0.25 / Math.max(Math.abs(k1.u), 1e-12) : Infinity, tmax);
    const k2 = deriv(axpy(S, dt / 2, k1)), k3 = deriv(axpy(S, dt / 2, k2)), k4 = deriv(axpy(S, dt, k3));
    const comb = (key) => (key === "u"
      ? S.u + (dt / 6) * (k1.u + 2 * k2.u + 2 * k3.u + k4.u)
      : S[key].map((c, i) => c + (dt / 6) * (k1[key][i] + 2 * k2[key][i] + 2 * k3[key][i] + k4[key][i])));
    S = lubricated ? { n: unit(comb("n")), u: comb("u") } : { x: comb("x") };
    t += dt;
    const p = toPos(S);
    pts.push(p);
    rmin = Math.min(rmin, norm(p));
    if (Math.abs(p[0]) > box || Math.abs(p[2]) > box) break;
    // on the compression axis the gap closes forever (in infinite time): stop before it underflows
    if (lubricated && S.u < -600) break;
    // closed: back past the start (distance from r0 to the last segment), after having left it
    const away = norm(vsub(p, r0)) > 0.5;
    if (away && !far) { far = true; prev = p; continue; }
    if (far && prev && segmentDistance(r0, prev, p) < 0.02) { pts.push(r0); return { pts, rmin, closed: true, t }; }
    prev = p;
  }
  return { pts, rmin, closed: false, t };
}
