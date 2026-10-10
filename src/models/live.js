// Interactive sandbox (page 14): spheres in the x-z plane pushed by external forces (a spring to
// the pointer, gravity) and a short-range repulsion, moving with the chosen hydrodynamic model.
//   "none": no hydrodynamic interaction, U_p = F_p / (6 pi mu a)
//   "rpy":  RPY mobility (model A), U = M F
//   "sd":   SD with forces and torques (model B, eq. 14.2), freely rotating spheres, R U = F
// The disturbance flow drawn on the page is the far-field estimate from the forces:
// each sphere contributes (1 + a^2/6 lap) G . F_p, which is exact for a single sphere.
import { solve, norm, vsub } from "../core/linalg.js";
import { mobilityVelocity, flatten, unflatten, rpyPairUnequal } from "../physics/rpy.js";
import { sdTranslationalFT } from "./sdFT.js";
import { heunStep } from "../physics/integrate.js";

export const MODELS = ["none", "rpy", "sd"];

/** Sphere velocities for forces F (one 3-vector per sphere). */
export function velocities(X, F, model, { a = 1, mu = 1 } = {}) {
  if (model === "none") return mobilityVelocity(X, F, { a, mu, model: "self" });
  if (model === "rpy") return mobilityVelocity(X, F, { a, mu, model: "rpy" });
  if (model === "sd") return unflatten(solve(sdTranslationalFT(X, { a, mu }), flatten(F)));
  throw new RangeError(`unknown model ${model}`);
}

/** Smallest surface gap between two spheres (Infinity for fewer than two). */
export function minGap(X, a = 1) {
  let g = Infinity;
  for (let p = 0; p < X.length; p++) for (let q = p + 1; q < X.length; q++) g = Math.min(g, norm(vsub(X[p], X[q])) - 2 * a);
  return g;
}

/**
 * Short-range repulsion F0 * strength * exp(-h / range) along the line of centres (h: gap).
 * Not hydrodynamics: it stands in for the surface forces that keep real particles apart.
 */
export function repulsion(X, { a = 1, F0 = 1, strength = 20, range = 0.02 } = {}) {
  const F = X.map(() => [0, 0, 0]);
  for (let p = 0; p < X.length; p++)
    for (let q = p + 1; q < X.length; q++) {
      const d = vsub(X[p], X[q]), r = norm(d), h = r - 2 * a;
      if (h > 20 * range) continue;
      const f = F0 * strength * Math.exp(-Math.max(h, 0) / range) / r;
      for (let k = 0; k < 3; k++) { F[p][k] += f * d[k]; F[q][k] -= f * d[k]; }
    }
  return F;
}

/** External forces plus the repulsion, as used for the motion (and for the drawn flow). */
export function totalForces(X, forceFn, { a = 1, F0 = 1 } = {}) {
  const Fe = forceFn(X), Fr = repulsion(X, { a, F0 });
  return Fe.map((f, p) => f.map((v, k) => v + Fr[p][k]));
}

/**
 * Advance by dt (or less: each substep is cut so that no sphere moves more than a fifth of the
 * smallest gap, and at most maxSub substeps are taken). forceFn: X -> external forces.
 * U0, if given, are the velocities at X for the same forces (saves one solve).
 * Returns { X, t } with the simulated time actually advanced.
 */
export function advance(X, forceFn, model, dt, { a = 1, mu = 1, maxSub = 20, F0 = 1, U0 = null } = {}) {
  const vel = (Y) => velocities(Y, totalForces(Y, forceFn, { a, F0 }), model, { a, mu });
  let t = 0, U = U0;
  for (let n = 0; n < maxSub && t < dt; n++) {
    U ??= vel(X);
    const vmax = Math.max(1e-12, ...U.map((u) => norm(u)));
    let step = Math.min(dt - t, (0.2 * Math.min(Math.max(minGap(X, a), 1e-4), 1)) / vmax);
    let next = null;
    // a predictor that overlaps (or a solve that fails there) halves the step
    for (let tries = 0; tries < 8; tries++) {
      try {
        const Y = heunStep(X, (Z) => (Z === X ? U : vel(Z)), step);
        if (minGap(Y, a) > 0) { next = Y; break; }
      } catch { /* overlapping predictor: try a smaller step */ }
      step /= 2;
    }
    if (!next) break; // keep the last separated state
    X = next;
    t += step;
    U = null;
  }
  return { X, t };
}

/** Far-field disturbance velocity at x from forces F on spheres X (null inside a sphere). */
export function flowAt(x, X, F, { a = 1, mu = 1 } = {}) {
  const u = [0, 0, 0];
  for (let p = 0; p < X.length; p++) {
    const r = vsub(x, X[p]);
    if (norm(r) <= a) return null;
    const G = rpyPairUnequal(r, a, 0, { mu }); // (1 + a^2/6 lap) G
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) u[i] += G[i][j] * F[p][j];
  }
  return u;
}
