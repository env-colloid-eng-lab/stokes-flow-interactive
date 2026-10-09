// Time stepping for overdamped (inertialess) particle motion, and Brownian increments.
import { applyNoise, noiseFactor } from "./brownian.js";

const axpy = (X, dt, U) => X.map((x, p) => x.map((xi, k) => xi + dt * U[p][k]));

// Heun (explicit trapezoid). vfun: X -> U, deterministic. No stochastic interpretation implied.
export function heunStep(X, vfun, dt) {
  if (!(dt > 0)) throw new RangeError("dt > 0 required");
  const k1 = vfun(X);
  const k2 = vfun(axpy(X, dt, k1));
  return X.map((x, p) => x.map((xi, k) => xi + (dt / 2) * (k1[p][k] + k2[p][k])));
}

// Explicit midpoint, as used by TwoSphereLab.pair_motion.
export function midpointStep(X, vfun, dt) {
  const k1 = vfun(X);
  return axpy(X, dt, vfun(axpy(X, dt / 2, k1)));
}

// Frozen-configuration Brownian displacement with covariance 2 kBT M dt.
// Uses the Cholesky factor; an element-wise square root of M would be wrong.
export function brownianIncrement(M, kBT, dt, randn) {
  if (!(kBT >= 0 && dt > 0)) throw new RangeError("kBT >= 0, dt > 0 required");
  return applyNoise(noiseFactor(M, "cholesky"), kBT, dt, randn);
}
