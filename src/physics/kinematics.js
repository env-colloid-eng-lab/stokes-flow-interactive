// Tracers in the simple shear u_x = gamma_dot z, run forward and then backward in time,
// optionally with Brownian steps of diffusivity D. Without diffusion the history is undone exactly.
import { mulberry32, normalSampler } from "../core/random.js";

// Returns frames[k] = [[x, z], ...] for k = 0 .. 2 * steps (forward steps, then the same number back).
export function shearReversal(points, { strain = 10, steps = 200, D = 0, seed = 1 } = {}) {
  const dt = strain / steps; // gamma_dot = 1
  const randn = normalSampler(mulberry32(seed));
  const amp = Math.sqrt(2 * D * dt);
  let cur = points.map((p) => p.slice());
  const frames = [cur];
  for (let k = 0; k < 2 * steps; k++) {
    const sign = k < steps ? 1 : -1;
    // exact advection over dt (u_x depends only on z), then the random step
    cur = cur.map(([x, z]) => {
      const xn = x + sign * z * dt;
      return D > 0 ? [xn + amp * randn(), z + amp * randn()] : [xn, z];
    });
    frames.push(cur);
  }
  return frames;
}

// Root-mean-square distance between two point sets of equal length.
export const rmsDistance = (A, B) => Math.sqrt(A.reduce((s, p, i) => s + (p[0] - B[i][0]) ** 2 + (p[1] - B[i][1]) ** 2, 0) / A.length);
