// Head-on approach of two equal spheres under opposite forces +F, -F along the line of centres.
// The gap obeys dh/dt = -V(h), V = 2F / (6 pi mu a (X11 - X12)).

/**
 * rel(h) -> { value: X11 - X12, ok: boolean }   (ok = false when the model is not trustworthy at h)
 * Options:
 *   contactBelow  stop and report "contact" once h falls below this value (models that do reach h = 0);
 *                 leave undefined for models whose gap only decays (e.g. with lubrication)
 * Returns { ts, hs, stop } with stop in { null, "contact", "series" }.
 */
export function approachTrajectory(rel, { h0, Tmax, F = 1, a = 1, mu = 1, contactBelow } = {}) {
  const unit = 6 * Math.PI * mu * a;
  const ts = [0], hs = [h0];
  let t = 0, hh = h0, stop = null;
  // velocity at h, or null when the model cannot be used there
  const V = (h) => {
    if (!(h > 0)) return null;
    const r = rel(h);
    return r.ok ? (2 * F) / (unit * r.value) : null;
  };
  const halt = (h) => (h > 0 ? "series" : "contact");
  while (t < Tmax) {
    const v1 = V(hh);
    if (v1 === null) { stop = halt(hh); break; }
    const dt = Math.min(0.5, (0.02 * hh) / v1, Tmax - t);
    const h2 = hh - (dt / 2) * v1, v2 = V(h2);
    if (v2 === null) { stop = halt(h2); break; }
    const h3 = hh - (dt / 2) * v2, v3 = V(h3);
    if (v3 === null) { stop = halt(h3); break; }
    const h4 = hh - dt * v3, v4 = V(h4);
    if (v4 === null) { stop = halt(h4); break; }
    hh -= (dt / 6) * (v1 + 2 * v2 + 2 * v3 + v4);
    t += dt;
    if (contactBelow !== undefined && hh < contactBelow) {
      ts.push(t); hs.push(Math.max(hh, contactBelow)); stop = "contact"; break;
    }
    if (!(hh > 0)) { stop = "contact"; break; }
    ts.push(t); hs.push(hh);
  }
  return { ts, hs, stop };
}
