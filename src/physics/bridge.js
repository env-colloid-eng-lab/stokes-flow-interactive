// Low-order reflections -> JO coefficients f0..f4 (BridgeLab.jl, section 13.3).
// Inverting the far-field mobility POWER SERIES gives the resistance series.
// Removing the stresslet term changes f4. This is not an all-order generator.
import { Rational, Q } from "../core/rational.js";

const R0 = () => [[Q(0), Q(0)], [Q(0), Q(0)]];
const madd = (A, B) => A.map((row, i) => row.map((a, j) => a.add(B[i][j])));
const mmul = (A, B) => A.map((row, i) => [0, 1].map((j) => row[0].mul(B[0][j]).add(row[1].mul(B[1][j]))));

// lambda: radius ratio as a Rational (exact).
export function reflectionCoefficients(lambda = Q(1), { stresslet = true } = {}) {
  lambda = Rational.of(lambda);
  if (!(lambda.num > 0n)) throw new RangeError("positive radius ratio required");
  const h = Q(1, 2);
  const M = [0, 1, 2, 3, 4].map(R0);
  M[0] = [[Q(1), Q(0)], [Q(0), Q(1).div(lambda)]];
  M[1] = [[Q(0), Q(3).mul(h)], [Q(3).mul(h), Q(0)]];
  const c3 = Q(1).add(lambda.mul(lambda)).mul(h).neg();
  M[3] = [[Q(0), c3], [c3, Q(0)]];
  if (stresslet) {
    const m = Q(-15, 4);
    M[4] = [[m.mul(lambda.mul(lambda).mul(lambda)), Q(0)], [Q(0), m]];
  }
  const R = [0, 1, 2, 3, 4].map(R0);
  R[0] = [[Q(1), Q(0)], [Q(0), lambda]];
  for (let k = 1; k <= 4; k++) {
    let product = R0();
    for (let j = 1; j <= k; j++) product = madd(product, mmul(M[j], R[k - j]));
    R[k] = mmul(R[0], product).map((row) => row.map((x) => x.neg()));
  }
  return [0, 1, 2, 3, 4].map((k) =>
    k % 2 === 0 ? Q(2 ** k).mul(R[k][0][0]) : Q(2 ** k).mul(R[k][0][1]).neg());
}

// Far-field axial mobility keeping the r^-4 self correction. Positive definiteness of this
// truncation is NOT guaranteed near contact.
export function reflectedMobility(a1, a2, r, { mu = 1, stresslet = true } = {}) {
  if (!(Math.min(a1, a2, mu) > 0 && r > a1 + a2)) throw new RangeError("separated positive spheres required");
  const c = (1 - (a1 * a1 + a2 * a2) / (3 * r * r)) / (4 * Math.PI * mu * r);
  let m1 = 1 / (6 * Math.PI * mu * a1), m2 = 1 / (6 * Math.PI * mu * a2);
  if (stresslet) {
    m1 -= (5 * a2 ** 3) / (8 * Math.PI * mu * r ** 4);
    m2 -= (5 * a1 ** 3) / (8 * Math.PI * mu * r ** 4);
  }
  return [[m1, c], [c, m2]];
}
