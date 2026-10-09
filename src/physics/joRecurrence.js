// Jeffrey & Onishi (1984) recurrences for all five families of two-sphere resistance functions.
//   XA: axisymmetric translation, eqs. (3.6)-(3.9)
//   YA, YB: transverse translation, eqs. (4.6)-(4.11) with P = V = delta_1n, Q = 0
//   YC: transverse rotation, the same recurrences with P = V = 0, Q = delta_1n, eqs. (7.3)-(7.5)
//   XC: axisymmetric rotation, eq. (6.5)
// f_k(lambda) are the coefficients of the far-field series (3.13)-(3.14), (4.13)-(4.14), (5.3)-(5.4),
// (6.7)-(6.8), (7.7)-(7.8). The exact back-end reproduces the polynomials printed in the paper;
// the Float64 back-end agrees with it to rounding and is used to generate many terms.
import { Rational, Q } from "../core/rational.js";

const exactOps = {
  zero: Q(0), one: Q(1), frac: (n, d) => Q(n, d), int: (n) => new Rational(BigInt(n)),
  add: (a, b) => a.add(b), sub: (a, b) => a.sub(b), mul: (a, b) => a.mul(b),
};
const floatOps = {
  zero: 0, one: 1, frac: (n, d) => n / d, int: (n) => Number(n),
  add: (a, b) => a + b, sub: (a, b) => a - b, mul: (a, b) => a * b,
};

function binomial(n, k) {
  let r = 1n;
  for (let i = 1n; i <= BigInt(k); i++) r = (r * (BigInt(n) - BigInt(k) + i)) / i;
  return r;
}

// Memoised recursive coefficient tables. `init` selects the starting values of the transverse problem.
function tables(T) {
  const memo = () => new Map();
  const key = (n, p, q) => (n * 1024 + p) * 1024 + q;
  const bin = memo(), choose = (n, k) => {
    const kk = n * 1024 + k;
    if (!bin.has(kk)) bin.set(kk, T.int(binomial(n, k)));
    return bin.get(kk);
  };
  const guard = (n, p, q) => n < 1 || p < 0 || q < 0;

  // X^A: P, V
  const xaP = memo(), xaV = memo();
  function PA(n, p, q) {
    if (guard(n, p, q)) return T.zero;
    if (p === 0 && q === 0) return n === 1 ? T.one : T.zero;
    const k = key(n, p, q); if (xaP.has(k)) return xaP.get(k);
    let v = T.zero;
    for (let s = 1; s <= q; s++) {
      const c1 = T.frac(n * (2 * n + 1) * (2 * n * s - n - s + 2), 2 * (n + 1) * (2 * s - 1) * (n + s));
      const c2 = T.frac(n * (2 * n - 1), 2 * (n + 1));
      const c3 = T.frac(n * (4 * n * n - 1), 2 * (n + 1) * (2 * s + 1));
      const inner = T.sub(T.sub(T.mul(c1, PA(s, q - s, p - n + 1)), T.mul(c2, PA(s, q - s, p - n - 1))), T.mul(c3, VA(s, q - s - 2, p - n + 1)));
      v = T.add(v, T.mul(choose(n + s, n), inner));
    }
    xaP.set(k, v); return v;
  }
  function VA(n, p, q) {
    if (guard(n, p, q)) return T.zero;
    if (p === 0 && q === 0) return n === 1 ? T.one : T.zero;
    const k = key(n, p, q); if (xaV.has(k)) return xaV.get(k);
    let v = PA(n, p, q);
    const fac = T.frac(2 * n, (n + 1) * (2 * n + 3));
    for (let s = 1; s <= q; s++) v = T.sub(v, T.mul(T.mul(fac, choose(n + s, n)), PA(s, q - s, p - n - 1)));
    xaV.set(k, v); return v;
  }

  // transverse problems: P, V, Q with starting values chosen by `which` ("A" or "C")
  function transverse(which) {
    const cP = memo(), cV = memo(), cQ = memo();
    const start = (n, w) => (n === 1 && (which === "A" ? w !== "Q" : w === "Q") ? T.one : T.zero);
    function P(n, p, q) {
      if (guard(n, p, q)) return T.zero;
      if (p === 0 && q === 0) return start(n, "P");
      const k = key(n, p, q); if (cP.has(k)) return cP.get(k);
      let v = T.zero;
      for (let s = 1; s <= q; s++) {
        const c1 = T.frac((2 * n + 1) * (3 * (n + s) - (n * s + 1) * (2 * n * s - s - n + 2)), 2 * (n + 1) * s * (n + s) * (2 * s - 1));
        const c2 = T.frac(n * (2 * n - 1), 2 * (n + 1));
        const c3 = T.frac(n * (4 * n * n - 1), 2 * (n + 1) * (2 * s + 1));
        const c4 = T.frac(2 * (4 * n * n - 1), 3 * (n + 1));
        let inner = T.add(T.mul(c1, P(s, q - s, p - n + 1)), T.mul(c2, P(s, q - s, p - n - 1)));
        inner = T.sub(T.add(inner, T.mul(c3, V(s, q - s - 2, p - n + 1))), T.mul(c4, Qf(s, q - s - 1, p - n + 1)));
        v = T.add(v, T.mul(choose(n + s, n + 1), inner));
      }
      cP.set(k, v); return v;
    }
    function V(n, p, q) {
      if (guard(n, p, q)) return T.zero;
      if (p === 0 && q === 0) return start(n, "V");
      const k = key(n, p, q); if (cV.has(k)) return cV.get(k);
      let v = P(n, p, q);
      const fac = T.frac(2 * n, (n + 1) * (2 * n + 3));
      for (let s = 1; s <= q; s++) v = T.add(v, T.mul(T.mul(fac, choose(n + s, n + 1)), P(s, q - s, p - n - 1)));
      cV.set(k, v); return v;
    }
    function Qf(n, p, q) {
      if (guard(n, p, q)) return T.zero;
      if (p === 0 && q === 0) return start(n, "Q");
      const k = key(n, p, q); if (cQ.has(k)) return cQ.get(k);
      let v = T.zero;
      for (let s = 1; s <= q; s++) {
        const inner = T.sub(T.mul(T.frac(s, n + 1), Qf(s, q - s - 1, p - n)), T.mul(T.frac(3, 2 * n * s * (n + 1)), P(s, q - s, p - n)));
        v = T.add(v, T.mul(choose(n + s, n + 1), inner));
      }
      cQ.set(k, v); return v;
    }
    return { P, Q: Qf };
  }

  // X^C: Q only
  const xcQ = memo();
  function QC(n, p, q) {
    if (guard(n, p, q)) return T.zero;
    if (p === 0 && q === 0) return n === 1 ? T.one : T.zero;
    const k = key(n, p, q); if (xcQ.has(k)) return xcQ.get(k);
    let v = T.zero;
    for (let s = 1; s <= q; s++) v = T.add(v, T.mul(T.mul(choose(n + s, n), T.frac(s, n + 1)), QC(s, q - s - 1, p - n)));
    xcQ.set(k, v); return v;
  }
  return { PA, transA: transverse("A"), transC: transverse("C"), QC };
}

/**
 * f_k(lambda) for k = 0..K as polynomial coefficients (increasing powers of lambda).
 * Returns { XA, YA, YB, XC, YC }, each an array of K+1 coefficient arrays.
 */
export function joAllPolynomials(K, { exact = false } = {}) {
  if (!(Number.isInteger(K) && K >= 0)) throw new RangeError("K must be a nonnegative integer");
  const T = exact ? exactOps : floatOps;
  const t = tables(T);
  const poly = (X, k) => Array.from({ length: k + 1 }, (_, q) => X(1, k - q, q));
  const scale = (c, f) => c.map((x) => T.mul(x, f));
  const timesLambda = (c) => [T.zero, ...c];
  const out = { XA: [], YA: [], YB: [], XC: [], YC: [] };
  for (let k = 0; k <= K; k++) {
    const two = T.int(2n ** BigInt(k));
    out.XA.push(scale(poly(t.PA, k), two));                          // (3.15)
    out.YA.push(scale(poly(t.transA.P, k), two));                    // same form for (4.13)-(4.14)
    out.YB.push(scale(poly(t.transA.Q, k), T.mul(two, T.int(2))));   // from (5.2)
    const xc = scale(poly(t.QC, k), two), yc = scale(poly(t.transC.Q, k), two);
    out.XC.push(k % 2 ? timesLambda(xc) : xc);                       // from (6.6)
    out.YC.push(k % 2 ? timesLambda(yc) : yc);                       // from (7.6)
  }
  return out;
}

// f_k at a given lambda (Float64), k = 0..K.
export function joAllValues(K, lambda = 1) {
  const P = joAllPolynomials(K);
  const ev = (c) => c.reduceRight((acc, x) => acc * lambda + x, 0);
  return Object.fromEntries(Object.entries(P).map(([name, rows]) => [name, rows.map(ev)]));
}
