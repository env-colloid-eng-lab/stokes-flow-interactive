// Exact rational arithmetic on BigInt, enough for generating series coefficients.

const babs = (x) => (x < 0n ? -x : x);
function gcd(a, b) {
  a = babs(a); b = babs(b);
  while (b) [a, b] = [b, a % b];
  return a;
}

export class Rational {
  constructor(num, den = 1n) {
    num = BigInt(num); den = BigInt(den);
    if (den === 0n) throw new RangeError("zero denominator");
    if (den < 0n) { num = -num; den = -den; }
    const g = gcd(num, den) || 1n;
    this.num = num / g;
    this.den = den / g;
  }
  static of(x) { return x instanceof Rational ? x : new Rational(BigInt(x)); }
  add(o) { o = Rational.of(o); return new Rational(this.num * o.den + o.num * this.den, this.den * o.den); }
  sub(o) { o = Rational.of(o); return new Rational(this.num * o.den - o.num * this.den, this.den * o.den); }
  mul(o) { o = Rational.of(o); return new Rational(this.num * o.num, this.den * o.den); }
  div(o) { o = Rational.of(o); return new Rational(this.num * o.den, this.den * o.num); }
  neg() { return new Rational(-this.num, this.den); }
  isZero() { return this.num === 0n; }
  equals(o) { o = Rational.of(o); return this.num === o.num && this.den === o.den; }
  toNumber() {
    // Keep the leading ~64 bits of numerator and denominator separately, then rescale,
    // so that a small numerator or denominator is never truncated to zero.
    const bits = (x) => babs(x).toString(2).length;
    const sn = Math.max(0, bits(this.num) - 64), sd = Math.max(0, bits(this.den) - 64);
    const q = Number(this.num >> BigInt(sn)) / Number(this.den >> BigInt(sd));
    const e = sn - sd;
    // apply 2^e in two halves to avoid intermediate overflow / underflow
    return q * 2 ** Math.trunc(e / 2) * 2 ** (e - Math.trunc(e / 2));
  }
  toString() { return this.den === 1n ? `${this.num}` : `${this.num}/${this.den}`; }
}

export const Q = (n, d = 1) => new Rational(BigInt(n), BigInt(d));
