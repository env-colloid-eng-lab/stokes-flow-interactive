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
    // Scale to keep precision when num and den are both huge.
    const shift = BigInt(Math.max(0, Math.max(this.num.toString(2).length, this.den.toString(2).length) - 1000));
    return Number(this.num >> shift) / Number(this.den >> shift);
  }
  toString() { return this.den === 1n ? `${this.num}` : `${this.num}/${this.den}`; }
}

export const Q = (n, d = 1) => new Rational(BigInt(n), BigInt(d));
