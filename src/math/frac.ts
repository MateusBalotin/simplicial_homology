// src/math/frac.ts

// Frações exatas sobre Q usando BigInt
export class Frac {
  num: bigint;
  den: bigint;

  constructor(num: bigint | number = 0n, den: bigint | number = 1n) {
    let N = typeof num === "number" ? BigInt(num) : (num as bigint);
    let D = typeof den === "number" ? BigInt(den) : (den as bigint);
    if (D === 0n) throw new Error("Denominador zero");
    if (D < 0n) {
      N = -N;
      D = -D;
    }
    const g = Frac._gcd(Frac._abs(N), D);
    this.num = N / g;
    this.den = D / g;
  }

  static from(x: Frac | number | bigint) {
    return x instanceof Frac ? x : new Frac(x, 1n);
  }
  static zero() {
    return new Frac(0n, 1n);
  }
  static one() {
    return new Frac(1n, 1n);
  }
  static _abs(a: bigint) {
    return a < 0n ? -a : a;
  }
  static _gcd(a: bigint, b: bigint) {
    a = Frac._abs(a);
    b = Frac._abs(b);
    while (b !== 0n) {
      const t = a % b;
      a = b;
      b = t;
    }
    return a;
  }

  add(b: Frac | number | bigint) {
    const B = Frac.from(b);
    return new Frac(this.num * B.den + B.num * this.den, this.den * B.den);
  }
  sub(b: Frac | number | bigint) {
    const B = Frac.from(b);
    return new Frac(this.num * B.den - B.num * this.den, this.den * B.den);
  }
  mul(b: Frac | number | bigint) {
    const B = Frac.from(b);
    return new Frac(this.num * B.num, this.den * B.den);
  }
  div(b: Frac | number | bigint) {
    const B = Frac.from(b);
    if (B.num === 0n) throw new Error("/0");
    return new Frac(this.num * B.den, this.den * B.num);
  }
  neg() {
    return new Frac(-this.num, this.den);
  }
  isZero() {
    return this.num === 0n;
  }
  toString() {
    return this.den === 1n ? this.num.toString() : `${this.num}/${this.den}`;
  }
  static abs(x: bigint): bigint {
    return x < 0n ? -x : x;
  }

  static gcd(a: bigint, b: bigint): bigint {
    let A = a < 0n ? -a : a;
    let B = b < 0n ? -b : b;

    while (B !== 0n) {
      const t = A % B;
      A = B;
      B = t;
    }
    return A;
    // ou, se preferir usar o _gcd:
    // return Frac._gcd(Frac.abs(a), Frac.abs(b));
  }
}
