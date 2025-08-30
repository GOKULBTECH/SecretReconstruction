#!/usr/bin/env node
// reconstruct.js
// Shamir-style interpolation with exact BigInt fractions (no external libs).
// Works with: node reconstruct.js testcase.json  OR  node reconstruct.js < testcase.json

console.log("Script started"); // Debugging line

const fs = require("fs");

// Make sure we got a filename
if (process.argv.length < 3) {
  console.error("Usage: node reconstruct.js <input.json>");
  process.exit(1);
}

// Read file
const filename = process.argv[2];
const data = fs.readFileSync(filename, "utf8");
console.log("File contents:", data); // Debugging line

// Parse JSON
let input;
try {
  input = JSON.parse(data);
  console.log("Parsed JSON:", input); // Debugging line
} catch (err) {
  console.error("Error parsing JSON:", err.message);
  process.exit(1);
}

// TODO: your reconstruction logic here
console.log("Final output goes here...");


// --------- Read Input ----------
let raw = "";
if (process.argv.length > 2) {
  // If a filename is provided, read from it
  raw = fs.readFileSync(process.argv[2], "utf8").trim();
} else {
  // Otherwise, read from stdin
  raw = fs.readFileSync(0, "utf8").trim();
}

if (!raw) {
  console.error("❌ No JSON input found. Usage:");
  console.error("   node reconstruct.js testcase.json");
  console.error("   OR");
  console.error("   node reconstruct.js < testcase.json");
  process.exit(1);
}

const obj = JSON.parse(raw);

// --------- Utilities: BigInt GCD & Fraction ----------
function bigAbs(a) { return a < 0n ? -a : a; }
function bigGcd(a, b) {
  a = bigAbs(a); b = bigAbs(b);
  while (b) { const t = a % b; a = b; b = t; }
  return a;
}
class Frac {
  constructor(num, den = 1n) {
    if (den === 0n) throw new Error('Division by zero');
    // normalize sign
    if (den < 0n) { num = -num; den = -den; }
    const g = bigGcd(num, den);
    this.n = num / g;
    this.d = den / g;
  }
  static fromBig(n) { return new Frac(n, 1n); }
  add(o) { return new Frac(this.n * o.d + o.n * this.d, this.d * o.d); }
  sub(o) { return new Frac(this.n * o.d - o.n * this.d, this.d * o.d); }
  mul(o) { return new Frac(this.n * o.n, this.d * o.d); }
  div(o) {
    if (o.n === 0n) throw new Error('Division by zero');
    return new Frac(this.n * o.d, this.d * o.n);
  }
  eqBig(b) { return this.d === 1n && this.n === b; }
  toBig() {
    if (this.d !== 1n) throw new Error('Not an integer result');
    return this.n;
  }
}

// --------- Base conversion to BigInt ----------
function charVal(c) {
  if (c >= '0' && c <= '9') return BigInt(c.charCodeAt(0) - 48);
  const lower = c.toLowerCase();
  if (lower >= 'a' && lower <= 'z') return BigInt(10 + (lower.charCodeAt(0) - 97));
  throw new Error(`Invalid digit '${c}'`);
}
function parseInBaseToBigInt(str, base) {
  const B = BigInt(base);
  let v = 0n;
  for (const ch of str) {
    const d = charVal(ch);
    if (d >= B) throw new Error(`Digit '${ch}' not valid for base ${base}`);
    v = v * B + d;
  }
  return v;
}

// --------- Lagrange evaluation f(x0) from a subset of k points ----------
function lagrangeEvalAtX0(points, x0) {
  // points: array of {x: BigInt, y: BigInt}. x0: BigInt
  let acc = new Frac(0n, 1n);
  for (let i = 0; i < points.length; i++) {
    const xi = points[i].x, yi = points[i].y;
    let num = new Frac(1n, 1n);
    let den = new Frac(1n, 1n);
    for (let j = 0; j < points.length; j++) {
      if (i === j) continue;
      num = num.mul(new Frac(x0 - points[j].x, 1n));
      den = den.mul(new Frac(xi - points[j].x, 1n));
    }
    const Li = num.div(den);
    acc = acc.add(Frac.fromBig(yi).mul(Li));
  }
  return acc;
}

// --------- Generate k-combinations ----------
function* kCombinations(arr, k, start = 0, prefix = []) {
  if (prefix.length === k) { yield prefix.slice(); return; }
  for (let i = start; i < arr.length; i++) {
    prefix.push(arr[i]);
    yield* kCombinations(arr, k, i + 1, prefix);
    prefix.pop();
  }
}

// --------- Find best polynomial model ----------
function findBestModel(points, k) {
  let best = null;
  for (const subset of kCombinations(points, k)) {
    const mismatches = [];
    for (const p of points) {
      const yFrac = lagrangeEvalAtX0(subset, p.x);
      try {
        const yPred = yFrac.toBig();
        if (yPred !== p.y) mismatches.push({ x: p.x, given: p.y, expected: yPred });
      } catch {
        mismatches.push({ x: p.x, given: p.y, expected: null });
      }
    }
    const score = points.length - mismatches.length;
    if (!best || score > best.score) {
      const secretFrac = lagrangeEvalAtX0(subset, 0n);
      best = { subset, mismatches, score, secretFrac };
      if (mismatches.length === 0) break; // perfect fit
    }
  }
  return best;
}

// --------- Main ----------
function main() {
  const n = obj.keys.n;
  const k = obj.keys.k;

  const points = [];
  for (const key of Object.keys(obj)) {
    if (key === 'keys') continue;
    const x = BigInt(key);
    const base = parseInt(obj[key].base, 10);
    const y = parseInBaseToBigInt(obj[key].value, base);
    points.push({ x, y });
  }
  points.sort((a, b) => (a.x < b.x ? -1 : a.x > b.x ? 1 : 0));

  const best = findBestModel(points, k);
  if (!best) { console.error('No valid polynomial found'); process.exit(1); }

  let secretStr;
  try {
    secretStr = best.secretFrac.toBig().toString();
  } catch {
    secretStr = `${best.secretFrac.n.toString()}/${best.secretFrac.d.toString()}`;
  }

  const wrong = best.mismatches.map(m => ({
    index: m.x.toString(),
    given: m.given.toString(),
    expected: m.expected === null ? 'non-integer (inconsistent)' : m.expected.toString()
  }));

  console.log(JSON.stringify({
    n, k,
    secret: secretStr,
    wrongShares: wrong
  }, null, 2));
}

if (require.main === module) main();
