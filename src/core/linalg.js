// Small dense linear algebra for teaching-scale problems (up to a few hundred unknowns).
// Matrices are arrays of rows (number[][]); vectors are number[].

export function zeros(m, n = m) {
  return Array.from({ length: m }, () => new Array(n).fill(0));
}

export function eye(n) {
  const A = zeros(n);
  for (let i = 0; i < n; i++) A[i][i] = 1;
  return A;
}

export const clone = (A) => A.map((row) => row.slice());

export function transpose(A) {
  const m = A.length, n = A[0].length;
  const T = zeros(n, m);
  for (let i = 0; i < m; i++) for (let j = 0; j < n; j++) T[j][i] = A[i][j];
  return T;
}

export function matmul(A, B) {
  const m = A.length, k = B.length, n = B[0].length;
  if (A[0].length !== k) throw new RangeError("matmul: inner dimensions differ");
  const C = zeros(m, n);
  for (let i = 0; i < m; i++) {
    const Ai = A[i], Ci = C[i];
    for (let p = 0; p < k; p++) {
      const a = Ai[p];
      if (a === 0) continue;
      const Bp = B[p];
      for (let j = 0; j < n; j++) Ci[j] += a * Bp[j];
    }
  }
  return C;
}

export function matvec(A, x) {
  if (A[0].length !== x.length) throw new RangeError("matvec: dimensions differ");
  return A.map((row) => row.reduce((s, a, j) => s + a * x[j], 0));
}

export const add = (A, B) => A.map((row, i) => row.map((a, j) => a + B[i][j]));
export const sub = (A, B) => A.map((row, i) => row.map((a, j) => a - B[i][j]));
export const scale = (A, c) => A.map((row) => row.map((a) => a * c));
export const outer = (u, v) => u.map((ui) => v.map((vj) => ui * vj));

export const dot = (u, v) => u.reduce((s, ui, i) => s + ui * v[i], 0);
export const norm = (v) => Math.sqrt(dot(v, v));
export const vadd = (u, v) => u.map((ui, i) => ui + v[i]);
export const vsub = (u, v) => u.map((ui, i) => ui - v[i]);
export const vscale = (u, c) => u.map((ui) => ui * c);
export const cross = (u, v) => [
  u[1] * v[2] - u[2] * v[1],
  u[2] * v[0] - u[0] * v[2],
  u[0] * v[1] - u[1] * v[0],
];

// Frobenius norm.
export const fnorm = (A) => Math.sqrt(A.reduce((s, row) => s + row.reduce((t, a) => t + a * a, 0), 0));

// ||A - A^T||_F / ||A||_F. Check this BEFORE symmetrising anything.
export function asymmetry(A) {
  const nA = fnorm(A);
  return nA === 0 ? 0 : fnorm(sub(A, transpose(A))) / nA;
}

export function getBlock(A, r0, c0, m, n = m) {
  const B = zeros(m, n);
  for (let i = 0; i < m; i++) for (let j = 0; j < n; j++) B[i][j] = A[r0 + i][c0 + j];
  return B;
}

export function setBlock(A, r0, c0, B) {
  for (let i = 0; i < B.length; i++) for (let j = 0; j < B[0].length; j++) A[r0 + i][c0 + j] = B[i][j];
  return A;
}

// LU with partial pivoting. Returns { LU, piv } or throws for a singular matrix.
export function lu(A) {
  const n = A.length;
  const LU = clone(A);
  const piv = Array.from({ length: n }, (_, i) => i);
  for (let k = 0; k < n; k++) {
    let p = k, big = Math.abs(LU[k][k]);
    for (let i = k + 1; i < n; i++) {
      const v = Math.abs(LU[i][k]);
      if (v > big) { big = v; p = i; }
    }
    if (big === 0) throw new RangeError("lu: matrix is singular");
    if (p !== k) {
      [LU[p], LU[k]] = [LU[k], LU[p]];
      [piv[p], piv[k]] = [piv[k], piv[p]];
    }
    for (let i = k + 1; i < n; i++) {
      const f = (LU[i][k] /= LU[k][k]);
      if (f === 0) continue;
      for (let j = k + 1; j < n; j++) LU[i][j] -= f * LU[k][j];
    }
  }
  return { LU, piv };
}

function luSolveVec({ LU, piv }, b) {
  const n = LU.length;
  const x = piv.map((p) => b[p]);
  for (let i = 0; i < n; i++) for (let j = 0; j < i; j++) x[i] -= LU[i][j] * x[j];
  for (let i = n - 1; i >= 0; i--) {
    for (let j = i + 1; j < n; j++) x[i] -= LU[i][j] * x[j];
    x[i] /= LU[i][i];
  }
  return x;
}

// Solve A x = b (b a vector) or A X = B (B a matrix), like Julia's A \ b.
export function solve(A, b) {
  const f = lu(A);
  if (!Array.isArray(b[0])) return luSolveVec(f, b);
  const cols = transpose(b).map((col) => luSolveVec(f, col));
  return transpose(cols);
}

export const inv = (A) => solve(A, eye(A.length));

// Cholesky A = L L^T for symmetric positive definite A. Only the lower triangle is read.
export function cholesky(A) {
  const n = A.length;
  const L = zeros(n);
  for (let j = 0; j < n; j++) {
    let d = A[j][j];
    for (let k = 0; k < j; k++) d -= L[j][k] * L[j][k];
    if (!(d > 0)) throw new RangeError("cholesky: matrix is not positive definite");
    L[j][j] = Math.sqrt(d);
    for (let i = j + 1; i < n; i++) {
      let s = A[i][j];
      for (let k = 0; k < j; k++) s -= L[i][k] * L[j][k];
      L[i][j] = s / L[j][j];
    }
  }
  return L;
}

// Solve A X = B given the Cholesky factor L of A (B a vector or a matrix).
export function choleskySolve(L, B) {
  const n = L.length;
  const one = (b) => {
    const y = b.slice();
    for (let i = 0; i < n; i++) { for (let k = 0; k < i; k++) y[i] -= L[i][k] * y[k]; y[i] /= L[i][i]; }
    for (let i = n - 1; i >= 0; i--) { for (let k = i + 1; k < n; k++) y[i] -= L[k][i] * y[k]; y[i] /= L[i][i]; }
    return y;
  };
  return Array.isArray(B[0]) ? transpose(transpose(B).map(one)) : one(B);
}

export function isPosDef(A) {
  try { cholesky(A); return true; } catch { return false; }
}

// Least squares min ||A X - B|| by Householder QR (A is m x n, m >= n, full column rank).
export function lstsq(A, B) {
  const m = A.length, n = A[0].length;
  if (m < n) throw new RangeError("lstsq: need m >= n");
  const vecRhs = !Array.isArray(B[0]);
  const R = clone(A);
  const Y = vecRhs ? B.map((b) => [b]) : clone(B);
  const k = Y[0].length;
  for (let j = 0; j < n; j++) {
    let s = 0;
    for (let i = j; i < m; i++) s += R[i][j] * R[i][j];
    const alpha = R[j][j] > 0 ? -Math.sqrt(s) : Math.sqrt(s);
    if (alpha === 0) throw new RangeError("lstsq: rank deficient");
    const v = new Array(m - j);
    for (let i = j; i < m; i++) v[i - j] = R[i][j];
    v[0] -= alpha;
    const vv = v.reduce((t, x) => t + x * x, 0);
    if (vv > 0) {
      for (let c = j; c < n; c++) {
        let t = 0;
        for (let i = j; i < m; i++) t += v[i - j] * R[i][c];
        t = (2 * t) / vv;
        for (let i = j; i < m; i++) R[i][c] -= t * v[i - j];
      }
      for (let c = 0; c < k; c++) {
        let t = 0;
        for (let i = j; i < m; i++) t += v[i - j] * Y[i][c];
        t = (2 * t) / vv;
        for (let i = j; i < m; i++) Y[i][c] -= t * v[i - j];
      }
    }
  }
  const X = zeros(n, k);
  for (let c = 0; c < k; c++) {
    for (let i = n - 1; i >= 0; i--) {
      let s = Y[i][c];
      for (let j = i + 1; j < n; j++) s -= R[i][j] * X[j][c];
      X[i][c] = s / R[i][i];
    }
  }
  return vecRhs ? X.map((row) => row[0]) : X;
}

// Symmetric eigenproblem by cyclic Jacobi rotations. Reads A as symmetric.
// Returns eigenvalues in ascending order and eigenvectors as COLUMNS of `vectors`.
export function symEig(A, { tol = 1e-14, maxSweeps = 100 } = {}) {
  const n = A.length;
  const S = clone(A);
  const V = eye(n);
  const scaleA = fnorm(A) || 1;
  for (let sweep = 0; sweep < maxSweeps; sweep++) {
    let off = 0;
    for (let p = 0; p < n; p++) for (let q = p + 1; q < n; q++) off += S[p][q] * S[p][q];
    if (Math.sqrt(off) <= tol * scaleA) break;
    for (let p = 0; p < n; p++) {
      for (let q = p + 1; q < n; q++) {
        const apq = S[p][q];
        if (apq === 0) continue;
        const theta = (S[q][q] - S[p][p]) / (2 * apq);
        const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
        const c = 1 / Math.sqrt(t * t + 1), s = t * c;
        for (let k = 0; k < n; k++) {
          const skp = S[k][p], skq = S[k][q];
          S[k][p] = c * skp - s * skq;
          S[k][q] = s * skp + c * skq;
        }
        for (let k = 0; k < n; k++) {
          const spk = S[p][k], sqk = S[q][k];
          S[p][k] = c * spk - s * sqk;
          S[q][k] = s * spk + c * sqk;
        }
        for (let k = 0; k < n; k++) {
          const vkp = V[k][p], vkq = V[k][q];
          V[k][p] = c * vkp - s * vkq;
          V[k][q] = s * vkp + c * vkq;
        }
      }
    }
  }
  const order = Array.from({ length: n }, (_, i) => i).sort((i, j) => S[i][i] - S[j][j]);
  return {
    values: order.map((i) => S[i][i]),
    vectors: V.map((row) => order.map((i) => row[i])),
  };
}

export const eigvalsSym = (A) => symEig(A).values;

// Singular values (descending) by one-sided Jacobi rotations on the columns of A (m >= n).
// Works on A directly, so the accuracy is not lost by squaring the condition number.
export function singularValues(A, { tol = 1e-15, maxSweeps = 60 } = {}) {
  const U = transpose(A); // rows of U are the columns of A
  const n = U.length;
  for (let sweep = 0; sweep < maxSweeps; sweep++) {
    let rotated = false;
    for (let p = 0; p < n - 1; p++)
      for (let q = p + 1; q < n; q++) {
        const up = U[p], uq = U[q];
        let alpha = 0, beta = 0, gamma = 0;
        for (let i = 0; i < up.length; i++) { alpha += up[i] * up[i]; beta += uq[i] * uq[i]; gamma += up[i] * uq[i]; }
        if (Math.abs(gamma) <= tol * Math.sqrt(alpha * beta)) continue;
        rotated = true;
        const zeta = (beta - alpha) / (2 * gamma);
        const t = Math.sign(zeta || 1) / (Math.abs(zeta) + Math.sqrt(1 + zeta * zeta));
        const c = 1 / Math.sqrt(1 + t * t), s = c * t;
        for (let i = 0; i < up.length; i++) {
          const x = up[i], y = uq[i];
          up[i] = c * x - s * y;
          uq[i] = s * x + c * y;
        }
      }
    if (!rotated) break;
  }
  return U.map(norm).sort((x, y) => y - x);
}

// 2-norm condition number sigma_max / sigma_min (as Julia's cond).
export function cond(A) {
  const sv = singularValues(A.length >= A[0].length ? A : transpose(A));
  return sv[0] / sv[sv.length - 1];
}
