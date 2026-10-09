// Symmetric traceless (STF) 3x3 tensors as 5-vectors. Port of StokesLab.jl (stf_basis, stf5, unstf5).
// The basis is orthonormal under A:B = sum_ij A_ij B_ij, so S:E = s . e.

export function stfBasis() {
  const z = () => [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  const B = [z(), z(), z(), z(), z()];
  const r2 = Math.SQRT1_2, r6 = 1 / Math.sqrt(6);
  B[0][0][0] = r2; B[0][1][1] = -r2;                      // diag(1,-1,0)/sqrt2
  B[1][0][0] = r6; B[1][1][1] = r6; B[1][2][2] = -2 * r6;  // diag(1,1,-2)/sqrt6
  B[2][0][1] = B[2][1][0] = r2;                            // xy
  B[3][0][2] = B[3][2][0] = r2;                            // xz
  B[4][1][2] = B[4][2][1] = r2;                            // yz
  return B;
}

const BASIS = stfBasis();
export const ddot = (A, B) => A.reduce((s, row, i) => s + row.reduce((t, a, j) => t + a * B[i][j], 0), 0);

// Symmetric traceless part (A + A^T)/2 - tr(A)/3 I.
export function stf(A) {
  const tr = (A[0][0] + A[1][1] + A[2][2]) / 3;
  return A.map((row, i) => row.map((a, j) => (a + A[j][i]) / 2 - (i === j ? tr : 0)));
}

export const stf5 = (A) => BASIS.map((B) => ddot(A, B));
export const unstf5 = (v) => [0, 1, 2].map((i) => [0, 1, 2].map((j) => v.reduce((s, vk, k) => s + vk * BASIS[k][i][j], 0)));

// The "naive" 5 numbers (Exx, Eyy, Exy, Exz, Eyz): their plain dot product is NOT E:S.
export const naive5 = (A) => [A[0][0], A[1][1], A[0][1], A[0][2], A[1][2]];
