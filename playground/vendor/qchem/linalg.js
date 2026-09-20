// Dense symmetric linear algebra (row-major Float64Array), enough for SCF.
(function (App) {
  "use strict";

  function zeros(n) { return new Float64Array(n); }

  function matmul(A, B, n) {
    var C = new Float64Array(n * n);
    for (var i = 0; i < n; i++) {
      for (var k = 0; k < n; k++) {
        var a = A[i * n + k];
        if (a === 0) continue;
        for (var j = 0; j < n; j++) C[i * n + j] += a * B[k * n + j];
      }
    }
    return C;
  }

  function transpose(A, n) {
    var T = new Float64Array(n * n);
    for (var i = 0; i < n; i++) for (var j = 0; j < n; j++) T[j * n + i] = A[i * n + j];
    return T;
  }

  function trace2(A, B, n) {
    var s = 0;
    for (var i = 0; i < n; i++) for (var j = 0; j < n; j++) s += A[i * n + j] * B[j * n + i];
    return s;
  }

  // Cyclic Jacobi eigensolver for a symmetric matrix.
  // Returns { values: Float64Array (ascending), vectors: columns in row-major }.
  function eighSym(Ain, n) {
    var A = new Float64Array(Ain);
    var V = new Float64Array(n * n);
    for (var i = 0; i < n; i++) V[i * n + i] = 1;

    for (var sweep = 0; sweep < 100; sweep++) {
      var off = 0;
      for (var p = 0; p < n; p++) for (var q = p + 1; q < n; q++) off += A[p * n + q] * A[p * n + q];
      if (Math.sqrt(off) < 1e-12) break;

      for (p = 0; p < n; p++) {
        for (q = p + 1; q < n; q++) {
          var apq = A[p * n + q];
          if (Math.abs(apq) < 1e-15) continue;
          var theta = (A[q * n + q] - A[p * n + p]) / (2 * apq);
          var t = Math.sign(theta) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
          if (theta === 0) t = 1;
          var c = 1 / Math.sqrt(t * t + 1), s = t * c;
          for (var k = 0; k < n; k++) {
            var akp = A[k * n + p], akq = A[k * n + q];
            A[k * n + p] = c * akp - s * akq;
            A[k * n + q] = s * akp + c * akq;
          }
          for (k = 0; k < n; k++) {
            var apk = A[p * n + k], aqk = A[q * n + k];
            A[p * n + k] = c * apk - s * aqk;
            A[q * n + k] = s * apk + c * aqk;
          }
          for (k = 0; k < n; k++) {
            var vkp = V[k * n + p], vkq = V[k * n + q];
            V[k * n + p] = c * vkp - s * vkq;
            V[k * n + q] = s * vkp + c * vkq;
          }
        }
      }
    }

    var order = [];
    for (i = 0; i < n; i++) order.push(i);
    order.sort(function (a, b) { return A[a * n + a] - A[b * n + b]; });

    var values = new Float64Array(n);
    var vectors = new Float64Array(n * n);
    for (var col = 0; col < n; col++) {
      var src = order[col];
      values[col] = A[src * n + src];
      for (i = 0; i < n; i++) vectors[i * n + col] = V[i * n + src];
    }
    return { values: values, vectors: vectors };
  }

  // S^(-1/2) via eigendecomposition; guards against linear dependence.
  function invSqrtSym(S, n) {
    var eig = eighSym(S, n);
    var X = new Float64Array(n * n);
    for (var k = 0; k < n; k++) {
      var lam = eig.values[k];
      if (lam < 1e-10) throw new Error("Overlap matrix is near-singular (linear dependence)");
      var w = 1 / Math.sqrt(lam);
      for (var i = 0; i < n; i++) {
        var vik = eig.vectors[i * n + k] * w;
        for (var j = 0; j < n; j++) X[i * n + j] += vik * eig.vectors[j * n + k];
      }
    }
    return X;
  }

  // Solve small dense linear system (Gauss with partial pivoting) - used by DIIS.
  function solveLin(Ain, bin) {
    var n = bin.length;
    var A = [], b = bin.slice();
    for (var i = 0; i < n; i++) A.push(Ain[i].slice());
    for (var col = 0; col < n; col++) {
      var piv = col;
      for (i = col + 1; i < n; i++) if (Math.abs(A[i][col]) > Math.abs(A[piv][col])) piv = i;
      if (Math.abs(A[piv][col]) < 1e-14) return null;
      var tmp = A[col]; A[col] = A[piv]; A[piv] = tmp;
      var tb = b[col]; b[col] = b[piv]; b[piv] = tb;
      for (i = col + 1; i < n; i++) {
        var f = A[i][col] / A[col][col];
        for (var j = col; j < n; j++) A[i][j] -= f * A[col][j];
        b[i] -= f * b[col];
      }
    }
    var x = new Array(n).fill(0);
    for (i = n - 1; i >= 0; i--) {
      var s = b[i];
      for (j = i + 1; j < n; j++) s -= A[i][j] * x[j];
      x[i] = s / A[i][i];
    }
    return x;
  }

  App.linalg = {
    zeros: zeros, matmul: matmul, transpose: transpose,
    trace2: trace2, eighSym: eighSym, invSqrtSym: invSqrtSym, solveLin: solveLin
  };
})(typeof globalThis.App === "object" ? globalThis.App : (globalThis.App = {}));
