// Two-electron repulsion integrals (ij|kl) with 8-fold symmetric packed storage
// and Schwarz screening. Same MD machinery as integrals.js.
(function (App) {
  "use strict";

  var TWO_PI_POW = 2 * Math.pow(Math.PI, 2.5);
  var SCREEN_TOL = 1e-10;

  function pairIndex(i, j) { return i >= j ? i * (i + 1) / 2 + j : j * (j + 1) / 2 + i; }

  // Precompute primitive-pair data for every contracted pair i >= j.
  function buildPairs(basis) {
    var nb = basis.length, pairs = [];
    for (var i = 0; i < nb; i++) {
      for (var j = 0; j <= i; j++) {
        var fi = basis[i], fj = basis[j];
        var AB = [fi.center[0] - fj.center[0], fi.center[1] - fj.center[1], fi.center[2] - fj.center[2]];
        var L = [fi.l + fj.l, fi.m + fj.m, fi.n + fj.n];
        var prims = [];
        for (var pi = 0; pi < fi.exps.length; pi++) {
          for (var pj = 0; pj < fj.exps.length; pj++) {
            var a = fi.exps[pi], b = fj.exps[pj], p = a + b;
            var Ex = App.integrals.makeE(a, b, AB[0]);
            var Ey = App.integrals.makeE(a, b, AB[1]);
            var Ez = App.integrals.makeE(a, b, AB[2]);
            var ex = [], ey = [], ez = [];
            for (var t = 0; t <= L[0]; t++) ex.push(Ex(fi.l, fj.l, t));
            for (t = 0; t <= L[1]; t++) ey.push(Ey(fi.m, fj.m, t));
            for (t = 0; t <= L[2]; t++) ez.push(Ez(fi.n, fj.n, t));
            prims.push({
              p: p,
              P: [(a * fi.center[0] + b * fj.center[0]) / p,
                  (a * fi.center[1] + b * fj.center[1]) / p,
                  (a * fi.center[2] + b * fj.center[2]) / p],
              c: fi.coefs[pi] * fj.coefs[pj],
              ex: ex, ey: ey, ez: ez
            });
          }
        }
        pairs.push({ prims: prims, L: L, Q: 0 });
      }
    }
    return pairs;
  }

  // Contracted quartet (bra|ket) from precomputed pair data.
  function quartet(bra, ket) {
    var tmax = bra.L[0] + bra.L[1] + bra.L[2] + ket.L[0] + ket.L[1] + ket.L[2];
    var g = 0;
    for (var bi = 0; bi < bra.prims.length; bi++) {
      var B = bra.prims[bi];
      for (var ki = 0; ki < ket.prims.length; ki++) {
        var K = ket.prims[ki];
        var p = B.p, q = K.p;
        var alpha = p * q / (p + q);
        var R = App.integrals.makeR(alpha, B.P[0] - K.P[0], B.P[1] - K.P[1], B.P[2] - K.P[2], tmax);
        var sum = 0;
        for (var t = 0; t <= bra.L[0]; t++) {
          if (B.ex[t] === 0) continue;
          for (var u = 0; u <= bra.L[1]; u++) {
            if (B.ey[u] === 0) continue;
            for (var v = 0; v <= bra.L[2]; v++) {
              var eb = B.ex[t] * B.ey[u] * B.ez[v];
              if (eb === 0) continue;
              for (var tt = 0; tt <= ket.L[0]; tt++) {
                if (K.ex[tt] === 0) continue;
                for (var uu = 0; uu <= ket.L[1]; uu++) {
                  if (K.ey[uu] === 0) continue;
                  for (var vv = 0; vv <= ket.L[2]; vv++) {
                    var ek = K.ex[tt] * K.ey[uu] * K.ez[vv];
                    if (ek === 0) continue;
                    var sign = ((tt + uu + vv) % 2 === 0) ? 1 : -1;
                    sum += eb * ek * sign * R(t + tt, u + uu, v + vv, 0);
                  }
                }
              }
            }
          }
        }
        g += B.c * K.c * TWO_PI_POW / (p * q * Math.sqrt(p + q)) * sum;
      }
    }
    return g;
  }

  // Full ERI tensor in 8-fold packed storage.
  function computeERI(basis, onProgress) {
    var nb = basis.length;
    var pairs = buildPairs(basis);
    var np = pairs.length;

    for (var ij = 0; ij < np; ij++) pairs[ij].Q = Math.sqrt(Math.abs(quartet(pairs[ij], pairs[ij])));

    var data = new Float64Array(np * (np + 1) / 2);
    var done = 0, total = np * (np + 1) / 2;
    for (ij = 0; ij < np; ij++) {
      for (var kl = 0; kl <= ij; kl++) {
        if (pairs[ij].Q * pairs[kl].Q > SCREEN_TOL) {
          data[ij * (ij + 1) / 2 + kl] = quartet(pairs[ij], pairs[kl]);
        }
        done++;
      }
      if (onProgress && (ij & 31) === 0) onProgress(done / total);
    }
    return { data: data, nb: nb };
  }

  function getERI(eri, i, j, k, l) {
    var ij = pairIndex(i, j), kl = pairIndex(k, l);
    return ij >= kl ? eri.data[ij * (ij + 1) / 2 + kl] : eri.data[kl * (kl + 1) / 2 + ij];
  }

  App.eri = { computeERI: computeERI, get: getERI, pairIndex: pairIndex, buildPairs: buildPairs };
})(typeof globalThis.App === "object" ? globalThis.App : (globalThis.App = {}));
