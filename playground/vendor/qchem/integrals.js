// One-electron integrals over contracted cartesian Gaussians (s, p):
// overlap S, kinetic T, nuclear attraction V. McMurchie-Davidson scheme.
// Reference: Helgaker, Jorgensen, Olsen, "Molecular Electronic-Structure Theory", ch. 9.
(function (App) {
  "use strict";

  // Boys function F_m(x) for m = 0..mmax, downward recursion from a series value.
  function boys(mmax, x) {
    var F = new Float64Array(mmax + 1);
    var m;
    if (x < 1e-12) {
      for (m = 0; m <= mmax; m++) F[m] = 1 / (2 * m + 1);
      return F;
    }
    if (x > 35) {
      // asymptotic: F_0 = sqrt(pi/x)/2, upward via F_{m+1} = ((2m+1)F_m - e^-x)/(2x); e^-x negligible
      F[0] = 0.5 * Math.sqrt(Math.PI / x);
      for (m = 0; m < mmax; m++) F[m + 1] = (2 * m + 1) * F[m] / (2 * x);
      return F;
    }
    // series at m = mmax, then stable downward recursion
    var term = 1 / (2 * mmax + 1), sum = term, k = 1;
    while (term > 1e-17 * sum || k < 5) {
      term *= 2 * x / (2 * mmax + 2 * k + 1);
      sum += term;
      if (++k > 200) break;
    }
    var ex = Math.exp(-x);
    F[mmax] = sum * ex;
    for (m = mmax - 1; m >= 0; m--) F[m] = (2 * x * F[m + 1] + ex) / (2 * m + 1);
    return F;
  }

  // Hermite expansion coefficients E_t^{i,j} for one dimension.
  // Returns table Et(i, j, t) flattened into a closure-backed function.
  function makeE(a, b, XAB) {
    var p = a + b;
    var XPA = -b * XAB / p, XPB = a * XAB / p;
    var K = Math.exp(-a * b * XAB * XAB / p);
    var memo = {};
    function E(i, j, t) {
      if (t < 0 || t > i + j) return 0;
      if (i === 0 && j === 0) return t === 0 ? K : 0;
      var key = i * 100 + j * 10 + t;
      if (key in memo) return memo[key];
      var v;
      if (i > 0) {
        v = E(i - 1, j, t - 1) / (2 * p) + XPA * E(i - 1, j, t) + (t + 1) * E(i - 1, j, t + 1);
      } else {
        v = E(i, j - 1, t - 1) / (2 * p) + XPB * E(i, j - 1, t) + (t + 1) * E(i, j - 1, t + 1);
      }
      memo[key] = v;
      return v;
    }
    return E;
  }

  // Hermite Coulomb integrals R_{tuv} at order n=0, built on Boys values.
  function makeR(alpha, X, Y, Z, tmax) {
    var R2 = X * X + Y * Y + Z * Z;
    var F = boys(tmax, alpha * R2);
    var base = new Float64Array(tmax + 1);
    for (var n = 0; n <= tmax; n++) base[n] = Math.pow(-2 * alpha, n) * F[n];
    var memo = {};
    function R(t, u, v, n) {
      if (t < 0 || u < 0 || v < 0) return 0;
      if (t === 0 && u === 0 && v === 0) return base[n];
      var key = ((t * 16 + u) * 16 + v) * 16 + n;
      if (key in memo) return memo[key];
      var val;
      if (t > 0) val = (t - 1) * R(t - 2, u, v, n + 1) + X * R(t - 1, u, v, n + 1);
      else if (u > 0) val = (u - 1) * R(t, u - 2, v, n + 1) + Y * R(t, u - 1, v, n + 1);
      else val = (v - 1) * R(t, u, v - 2, n + 1) + Z * R(t, u, v - 1, n + 1);
      memo[key] = val;
      return val;
    }
    return R;
  }

  // S, T, V matrices for the whole basis. atoms: [{ Z, xyz }] in bohr.
  function oneElectron(basis, atoms) {
    var nb = basis.length;
    var S = new Float64Array(nb * nb);
    var T = new Float64Array(nb * nb);
    var V = new Float64Array(nb * nb);

    for (var i = 0; i < nb; i++) {
      var fi = basis[i];
      for (var j = i; j < nb; j++) {
        var fj = basis[j];
        var AB = [fi.center[0] - fj.center[0], fi.center[1] - fj.center[1], fi.center[2] - fj.center[2]];
        var li = [fi.l, fi.m, fi.n], lj = [fj.l, fj.m, fj.n];
        var s = 0, t = 0, v = 0;

        for (var pi = 0; pi < fi.exps.length; pi++) {
          for (var pj = 0; pj < fj.exps.length; pj++) {
            var a = fi.exps[pi], b = fj.exps[pj];
            var c = fi.coefs[pi] * fj.coefs[pj];
            var p = a + b;
            var P = [
              (a * fi.center[0] + b * fj.center[0]) / p,
              (a * fi.center[1] + b * fj.center[1]) / p,
              (a * fi.center[2] + b * fj.center[2]) / p
            ];
            var Ex = makeE(a, b, AB[0]), Ey = makeE(a, b, AB[1]), Ez = makeE(a, b, AB[2]);
            var pref = Math.pow(Math.PI / p, 1.5);

            // overlap per dimension
            var sx = Ex(li[0], lj[0], 0), sy = Ey(li[1], lj[1], 0), sz = Ez(li[2], lj[2], 0);
            s += c * pref * sx * sy * sz;

            // kinetic: T_j-relations per dimension
            function kin1d(E, ia, ja) {
              return -2 * b * b * E(ia, ja + 2, 0) + b * (2 * ja + 1) * E(ia, ja, 0)
                - 0.5 * ja * (ja - 1) * E(ia, ja - 2, 0);
            }
            var tx = kin1d(Ex, li[0], lj[0]), ty = kin1d(Ey, li[1], lj[1]), tz = kin1d(Ez, li[2], lj[2]);
            t += c * pref * (tx * sy * sz + sx * ty * sz + sx * sy * tz);

            // nuclear attraction over all nuclei
            var tmax = li[0] + lj[0] + li[1] + lj[1] + li[2] + lj[2];
            for (var ai = 0; ai < atoms.length; ai++) {
              var C = atoms[ai].xyz;
              var R = makeR(p, P[0] - C[0], P[1] - C[1], P[2] - C[2], tmax);
              var sumR = 0;
              for (var tt = 0; tt <= li[0] + lj[0]; tt++) {
                var ex = Ex(li[0], lj[0], tt);
                if (ex === 0) continue;
                for (var uu = 0; uu <= li[1] + lj[1]; uu++) {
                  var ey = Ey(li[1], lj[1], uu);
                  if (ey === 0) continue;
                  for (var vv = 0; vv <= li[2] + lj[2]; vv++) {
                    var ez = Ez(li[2], lj[2], vv);
                    if (ez === 0) continue;
                    sumR += ex * ey * ez * R(tt, uu, vv, 0);
                  }
                }
              }
              v += -atoms[ai].Z * c * (2 * Math.PI / p) * sumR;
            }
          }
        }
        S[i * nb + j] = S[j * nb + i] = s;
        T[i * nb + j] = T[j * nb + i] = t;
        V[i * nb + j] = V[j * nb + i] = v;
      }
    }
    return { S: S, T: T, V: V };
  }

  App.integrals = { boys: boys, makeE: makeE, makeR: makeR, oneElectron: oneElectron };
})(typeof globalThis.App === "object" ? globalThis.App : (globalThis.App = {}));
