// Molecular properties from the converged density.
// Dipole: <i|r|j> via the l-raising identity x·G(l) = G(l+1) + Bx·G(l),
// so only MD overlap coefficients are needed. Mayer bond orders from (DS).
(function (App) {
  "use strict";

  var DEBYE_PER_AU = 2.5417464;

  // Three nb x nb matrices <i|x|j>, <i|y|j>, <i|z|j> (absolute coordinates).
  function dipoleMatrices(basis) {
    var nb = basis.length;
    var M = [new Float64Array(nb * nb), new Float64Array(nb * nb), new Float64Array(nb * nb)];
    for (var i = 0; i < nb; i++) {
      var fi = basis[i];
      for (var j = i; j < nb; j++) {
        var fj = basis[j];
        var AB = [fi.center[0] - fj.center[0], fi.center[1] - fj.center[1], fi.center[2] - fj.center[2]];
        var li = [fi.l, fi.m, fi.n], lj = [fj.l, fj.m, fj.n];
        var acc = [0, 0, 0];
        for (var pi = 0; pi < fi.exps.length; pi++) {
          for (var pj = 0; pj < fj.exps.length; pj++) {
            var a = fi.exps[pi], b = fj.exps[pj], p = a + b;
            var c = fi.coefs[pi] * fj.coefs[pj];
            var pref = c * Math.pow(Math.PI / p, 1.5);
            var E = [App.integrals.makeE(a, b, AB[0]), App.integrals.makeE(a, b, AB[1]), App.integrals.makeE(a, b, AB[2])];
            var s0 = [E[0](li[0], lj[0], 0), E[1](li[1], lj[1], 0), E[2](li[2], lj[2], 0)];
            for (var d = 0; d < 3; d++) {
              var s1 = E[d](li[d], lj[d] + 1, 0) + fj.center[d] * s0[d];
              acc[d] += pref * s1 * s0[(d + 1) % 3] * s0[(d + 2) % 3];
            }
          }
        }
        for (d = 0; d < 3; d++) { M[d][i * nb + j] = acc[d]; M[d][j * nb + i] = acc[d]; }
      }
    }
    return M;
  }

  // Mayer bond order between atom blocks of (DS); spin term covers UHF.
  function mayerBonds(atoms, basis, S, D, Ds) {
    var nb = basis.length;
    var DtS = App.linalg.matmul(D, S, nb);
    var DsS = Ds ? App.linalg.matmul(Ds, S, nb) : null;
    var bonds = [];
    for (var a = 0; a < atoms.length; a++) {
      for (var b = a + 1; b < atoms.length; b++) {
        var B = 0;
        for (var i = 0; i < nb; i++) {
          if (basis[i].atom !== a) continue;
          for (var j = 0; j < nb; j++) {
            if (basis[j].atom !== b) continue;
            B += DtS[i * nb + j] * DtS[j * nb + i];
            if (DsS) B += DsS[i * nb + j] * DsS[j * nb + i];
          }
        }
        if (B > 0.05) bonds.push({ a: a, b: b, order: B });
      }
    }
    return bonds;
  }

  // scf must carry D (total density) and, for UHF, Ds (spin density).
  function compute(atoms, basis, S, scf) {
    var nb = scf.nb, D = scf.D;
    var M = dipoleMatrices(basis);

    var mu = [0, 0, 0], O = [0, 0, 0], Zsum = 0;
    atoms.forEach(function (at) {
      Zsum += at.Z;
      for (var d = 0; d < 3; d++) { mu[d] += at.Z * at.xyz[d]; O[d] += at.Z * at.xyz[d]; }
    });
    for (var d = 0; d < 3; d++) O[d] /= Zsum;
    for (d = 0; d < 3; d++) {
      var el = 0;
      for (var i = 0; i < nb * nb; i++) el += D[i] * M[d][i];
      mu[d] -= el;
    }
    // ions: report about the center of nuclear charge, not the lab origin
    var q = Zsum - scf.nelec;
    for (d = 0; d < 3; d++) mu[d] -= q * O[d];

    var au = Math.hypot(mu[0], mu[1], mu[2]);
    return {
      dipole: { vec: mu, au: au, debye: au * DEBYE_PER_AU },
      mayer: mayerBonds(atoms, basis, S, D, scf.Ds || null)
    };
  }

  App.props = { compute: compute, dipoleMatrices: dipoleMatrices, DEBYE_PER_AU: DEBYE_PER_AU };
})(typeof globalThis.App === "object" ? globalThis.App : (globalThis.App = {}));
