// Restricted Hartree-Fock with DIIS convergence acceleration.
// Closed-shell only: every occupied field mode holds two quanta (spins up/down).
(function (App) {
  "use strict";

  var LA;

  function nuclearRepulsion(atoms) {
    var e = 0;
    for (var a = 0; a < atoms.length; a++) {
      for (var b = a + 1; b < atoms.length; b++) {
        var dx = atoms[a].xyz[0] - atoms[b].xyz[0];
        var dy = atoms[a].xyz[1] - atoms[b].xyz[1];
        var dz = atoms[a].xyz[2] - atoms[b].xyz[2];
        e += atoms[a].Z * atoms[b].Z / Math.sqrt(dx * dx + dy * dy + dz * dz);
      }
    }
    return e;
  }

  // F = H + J - K/2 built from density and packed ERIs; also returns EJ, EK.
  function buildFock(H, D, eri, nb) {
    var J = new Float64Array(nb * nb), K = new Float64Array(nb * nb);
    for (var m = 0; m < nb; m++) {
      for (var n = m; n < nb; n++) {
        var j = 0, k = 0;
        for (var l = 0; l < nb; l++) {
          for (var s = 0; s < nb; s++) {
            var d = D[l * nb + s];
            if (d === 0) continue;
            j += d * App.eri.get(eri, m, n, l, s);
            k += d * App.eri.get(eri, m, l, n, s);
          }
        }
        J[m * nb + n] = J[n * nb + m] = j;
        K[m * nb + n] = K[n * nb + m] = k;
      }
    }
    var F = new Float64Array(nb * nb);
    for (var i = 0; i < nb * nb; i++) F[i] = H[i] + J[i] - 0.5 * K[i];
    return {
      F: F,
      EJ: 0.5 * LA.trace2(D, J, nb),
      EK: -0.25 * LA.trace2(D, K, nb)
    };
  }

  function density(C, nb, nocc) {
    var D = new Float64Array(nb * nb);
    for (var i = 0; i < nb; i++) {
      for (var j = 0; j < nb; j++) {
        var d = 0;
        for (var o = 0; o < nocc; o++) d += C[i * nb + o] * C[j * nb + o];
        D[i * nb + j] = 2 * d;
      }
    }
    return D;
  }

  function solveFock(F, X, nb) {
    var Ft = LA.matmul(LA.matmul(LA.transpose(X, nb), F, nb), X, nb);
    var eig = LA.eighSym(Ft, nb);
    return { eps: eig.values, C: LA.matmul(X, eig.vectors, nb) };
  }

  // atoms: [{Z, xyz(bohr)}]; charge: total molecular charge.
  function runRHF(atoms, basis, ints, eri, charge) {
    LA = App.linalg;
    var nb = basis.length;
    var nelec = atoms.reduce(function (s, a) { return s + a.Z; }, 0) - (charge || 0);
    if (nelec % 2 !== 0) throw new Error(App.tr("err.rhf.open"));
    var nocc = nelec / 2;
    if (nocc > nb) throw new Error(App.tr("err.rhf.basis"));

    var H = new Float64Array(nb * nb);
    for (var i = 0; i < nb * nb; i++) H[i] = ints.T[i] + ints.V[i];
    var X = LA.invSqrtSym(ints.S, nb);
    var Enuc = nuclearRepulsion(atoms);

    // GWH guess: core guess alone often converges N2-like systems to a wrong state
    var F0 = new Float64Array(nb * nb);
    for (i = 0; i < nb; i++) {
      for (var j0 = 0; j0 < nb; j0++) {
        F0[i * nb + j0] = i === j0
          ? H[i * nb + i]
          : 0.875 * ints.S[i * nb + j0] * (H[i * nb + i] + H[j0 * nb + j0]);
      }
    }
    var sol = solveFock(F0, X, nb);
    var D = density(sol.C, nb, nocc);
    var Eold = 0, result = null;
    var diisF = [], diisE = [];

    for (var iter = 1; iter <= 200; iter++) {
      var fock = buildFock(H, D, eri, nb);
      var F = fock.F;

      // DIIS error e = X^T (FDS - SDF) X
      var FDS = LA.matmul(LA.matmul(F, D, nb), ints.S, nb);
      var SDF = LA.matmul(LA.matmul(ints.S, D, nb), F, nb);
      var err = new Float64Array(nb * nb);
      for (i = 0; i < nb * nb; i++) err[i] = FDS[i] - SDF[i];
      err = LA.matmul(LA.matmul(LA.transpose(X, nb), err, nb), X, nb);

      diisF.push(F); diisE.push(err);
      if (diisF.length > 8) { diisF.shift(); diisE.shift(); }
      var nd = diisF.length;
      if (iter <= 3) {
        // damped start keeps early iterations from overshooting into excited states
        nd = 0;
      } else if (nd >= 2) {
        var B = [], rhs = [];
        for (i = 0; i < nd; i++) {
          B.push([]);
          for (var j = 0; j < nd; j++) B[i].push(LA.trace2(diisE[i], LA.transpose(diisE[j], nb), nb));
          B[i].push(-1); rhs.push(0);
        }
        B.push(new Array(nd).fill(-1).concat([0])); rhs.push(-1);
        var cs = LA.solveLin(B, rhs);
        if (cs) {
          var Fd = new Float64Array(nb * nb);
          for (i = 0; i < nd; i++) for (j = 0; j < nb * nb; j++) Fd[j] += cs[i] * diisF[i][j];
          F = Fd;
        }
      }

      sol = solveFock(F, X, nb);
      var Dnew = density(sol.C, nb, nocc);
      var Eelec = 0.5 * LA.trace2(D, H, nb) + 0.5 * LA.trace2(D, fock.F, nb);

      var rms = 0;
      for (i = 0; i < nb * nb; i++) { var d = Dnew[i] - D[i]; rms += d * d; }
      rms = Math.sqrt(rms / (nb * nb));
      if (iter <= 3) {
        for (i = 0; i < nb * nb; i++) D[i] = 0.5 * D[i] + 0.5 * Dnew[i];
      } else {
        D = Dnew;
      }

      if (Math.abs(Eelec - Eold) < 1e-9 && rms < 1e-7) {
        result = {
          converged: true, iterations: iter,
          E: Eelec + Enuc, Eelec: Eelec, Enuc: Enuc,
          ET: LA.trace2(D, ints.T, nb), EVne: LA.trace2(D, ints.V, nb),
          EJ: fock.EJ, EK: fock.EK,
          eps: Array.from(sol.eps), C: Array.from(sol.C),
          D: Array.from(D),
          nocc: nocc, nelec: nelec, nb: nb
        };
        break;
      }
      Eold = Eelec;
    }

    if (!result) throw new Error(App.tr("err.rhf.conv"));

    // Mulliken charges: q_A = Z_A - sum_{mu in A} (DS)_mumu
    var DS = LA.matmul(D, ints.S, nb);
    var pop = atoms.map(function (a) { return a.Z; });
    for (i = 0; i < nb; i++) pop[basis[i].atom] -= DS[i * nb + i];
    result.mulliken = pop;
    return result;
  }

  App.scf = { runRHF: runRHF, nuclearRepulsion: nuclearRepulsion };
})(typeof globalThis.App === "object" ? globalThis.App : (globalThis.App = {}));
