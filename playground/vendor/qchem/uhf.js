// Unrestricted Hartree-Fock: separate alpha/beta modes for open shells
// (radicals, triplet O2). Same GWH start, damping and DIIS as the RHF path;
// additionally reports <S^2> so spin contamination is visible, not hidden.
(function (App) {
  "use strict";

  var LA;
  var UHF_MAX_ITER = 600;

  function density1(C, nb, nocc) {
    var D = new Float64Array(nb * nb);
    for (var i = 0; i < nb; i++) {
      for (var j = 0; j < nb; j++) {
        var d = 0;
        for (var o = 0; o < nocc; o++) d += C[i * nb + o] * C[j * nb + o];
        D[i * nb + j] = d;
      }
    }
    return D;
  }

  // One quadruple loop builds J(Dt) and both exchange matrices.
  function buildFocks(H, Da, Db, eri, nb) {
    var J = new Float64Array(nb * nb), Ka = new Float64Array(nb * nb), Kb = new Float64Array(nb * nb);
    for (var m = 0; m < nb; m++) {
      for (var n = m; n < nb; n++) {
        var j = 0, ka = 0, kb = 0;
        for (var l = 0; l < nb; l++) {
          for (var s = 0; s < nb; s++) {
            var dt = Da[l * nb + s] + Db[l * nb + s];
            var gmnls = App.eri.get(eri, m, n, l, s);
            var gmlns = App.eri.get(eri, m, l, n, s);
            j += dt * gmnls;
            ka += Da[l * nb + s] * gmlns;
            kb += Db[l * nb + s] * gmlns;
          }
        }
        J[m * nb + n] = J[n * nb + m] = j;
        Ka[m * nb + n] = Ka[n * nb + m] = ka;
        Kb[m * nb + n] = Kb[n * nb + m] = kb;
      }
    }
    var Fa = new Float64Array(nb * nb), Fb = new Float64Array(nb * nb);
    for (var i = 0; i < nb * nb; i++) {
      Fa[i] = H[i] + J[i] - Ka[i];
      Fb[i] = H[i] + J[i] - Kb[i];
    }
    return { Fa: Fa, Fb: Fb, J: J, Ka: Ka, Kb: Kb };
  }

  function solveFock(F, X, nb) {
    var Ft = LA.matmul(LA.matmul(LA.transpose(X, nb), F, nb), X, nb);
    var eig = LA.eighSym(Ft, nb);
    return { eps: eig.values, C: LA.matmul(X, eig.vectors, nb) };
  }

  function diisError(F, D, S, X, nb) {
    var FDS = LA.matmul(LA.matmul(F, D, nb), S, nb);
    var SDF = LA.matmul(LA.matmul(S, D, nb), F, nb);
    var e = new Float64Array(nb * nb);
    for (var i = 0; i < nb * nb; i++) e[i] = FDS[i] - SDF[i];
    return LA.matmul(LA.matmul(LA.transpose(X, nb), e, nb), X, nb);
  }

  function runUHF(atoms, basis, ints, eri, charge, mult) {
    LA = App.linalg;
    var nb = basis.length;
    var nelec = atoms.reduce(function (s, a) { return s + a.Z; }, 0) - (charge || 0);
    var nun = (mult || 1) - 1; // unpaired electrons
    if ((nelec - nun) % 2 !== 0 || nun > nelec || nun < 0)
      throw new Error(App.tr("err.uhf.mult", { mult: mult, n: nelec }));
    var na = (nelec + nun) / 2, nbe = nelec - na;
    if (na > nb) throw new Error(App.tr("err.uhf.basis"));

    var H = new Float64Array(nb * nb);
    for (var i = 0; i < nb * nb; i++) H[i] = ints.T[i] + ints.V[i];
    var X = LA.invSqrtSym(ints.S, nb);
    var Enuc = App.scf.nuclearRepulsion(atoms);

    var F0 = new Float64Array(nb * nb);
    for (i = 0; i < nb; i++) {
      for (var j0 = 0; j0 < nb; j0++) {
        F0[i * nb + j0] = i === j0
          ? H[i * nb + i]
          : 0.875 * ints.S[i * nb + j0] * (H[i * nb + i] + H[j0 * nb + j0]);
      }
    }
    var sol0 = solveFock(F0, X, nb);
    var solA = sol0, solB = sol0;
    var Da = density1(sol0.C, nb, na), Db = density1(sol0.C, nb, nbe);
    var Eold = 0, result = null;
    var dFa = [], dFb = [], dE = [];
    var prevRms = Infinity, prevErrNorm = Infinity;

    for (var iter = 1; iter <= UHF_MAX_ITER; iter++) {
      var fk = buildFocks(H, Da, Db, eri, nb);
      var Fa = fk.Fa, Fb = fk.Fb;
      var Dt = new Float64Array(nb * nb);
      for (i = 0; i < nb * nb; i++) Dt[i] = Da[i] + Db[i];

      var ea = diisError(fk.Fa, Da, ints.S, X, nb);
      var eb = diisError(fk.Fb, Db, ints.S, X, nb);
      var err = new Float64Array(2 * nb * nb);
      err.set(ea); err.set(eb, nb * nb);
      var errNorm = 0;
      for (var k = 0; k < err.length; k++) errNorm += err[k] * err[k];
      if (iter > 8 && isFinite(prevErrNorm) && errNorm > prevErrNorm * 1.4) {
        dFa.length = 0; dFb.length = 0; dE.length = 0;
      }
      prevErrNorm = errNorm;
      dFa.push(fk.Fa); dFb.push(fk.Fb); dE.push(err);
      if (dFa.length > 8) { dFa.shift(); dFb.shift(); dE.shift(); }
      var nd = iter <= 5 ? 0 : dFa.length; // damped start, DIIS after
      if (nd >= 2) {
        var B = [], rhs = [];
        for (i = 0; i < nd; i++) {
          B.push([]);
          for (var j = 0; j < nd; j++) {
            var dot = 0;
            for (var k = 0; k < 2 * nb * nb; k++) dot += dE[i][k] * dE[j][k];
            B[i].push(dot);
          }
          B[i].push(-1); rhs.push(0);
        }
        B.push(new Array(nd).fill(-1).concat([0])); rhs.push(-1);
        var cs = LA.solveLin(B, rhs);
        if (cs) {
          Fa = new Float64Array(nb * nb); Fb = new Float64Array(nb * nb);
          for (i = 0; i < nd; i++) {
            for (j = 0; j < nb * nb; j++) { Fa[j] += cs[i] * dFa[i][j]; Fb[j] += cs[i] * dFb[i][j]; }
          }
        }
      }

      solA = solveFock(Fa, X, nb);
      solB = solveFock(Fb, X, nb);
      var DaN = density1(solA.C, nb, na), DbN = density1(solB.C, nb, nbe);
      var Eelec = 0;
      for (i = 0; i < nb * nb; i++) Eelec += 0.5 * (Dt[i] * H[i] + Da[i] * fk.Fa[i] + Db[i] * fk.Fb[i]);

      var rms = 0;
      for (i = 0; i < nb * nb; i++) {
        var d = (DaN[i] + DbN[i]) - Dt[i];
        rms += d * d;
      }
      rms = Math.sqrt(rms / (nb * nb));
      var damping = iter <= 8 ? 0.5 : 0;
      if (iter > 8 && rms > prevRms * 0.98) damping = 0.2;
      if (damping > 0) {
        for (i = 0; i < nb * nb; i++) {
          Da[i] = damping * Da[i] + (1 - damping) * DaN[i];
          Db[i] = damping * Db[i] + (1 - damping) * DbN[i];
        }
      } else { Da = DaN; Db = DbN; }
      prevRms = rms;

      if (Math.abs(Eelec - Eold) < 1e-9 && rms < 1e-7) {
        var Ds = new Float64Array(nb * nb);
        for (i = 0; i < nb * nb; i++) { Dt[i] = Da[i] + Db[i]; Ds[i] = Da[i] - Db[i]; }
        // <S^2> = Sz(Sz+1) + Nb - sum |<phi_i^a|phi_j^b>|^2 over occupied pairs
        var O = LA.matmul(LA.matmul(LA.transpose(solA.C, nb), ints.S, nb), solB.C, nb);
        var ss = 0;
        for (i = 0; i < na; i++) for (j = 0; j < nbe; j++) ss += O[i * nb + j] * O[i * nb + j];
        var sz = nun / 2;
        result = {
          converged: true, iterations: iter, uhf: true, mult: mult,
          E: Eelec + Enuc, Eelec: Eelec, Enuc: Enuc,
          ET: LA.trace2(Dt, ints.T, nb), EVne: LA.trace2(Dt, ints.V, nb),
          EJ: 0.5 * LA.trace2(Dt, fk.J, nb),
          EK: -0.5 * (LA.trace2(Da, fk.Ka, nb) + LA.trace2(Db, fk.Kb, nb)),
          eps: Array.from(solA.eps), C: Array.from(solA.C), nocc: na,
          epsB: Array.from(solB.eps), CB: Array.from(solB.C), noccB: nbe,
          S2: sz * (sz + 1) + nbe - ss, S2exact: sz * (sz + 1),
          D: Array.from(Dt), Ds: Array.from(Ds),
          nelec: nelec, nb: nb
        };
        break;
      }
      Eold = Eelec;
    }

    if (!result) throw new Error(App.tr("err.uhf.conv"));

    var DS = LA.matmul(result.D, ints.S, nb);
    var DsS = LA.matmul(result.Ds, ints.S, nb);
    var pop = atoms.map(function (a) { return a.Z; });
    var spin = atoms.map(function () { return 0; });
    for (i = 0; i < nb; i++) {
      pop[basis[i].atom] -= DS[i * nb + i];
      spin[basis[i].atom] += DsS[i * nb + i];
    }
    result.mulliken = pop;
    result.spinPop = spin;
    return result;
  }

  App.uhf = { runUHF: runUHF };
})(typeof globalThis.App === "object" ? globalThis.App : (globalThis.App = {}));
