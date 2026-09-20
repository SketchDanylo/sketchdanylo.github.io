// STO-3G basis set data for H..Ne.
// Source: Basis Set Exchange (www.basissetexchange.org), revision "Data from Gaussian09".
// Each AO is a contraction of 3 primitive Gaussians; 2s/2p share exponents (sp shells).
(function (App) {
  "use strict";

  var C1S = [0.1543289673, 0.5353281423, 0.4446345422];
  var C2S = [-0.09996722919, 0.3995128261, 0.7001154689];
  var C2P = [0.1559162750, 0.6076837186, 0.3919573931];

  // Per element: list of shells { l, exps, coefs }
  var DATA = {
    1:  [{ l: 0, e: [3.425250914, 0.6239137298, 0.1688554040], c: C1S }],
    2:  [{ l: 0, e: [6.362421394, 1.158922999, 0.3136497915], c: C1S }],
    3:  [{ l: 0, e: [16.11957475, 2.936200663, 0.7946504870], c: C1S },
         { l: 0, e: [0.6362897469, 0.1478600533, 0.04808867840], c: C2S },
         { l: 1, e: [0.6362897469, 0.1478600533, 0.04808867840], c: C2P }],
    4:  [{ l: 0, e: [30.16787069, 5.495115306, 1.487192653], c: C1S },
         { l: 0, e: [1.314833110, 0.3055389383, 0.09937074560], c: C2S },
         { l: 1, e: [1.314833110, 0.3055389383, 0.09937074560], c: C2P }],
    5:  [{ l: 0, e: [48.79111318, 8.887362172, 2.405267040], c: C1S },
         { l: 0, e: [2.236956142, 0.5198204999, 0.1690617600], c: C2S },
         { l: 1, e: [2.236956142, 0.5198204999, 0.1690617600], c: C2P }],
    6:  [{ l: 0, e: [71.61683735, 13.04509632, 3.530512160], c: C1S },
         { l: 0, e: [2.941249355, 0.6834830964, 0.2222899159], c: C2S },
         { l: 1, e: [2.941249355, 0.6834830964, 0.2222899159], c: C2P }],
    7:  [{ l: 0, e: [99.10616896, 18.05231239, 4.885660238], c: C1S },
         { l: 0, e: [3.780455879, 0.8784966449, 0.2857143744], c: C2S },
         { l: 1, e: [3.780455879, 0.8784966449, 0.2857143744], c: C2P }],
    8:  [{ l: 0, e: [130.7093214, 23.80886605, 6.443608313], c: C1S },
         { l: 0, e: [5.033151319, 1.169596125, 0.3803889600], c: C2S },
         { l: 1, e: [5.033151319, 1.169596125, 0.3803889600], c: C2P }],
    9:  [{ l: 0, e: [166.6791340, 30.36081233, 8.216820672], c: C1S },
         { l: 0, e: [6.464803249, 1.502281245, 0.4885884864], c: C2S },
         { l: 1, e: [6.464803249, 1.502281245, 0.4885884864], c: C2P }],
    10: [{ l: 0, e: [207.0156070, 37.70815124, 10.20529731], c: C1S },
         { l: 0, e: [8.246315120, 1.916266291, 0.6232292721], c: C2S },
         { l: 1, e: [8.246315120, 1.916266291, 0.6232292721], c: C2P }]
  };

  var SYMBOLS = ["", "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne"];

  function doubleFactorial(n) {
    var r = 1;
    for (var k = n; k > 1; k -= 2) r *= k;
    return r;
  }

  // Norm of a primitive cartesian Gaussian x^l y^m z^n exp(-a r^2)
  function primNorm(a, l, m, n) {
    var L = l + m + n;
    return Math.pow(2 * a / Math.PI, 0.75) * Math.pow(4 * a, L / 2) /
      Math.sqrt(doubleFactorial(2 * l - 1) * doubleFactorial(2 * m - 1) * doubleFactorial(2 * n - 1));
  }

  // Cartesian components for angular momentum l (s; px,py,pz; 6 cartesian d)
  var CART = {
    0: [[0, 0, 0]],
    1: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
    2: [[2, 0, 0], [0, 2, 0], [0, 0, 2], [1, 1, 0], [1, 0, 1], [0, 1, 1]]
  };
  var SUF = { 1: ["px", "py", "pz"], 2: ["dxx", "dyy", "dzz", "dxy", "dxz", "dyz"] };

  // shell display name: k-th shell of given l on an atom ("2s'" = outer valence)
  function shellName(Z, l, k) {
    if (l === 1) return "2p" + (k > 1 ? "'".repeat(k - 1) : "");
    if (l === 2) return "3d";
    if (Z <= 2) return k === 1 ? "1s" : "1s" + "'".repeat(k - 1);
    return k === 1 ? "1s" : "2s" + (k > 2 ? "'".repeat(k - 2) : "");
  }

  // Build the list of contracted basis functions for a molecule.
  // atoms: [{ Z, xyz: [bohr] }]. Returns plain objects (JSON-able for the worker).
  function buildBasis(atoms, basisName) {
    var table = App.BASIS_TABLES[basisName || "STO-3G"];
    if (!table) throw new Error("Unknown basis: " + basisName);
    var fns = [];
    atoms.forEach(function (atom, ai) {
      var shells = table[atom.Z];
      if (!shells) throw new Error("Element Z=" + atom.Z + " not supported (H..Ne only)");
      var lCount = {};
      shells.forEach(function (sh) {
        lCount[sh.l] = (lCount[sh.l] || 0) + 1;
        var sn = shellName(atom.Z, sh.l, lCount[sh.l]);
        CART[sh.l].forEach(function (lmn, comp) {
          var l = lmn[0], m = lmn[1], n = lmn[2];
          // primitive-normalized contraction coefficients
          var coefs = sh.c.map(function (c, i) { return c * primNorm(sh.e[i], l, m, n); });
          // normalize the contracted function: <phi|phi> = 1
          var s = 0;
          for (var i = 0; i < sh.e.length; i++) {
            for (var j = 0; j < sh.e.length; j++) {
              var p = sh.e[i] + sh.e[j];
              var pref = Math.pow(Math.PI / p, 1.5) / Math.pow(2 * p, l + m + n) *
                doubleFactorial(2 * l - 1) * doubleFactorial(2 * m - 1) * doubleFactorial(2 * n - 1);
              s += coefs[i] * coefs[j] * pref;
            }
          }
          var nrm = 1 / Math.sqrt(s);
          var lbl = sh.l === 0 ? sn : sn.replace(sh.l === 1 ? "p" : "d", SUF[sh.l][comp]);
          fns.push({
            atom: ai,
            center: atom.xyz.slice(),
            l: l, m: m, n: n,
            exps: sh.e.slice(),
            coefs: coefs.map(function (c) { return c * nrm; }),
            label: SYMBOLS[atom.Z] + (ai + 1) + " " + lbl
          });
        });
      });
    });
    return fns;
  }

  App.SYMBOLS = SYMBOLS;
  App.BASIS_TABLES = { "STO-3G": DATA };
  App.buildBasis = buildBasis;
})(typeof globalThis.App === "object" ? globalThis.App : (globalThis.App = {}));
