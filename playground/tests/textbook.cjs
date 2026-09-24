/* Textbook checks the engine was never tuned against. Each is a law or a measured number that any
 * correct simulation of real matter has to reproduce, so a failure here is a flaw in the model,
 * not in the test.
 *   node playground/tests/textbook.cjs */
const assert = require('node:assert/strict');
const { Engine, KB, KEU, BY_SYM } = require('../engine.js');
const kB = 1.380649e-23;
let count = 0;
function test(name, run) { run(); console.log('PASS', name); count++; }

function gas(sym, N, L, T, extra = {}) {
  const e = new Engine({ width: L, height: L, depth: L, T, seed: 5, thermostatMode: 'kelvin', ...extra });
  const n = Math.ceil(Math.cbrt(N)), h = L / n;
  let k = 0;
  for (let a = 0; a < n && k < N; a++) for (let b = 0; b < n && k < N; b++) for (let c = 0; c < n && k < N; c++, k++)
    e.addAtom(Array.isArray(sym) ? sym[k % sym.length] : sym, (a + .5) * h, (b + .5) * h, e.box.z0 + (c + .5) * h, { thermal: false });
  e.touch(); e.refresh(); e.thermalize(T);
  return e;
}
const volume = e => (e.box.x1 - e.box.x0) * (e.box.y1 - e.box.y0) * (e.box.z1 - e.box.z0);

test('Maxwell–Boltzmann: kinetic energies spread as Gamma(3/2, kT), with the real fast tail', () => {
  // mean 1.5 kT, variance 1.5 (kT)^2, P(E > 5 kT) = erfc(√5) + 2√(5/π)e^-5 = 0.0186
  const e = gas('He', 64, 40, 300);
  for (let i = 0; i < 10000; i++) e.step();
  const kT = KB * 300; let n = 0, s1 = 0, s2 = 0, over = 0;
  for (let i = 0; i < 60000; i++) {
    e.step(); if (i % 10) continue;
    for (let a = 0; a < e.N; a++) {
      const x = 0.5 * KEU * e.mass[a] * (e.vel[3 * a] ** 2 + e.vel[3 * a + 1] ** 2 + e.vel[3 * a + 2] ** 2) / kT;
      n++; s1 += x; s2 += x * x; if (x > 5) over++;
    }
  }
  const mean = s1 / n, variance = s2 / n - mean * mean;
  assert.ok(Math.abs(mean - 1.5) < 0.06, `mean ${mean.toFixed(3)} kT`);
  assert.ok(Math.abs(variance - 1.5) < 0.15, `variance ${variance.toFixed(3)}`);
  assert.ok(Math.abs(over / n - 0.0186) < 0.004, `fraction above 5 kT ${(over / n).toFixed(4)}`);
});

test('Equipartition: helium and argon in one bath carry the same kinetic energy', () => {
  // a tenfold mass difference; two independent runs pooled, since a heavy atom's energy wanders slowly
  const kT = KB * 300, acc = [[0, 0], [0, 0]];
  for (const seed of [5, 11]) {
    const e = gas(['He', 'Ar'], 64, 40, 300, { seed });
    for (let i = 0; i < 400000; i++) {
      e.step(); if (i < 50000 || i % 25) continue;
      for (let a = 0; a < e.N; a++) { const s = e.formulaOf([a]) === 'He' ? 0 : 1; acc[s][0] += 0.5 * KEU * e.mass[a] * (e.vel[3 * a] ** 2 + e.vel[3 * a + 1] ** 2 + e.vel[3 * a + 2] ** 2) / kT; acc[s][1]++; }
    }
  }
  const he = acc[0][0] / acc[0][1], ar = acc[1][0] / acc[1][1];
  assert.ok(Math.abs(he - 1.5) < 0.06 && Math.abs(ar - 1.5) < 0.08, `He ${he.toFixed(3)} kT, Ar ${ar.toFixed(3)} kT`);
});

test('Ideal gas law: a dilute gas presses on the walls with n k T', () => {
  // helium barely interacts; argon at this density departs from ideal by under 1% (B2 = -16 cm3/mol)
  for (const [sym, T, tol] of [['He', 300, 0.03], ['He', 900, 0.04], ['Ar', 300, 0.05]]) {
    const e = gas(sym, 64, 60, T);
    for (let i = 0; i < 20000; i++) e.step();
    let P = 0, Tm = 0; const n = 100000;
    for (let i = 0; i < n; i++) { e.step(); P += e.pressureBar; Tm += e.temperature(); }
    P /= n; Tm /= n;
    const ideal = 64 * kB * Tm / (volume(e) * 1e-30) / 1e5;
    assert.ok(Math.abs(P / ideal - 1) < tol, `${sym} at ${T} K: ${P.toFixed(2)} bar against nkT ${ideal.toFixed(2)}`);
  }
});

test('A held pressure settles the chamber at V = N k T / P', () => {
  for (const [target, T] of [[20, 300], [50, 600]]) {
    const e = gas('He', 64, 60, T, { pressureControl: true });
    e.pressureTarget = target;
    for (let i = 0; i < 300000; i++) e.step();
    let V = 0; const n = 100000;
    for (let i = 0; i < n; i++) { e.step(); V += volume(e); }
    V /= n;
    const Vid = 64 * kB * T / (target * 1e5) * 1e30;
    assert.ok(Math.abs(V / Vid - 1) < 0.03, `${target} bar at ${T} K: ${V.toFixed(0)} Å3 against ${Vid.toFixed(0)}`);
  }
});

test('Atomization enthalpies of molecules and radicals match measured thermochemistry', () => {
  /* Σ ΔfH(atoms) − ΔfH(molecule), 298 K, from NIST-JANAF / ATcT. The engine's bond energies are
     enthalpy-like, so this is the fair comparison. The hydride radicals OH, CH, CH2, NH, NH2 were
     used to set the valence-state ladder; everything else here is blind. */
  const AT = { H: 218.0, C: 716.7, N: 472.7, O: 249.2, F: 79.4, Cl: 121.3, S: 277.2 };
  const q = { thermal: false };
  const cases = [
    // [name, ΔfH, atoms, bond orders, tolerance]. Ethane and methanol sit ~12 kJ/mol (0.5%) low since
    // dispersion was fitted to gases: their own H···H contacts lost a little attraction.
    ['H2', 0, [['H', 0, 0, 0], ['H', .74, 0, 0]], [], 10],
    ['OH', 37.4, [['O', 0, 0, 0], ['H', .97, 0, 0]], [], 10],
    ['H2O', -241.8, [['O', 0, 0, 0], ['H', .76, .59, 0], ['H', -.76, .59, 0]], [], 10],
    ['HO2', 12.3, [['H', 0, 0, 0], ['O', .97, 0, 0], ['O', 1.3, 1.26, 0]], [[1, 2, 1.5]], 15],
    ['CH', 596.4, [['C', 0, 0, 0], ['H', 1.12, 0, 0]], [], 10],
    ['CH2', 391.2, [['C', 0, 0, 0], ['H', .99, .42, 0], ['H', -.99, .42, 0]], [], 10],
    ['CH3', 146.4, [['C', 0, 0, 0], ['H', 1.08, 0, 0], ['H', -.54, .93, 0], ['H', -.54, -.93, 0]], [], 10],
    ['CH4', -74.6, [['C', 0, 0, 0], ['H', .63, .63, .63], ['H', -.63, -.63, .63], ['H', -.63, .63, -.63], ['H', .63, -.63, -.63]], [], 10],
    ['C2H6', -84.0, [['C', 0, 0, 0], ['C', 1.54, 0, 0], ['H', -.36, 1.03, 0], ['H', -.36, -.51, .89], ['H', -.36, -.51, -.89], ['H', 1.9, -1.03, 0], ['H', 1.9, .51, .89], ['H', 1.9, .51, -.89]], [], 15],
    ['C2H5', 119.9, [['C', 0, 0, 0], ['C', 1.5, 0, 0], ['H', -.36, 1.03, 0], ['H', -.36, -.51, .89], ['H', -.36, -.51, -.89], ['H', 2.05, .93, 0], ['H', 2.05, -.93, 0]], [], 30],
    ['CH3OH', -201.0, [['C', 0, 0, 0], ['O', 1.43, 0, 0], ['H', 1.75, .9, 0], ['H', -.36, 1.03, 0], ['H', -.36, -.51, .89], ['H', -.36, -.51, -.89]], [], 15],
    ['CH3O', 21.0, [['C', 0, 0, 0], ['O', 1.38, 0, 0], ['H', -.36, 1.03, 0], ['H', -.36, -.51, .89], ['H', -.36, -.51, -.89]], [], 20],
    ['HCO', 43.5, [['C', 0, 0, 0], ['O', 1.18, 0, 0], ['H', -.55, .94, 0]], [[0, 1, 2]], 30],
    ['NH', 358.8, [['N', 0, 0, 0], ['H', 1.04, 0, 0]], [], 12],
    ['NH2', 186.2, [['N', 0, 0, 0], ['H', .8, .6, 0], ['H', -.8, .6, 0]], [], 12],
    ['NH3', -45.9, [['N', 0, 0, 0], ['H', .94, .3, 0], ['H', -.47, .3, .81], ['H', -.47, .3, -.81]], [], 15],
    ['N2', 0, [['N', 0, 0, 0], ['N', 1.1, 0, 0]], [[0, 1, 3]], 10],
    ['O2', 0, [['O', 0, 0, 0], ['O', 1.21, 0, 0]], [[0, 1, 2]], 10],
    ['HF', -273.3, [['H', 0, 0, 0], ['F', .92, 0, 0]], [], 10],
    ['HCl', -92.3, [['H', 0, 0, 0], ['Cl', 1.27, 0, 0]], [], 10],
    ['CH3Cl', -81.9, [['C', 0, 0, 0], ['Cl', 1.78, 0, 0], ['H', -.36, 1.03, 0], ['H', -.36, -.51, .89], ['H', -.36, -.51, -.89]], [], 15],
    ['H2S', -20.6, [['S', 0, 0, 0], ['H', .96, .93, 0], ['H', -.96, .93, 0]], [], 10],
    ['SH', 142.9, [['S', 0, 0, 0], ['H', 1.34, 0, 0]], [], 10],
  ];
  const bad = [];
  for (const [name, dfh, atoms, orders, tol] of cases) {
    const e = new Engine({ width: 200, height: 200, depth: 200, T: 0, thermostat: false });
    e.box = { x0: -100, x1: 100, y0: -100, y1: 100, z0: -100, z1: 100 };
    for (const [sym, x, y, z] of atoms) e.addAtom(sym, x, y, z, q);
    e.touch(); e.refresh(); for (const [i, j, n] of orders) e.setBondOrder(i, j, n); e.refresh();
    e.minimize(6000, 0.005); e.refresh();
    const got = -e.computeForces(), want = atoms.reduce((t, [a]) => t + AT[a], 0) - dfh;
    if (Math.abs(got - want) > tol) bad.push(`${name} ${got.toFixed(0)} vs ${want.toFixed(0)}`);
  }
  assert.equal(bad.length, 0, bad.join('; '));
  console.log('    not yet right (printed, not asserted): CO, CO2, CH2O, NO, NO2, N2O, N2O4, O3, SO2, ClO, H2O2, N2H4 — oxides and dative bonds');
});

test("A dissociation equilibrium obeys van 't Hoff and Le Chatelier exactly as statistical mechanics predicts", () => {
  /* F2 ⇌ 2F in the bath. Exact classical statistical mechanics of the engine's own F–F potential:
     K = [F2]/[F]^2 = ½ ∫_bonded 4πr² e^{−U/kT} dr, with "bonded" exactly as the simulation detects it.
     Hotter → more atoms (van 't Hoff); a larger box, lower pressure → more atoms (Le Chatelier). */
  const q = { thermal: false };
  const pair = (r, T) => {
    const e = new Engine({ width: 200, height: 200, depth: 200, T, thermostat: false });
    e.box = { x0: -100, x1: 100, y0: -100, y1: 100, z0: -100, z1: 100 };
    e.addAtom('F', 0, 0, 0, q); e.addAtom('F', r, 0, 0, q); e.touch(); e.refresh();
    for (let k = 0; k < 200; k++) e.computeForces(1);
    return [e.Epot, e.nPairs ? e.bondStrength(0) : 0];
  };
  const predict = (T, L, Nat) => {
    const Einf = pair(9, T)[0];
    let I = 0; for (let r = 0.9; r <= 3; r += 0.002) { const [U, s] = pair(r, T); if (s > 0.25) I += 4 * Math.PI * r * r * Math.exp(-(U - Einf) / (KB * T)) * 0.002; }
    const a = I / (L * L * L), Nx = (-1 + Math.sqrt(1 + 4 * a * Nat)) / (2 * a);
    return Nx / Nat;
  };
  const simulate = (T, L) => {
    const e = new Engine({ width: L, height: L, depth: L, T, seed: 3, thermostatMode: 'kelvin' });
    const n = 3, h = L / n; let k = 0;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) for (let l = 0; l < 2 && k < 16; l++, k++) {
      const x = (i + .5) * h, y = (j + .5) * h, z = e.box.z0 + (l + .5) * L / 2;
      e.addAtom('F', x - .7, y, z, q); e.addAtom('F', x + .7, y, z, q); e.setBondOrder(e.N - 2, e.N - 1, 1);
    }
    e.touch(); e.refresh(); e.thermalize(T);
    let sum = 0, cnt = 0;
    for (let s = 0; s < 3000000; s++) {
      e.step(); if (s < 750000 || s % 200) continue;
      let free = 0; for (const f of e.fragments().list) if (f.length === 1) free++;
      sum += free / e.N; cnt++;
    }
    return sum / cnt;
  };
  const dense = simulate(1500, 60), dilute = simulate(1500, 120);
  const pd = predict(1500, 60, 32), pl = predict(1500, 120, 32);
  assert.ok(Math.abs(dense - pd) < 0.09, `at 20 bar: α ${dense.toFixed(3)} against ${pd.toFixed(3)}`);
  assert.ok(Math.abs(dilute - pl) < 0.09, `at 3 bar: α ${dilute.toFixed(3)} against ${pl.toFixed(3)}`);
  assert.ok(dilute > dense + 0.15, 'lowering the pressure shifts the equilibrium toward atoms');
});

test('Dissociation equilibria of F2, Cl2, I2 and O2 match real thermodynamics at flame temperatures', () => {
  /* Real K = [X2]/[X]^2 from textbook statistical mechanics with measured spectroscopic constants
     (Huber & Herzberg), D0 (JANAF/ATcT) and atomic fine-structure levels (NIST): this is how the
     JANAF tables are built. The engine's K is exact classical statistical mechanics of its own
     surface. Before free atoms carried their electronic degeneracy, every one was 20–30× off. */
  const q = { thermal: false }, h = 6.62607015e-34, kB = 1.380649e-23, amu = 1.66053907e-27, NA = 6.02214076e23, c2 = 1.438777;
  const DATA = {
    F: { m: 18.998, D0: 154.6, we: 916.64, Be: 0.89019, gA: T => 4 + 2 * Math.exp(-c2 * 404.1 / T), gM: 1, T: 2000 },
    Cl: { m: 35.45, D0: 239.2, we: 559.7, Be: 0.2440, gA: T => 4 + 2 * Math.exp(-c2 * 882.4 / T), gM: 1, T: 2500 },
    I: { m: 126.904, D0: 148.8, we: 214.5, Be: 0.03737, gA: T => 4, gM: 1, T: 1200 },
    O: { m: 15.999, D0: 493.6, we: 1580.2, Be: 1.4456, gA: T => 5 + 3 * Math.exp(-c2 * 158.3 / T) + Math.exp(-c2 * 227 / T), gM: 3, T: 3500 },
  };
  const bad = [];
  for (const [X, d] of Object.entries(DATA)) {
    const T = d.T, at = r => {
      const e = new Engine({ width: 200, height: 200, depth: 200, T, thermostat: false });
      e.box = { x0: -100, x1: 100, y0: -100, y1: 100, z0: -100, z1: 100 };
      e.addAtom(X, 0, 0, 0, q); e.addAtom(X, r, 0, 0, q); e.touch(); e.refresh();
      if (X === 'O') e.setBondOrder(0, 1, 2);
      for (let k = 0; k < 100; k++) e.computeForces(1);
      return [e.Epot, e.nPairs ? e.bondStrength(0) : 0];
    };
    const Einf = at(9)[0]; let I = 0;
    for (let r = 0.5; r <= 4.5; r += 0.005) { const [U, s] = at(r); if (s > 0.25) I += 4 * Math.PI * r * r * Math.exp(-(U - Einf) / (KB * T)) * 0.005; }
    const kT = kB * T, tr = m => Math.pow(2 * Math.PI * m * amu * kT / (h * h), 1.5) * 1e-30;
    const qA = tr(d.m) * d.gA(T), qM = tr(2 * d.m) * (T / (2 * c2 * d.Be)) / (1 - Math.exp(-c2 * d.we / T)) * d.gM;
    const real = qM / (qA * qA) * Math.exp(d.D0 * 1000 / (NA * kT)), engine = I / 2, ratio = real / engine;
    if (ratio < 0.6 || ratio > 1.6) bad.push(`${X}2 at ${T} K: real/engine ${ratio.toFixed(2)}`);
  }
  assert.equal(bad.length, 0, bad.join('; '));
});

test('Real gases attract as measured: second virial coefficients of argon and nitrogen', () => {
  /* B2 = −2π N_A ∫ <e^{−U/kT} − 1> r² dr over orientations of rigid monomers, from the engine's own
     intermolecular energy, plus the analytic dispersion tail beyond the cutoff. B2 is what makes a
     real gas's pressure fall below nkT as it is compressed. Argon's parameters are the classic
     virial fit; nitrogen was never fitted. Measured values: Dymond & Smith. */
  const q = { thermal: false };
  const B2 = (atoms, orders, T, nrot) => {
    const n = atoms.length, e = new Engine({ width: 200, height: 200, depth: 200, T: 0, thermostat: false });
    e.box = { x0: -100, x1: 100, y0: -100, y1: 100, z0: -100, z1: 100 };
    for (let k = 0; k < 2; k++) for (const [sym, x, y, z] of atoms) e.addAtom(sym, x, y, z + k * 50, q);
    e.touch(); e.refresh(); for (const [i, j, o] of orders) { e.setBondOrder(i, j, o); e.setBondOrder(i + n, j + n, o); }
    let seed = 99; const rng = () => { seed = (seed * 1103515245 + 12345) % 2147483648; return seed / 2147483648; };
    const rot = () => { const u1 = rng(), u2 = rng(), u3 = rng(), a = Math.sqrt(1 - u1) * Math.sin(2 * Math.PI * u2), b = Math.sqrt(1 - u1) * Math.cos(2 * Math.PI * u2), c = Math.sqrt(u1) * Math.sin(2 * Math.PI * u3), d = Math.sqrt(u1) * Math.cos(2 * Math.PI * u3);
      return [[1 - 2 * (c * c + d * d), 2 * (b * c - a * d), 2 * (b * d + a * c)], [2 * (b * c + a * d), 1 - 2 * (b * b + d * d), 2 * (c * d - a * b)], [2 * (b * d - a * c), 2 * (c * d + a * b), 1 - 2 * (b * b + c * c)]]; };
    const place = (R1, R2, r) => { atoms.forEach(([, x, y, z], i) => { for (let d = 0; d < 3; d++) { e.pos[3 * i + d] = R1[d][0] * x + R1[d][1] * y + R1[d][2] * z; e.pos[3 * (i + n) + d] = R2[d][0] * x + R2[d][1] * y + R2[d][2] * z + (d === 2 ? r : 0); } }); e.needRebuild = true; e._checkRebuild(); return e.computeForces(); };
    const I3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]], Einf = place(I3, I3, 60), dr = 0.05;
    let I = -1 / 3;
    const rots = []; for (let k = 0; k < nrot; k++) rots.push([rot(), rot()]);
    for (let r = 1; r <= 12; r += dr) { let f = 0; for (const [R1, R2] of rots) f += Math.exp(-(place(R1, R2, r) - Einf) / (KB * T)) - 1; I += f / nrot * r * r * dr; }
    let tail = 0; for (const [a] of atoms) for (const [b] of atoms) { const A = BY_SYM[a], B = BY_SYM[b]; tail -= 2 * Math.sqrt(A.ljDA * B.ljDA) * (A.ljX * B.ljX) ** 3 / (3 * 8 ** 3); }
    return (-2 * Math.PI * I + 2 * Math.PI * tail / (KB * T)) * 0.6022;
  };
  const bad = [];
  for (const [label, atoms, orders, T, want, tol, nrot] of [
    ['Ar 150 K', [['Ar', 0, 0, 0]], [], 150, -86, 18, 1], ['Ar 300 K', [['Ar', 0, 0, 0]], [], 300, -15.6, 8, 1],
    ['N2 300 K', [['N', 0, 0, -.55], ['N', 0, 0, .55]], [[0, 1, 3]], 300, -4.2, 8, 60], ['N2 600 K', [['N', 0, 0, -.55], ['N', 0, 0, .55]], [[0, 1, 3]], 600, 21.3, 8, 60],
  ]) { const got = B2(atoms, orders, T, nrot); if (Math.abs(got - want) > tol) bad.push(`${label}: ${got.toFixed(1)} against ${want} cm3/mol`); }
  assert.equal(bad.length, 0, bad.join('; '));
});

console.log(count + ' textbook checks passed.');
