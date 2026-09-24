/* Textbook checks the engine was never tuned against. Each is a law or a measured number that any
 * correct simulation of real matter has to reproduce, so a failure here is a flaw in the model,
 * not in the test.
 *   node playground/tests/textbook.cjs */
const assert = require('node:assert/strict');
const { Engine, KB, KEU } = require('../engine.js');
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

test('Equipartition: helium and xenon in one bath carry the same kinetic energy', () => {
  const e = gas(['He', 'Xe'], 64, 40, 300);
  const kT = KB * 300, acc = [[0, 0], [0, 0]];
  for (let i = 0; i < 400000; i++) {
    e.step(); if (i < 50000 || i % 25) continue;
    for (let a = 0; a < e.N; a++) { const s = e.formulaOf([a]) === 'He' ? 0 : 1; acc[s][0] += 0.5 * KEU * e.mass[a] * (e.vel[3 * a] ** 2 + e.vel[3 * a + 1] ** 2 + e.vel[3 * a + 2] ** 2) / kT; acc[s][1]++; }
  }
  const he = acc[0][0] / acc[0][1], xe = acc[1][0] / acc[1][1];
  assert.ok(Math.abs(he - 1.5) < 0.08 && Math.abs(xe - 1.5) < 0.08, `He ${he.toFixed(3)} kT, Xe ${xe.toFixed(3)} kT`);
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
    // [name, ΔfH, atoms, bond orders, tolerance]
    ['H2', 0, [['H', 0, 0, 0], ['H', .74, 0, 0]], [], 10],
    ['OH', 37.4, [['O', 0, 0, 0], ['H', .97, 0, 0]], [], 10],
    ['H2O', -241.8, [['O', 0, 0, 0], ['H', .76, .59, 0], ['H', -.76, .59, 0]], [], 10],
    ['HO2', 12.3, [['H', 0, 0, 0], ['O', .97, 0, 0], ['O', 1.3, 1.26, 0]], [[1, 2, 1.5]], 15],
    ['CH', 596.4, [['C', 0, 0, 0], ['H', 1.12, 0, 0]], [], 10],
    ['CH2', 391.2, [['C', 0, 0, 0], ['H', .99, .42, 0], ['H', -.99, .42, 0]], [], 10],
    ['CH3', 146.4, [['C', 0, 0, 0], ['H', 1.08, 0, 0], ['H', -.54, .93, 0], ['H', -.54, -.93, 0]], [], 10],
    ['CH4', -74.6, [['C', 0, 0, 0], ['H', .63, .63, .63], ['H', -.63, -.63, .63], ['H', -.63, .63, -.63], ['H', .63, -.63, -.63]], [], 10],
    ['C2H6', -84.0, [['C', 0, 0, 0], ['C', 1.54, 0, 0], ['H', -.36, 1.03, 0], ['H', -.36, -.51, .89], ['H', -.36, -.51, -.89], ['H', 1.9, -1.03, 0], ['H', 1.9, .51, .89], ['H', 1.9, .51, -.89]], [], 10],
    ['C2H5', 119.9, [['C', 0, 0, 0], ['C', 1.5, 0, 0], ['H', -.36, 1.03, 0], ['H', -.36, -.51, .89], ['H', -.36, -.51, -.89], ['H', 2.05, .93, 0], ['H', 2.05, -.93, 0]], [], 30],
    ['CH3OH', -201.0, [['C', 0, 0, 0], ['O', 1.43, 0, 0], ['H', 1.75, .9, 0], ['H', -.36, 1.03, 0], ['H', -.36, -.51, .89], ['H', -.36, -.51, -.89]], [], 10],
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

console.log(count + ' textbook checks passed.');
