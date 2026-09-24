/* Shapes. A molecule relaxed from a rough sketch must settle where the real one sits, and a
 * radical must hold its real shape while it lives — the bond angles are what a person reads first.
 *
 * Radicals are the hard part. The angle rule counted bonds and lone pairs and nothing else, so an
 * unpaired electron took no room: CH2 came out straight instead of bent, NH2 at 120° instead of
 * water's 103°, and a CH3 made from that straight CH2 got its third hydrogen at 90° and stayed a T
 * for as long as it lived. */
const assert = require('node:assert/strict');
const { Engine, TUNE } = require('../engine.js');
let count = 0;
function test(name, run) { run(); console.log('PASS', name); count++; }
const q = { thermal: false };
function relax(atoms, bonds) {
  const e = new Engine({ width: 60, height: 60, depth: 60, T: 0, thermostat: false });
  for (const [s, x, y, z] of atoms) e.addAtom(s, 30 + x, 30 + y, z || 0, q);
  e.touch(); e.refresh();
  if (bonds) for (const [i, j, n] of bonds) e.setBondOrder(i, j, n);
  for (let k = 0; k < 4; k++) { for (let s = 0; s < 300; s++) e._updateBondOrders(1); e.minimize(5000, 0.01); e.refresh(); }
  return e;
}
const P = (e, i) => [e.pos[3 * i], e.pos[3 * i + 1], e.pos[3 * i + 2]];
const vec = (e, a, b) => P(e, a).map((v, k) => v - P(e, b)[k]);
const len = v => Math.hypot(...v);
const angle = (e, a, c, b) => { const u = vec(e, a, c), v = vec(e, b, c); return Math.acos((u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / len(u) / len(v)) * 180 / Math.PI; };
const d = 0.63;
const MOL = {
  CH4: [[['C', 0, 0, 0], ['H', d, d, d], ['H', -d, -d, d], ['H', -d, d, -d], ['H', d, -d, -d]]],
  CH3: [[['C', 0, 0, 0], ['H', 1.05, 0.1, 0.2], ['H', -0.5, 0.95, -0.1], ['H', -0.55, -0.9, 0.15]]],
  CH2: [[['C', 0, 0, 0], ['H', -0.9, -0.55, 0.05], ['H', 0.9, -0.55, -0.05]]],
  NH3: [[['N', 0, 0, 0], ['H', 1.0, 0, 0.3], ['H', -0.5, 0.87, 0.3], ['H', -0.5, -0.87, 0.3]]],
  NH2: [[['N', 0, 0, 0], ['H', 0.8, 0.6, 0], ['H', -0.8, 0.6, 0.05]]],
  H2O: [[['O', 0, 0, 0], ['H', 0.76, 0.59, 0], ['H', -0.76, 0.59, 0.05]]],
  C2H4: [[['C', 0, 0, 0], ['C', 1.34, 0, 0], ['H', -0.55, 0.93, 0], ['H', -0.55, -0.93, 0.05], ['H', 1.89, 0.93, 0], ['H', 1.89, -0.93, 0]], [[0, 1, 2]]],
  C2H2: [[['C', 0, 0, 0], ['C', 1.2, 0, 0], ['H', -1.06, 0.02, 0], ['H', 2.26, -0.02, 0]], [[0, 1, 3]]],
  C3H6: [[['C', 0, 0.87, 0], ['C', -0.755, -0.435, 0], ['C', 0.755, -0.435, 0], ['H', 0, 1.46, 0.9], ['H', 0, 1.46, -0.9], ['H', -1.27, -0.73, 0.9], ['H', -1.27, -0.73, -0.9], ['H', 1.27, -0.73, 0.9], ['H', 1.27, -0.73, -0.9]]],
  C6H6: [[['C', 1.4, 0, 0], ['C', 0.7, 1.2124, 0], ['C', -0.7, 1.2124, 0], ['C', -1.4, 0, 0], ['C', -0.7, -1.2124, 0], ['C', 0.7, -1.2124, 0], ['H', 2.48, 0, 0], ['H', 1.24, 2.1477, 0], ['H', -1.24, 2.1477, 0], ['H', -2.48, 0, 0], ['H', -1.24, -2.1477, 0], ['H', 1.24, -2.1477, 0]], [[0, 1, 1.5], [1, 2, 1.5], [2, 3, 1.5], [3, 4, 1.5], [4, 5, 1.5], [5, 0, 1.5]]],
  CH2O: [[['C', 0, 0, 0], ['O', 1.21, 0, 0], ['H', -0.57, 0.94, 0], ['H', -0.57, -0.94, 0.05]], [[0, 1, 2]]],
  CO2: [[['C', 0, 0, 0], ['O', 1.16, 0, 0], ['O', -1.16, 0.05, 0]], [[0, 1, 2], [0, 2, 2]]],
};
const near = (got, want, tol, what) => assert.ok(Math.abs(got - want) <= tol, `${what}: ${got.toFixed(2)} against ${want}`);

test('Closed-shell molecules settle at their real angles and lengths', () => {
  let e = relax(...MOL.CH4); near(angle(e, 1, 0, 2), 109.5, 0.5, 'CH4 angle'); near(len(vec(e, 1, 0)), 1.087, 0.01, 'C–H');
  e = relax(...MOL.NH3); near(angle(e, 1, 0, 2), 106.7, 1.5, 'NH3 angle');
  e = relax(...MOL.H2O); near(angle(e, 1, 0, 2), 104.5, 1, 'H2O angle');
  e = relax(...MOL.C2H4); near(angle(e, 1, 0, 2) + angle(e, 1, 0, 3) + angle(e, 2, 0, 3), 360, 1, 'ethylene carbon is flat'); near(len(vec(e, 0, 1)), 1.339, 0.02, 'C=C');
  e = relax(...MOL.C2H2); near(angle(e, 1, 0, 2), 180, 1, 'acetylene is straight'); near(len(vec(e, 0, 1)), 1.20, 0.02, 'C≡C');
  e = relax(...MOL.C3H6); near(angle(e, 1, 0, 2), 60, 1, 'cyclopropane ring');
  e = relax(...MOL.C6H6); near(angle(e, 1, 0, 5), 120, 1, 'benzene ring');
  e = relax(...MOL.CO2); near(angle(e, 1, 0, 2), 180, 1, 'CO2 is straight');
});

test('Radicals take their real shapes: CH3 flat, CH2 bent, NH2 like water', () => {
  let e = relax(...MOL.CH3);
  near(angle(e, 1, 0, 2) + angle(e, 1, 0, 3) + angle(e, 2, 0, 3), 360, 2, 'CH3 is flat (angles sum)');
  e = relax(...MOL.CH2); const ch2 = angle(e, 1, 0, 2);
  assert.ok(ch2 > 115 && ch2 < 140, `CH2 is bent between singlet 102° and triplet 134°, not straight: ${ch2.toFixed(1)}`);
  e = relax(...MOL.NH2); near(angle(e, 1, 0, 2), 103.4, 3, 'NH2 angle');
});

test('A CH3 made from CH2 + H is trigonal while it lives, not a T', () => {
  let snaps = 0, tee = 0;
  for (let s = 0; s < 4; s++) {
    const e = new Engine({ width: 20, height: 16, depth: 11, T: 300, seed: 40 + s * 17, thermostatMode: 'kelvin' });
    e.addAtom('C', 6, 8, 0, q); e.addAtom('H', 5.01, 7.58, 0, q); e.addAtom('H', 6.99, 7.58, 0, q);
    e.touch(); e.minimize(1500, 0.05); e.refresh();
    e.addAtom('H', 13, 8.4, 0, q); e.touch(); e.refresh(); e.thermalize(300);
    for (let i = 0; i < 25000; i++) {
      e.step();
      if (i % 100) continue;
      const ch3 = e.fragments().list.find(f => f.length === 4); if (!ch3) continue;
      const c = ch3.find(k => e.type[k] === e.type[0]), hs = ch3.filter(k => k !== c);
      const widest = Math.max(angle(e, hs[0], c, hs[1]), angle(e, hs[0], c, hs[2]), angle(e, hs[1], c, hs[2]));
      snaps++; if (widest > 160) tee++;
    }
  }
  assert.ok(snaps > 200, 'CH3 formed');
  assert.ok(tee / snaps < 0.02, `T-shaped ${(100 * tee / snaps).toFixed(0)}% of the time`);
});

test('Insertion leaves every settled molecule as it was', () => {
  for (const name of Object.keys(MOL)) {
    TUNE.insert = 1; const on = relax(...MOL[name]); on.computeForces();
    // at the settled shape nothing is engaged at all, so the energy there is the same either way
    assert.equal(on.nIns, 0, `${name}: insertion engaged at rest`);
    const E1 = on.Epot; TUNE.insert = 0; on.computeForces(); const E0 = on.Epot;
    assert.ok(Math.abs(E1 - E0) < 1e-9, `${name}: ${E1} vs ${E0}`);
    // and relaxing from the same sketch lands on the same minimum, to the minimiser's tolerance
    const off = relax(...MOL[name]); TUNE.insert = 1;
    assert.ok(Math.abs(off.Epot - E1) < 1e-3, `${name}: minimum moved by ${(off.Epot - E1).toExponential(2)} kJ/mol`);
    for (let k = 0; k < 3 * on.N; k++) assert.ok(Math.abs(on.pos[k] - off.pos[k]) < 0.01, `${name}: an atom moved`);
  }
});

console.log(count + ' geometry checks passed.');
