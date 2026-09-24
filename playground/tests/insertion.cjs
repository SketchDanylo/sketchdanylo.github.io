/* Carbene insertion: CH2 + H2 -> CH4.
 *
 * Methylene with two free valences can take a sigma bond's own electron pair into its empty
 * orbital and bond to both ends at once, which is why singlet CH2 inserts into H2 with no barrier.
 * Before this the engine could not take a hydrogen back out of an H2 molecule at all: the carbon
 * stalled at 2.2 A even at 5000 K, and CH2 + H2 never made methane by any route. These lock in
 * that it does now, that the forces stay exact, and that nothing which should not react does. */
const assert = require('node:assert/strict');
const { Engine, TUNE } = require('../engine.js');
let count = 0;
function test(name, run) { run(); console.log('PASS', name); count++; }
const q = { thermal: false };
const products = e => {
  const c = {};
  for (const f of e.fragments().list) { const s = e.formulaOf(f); c[s] = (c[s] || 0) + 1; }
  return Object.keys(c).sort().map(k => (c[k] > 1 ? c[k] + ' ' : '') + k).join(' + ');
};
const CH2 = (e, x, y) => { e.addAtom('C', x, y, 0, q); e.addAtom('H', x - .99, y - .42, 0, q); e.addAtom('H', x + .99, y - .42, 0, q); };
const CH3 = (e, x, y) => { e.addAtom('C', x, y, 0, q); e.addAtom('H', x + 1.08, y, 0, q); e.addAtom('H', x - .54, y + .93, 0, q); e.addAtom('H', x - .54, y - .93, 0, q); };
const CH4 = (e, x, y) => { const d = .63; e.addAtom('C', x, y, 0, q); for (const [a, b, c] of [[d, d, d], [-d, -d, d], [-d, d, -d], [d, -d, -d]]) e.addAtom('H', x + a, y + b, c, q); };
const H2 = (e, x, y) => { e.addAtom('H', x, y, 0, q); e.addAtom('H', x + .74, y, 0, q); e.setBondOrder(e.N - 2, e.N - 1, 1); };
const CO = (e, x, y) => { e.addAtom('C', x, y, 0, q); e.addAtom('O', x + 1.13, y, 0, q); e.setBondOrder(e.N - 2, e.N - 1, 3); };
const CO2 = (e, x, y) => { e.addAtom('C', x, y, 0, q); e.addAtom('O', x + 1.16, y, 0, q); e.addAtom('O', x - 1.16, y, 0, q); e.setBondOrder(e.N - 3, e.N - 2, 2); e.setBondOrder(e.N - 3, e.N - 1, 2); };
function runs(build, seeds, steps = 60000, opts = {}) {
  const t = {};
  for (let s = 0; s < seeds; s++) {
    const e = new Engine({ width: 20, height: 16, depth: 11, T: 300, seed: 40 + s * 17, thermostatMode: 'kelvin', ...opts });
    build(e); e.touch(); e.minimize(1500, 0.05); e.refresh(); e.thermalize(300);
    for (let i = 0; i < steps; i++) e.step();
    const p = products(e); t[p] = (t[p] || 0) + 1;
  }
  return t;
}

test('Forces are the exact gradient of the energy through an insertion', () => {
  const geoms = [2.3, 2.0, 1.7, 1.4].map(d => [['C', 0, 0, 0], ['H', -0.9, -0.75, 0.05], ['H', 0.9, -0.75, -0.04], ['H', -0.37, d, 0.02], ['H', 0.37, d + 0.03, -0.03]]);
  geoms.push([['Si', 0, 0, 0], ['H', -1.2, -0.9, 0], ['H', 1.2, -0.9, 0.1], ['H', -0.37, 2.0, 0], ['H', 0.37, 2.02, 0.05]]);
  for (const atoms of geoms) {
    const e = new Engine({ width: 200, height: 200, depth: 200, T: 0, thermostat: false });
    e.box = { x0: -100, x1: 100, y0: -100, y1: 100, z0: -100, z1: 100 };
    for (const [s, x, y, z] of atoms) e.addAtom(s, x, y, z, q);
    e.touch(); e.refresh(); e.setBondOrder(3, 4, 1); e.refresh();
    e.computeForces();
    assert.ok(e.nIns > 0, 'the insertion is engaged at this geometry');
    const F = e.frc.slice(0, 3 * e.N), h = 1e-5;
    let worst = 0, fmax = 0;
    for (let k = 0; k < 3 * e.N; k++) {
      const x0 = e.pos[k];
      e.pos[k] = x0 + h; e.needRebuild = true; e._checkRebuild(); const Ep = e.computeForces();
      e.pos[k] = x0 - h; e.needRebuild = true; e._checkRebuild(); const Em = e.computeForces();
      e.pos[k] = x0;
      worst = Math.max(worst, Math.abs(-(Ep - Em) / (2 * h) - F[k])); fmax = Math.max(fmax, Math.abs(F[k]));
    }
    assert.ok(worst < 1e-4 * fmax, `force error ${worst} against ${fmax}`);
  }
});

test('Settled molecules are untouched: switching insertion off changes no energy', () => {
  const build = [CH4, CO2, CO, (e, x, y) => { CH2(e, x, y); }, H2];
  for (const b of build) {
    const at = on => {
      TUNE.insert = on;
      const e = new Engine({ width: 60, height: 60, depth: 60, T: 0, thermostat: false });
      b(e, 30, 30); e.touch(); e.minimize(2000, 0.02); e.refresh();
      return e.Epot;
    };
    const on = at(1), off = at(0);
    TUNE.insert = 1;
    assert.ok(Math.abs(on - off) < 1e-6, `${on} vs ${off}`);
  }
});

test('CH2 + H2 makes methane at room temperature', () => {
  const t = runs(e => { CH2(e, 6, 8); H2(e, 13, 8.4); }, 8);
  assert.ok((t['CH4'] || 0) >= 7, JSON.stringify(t));
});

test('...and without insertion it never did', () => {
  TUNE.insert = 0;
  const t = runs(e => { CH2(e, 6, 8); H2(e, 13, 8.4); }, 6);
  TUNE.insert = 1;
  assert.equal(t['CH2 + H2'], 6, JSON.stringify(t));
});

test('A bare carbon atom and CH insert too: C + 2 H2 makes methane by way of CH2', () => {
  assert.ok((runs(e => { e.addAtom('C', 10, 8, 0, q); H2(e, 4, 5); H2(e, 15, 11); }, 6)['CH4'] || 0) >= 5);
  assert.ok((runs(e => { e.addAtom('C', 8, 8, 0, q); e.addAtom('H', 8, 6.88, 0, q); H2(e, 14, 8.4); }, 6)['CH3'] || 0) >= 5);
});

test('Carbide cannot: C4- has a closed shell and nothing to share, and making CH4 from H2 would strand four electrons', () => {
  assert.equal(runs(e => { e.addAtom('C', 10, 8, 0, { thermal: false, charge: -4 }); H2(e, 4, 5); H2(e, 15, 11); }, 3)['C4− + 2 H2'], 3);
});

test('A radical cannot insert: CH3 + H2 and H + H2 are left alone', () => {
  assert.equal(runs(e => { CH3(e, 6, 8); H2(e, 13, 8.4); }, 6)['CH3 + H2'], 6);
  assert.equal(runs(e => { e.addAtom('H', 6, 8, 0, q); H2(e, 13, 8.4); }, 6)['H + H2'], 6);
});

test('Carbon beside a multiple bond is not a carbene: CO + H2 and CO2 + H2 are left alone', () => {
  assert.equal(runs(e => { CO(e, 6, 8); H2(e, 13, 8.4); }, 6)['CO + H2'], 6);
  assert.equal(runs(e => { CO2(e, 6, 8); H2(e, 13, 8.4); }, 6)['CO2 + H2'], 6);
});

test('Saturated molecules stay what they are', () => {
  assert.equal(runs(e => { CH4(e, 6, 8); H2(e, 13, 8.4); }, 6)['CH4 + H2'], 6);
  assert.equal(runs(e => { CH4(e, 6, 8); CH4(e, 13, 8); }, 6)['2 CH4'], 6);
});

console.log(count + ' insertion checks passed.');
