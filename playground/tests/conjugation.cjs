const assert = require('node:assert/strict');
const { Engine } = require('../engine.js');
let count = 0;
function test(name, run) { run(); console.log('PASS', name); count++; }
const q = { thermal: false };

function relaxed(atoms, orders = []) {
  const e = new Engine({ width: 200, height: 200, depth: 200, T: 0, thermostat: false });
  e.box = { x0: -100, x1: 100, y0: -100, y1: 100, z0: -100, z1: 100 };
  for (const [s, x, y, z] of atoms) e.addAtom(s, x, y, z || 0, q);
  e.touch(); e.refresh(); for (const [i, j, n] of orders) e.setBondOrder(i, j, n); e.refresh();
  for (let k = 0; k < 3; k++) e.minimize(4000, 0.02);
  e.refresh(); for (let k = 0; k < 300; k++) e.computeForces(1);
  return e;
}
const order = (e, i, j) => e.pN[e.pairMap.get(Math.min(i, j) * 1048576 + Math.max(i, j))];
const dist = (e, i, j) => Math.hypot(e.pos[3 * i] - e.pos[3 * j], e.pos[3 * i + 1] - e.pos[3 * j + 1], e.pos[3 * i + 2] - e.pos[3 * j + 2]);

const TRANS = [['C', -1.8, -.4, 0], ['C', -.7, .35, 0], ['C', .7, -.35, 0], ['C', 1.8, .4, 0], ['H', -2.75, .1, 0], ['H', -1.8, -1.48, 0], ['H', -.7, 1.43, 0], ['H', .7, -1.43, 0], ['H', 2.75, -.1, 0], ['H', 1.8, 1.48, 0]];
const CIS = [['C', -1.45, 1.0, 0], ['C', -.73, -.1, 0], ['C', .73, -.1, 0], ['C', 1.45, 1.0, 0], ['H', -2.53, .95, 0], ['H', -1.0, 1.98, 0], ['H', -1.25, -1.05, 0], ['H', 1.25, -1.05, 0], ['H', 2.53, .95, 0], ['H', 1.0, 1.98, 0]];
const BENZENE = [['C', 1.4, 0, 0], ['C', .7, 1.2124, 0], ['C', -.7, 1.2124, 0], ['C', -1.4, 0, 0], ['C', -.7, -1.2124, 0], ['C', .7, -1.2124, 0], ['H', 2.48, 0, 0], ['H', 1.24, 2.1477, 0], ['H', -1.24, 2.1477, 0], ['H', -2.48, 0, 0], ['H', -1.24, -2.1477, 0], ['H', 1.24, -2.1477, 0]];

test('Butadiene alternates: C=C 1.34 Å and C–C 1.47 Å, not three equal bonds', () => {
  const e = relaxed(TRANS, [[0, 1, 2], [2, 3, 2]]);
  assert.ok(order(e, 0, 1) > 1.75 && order(e, 1, 2) < 1.25, `orders ${order(e, 0, 1).toFixed(2)} / ${order(e, 1, 2).toFixed(2)}`);
  assert.ok(Math.abs(dist(e, 0, 1) - 1.338) < 0.03 && Math.abs(dist(e, 1, 2) - 1.467) < 0.05, `lengths ${dist(e, 0, 1).toFixed(3)} / ${dist(e, 1, 2).toFixed(3)}`);
  assert.equal(e.isRadical([...Array(e.N).keys()]), false);
});

test('s-trans butadiene is more stable than s-cis (exp. ≈ 12 kJ/mol)', () => {
  const d = relaxed(TRANS, [[0, 1, 2], [2, 3, 2]]).Epot - relaxed(CIS, [[0, 1, 2], [2, 3, 2]]).Epot;
  assert.ok(d < -4 && d > -20, `s-trans − s-cis ${d.toFixed(1)} kJ/mol`);
});

test('Benzene stays aromatic: six equal bonds of order 1.5', () => {
  const e = relaxed(BENZENE);
  for (let k = 0; k < 6; k++) assert.ok(Math.abs(order(e, k, (k + 1) % 6) - 1.5) < 0.01, `bond ${k}: ${order(e, k, (k + 1) % 6).toFixed(3)}`);
});

test('Butadiene sits unreacted at room temperature: no radicals, no lost hydrogens', () => {
  const L = 22, e = new Engine({ width: L, height: L, depth: L, T: 300, seed: 3, thermostatMode: 'kelvin' });
  const m = relaxed(TRANS, [[0, 1, 2], [2, 3, 2]]), g = [];
  for (let i = 0; i < m.N; i++) g.push([TRANS[i][0], m.pos[3 * i], m.pos[3 * i + 1], m.pos[3 * i + 2]]);
  for (const [cx, cy, cz] of [[6, 6, 6], [16, 6, 16], [6, 16, 16], [16, 16, 6]]) {
    const b = e.N;
    for (const [s, x, y, z] of g) e.addAtom(s, cx + x, cy + y, e.box.z0 + cz + z, q);
    e.setBondOrder(b, b + 1, 2); e.setBondOrder(b + 2, b + 3, 2);
  }
  e.touch(); e.refresh(); e.minimize(300, 1); e.refresh(); e.thermalize(300);
  for (let i = 0; i < 5000; i++) e.step();
  const frags = e.fragments().list;
  assert.equal(frags.length, 4, frags.map(f => e.formulaOf(f)).join(', '));
  for (const f of frags) { assert.equal(e.formulaOf(f), 'C4H6'); assert.equal(e.isRadical(f), false); }
});

console.log(count + ' conjugation checks passed.');
