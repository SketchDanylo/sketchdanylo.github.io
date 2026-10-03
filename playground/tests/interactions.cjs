const assert = require('node:assert/strict');
const { Engine } = require('../engine.js');
const { geo } = require('./reactions.cjs');
let count = 0;
function test(name, run) { run(); console.log('PASS', name); count++; }
const q = { thermal: false };
const ledger = e => e.Epot + e.kinetic() + e.voidHeat - e.servoWorkTotal - (e.kelvinWork || 0);
const frags = e => e.fragments().list.map(f => e.formulaOf(f)).sort().join(' + ');
const outside = e => {
  const b = e.box; let n = 0;
  for (let i = 0; i < e.N; i++) {
    const x = e.pos[3 * i], y = e.pos[3 * i + 1], z = e.pos[3 * i + 2];
    if (x < b.x0 - 0.1 || x > b.x1 + 0.1 || y < b.y0 - 0.1 || y > b.y1 + 0.1 || z < b.z0 - 0.1 || z > b.z1 + 0.1) n++;
  }
  return n;
};
function place(e, name, cx, cy, cz, v, flip) {
  const base = e.N;
  for (const [s, x, y, z] of geo(name)) {
    const Y = flip ? (z || 0) : y, Z = flip ? y : (z || 0);
    e.addAtom(s, cx + x, cy + Y, cz + Z, v ? { v } : q);
  }
  return base;
}
function dragTo(e, h, tx, ty, fs, each) {
  e.tweezer = { i: h, x: tx, y: ty };
  for (let s = 0; s < fs / e.dt; s++) { e.step(); if (each) each(); }
  e.tweezer = null;
}

test('A molecule thrown into a velocity-void wall stops there, and every joule it loses is in voidHeat', () => {
  const e = new Engine({ width: 30, height: 20, depth: 12, T: 0, thermostat: false, voidVelocity: true });
  place(e, 'CH4', 22, 10, e.box.z0 + 6, [0.03, 0, 0]);
  e.touch(); e.refresh(); for (let i = 0; i < 300; i++) e.computeForces(1);
  for (let k = 0; k < 4; k++) e.step();
  const L0 = ledger(e); let worst = 0, out = 0;
  for (let i = 0; i < 4000; i++) { e.step(); worst = Math.max(worst, Math.abs(ledger(e) - L0)); out = Math.max(out, outside(e)); }
  assert.ok(worst < 3, `energy not accounted for: ${worst.toFixed(2)} kJ/mol`);
  assert.equal(out, 0);
  assert.equal(frags(e), 'CH4');
  assert.ok(e.voidHeat > 50, `the wall took only ${e.voidHeat.toFixed(1)} kJ/mol`);
});

test('A reacting hot gas in a velocity-void chamber keeps its books', () => {
  const e = new Engine({ width: 16, height: 16, depth: 12, T: 3000, seed: 5, thermostat: false, voidVelocity: true });
  for (let k = 0; k < 6; k++) place(e, 'H2', 3 + (k % 3) * 5, 4 + Math.floor(k / 3) * 8, e.box.z0 + 4);
  for (let k = 0; k < 3; k++) place(e, 'O2', 3 + k * 5, 8, e.box.z0 + 8);
  e.touch(); e.refresh(); e.thermalize(3000); for (let k = 0; k < 60; k++) e.step();
  const L0 = ledger(e); let worst = 0, out = 0;
  for (let i = 0; i < 6000; i++) { e.step(); worst = Math.max(worst, Math.abs(ledger(e) - L0)); out = Math.max(out, outside(e)); }
  assert.ok(worst < 10, `energy not accounted for: ${worst.toFixed(2)} kJ/mol`);
  assert.equal(out, 0);
});

test('Dragging a molecule hard into a velocity-void wall presses it there without tearing it or losing track', () => {
  const e = new Engine({ width: 30, height: 20, depth: 12, T: 0, thermostat: false, voidVelocity: true });
  const b = place(e, 'C2H6', 12, 10, e.box.z0 + 6);
  e.touch(); e.refresh(); e.minimize(300, .5); e.refresh();
  const L0 = ledger(e); let worst = 0, out = 0;
  dragTo(e, b, 40, 10, 3000, () => { worst = Math.max(worst, Math.abs(ledger(e) - L0)); out = Math.max(out, outside(e)); });
  for (let i = 0; i < 1000; i++) { e.step(); worst = Math.max(worst, Math.abs(ledger(e) - L0)); }
  let xmax = -Infinity; for (let i = 0; i < e.N; i++) xmax = Math.max(xmax, e.pos[3 * i]);
  assert.ok(worst < 5, `energy not accounted for: ${worst.toFixed(2)} kJ/mol`);
  assert.equal(out, 0);
  assert.equal(frags(e), 'C2H6');
  assert.ok(xmax > 29, `it never reached the wall: ${xmax.toFixed(2)}`);
});

test('A radical dragged onto a cold double bond adds to it, the hand paying for the π bond giving way', () => {
  const e = new Engine({ width: 30, height: 20, depth: 12, T: 0, thermostat: false });
  place(e, 'C2H4', 15, 10, e.box.z0 + 6, null, true);
  const cl = e.addAtom('Cl', 15, 15, e.box.z0 + 6, q);
  e.touch(); e.refresh(); e.setBondOrder(0, 1, 2); e.refresh();
  for (let k = 0; k < 100; k++) e.step();
  const L0 = ledger(e); let worst = 0;
  e.tweezer = { i: cl, x: e.pos[3 * cl], y: e.pos[3 * cl + 1] };
  for (let s = 0; s < 2500 / e.dt; s++) {
    e.tweezer.x = e.pos[3 * cl] + (15 - e.pos[3 * cl]) * .02;
    e.tweezer.y = e.pos[3 * cl + 1] + (12 - e.pos[3 * cl + 1]) * .02;
    e.step(); worst = Math.max(worst, Math.abs(ledger(e) - L0));
  }
  e.tweezer = null;
  for (let i = 0; i < 2000; i++) { e.step(); worst = Math.max(worst, Math.abs(ledger(e) - L0)); }
  assert.equal(frags(e), 'C2H4Cl');
  assert.ok(worst < 25, `energy not accounted for: ${worst.toFixed(2)} kJ/mol`);
});

test('A molecule dragged through a bath of others arrives, and nothing it passes is damaged', () => {
  const e = new Engine({ width: 34, height: 24, depth: 12, T: 300, seed: 3, thermostatMode: 'kelvin' });
  const first = place(e, 'C6H6', 6, 6, e.box.z0 + 6);
  for (let k = 0; k < 5; k++) place(e, 'C6H6', 5 + k * 6, 18, e.box.z0 + 6);
  e.touch(); e.refresh(); e.minimize(300, 1); e.refresh(); e.thermalize(300);
  for (let k = 0; k < 500; k++) e.step();
  dragTo(e, first, 28, 6, 3000);
  let cx = 0, cy = 0, M = 0;
  for (let i = first; i < first + 12; i++) { cx += e.pos[3 * i] * e.mass[i]; cy += e.pos[3 * i + 1] * e.mass[i]; M += e.mass[i]; }
  assert.ok(Math.hypot(cx / M - 28, cy / M - 6) < 2, 'the dragged benzene did not arrive');
  assert.ok(e.fragments().list.every(f => e.formulaOf(f) === 'C6H6'), frags(e));
});

test('A saved scene comes back with its double and aromatic bonds, not as single bonds', () => {
  const e = new Engine({ width: 24, height: 20, depth: 12, T: 300, seed: 1, thermostatMode: 'kelvin' });
  place(e, 'C6H6', 8, 10, e.box.z0 + 6); place(e, 'O2', 18, 10, e.box.z0 + 6);
  e.touch(); e.refresh(); e.setBondOrder(12, 13, 2); e.refresh(); e.thermalize(300);
  for (let i = 0; i < 500; i++) e.step();
  const json = JSON.parse(JSON.stringify(e.toJSON()));
  const f = new Engine({ width: 24, height: 20, depth: 12, T: 300, seed: 1, thermostatMode: 'kelvin' }); f.box = { ...json.box };
  const index = json.atoms.map(a => f.addAtom(a[0], a[1], a[2], a[3], { v: [a[4], a[5], a[6]], charge: a[7], V: a[8] }));
  f.touch(); f.refresh(); for (const [i, j, o] of json.bonds) f.setBondOrder(index[i], index[j], o); f.refresh();
  assert.equal(json.bonds.length, 7);
  let worst = 0;
  for (let p = 0; p < e.nPairs; p++) { const q = f.pairMap.get(e.pI[p] * 1048576 + e.pJ[p]); if (q !== undefined) worst = Math.max(worst, Math.abs(e.pN[p] - f.pN[q])); }
  assert.ok(worst < 1e-3, `bond orders changed by ${worst}`);
  assert.ok(Math.abs(e.computeForces() - f.computeForces()) < 0.5);
});

console.log(count + ' interaction checks passed.');
