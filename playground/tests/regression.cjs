const assert = require('node:assert/strict');
const { Engine, BY_SYM } = require('../engine.js');
const { stepsPerSecond, validMolecule, thermalTranslation } = require('../protocol.js');
let passed = 0;
function test(name, fn) { fn(); passed++; console.log('PASS', name); }
function engine(T = 300) {
  const e = new Engine({ T, width: 50, height: 50, depth: 50 });
  e.box = { x0: -25, x1: 25, y0: -25, y1: 25, z0: -25, z1: 25 }; return e;
}
function close(a, b, tolerance = 1e-10) { assert.ok(Math.abs(a - b) < tolerance, `${a} ≠ ${b}`); }
test('Playback scale uses femtoseconds, not nanoseconds', () => {
  close(stepsPerSecond(.1), 1); close(stepsPerSecond(1), 20000); close(stepsPerSecond(5), 100000);
  close(stepsPerSecond(1) * 1e-15, 20e-12, 1e-25);
  assert.equal(stepsPerSecond(NaN), 0);
});
test('Force queries do not evolve bond order or energy', () => {
  const e = engine(0); e.addAtom('O', 0, 0, 0); e.addAtom('O', 1.3, 0, 0); e.refresh();
  e.setBondOrder(0, 1, 1.6); const E = e.computeForces(), n = e.pN[0];
  for (let i = 0; i < 50; i++) close(e.computeForces(), E);
  assert.equal(e.pN[0], n);
});
test('Warm imported molecules receive translational motion at the molecular mass', () => {
  const waterMass = 18.015, T = 298.15, v = thermalTranslation(waterMass, T, () => 1);
  close(waterMass * 1e4 * v[0] ** 2, Engine.UNITS.KB * T);
  assert.ok(v.every(x => x > 0));
  const heavy = thermalTranslation(waterMass * 4, T, () => 1);
  close(heavy[0], v[0] / 2);
  assert.deepEqual(thermalTranslation(waterMass, 0, () => { throw new Error('No RNG required at 0 K'); }), [0, 0, 0]);
});
test('Polarization conserves total formal charge', () => {
  const e = engine(0);
  e.addAtom('O', 0, 0, 0); e.addAtom('H', .96, 0, 0); e.addAtom('H', -.24, .93, 0); e.addAtom('Na', 3.5, 0, 0, { charge: 1 });
  e.refresh(); close(e.q.slice(0, e.N).reduce((a, b) => a + b, 0), 1);
  for (let n = 0; n < 80; n++) { e.step(); close(e.q.slice(0, e.N).reduce((a, b) => a + b, 0), 1); }
});
test('Forces include changing polar hydrogen radius', () => {
  const e = engine(0); e.thermostat = false;
  e.addAtom('H', 0, 0, 0); e.addAtom('O', 1.28, .1, .2); e.addAtom('Ar', -.3, 2.7, .2); e.refresh();
  assert.ok(e.ljq[0] !== 0, 'test must exercise the variable-radius region');
  const F = e.frc.slice(0, 9), h = 1e-5;
  for (let k = 0; k < 9; k++) {
    const x = e.pos[k]; e.pos[k] = x + h; const plus = e.computeForces();
    e.pos[k] = x - h; const minus = e.computeForces(); e.pos[k] = x;
    close(F[k], -(plus - minus) / (2 * h), 1e-3);
  }
});
test('Hot trajectory restores all state and replays across checkpoints', () => {
  const e = engine(3500);
  e.addAtom('O', 0, 0, 0); e.addAtom('H', .96, 0, 0); e.addAtom('H', -.24, .93, 0); e.refresh();
  for (let i = 0; i < 195; i++) e.step();
  const expected = e.snapshot();
  for (let i = 0; i < 140; i++) e.step();
  assert.equal(e.stepBack(140), true);
  const actual = e.snapshot();
  for (const key of ['pos', 'vel', 'frc', 'pN', 'rng', 'nextSub', 'Epot', 'pressureEMA', 'stepCount']) assert.deepEqual(actual[key], expected[key], key);
  const branch = new Engine(); branch.restore(expected);
  for (let i = 0; i < 40; i++) { e.step(); branch.step(); }
  assert.deepEqual(e.pos.slice(0, 9), branch.pos.slice(0, 9));
});
test('Pinned atoms stay motionless through heating and momentum removal', () => {
  const e = engine(); e.addAtom('Ar', 0, 0, 0); e.addAtom('Ar', 5, 0, 0); e.pinned[0] = 1;
  e.thermalize(1000); e.zeroMomentum();
  assert.deepEqual([...e.vel.slice(0, 3)], [0, 0, 0]);
  assert.throws(() => e.thermalize(-1), RangeError);
});
test('Numerical failures stop without advancing the clock', () => {
  const e = engine(); e.addAtom('Ar', 0, 0, 0, { v: [NaN, 0, 0] }); e.refresh();
  const before = e.pos.slice(0, 3); assert.throws(() => e.step(), /numerical failure/);
  assert.equal(e.time, 0); assert.deepEqual(e.pos.slice(0, 3), before);
});
test('Nomenclature schema rejects malformed topology and non-finite positions', () => {
  const valid = { format: 'chem-playground/molecule@1', id: 'test-water', units: 'angstrom', formula: 'H2O', atoms: [{ el: 'O', x: 0, y: 0 }, { el: 'H', x: .96, y: 0 }, { el: 'H', x: -.24, y: .93 }], bonds: [{ a: 0, b: 1, order: 1 }, { a: 0, b: 2, order: 1 }] };
  assert.ok(validMolecule(valid, BY_SYM));
  for (const mutate of [m => m.atoms[0].x = NaN, m => m.bonds[0].b = 80, m => m.units = 'nm', m => m.atoms[0].el = '__proto__', m => m.bonds.push(m.bonds[0]), m => m.atoms[0].charge = 1.5]) {
    const m = structuredClone(valid); mutate(m); assert.equal(validMolecule(m, BY_SYM), false);
  }
});
console.log(`${passed} regression checks passed.`);
