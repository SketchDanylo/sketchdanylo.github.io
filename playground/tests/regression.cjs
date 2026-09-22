const assert = require('node:assert/strict');
const { Engine, BY_SYM, TUNE } = require('../engine.js');
const { stepsPerSecond, validMolecule, thermalTranslation } = require('../protocol.js');
const near = (x, y, tol = 1e-9) => assert.ok(Math.abs(x - y) < tol, `${x} vs ${y}`);
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
test('Valence sharing leaves every atom at or under its valence untouched', () => {
  // equilibrium structures never over-claim, so sharing must be exactly inert there
  const build = share => {
    TUNE.share = share;
    const e = new Engine({ width: 40, height: 30, depth: 20, T: 0, wallT: 0, thermostat: false, seed: 4 });
    e.addAtom('C', 20, 15, 0, { thermal: false });
    for (const [dx, dy, dz] of [[.63, .63, .63], [-.63, -.63, .63], [-.63, .63, -.63], [.63, -.63, -.63]]) e.addAtom('H', 20 + dx, 15 + dy, dz, { thermal: false });
    e.refresh();
    for (let s = 0; s < 8000; s++) { e.step(); e.vel.fill(0); }
    e.refresh();
    return { E: e.Epot, r: Math.hypot(e.pos[3] - e.pos[0], e.pos[4] - e.pos[1], e.pos[5] - e.pos[2]) };
  };
  const off = build(0), on = build(0.5);
  TUNE.share = 0.5;
  near(on.E, off.E, 1e-9);
  near(on.r, off.r, 1e-12);
});
test('Valence sharing conserves total bonding through a handover', () => {
  // one hydrogen between its old partner and an incoming carbon: it cannot give a whole bond to
  // both, and what it does give must add up to about one bond rather than collapsing
  const crossing = share => {
    TUNE.share = share;
    const e = new Engine({ width: 40, height: 30, depth: 20, T: 0, wallT: 0, thermostat: false, seed: 4 });
    e.addAtom('C', 20, 15, 0, { thermal: false });
    e.addAtom('H', 21.08, 15, 0, { thermal: false });
    e.addAtom('H', 19.46, 15.94, 0, { thermal: false });
    e.addAtom('H', 20, 13.9, 0, { thermal: false });
    e.addAtom('H', 20, 13.16, 0, { thermal: false });
    e.refresh();
    for (let s = 0; s < 400; s++) e._updateBondOrders(1);
    e.computeForces();
    let cs = 0, hs = 0;
    for (const b of e.bonds(0.0001)) {
      const k = [b.i, b.j].sort().join('-');
      if (k === '0-3') cs = b.strength;
      if (k === '3-4') hs = b.strength;
    }
    return cs + hs;
  };
  const off = crossing(0), on = crossing(0.5);
  TUNE.share = 0.5;
  assert.ok(off < 0.7, `without sharing the handover collapses, total ${off.toFixed(2)}`);
  assert.ok(on > 0.75 && on < 1.25, `with sharing it is conserved, total ${on.toFixed(2)}`);
});
test('Forces stay the exact gradient of the energy with sharing active', () => {
  for (const share of [0, 0.5, 1]) {
    TUNE.share = share;
    const e = new Engine({ width: 40, height: 30, depth: 20, T: 0, wallT: 0, thermostat: false, seed: 2 });
    e.addAtom('C', 20, 15, 0, { thermal: false });
    e.addAtom('H', 21.08, 15, 0, { thermal: false });
    e.addAtom('H', 19.46, 15.94, 0, { thermal: false });
    e.addAtom('H', 20.1, 13.55, 0.2, { thermal: false });   // deliberately over-coordinated
    e.addAtom('H', 20.0, 12.82, -0.1, { thermal: false });
    e.refresh(); e.computeForces();
    const analytic = [...e.frc.slice(0, 3 * e.N)], h = 2e-5;
    let worst = 0;
    for (let k = 0; k < 3 * e.N; k++) {
      const p0 = e.pos[k];
      e.pos[k] = p0 + h; e.touch(); e.refresh(); e.computeForces(); const Ep = e.Epot;
      e.pos[k] = p0 - h; e.touch(); e.refresh(); e.computeForces(); const Em = e.Epot;
      e.pos[k] = p0; e.touch(); e.refresh(); e.computeForces();
      worst = Math.max(worst, Math.abs(-(Ep - Em) / (2 * h) - analytic[k]));
    }
    const scale = Math.max(...analytic.map(Math.abs));
    assert.ok(worst / scale < 1e-5, `share ${share}: relative gradient error ${(worst / scale).toExponential(2)}`);
  }
  TUNE.share = 0.5;
});
console.log(`${passed} regression checks passed.`);
