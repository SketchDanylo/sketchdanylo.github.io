const assert = require('node:assert/strict');
const { Engine, KB } = require('../engine.js');
const near = (a, b, tol = 1e-9) => assert.ok(Math.abs(a - b) < tol, `${a} vs ${b}`);
let count = 0;
function test(name, run) { run(); console.log('PASS', name); count++; }
function chamber(options = {}) {
  return new Engine({ width: 100, height: 100, depth: 100, T: 300, wallT: 300, ...options });
}
test('A setpoint heats the wall gradually and leaves the interior velocities unchanged', () => {
  const e = chamber();
  e.addAtom('Ar', 50, 50, 0, { v: [.01, -.02, .03] });
  const before = [...e.vel.slice(0, 3)];
  e.setTemperature(623.15);
  assert.equal(e.wallT, 300); assert.deepEqual([...e.vel.slice(0, 3)], before);
  for (let i = 0; i < 10000; i++) e._wallBath(1);
  near(e.wallT, 623.15 + (300 - 623.15) / Math.E, 1e-8);
  assert.deepEqual([...e.vel.slice(0, 3)], before);
  assert.equal(e.heatToSample, 0);
});
test('Wall response is independent of how elapsed simulation time is partitioned', () => {
  const a = chamber(), b = chamber(); a.setTemperature(500); b.setTemperature(500);
  for (let i = 0; i < 1000; i++) a._wallBath(1);
  b._wallBath(1000); near(a.wallT, b.wallT);
});
test('Boundary heat has an equal and opposite reservoir energy change', () => {
  const e = chamber(); e.addAtom('Ar', .1, 50, 0, { thermal: false });
  const E0 = e.kinetic() + e.wallT * e.wallCapacity;
  for (let i = 0; i < 1000; i++) e._wallBath(1);
  near(e.heatToSample, e.kinetic(), 1e-10);
  near(e.kinetic() + e.wallT * e.wallCapacity - E0, e.heaterWork, 1e-8);
  assert.ok(e.kinetic() > 0);
});
test('Stationary wall coupling reproduces Maxwell velocity variance at the wall', () => {
  const e = chamber({ wallT: 300, wallTarget: 300, wallCapacity: 1e14 });
  e.addAtom('Ar', 0, 50, 0, { thermal: false });
  let v = 0, v2 = 0, n = 0;
  for (let i = 0; i < 31000; i++) {
    e._wallBath(300);
    if (i < 1000) continue;
    for (let d = 0; d < 3; d++) { v += e.vel[d]; v2 += e.vel[d] ** 2; n++; }
  }
  const expected = KB * 300 / (e.mass[0] * 1e4);
  near(v2 / n / expected, 1, .025); near(v / n / Math.sqrt(expected), 0, .02);
});
test('Thermostat-off edits immediately set kinetic temperature, including from rest', () => {
  const e = chamber({ thermostat: false });
  for (let i = 0; i < 10; i++) e.addAtom('Ar', 4 * i, 30, 0, { thermal: false });
  e.pinned[0] = 1;
  e.setTemperature(273.15); near(e.temperature(), 273.15);
  e.setTemperature(1200); near(e.temperature(), 1200);
  assert.deepEqual([...e.vel.slice(0, 3)], [0, 0, 0]);
  e.setTemperature(0); assert.equal(e.kinetic(), 0);
  assert.throws(() => e.setTemperature(-1), RangeError);
  const empty = chamber({ thermostat: false }); empty.setTemperature(300); assert.equal(empty.kinetic(), 0);
});
test('Wall temperature target is constrained to 15–350 °C', () => {
  const e = chamber(); e.setTemperature(0); near(e.wallTarget, 288.15);
  e.setTemperature(5000); near(e.wallTarget, 623.15);
});
test('Wall history reproduces velocities, reservoir state and heat on replay', () => {
  const e = chamber(); e.addAtom('Ar', .5, 50, 0); e.refresh(); e.setTemperature(500);
  for (let i = 0; i < 180; i++) e.step();
  const expected = e.snapshot();
  for (let i = 0; i < 100; i++) e.step();
  assert.ok(e.stepBack(100)); const actual = e.snapshot();
  for (const key of ['pos', 'vel', 'rng', 'wallT', 'wallTarget', 'wallTau', 'heatToSample', 'heaterWork']) assert.deepEqual(actual[key], expected[key], key);
});
test('The Kelvin stat holds the setpoint exactly, every step and everywhere', () => {
  const e = chamber({ thermostatMode: 'kelvin', T: 450 });
  for (let i = 0; i < 20; i++) e.addAtom('Ar', 8 + (i % 5) * 20, 12 + ((i / 5) | 0) * 22, 0, { thermal: true });
  for (let s = 0; s < 2000; s++) { e.step(); near(e.temperature(), 450, 1e-9); }
});
test('It replaces exactly what a void wall removes', () => {
  const e = chamber({ thermostatMode: 'kelvin', T: 450, voidTemperature: true, voidVelocity: true });
  for (let i = 0; i < 20; i++) e.addAtom('Ar', 8 + (i % 5) * 20, 12 + ((i / 5) | 0) * 22, 0, { thermal: true });
  for (let s = 0; s < 4000; s++) e.step();
  near(e.temperature(), 450, 1e-9);
  assert.ok(e.voidHeat > 0, 'the wall took energy');
  assert.ok(e.kelvinWork > 0, 'and the stat put it back');
});
test('It starts a chamber from rest and can set it back to rest', () => {
  const e = chamber({ thermostatMode: 'kelvin', T: 300 });
  for (let i = 0; i < 6; i++) e.addAtom('Ar', 10 + i * 15, 40, 0, { thermal: false });
  assert.equal(e.kinetic(), 0);
  e.step();
  near(e.temperature(), 300, 1e-9);
  e.setTemperature(0);
  near(e.temperature(), 0, 1e-12);
  e.setTemperature(750);
  near(e.temperature(), 750, 1e-9);        // instant, by definition
});
test('It leaves pinned atoms fixed and scales only the total, not the distribution', () => {
  const e = chamber({ thermostatMode: 'kelvin', T: 500 });
  for (let i = 0; i < 6; i++) e.addAtom('Ar', 10 + i * 15, 40, 0, { thermal: true });
  e.addAtom('Ar', 50, 20, 0, { thermal: true }); e.pinned[6] = 1;
  const ratio = Math.hypot(e.vel[0], e.vel[1], e.vel[2]) / Math.hypot(e.vel[3], e.vel[4], e.vel[5]);
  e._kelvin();
  near(Math.hypot(e.vel[0], e.vel[1], e.vel[2]) / Math.hypot(e.vel[3], e.vel[4], e.vel[5]), ratio, 1e-12);
  for (let d = 18; d < 21; d++) assert.equal(e.vel[d], 0);
});
test('Kelvin trajectories replay bit-exactly and survive a scene round trip', () => {
  const e = chamber({ thermostatMode: 'kelvin', T: 500, voidTemperature: true });
  e.recording = true;
  for (let i = 0; i < 8; i++) e.addAtom('Ar', 10 + i * 11, 30 + (i % 3) * 14, 0, { thermal: true });
  for (let s = 0; s < 200; s++) e.step();
  const snap = e.snapshot(), trace = [];
  for (let i = 0; i < 80; i++) { e.step(); trace.push([e.pos[0], e.vel[0], e.kelvinWork]); }
  e.restore(snap);
  assert.equal(e.thermostatMode, 'kelvin');
  for (const want of trace) { e.step(); near(e.pos[0], want[0], 1e-12); near(e.vel[0], want[1], 1e-14); near(e.kelvinWork, want[2], 1e-9); }
  assert.equal(e.toJSON().thermostatMode, 'kelvin');
});
test('Kelvin startup from rest accounts for all injected energy',()=>{
  const e=chamber({thermostatMode:'kelvin',T:300});e.addAtom('Ar',50,50,0,{thermal:false});
  e._kelvin();near(e.kelvinWork,e.kinetic(),1e-10);
  e.T=0;e._kelvin();near(e.kelvinWork,0,1e-10);near(e.kinetic(),0);
});
test('Dragging one atom does not stop every other one', () => {
  // an atom held by the pointer is driven from outside; counting its motion as temperature made
  // a stat holding the total scale the whole chamber down to compensate
  for (const opts of [{ thermostatMode: 'kelvin' }, { thermostat: false, voidTemperature: true }]) {
    const e = chamber({ T: 400, ...opts });
    for (let i = 0; i < 20; i++) e.addAtom('Ar', 8 + (i % 5) * 18, 9 + ((i / 5) | 0) * 22, 0, { thermal: true });
    for (let s = 0; s < 2000; s++) e.step();
    const rest = () => { let k = 0; for (let i = 1; i < e.N; i++) { const j = 3 * i; k += e.mass[i] * (e.vel[j] ** 2 + e.vel[j + 1] ** 2 + e.vel[j + 2] ** 2); } return k; };
    const before = rest();
    e.tweezer = { i: 0, x: e.pos[0], y: e.pos[1], k: 30 };
    for (let s = 0; s < 3000; s++) {
      e.tweezer.x = 50 + 30 * Math.sin(s / 300); e.tweezer.y = 50 + 20 * Math.cos(s / 300);
      e.needForces = true; e.step();
    }
    assert.ok(rest() > 0.5 * before, `the rest of the chamber keeps moving: ${rest()} vs ${before}`);
    /* The reading is the fluid's, not the drag's: pinned exactly by the stat, and with only a
       void the chamber is free to cool, since a drag is no longer a heat source. */
    if (opts.thermostatMode === 'kelvin') near(e.temperature(), 400, 1e-9);
    else assert.ok(e.temperature() <= 400 + 1e-9, `never above its setting, got ${e.temperature()}`);
    e.tweezer = null;
  }
});
test('A driven atom is excluded from the temperature, and only while it is held', () => {
  const e = chamber({ thermostat: false, T: 300 });
  e.addAtom('Ar', 30, 30, 0, { thermal: false, v: [0.004, 0, 0] });
  e.addAtom('Ar', 50, 50, 0, { thermal: false, v: [0.004, 0, 0] });
  const both = e.temperature();
  e.tweezer = { i: 0, x: 30, y: 30, k: 30 };
  const one = e.temperature();
  assert.equal(e.dof(), 3, 'the held atom leaves the count');
  near(one, both, 1e-9);                          // identical atoms: same temperature, fewer of them
  e.vel[0] = 0.4;                                 // fling the held one
  near(e.temperature(), one, 1e-12, 'its motion never enters the reading');
  e.tweezer = null;
  assert.ok(e.temperature() > one * 10, 'and counts again the moment it is let go');
});
test('A spark survives the thermostat long enough to light something', () => {
  // a stat that holds a total would erase a spark on the step it landed, so ignition holds it off
  const build = () => {
    const e = chamber({ T: 300, thermostatMode: 'kelvin' });
    for (let m = 0; m < 3; m++) { e.addAtom('H', 10 + m * 14, 20, 0, { thermal: true }); e.addAtom('H', 10.74 + m * 14, 20, 0, { thermal: true }); e.setBondOrder(e.N - 2, e.N - 1, 1); }
    e.refresh();
    for (let s = 0; s < 2000; s++) e.step();
    return e;
  };
  const quiet = build(), lit = build();
  const K = lit.kinetic();
  lit.vel[0] -= 0.25; lit.vel[3] += 0.25;         // pull one bond apart, as the spark does
  lit.sparkHold = lit.time + 2000;
  lit.step();
  assert.ok(lit.kinetic() > K * 2, 'the stat leaves the spark alone while the window is open');
  while (lit.time < lit.sparkHold) lit.step();
  for (let s = 0; s < 3000; s++) lit.step();
  near(lit.temperature(), 300, 1e-6);             // and takes hold again once it closes
  // without the window the same kick is erased immediately
  quiet.vel[0] -= 0.25; quiet.vel[3] += 0.25;
  quiet.step();
  near(quiet.temperature(), 300, 1e-6);
});
test('The ignition window is part of the state and replays exactly', () => {
  const e = chamber({ T: 400, thermostatMode: 'kelvin' });
  e.recording = true;
  for (let i = 0; i < 6; i++) e.addAtom('Ar', 12 + i * 12, 30, 0, { thermal: true });
  for (let s = 0; s < 200; s++) e.step();
  e.vel[0] += 0.3; e.sparkHold = e.time + 500;
  const snap = e.snapshot(), trace = [];
  for (let i = 0; i < 300; i++) { e.step(); trace.push([e.pos[0], e.vel[0]]); }
  e.restore(snap);
  assert.equal(e.sparkHold, snap.sparkHold);
  for (const want of trace) { e.step(); near(e.pos[0], want[0], 1e-12); near(e.vel[0], want[1], 1e-14); }
});
test('Dragging a molecule does not stop the rest either', () => {
  // holding one atom drags its bonded partners along; that motion is the pointer's, not the
  // sample's, and counting it as heat let a stat scale everything else to a halt
  for (const opts of [{ thermostatMode: 'kelvin' }, { thermostat: false, voidTemperature: true }]) {
    const e = chamber({ T: 400, ...opts });
    e.addAtom('O', 20, 30, 0, { thermal: true });
    e.addAtom('H', 20.76, 30.59, 0, { thermal: true });
    e.addAtom('H', 19.24, 30.59, 0, { thermal: true });
    e.setBondOrder(0, 1, 1); e.setBondOrder(0, 2, 1);
    for (let i = 0; i < 12; i++) e.addAtom('Ar', 40 + (i % 4) * 8, 12 + ((i / 4) | 0) * 16, 0, { thermal: true });
    e.refresh();
    for (let s = 0; s < 2000; s++) e.step();
    const rest = () => { let k = 0; for (let i = 3; i < e.N; i++) { const j = 3 * i; k += e.mass[i] * (e.vel[j] ** 2 + e.vel[j + 1] ** 2 + e.vel[j + 2] ** 2); } return k; };
    const before = rest();
    e.tweezer = { i: 0, x: e.pos[0], y: e.pos[1], k: 30 };      // hold the oxygen
    assert.ok(e.driven(1) && e.driven(2), 'the whole molecule counts as dragged');
    assert.ok(!e.driven(3), 'and nothing else does');
    for (let s = 0; s < 3000; s++) {
      e.tweezer.x = 40 + 20 * Math.sin(s / 300); e.tweezer.y = 40 + 15 * Math.cos(s / 300);
      e.needForces = true; e.step();
    }
    assert.ok(rest() > 0.4 * before, `the chamber keeps moving: ${(rest() / before * 100).toFixed(0)}%`);
    e.tweezer = null;
  }
});
console.log(`${count} thermal-wall checks passed.`);
