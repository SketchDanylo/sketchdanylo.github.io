const assert = require('node:assert/strict');
const { Engine } = require('../engine.js');
const near = (a, b, tol = 1e-9) => assert.ok(Math.abs(a - b) < tol, `${a} vs ${b}`);
let count = 0;
function test(name, run) { run(); console.log('PASS', name); count++; }
function chamber(options = {}) {
  return new Engine({ width: 60, height: 60, depth: 60, T: 300, wallT: 300, thermostat: false, ...options });
}

test('Solid bounds push back only past the face; a forcefield catches an atom before it', () => {
  const solid = chamber(), field = chamber({ boundsMode: 'forcefield' });
  // just inside the face: solid is silent, the forcefield already pushes inward.
  for (const e of [solid, field]) e.addAtom('Ar', 59.6, 30, 0, { thermal: false });
  solid.computeForces(); field.computeForces();
  near(solid.Ewall, 0);
  assert.ok(field.frc[0] < -0.5, `forcefield pushes inward, got ${field.frc[0]}`);
  assert.ok(field.Ewall > 0);
});

test('Both bounds contain a fast atom, the forcefield over a longer distance', () => {
  const solid = chamber(), field = chamber({ boundsMode: 'forcefield' });
  let deepSolid = 0, deepField = 0;
  for (const [e, record] of [[solid, v => deepSolid = Math.max(deepSolid, v)], [field, v => deepField = Math.max(deepField, v)]]) {
    e.addAtom('Ar', 55, 30, 0, { thermal: false, v: [0.02, 0, 0] }); // ~ 2 km/s outward
    for (let i = 0; i < 4000; i++) { e.step(); record(e.pos[0] - e.box.x1); }
    assert.ok(e.pos[0] < e.box.x1, 'atom ends inside the chamber');
  }
  assert.ok(deepSolid < deepField, `forcefield yields further: ${deepSolid} vs ${deepField}`);
  assert.ok(deepField < 8, 'the forcefield still turns the atom around');
});

test('A contained atom keeps its speed when no void is selected', () => {
  const e = chamber({ boundsMode: 'forcefield' });
  e.addAtom('Ar', 55, 30, 0, { thermal: false, v: [0.02, 0, 0] });
  e.computeForces();
  const E0 = e.kinetic() + e.Epot;
  for (let i = 0; i < 4000; i++) e.step();
  const E1 = e.kinetic() + e.Epot;
  assert.ok(Math.abs(E1 - E0) < 0.02 * E0, `energy conserved: ${E1} vs ${E0}`);
  assert.equal(e.voidHeat, 0);
});

test('Void temperature drains the energy an atom carries past the face', () => {
  const e = chamber({ boundsMode: 'forcefield', voidTemperature: true });
  e.addAtom('Ar', 55, 30, 0, { thermal: false, v: [0.02, 0, 0] });
  const K0 = e.kinetic();
  for (let i = 0; i < 4000; i++) e.step();
  assert.ok(e.kinetic() < 0.5 * K0, `energy left with the void: ${e.kinetic()} vs ${K0}`);
  assert.ok(e.voidHeat > 0);
  assert.ok(e.pos[0] < e.box.x1, 'and the atom is still pushed back inside');
});

test('The void drain never acts inside the chamber', () => {
  const e = chamber({ voidTemperature: true });
  e.addAtom('Ar', 30, 30, 0, { thermal: false, v: [0.004, -0.002, 0.001] });
  const before = [...e.vel.slice(0, 3)];
  for (let i = 0; i < 500; i++) e.step();
  assert.equal(e.voidHeat, 0);
  for (let d = 0; d < 3; d++) near(e.vel[d], before[d], 1e-12);
});

test('Void pressure absorbs the impulse the gauge would report', () => {
  const open = chamber({ voidPressure: true }), closed = chamber();
  for (const e of [open, closed]) {
    e.addAtom('Ar', 59.9, 30, 0, { thermal: false, v: [0.01, 0, 0] });
    for (let i = 0; i < 200; i++) e.step();
  }
  assert.ok(closed.wallForce > 0, 'a closed chamber registers the collision');
  assert.equal(open.wallForce, 0);
  assert.equal(open.pressureBar, 0);
  assert.ok(open.voidForce > 0, 'the absorbed impulse is still tallied');
  near(open.pos[0], closed.pos[0], 1e-9); // the atom is pushed back identically
});

test('Boundary settings survive a checkpoint replay', () => {
  const e = chamber({ boundsMode: 'forcefield', voidTemperature: true, voidPressure: true });
  e.recording = true;
  e.addAtom('Ar', 55, 30, 0, { thermal: false, v: [0.015, 0.004, 0] });
  for (let i = 0; i < 300; i++) e.step();
  const s = e.snapshot();
  assert.equal(s.boundsMode, 'forcefield');
  const trace = [];
  for (let i = 0; i < 120; i++) { e.step(); trace.push(e.pos[0], e.vel[0], e.voidHeat); }
  e.restore(s);
  assert.equal(e.boundsMode, 'forcefield');
  assert.equal(e.voidTemperature, true);
  for (let i = 0; i < 120; i++) { e.step(); near(e.pos[0], trace[3 * i], 1e-12); near(e.vel[0], trace[3 * i + 1], 1e-14); near(e.voidHeat, trace[3 * i + 2], 1e-12); }
});

test('A serialized scene carries the boundary configuration', () => {
  const e = chamber({ boundsMode: 'forcefield', voidTemperature: true });
  e.addAtom('Ar', 30, 30, 0, { thermal: false });
  const scene = e.toJSON();
  assert.equal(scene.boundsMode, 'forcefield');
  assert.equal(scene.voidTemperature, true);
  assert.equal(scene.voidPressure, false);
});

test('A forcefield lets ordinary thermal atoms reach the void beyond the face', () => {
  const e = chamber({ boundsMode: 'forcefield', T: 300 });
  for (let i = 0; i < 24; i++) e.addAtom('Ar', 6 + (i % 6) * 9, 6 + ((i / 6) | 0) * 9, 0, { thermal: true });
  let out = 0;
  for (let s = 0; s < 8000; s++) { e.step(); for (let i = 0; i < e.N; i++) out = Math.max(out, e.pos[3 * i] - e.box.x1, e.box.x0 - e.pos[3 * i]); }
  assert.ok(out > 0.05, `atoms lean past the face at 300 K, deepest ${out}`);
});
console.log(`${count} boundary checks passed.`);
