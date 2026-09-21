const assert = require('node:assert/strict');
const { Engine } = require('../engine.js');
const near = (a, b, tol = 1e-9) => assert.ok(Math.abs(a - b) < tol, `${a} vs ${b}`);
let count = 0;
function test(name, run) { run(); console.log('PASS', name); count++; }
function chamber(options = {}) {
  return new Engine({ width: 60, height: 60, depth: 60, T: 300, wallT: 300, thermostat: false, ...options });
}

test('Solid bounds exert no soft force; a forcefield catches an atom before the face', () => {
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
    if(e.boundsMode==='solid') assert.equal(deepSolid,0);
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

test('Solid walls reflect all six faces, conserve energy and report momentum flux', () => {
  for(let axis=0;axis<3;axis++) for(const sign of [-1,1]) {
    const e=chamber();   // no void selected: the face gives back everything it receives
    const xyz=[30,30,0], v=[0,0,0]; xyz[axis]=(axis===2?0:30)+sign*29.99; v[axis]=sign*.02;
    e.addAtom('Ar',...xyz,{thermal:false,v});
    const K=e.kinetic(); e.step();
    near(e.kinetic(),K); assert.equal(Math.sign(e.vel[axis]),-sign);
    const area=6*60*60, force=2*e.mass[0]*.02*1e4/e.dt;
    near(e.wallForce,force); near(e.pressureBar,force/area*16605.39,1e-5);
    for(let k=0;k<5000;k++) {e.step(); assert.ok(e.pos[axis]>= (axis===2?-30:0) && e.pos[axis]<= (axis===2?30:60));}
    near(e.kinetic(),K);
  }
});
test('Corner and multiple-face drift crossings are folded without losing energy', () => {
  const e=chamber(); e.addAtom('Ar',30,30,0,{thermal:false,v:[.02,.03,.04]});
  const K=e.kinetic(); e.pos.set([181,-121,151]); e._wallImpulse=0; e._reflectWalls();
  [59,1,29].forEach((x,k)=>near(e.pos[k],x)); near(e.kinetic(),K);
  assert.ok(e._wallImpulse>0);
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
  assert.equal(scene.voidVelocity, false);
});

test('A forcefield lets ordinary thermal atoms reach the void beyond the face', () => {
  const e = chamber({ boundsMode: 'forcefield', T: 300 });
  for (let i = 0; i < 24; i++) e.addAtom('Ar', 6 + (i % 6) * 9, 6 + ((i / 6) | 0) * 9, 0, { thermal: true });
  let out = 0;
  for (let s = 0; s < 8000; s++) { e.step(); for (let i = 0; i < e.N; i++) out = Math.max(out, e.pos[3 * i] - e.box.x1, e.box.x0 - e.pos[3 * i]); }
  assert.ok(out > 0.05, `atoms lean past the face at 300 K, deepest ${out}`);
});
test('Solid collision pressure and velocities replay exactly',()=>{
  const e=chamber();e.recording=true;
  e.addAtom('Ar',59.8,59.7,29.6,{thermal:false,v:[.02,.03,.04]});
  for(let i=0;i<7;i++)e.step();const snap=e.snapshot(), trace=[];
  for(let i=0;i<100;i++){e.step();trace.push([...e.pos.slice(0,3),...e.vel.slice(0,3),e.pressureEMA]);}
  e.restore(snap);
  for(const expected of trace){e.step();assert.deepEqual([...e.pos.slice(0,3),...e.vel.slice(0,3),e.pressureEMA],expected);}
});
test('Time-averaged collision pressure agrees with ideal-gas momentum balance',()=>{
  const e=chamber();e.addAtom('Ar',30,30,0,{thermal:false,v:[.02,.02,.02]});
  let p=0;const n=6000;for(let i=0;i<n;i++){e.step();p+=e.pressureBar;}
  near(p/n,(2*e.kinetic()/3)/(60**3)*16605.39,1e-6);
});
test('Every void channel removes energy at both wall kinds, and none can add any', () => {
  for (const mode of ['solid', 'forcefield']) for (const ch of ['voidTemperature', 'voidPressure', 'voidVelocity']) {
    const e = chamber({ boundsMode: mode, [ch]: true, T: 600 });
    for (let i = 0; i < 18; i++) e.addAtom('Ar', 6 + (i % 6) * 10, 8 + ((i / 6) | 0) * 18, 0, { thermal: true });
    const K0 = e.kinetic();
    for (let s = 0; s < 6000; s++) e.step();
    assert.ok(e.voidHeat > 0, `${mode}/${ch} removed energy, got ${e.voidHeat}`);
    assert.ok(e.kinetic() < K0, `${mode}/${ch} cooled the sample`);
  }
});
test('A reflecting chamber keeps its energy; every channel is opt-in', () => {
  for (const mode of ['solid', 'forcefield']) {
    const e = chamber({ boundsMode: mode, T: 600 });
    for (let i = 0; i < 12; i++) e.addAtom('Ar', 8 + (i % 4) * 14, 10 + ((i / 4) | 0) * 18, 0, { thermal: true });
    e.computeForces();
    const E0 = e.kinetic() + e.Epot;
    for (let s = 0; s < 6000; s++) e.step();
    assert.equal(e.voidHeat, 0);
    assert.ok(Math.abs(e.kinetic() + e.Epot - E0) < 0.05 * E0, `${mode} conserves energy: ${e.kinetic() + e.Epot} vs ${E0}`);
  }
});
test('A pressure void reports no wall impulse while still containing the sample', () => {
  const e = chamber({ voidPressure: true, T: 600 });
  for (let i = 0; i < 12; i++) e.addAtom('Ar', 8 + (i % 4) * 14, 10 + ((i / 4) | 0) * 18, 0, { thermal: true });
  for (let s = 0; s < 4000; s++) e.step();
  near(e.wallForce, 0, 1e-12);
  near(e.pressureBar, 0, 1e-12);
  for (let i = 0; i < e.N; i++) for (let d = 0; d < 3; d++) {
    const lo = [e.box.x0, e.box.y0, e.box.z0][d], hi = [e.box.x1, e.box.y1, e.box.z1][d];
    assert.ok(e.pos[3 * i + d] >= lo - 1e-9 && e.pos[3 * i + d] <= hi + 1e-9, 'atoms stay inside');
  }
});
test('A temperature void takes the agitation and leaves the molecule whole', () => {
  // one N2 drifting into a wall: the bond must survive, the shaking must not
  const e = chamber({ voidTemperature: true, thermostat: false });
  e.addAtom('N', 57.0, 30, 0, { thermal: false, v: [0.012, 0, 0] });
  e.addAtom('N', 58.1, 30, 0, { thermal: false, v: [0.012, 0, 0] });
  e.setBondOrder(0, 1, 3);
  const bond = () => Math.abs(e.pos[0] - e.pos[3]);
  for (let s = 0; s < 4000; s++) e.step();
  assert.ok(bond() > 0.9 && bond() < 1.4, `the bond survived, length ${bond().toFixed(3)} A`);
  assert.ok(e.voidHeat > 0, 'and the wall took energy');
});
test('A velocity void stops what reaches the wall', () => {
  const e = chamber({ voidVelocity: true, thermostat: false });
  e.addAtom('Ar', 30, 30, 0, { thermal: false, v: [0.02, 0, 0] });
  for (let s = 0; s < 4000; s++) e.step();
  assert.ok(e.kinetic() < 0.02 * 0.5 * 39.948 * 0.02 * 0.02 * 1e4, `atom came to rest, K = ${e.kinetic()}`);
});
test('Void settings replay exactly and survive a scene round trip', () => {
  const e = chamber({ boundsMode: 'forcefield', voidTemperature: true, voidVelocity: true, T: 500 });
  e.recording = true;
  for (let i = 0; i < 8; i++) e.addAtom('Ar', 10 + i * 6, 20 + (i % 3) * 12, 0, { thermal: true });
  for (let s = 0; s < 300; s++) e.step();
  const snap = e.snapshot(), trace = [];
  for (let i = 0; i < 80; i++) { e.step(); trace.push([e.pos[0], e.vel[0], e.voidHeat]); }
  e.restore(snap);
  assert.equal(e.voidVelocity, true);
  for (const want of trace) { e.step(); near(e.pos[0], want[0], 1e-12); near(e.vel[0], want[1], 1e-14); near(e.voidHeat, want[2], 1e-10); }
  const json = e.toJSON();
  assert.equal(json.dampingVersion, 3);
  assert.equal(json.voidVelocity, true);
  assert.equal(json.voidPressure, false);
});
test('Pressure control moves the chamber toward the target and leaves bonds alone', () => {
  const e = chamber({ thermostat: false, pressureControl: true, pressureTarget: 4, pressureTau: 400 });
  for (let i = 0; i < 16; i++) e.addAtom('Ar', 8 + (i % 4) * 14, 10 + ((i / 4) | 0) * 15, 0, { thermal: true });
  const w0 = e.box.x1 - e.box.x0;
  for (let s = 0; s < 4000; s++) e.step();
  const w1 = e.box.x1 - e.box.x0;
  assert.ok(Math.abs(w1 - w0) > 1e-6, `the chamber moved: ${w0} -> ${w1}`);
  near((e.box.z1 - e.box.z0), 60, 1e-9);           // the slab depth is the user's, untouched
  const held = chamber({ thermostat: false });
  for (let s = 0; s < 500; s++) held.step();
  near(held.box.x1 - held.box.x0, 60, 1e-12);      // and passive chambers never move
});
test('Pressure control keeps a molecule rigid while the walls move', () => {
  const e = chamber({ thermostat: false, pressureControl: true, pressureTarget: 200, pressureTau: 200 });
  e.addAtom('N', 29.4, 30, 0, { thermal: false });
  e.addAtom('N', 30.5, 30, 0, { thermal: false });
  e.setBondOrder(0, 1, 3);
  const before = Math.abs(e.pos[0] - e.pos[3]);
  for (let s = 0; s < 2000; s++) e.step();
  near(Math.abs(e.pos[0] - e.pos[3]), before, 0.05);
});
console.log(`${count} boundary checks passed.`);
