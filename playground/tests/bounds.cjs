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
    // a solid face is a stiff force starting at the face: a fast argon goes a fraction of an ångström in
    if(e.boundsMode==='solid') assert.ok(deepSolid < 0.3, `solid face yielded ${deepSolid} Å`);
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

test('Void temperature radiates self-made heat away instantly, everywhere', () => {
  /* A chamber that really does heat itself: loose hydrogen and oxygen atoms, which pair off and
     release their bond energy into the gas. A settled mixture of H2 and O2 will not do — at 300 K
     it is metastable and stays that way, which is the whole point of needing a spark. */
  const build = opts => {
    const e = new Engine({ width: 26, height: 22, depth: 14, T: 300, wallT: 300, seed: 11, thermostat: false, ...opts });
    let x = 5, y = 5;
    for (let m = 0; m < 8; m++) { e.addAtom('H', x, y, 0, { thermal: true }); x += 6; if (x > 22) { x = 5; y += 6; } }
    for (let m = 0; m < 4; m++) { e.addAtom('O', x, y, 0, { thermal: true }); x += 6; if (x > 22) { x = 5; y += 6; } }
    return e;
  };
  const closed = build({}), open = build({ voidTemperature: true });
  let hottestClosed = 0, hottestOpen = 0;
  for (let s = 0; s < 40000; s++) {
    closed.step(); open.step();
    hottestClosed = Math.max(hottestClosed, closed.temperature());
    hottestOpen = Math.max(hottestOpen, open.temperature());
  }
  assert.ok(hottestClosed > 1000, `a closed chamber heats itself, reached ${hottestClosed.toFixed(0)} K`);
  assert.ok(hottestOpen <= 300 + 1e-6, `a voided one never exceeds its setting, reached ${hottestOpen.toFixed(2)} K`);
  assert.ok(open.voidHeat > 0);
});

test('Void temperature only ever removes: a cold chamber is left alone', () => {
  const e = chamber({ voidTemperature: true, T: 900, thermostat: false });
  e.addAtom('Ar', 30, 30, 0, { thermal: false, v: [0.004, -0.002, 0.001] });
  const before = [...e.vel.slice(0, 3)];
  for (let i = 0; i < 500; i++) e.step();
  assert.equal(e.voidHeat, 0);                       // well below 900 K: nothing to radiate
  for (let d = 0; d < 3; d++) near(e.vel[d], before[d], 1e-12);
});

test('A velocity void stops an atom dead at the wall and leaves it there', () => {
  const e = chamber({ voidVelocity: true, thermostat: false });
  e.addAtom('Ar', 30, 30, 0, { thermal: false, v: [0.02, 0.004, -0.003] });
  for (let i = 0; i < 6000; i++) e.step();
  for (let d = 0; d < 3; d++) assert.equal(e.vel[d], 0, `component ${d} is exactly zero, not damped`);
  const resting = [...e.pos.slice(0, 3)];
  for (let i = 0; i < 2000; i++) e.step();
  for (let d = 0; d < 3; d++) assert.equal(e.pos[d], resting[d], 'and it does not drift afterwards');
});

test('A pressure void silences the gauge and changes nothing about the motion', () => {
  const open = chamber({ voidPressure: true, thermostat: false });
  const closed = chamber({ thermostat: false });
  for (const e of [open, closed]) for (let i = 0; i < 18; i++) e.addAtom('Ar', 6 + (i % 6) * 10, 8 + ((i / 6) | 0) * 18, 0, { thermal: true });
  const width = open.box.x1 - open.box.x0;
  for (let i = 0; i < 4000; i++) { open.step(); closed.step(); }
  // a voided chamber is trajectory-identical to a reflecting one
  for (let k = 0; k < 3 * open.N; k++) { near(open.pos[k], closed.pos[k], 1e-12); near(open.vel[k], closed.vel[k], 1e-12); }
  near(open.temperature(), closed.temperature(), 1e-12);
  assert.equal(open.voidHeat, 0, 'voiding a reading takes no energy');
  near(open.wallForce, 0, 1e-12);
  near(open.pressureBar, 0, 1e-12);
  assert.ok(closed.pressureEMA > 0, 'while the same chamber otherwise reports pressure');
  near(open.box.x1 - open.box.x0, width, 1e-12);      // and the chamber never grew
});

test('The wall reports the fluid lying against it, not a heater setting', () => {
  const e = chamber({ voidTemperature: true, T: 600, thermostat: true });
  for (let i = 0; i < 24; i++) e.addAtom('Ar', 6 + (i % 6) * 9, 7 + ((i / 6) | 0) * 9, 0, { thermal: true });
  for (let s = 0; s < 4000; s++) e.step();
  assert.ok(e.wallContact > 0, 'atoms are in the boundary layer');
  near(e.wallT, e.wallMeasured, 1e-12);               // the wall is what the fluid is
  assert.ok(e.wallMeasured > 0 && e.wallMeasured < 1e5);
});

test('Solid walls turn an atom round at all six faces, keep its energy and report the momentum it delivered', () => {
  for(let axis=0;axis<3;axis++) for(const sign of [-1,1]) {
    const e=chamber();   // no void selected: the face gives back everything it receives
    const xyz=[30,30,0], v=[0,0,0]; xyz[axis]=(axis===2?0:30)+sign*29.9; v[axis]=sign*.02;
    e.addAtom('Ar',...xyz,{thermal:false,v});
    e.computeForces();
    const K=e.kinetic();
    let J=0;                                       // impulse the face delivered, from its own force
    for(let k=0;k<120;k++){ e.step(); J+=e.wallForce*e.dt; }
    assert.equal(Math.sign(e.vel[axis]),-sign,'turned round');
    near(e.kinetic(),K,1e-4*K);                    // off the face again, with the speed it arrived with
    near(J,2*e.mass[0]*.02*1e4,0.01*2*e.mass[0]*.02*1e4);
    for(let k=0;k<5000;k++) {e.step(); assert.ok(e.pos[axis]>= (axis===2?-30:0)-0.3 && e.pos[axis]<= (axis===2?30:60)+0.3);}
    near(e.kinetic()+e.Epot,K,1e-3*K);
  }
});
test('Corner and multiple-face drift crossings are folded without losing energy', () => {
  const e=chamber(); e.addAtom('Ar',30,30,0,{thermal:false,v:[.02,.03,.04]});
  const K=e.kinetic(); e.pos.set([181,-121,151]); e._wallImpulse=0; e._reflectWalls();
  [59,1,29].forEach((x,k)=>near(e.pos[k],x)); near(e.kinetic(),K);
  assert.ok(e._wallImpulse>0);
});

test('A solid wall bounces a molecule without making or losing energy', () => {
  /* Bouncing each atom on its own mirrored one across the wall while its partners stayed put:
     a cold water molecule thrown at a wall at thermal speed came away 110 kJ/mol hotter. Carrying
     the whole molecule back fixed that alone but shoved it into its neighbours in a crowd. The face
     is a stiff force now, and a molecule hitting it is knocked into vibration the way a real one
     is — energy moves from its travel into its bonds, and the total stays put. */
  for (const speed of [0.005, 0.01, 0.02, 0.05]) {
    const e = new Engine({ width: 30, height: 24, depth: 14, T: 0, seed: 7, thermostat: false, boundsMode: 'solid' });
    e.addAtom('O', 15, 12, 0, { thermal: false });
    e.addAtom('H', 15.76, 12.6, 0, { thermal: false });
    e.addAtom('H', 14.24, 12.6, 0, { thermal: false });
    e.touch(); e.minimize(1500, 0.02); e.refresh();
    for (let i = 0; i < 3; i++) { e.vel[3 * i] = speed; e.vel[3 * i + 1] = speed * 0.6; }
    const E0 = e.Epot + e.kinetic();
    for (let i = 0; i < 20000; i++) e.step();
    const drift = e.Epot + e.kinetic() - E0;
    // up to 2 km/s — a hot gas — the wall keeps the books to a couple of kJ/mol. At 5 km/s the
    // molecule arrives with some 300 kJ/mol, about 24 000 K of travel, is thrown into violent
    // vibration, and the trajectory turns chaotic: a few percent of that is honest there.
    const allowed = speed < 0.03 ? 2 : 15;
    assert.ok(Math.abs(drift) < allowed, `at ${speed} Å/fs the wall changed the energy by ${drift.toFixed(1)} kJ/mol`);
    assert.equal(e.fragments().list.map(f => e.formulaOf(f)).join('+'), 'H2O');
  }
});
test('A crowded chamber against solid walls keeps its energy', () => {
  const e = new Engine({ width: 20, height: 14, depth: 10, T: 300, seed: 3, thermostat: false, boundsMode: 'solid' });
  for (let i = 0; i < 5; i++) for (let j = 0; j < 3; j++) {
    const x = 2.5 + i * 3.8, y = 2.5 + j * 4.2;
    e.addAtom('O', x, y, 0); e.addAtom('H', x + 0.76, y + 0.6, 0); e.addAtom('H', x - 0.76, y + 0.6, 0);
  }
  e.touch(); e.minimize(800, 0.2); e.refresh(); e.thermalize(300);
  const E0 = e.Epot + e.kinetic();
  for (let s = 0; s < 20000; s++) e.step();
  // before this was +173 kJ/mol in 20 ps, and +4122 with the per-atom mirror before that
  assert.ok(Math.abs(e.Epot + e.kinetic() - E0) < 10, `drifted ${(e.Epot + e.kinetic() - E0).toFixed(1)} kJ/mol`);
});
test('A molecule crossing a corner comes back whole', () => {
  const e = new Engine({ width: 30, height: 24, depth: 14, T: 0, seed: 7, thermostat: false, boundsMode: 'solid' });
  e.addAtom('O', 15, 12, 0, { thermal: false });
  e.addAtom('H', 15.76, 12.6, 0, { thermal: false });
  e.addAtom('H', 14.24, 12.6, 0, { thermal: false });
  e.touch(); e.minimize(1500, 0.02); e.refresh();
  const d = (a, b) => Math.hypot(e.pos[3 * a] - e.pos[3 * b], e.pos[3 * a + 1] - e.pos[3 * b + 1], e.pos[3 * a + 2] - e.pos[3 * b + 2]);
  const was = [d(0, 1), d(0, 2), d(1, 2)];
  for (let i = 0; i < 3; i++) { e.pos[3 * i] += 20; e.pos[3 * i + 1] -= 16; e.vel[3 * i] = .02; e.vel[3 * i + 1] = -.02; }
  const K = e.kinetic();                          // out through a corner, and travelling further out
  e._reflectWalls();
  near(d(0, 1), was[0], 1e-12); near(d(0, 2), was[1], 1e-12); near(d(1, 2), was[2], 1e-12);
  near(e.kinetic(), K, 1e-10);
  for (let i = 0; i < 3; i++) { assert.ok(e.pos[3 * i] >= 0 && e.pos[3 * i] <= 30); assert.ok(e.pos[3 * i + 1] >= 0 && e.pos[3 * i + 1] <= 24); }
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
  // a contact lasts a few femtoseconds rather than none, so this agrees to a fraction of a percent, not to the digit
  const ideal=(2*e.kinetic()/3)/(60**3)*16605.39; near(p/n,ideal,0.005*ideal);
});
test('Every void channel removes energy at both wall kinds, and none can add any', () => {
  for (const mode of ['solid', 'forcefield']) for (const ch of ['voidTemperature', 'voidPressure', 'voidVelocity']) {
    const e = chamber({ boundsMode: mode, [ch]: true, T: 600 });
    for (let i = 0; i < 18; i++) e.addAtom('Ar', 6 + (i % 6) * 10, 8 + ((i / 6) | 0) * 18, 0, { thermal: true });
    // a sample above its setting, so the temperature channel has excess heat to radiate
    if (ch === 'voidTemperature') e.T = 100;
    const K0 = e.kinetic();
    for (let s = 0; s < 6000; s++) e.step();
    if (ch === 'voidPressure') {                     // a voided reading, not an energy sink
      assert.equal(e.voidHeat, 0, `${mode}/${ch} takes no energy`);
      near(e.pressureBar, 0, 1e-12);
    } else {
      assert.ok(e.voidHeat > 0, `${mode}/${ch} removed energy, got ${e.voidHeat}`);
      assert.ok(e.kinetic() < K0, `${mode}/${ch} cooled the sample`);
    }
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
test('Stopping at a wall reports m*v, and a pressure void deletes that reading', () => {
  // velocity: the atom stops, and the momentum it delivered is still wall stress
  {
    const e=chamber({voidVelocity:true});e.addAtom('Ar',59.999,30,0,{thermal:false,v:[.02,.005,0]});
    const K=e.kinetic(), impulse=e.mass[0]*.02*1e4;e.step();
    near(e.wallForce,impulse/e.dt,1e-9);near(e.kinetic()+e.voidHeat,K,1e-9);near(e.pos[0],60);
    assert.ok(e.pressureBar>0,'a wall that is struck reports the blow');
    const snap=e.snapshot();e.step();const expected=e.snapshot();e.restore(snap);e.step();
    near(e.voidHeat,expected.voidHeat);near(e.pressureEMA,expected.pressureEMA);
  }
  // pressure: the same blow lands and rebounds exactly as always; only the reading is voided
  {
    const open=chamber({voidPressure:true}), closed=chamber();
    for(const e of [open,closed]) e.addAtom('Ar',59.999,30,0,{thermal:false,v:[.02,.005,0]});
    const K=open.kinetic(); open.step(); closed.step();
    assert.deepEqual([...open.vel.slice(0,3)],[...closed.vel.slice(0,3)],'the bounce is untouched');
    near(open.pos[0],closed.pos[0],1e-12);
    // mid-contact some of the energy sits in the face's spring; none of it is taken away
    near(open.kinetic()+open.Ewall,K,1e-3*K); assert.equal(open.voidHeat,0);
    near(open.wallForce,0,1e-12); near(open.pressureBar,0,1e-12);
    assert.ok(closed.pressureBar>0,'while a closed chamber reports the blow');
    const snap=open.snapshot();open.step();const expected=open.snapshot();open.restore(snap);open.step();
    near(open.voidHeat,expected.voidHeat);near(open.pressureEMA,expected.pressureEMA);
  }
});
test('Absorbing corner impacts stop at the faces and count each normal impulse',()=>{
  const e=chamber({voidVelocity:true});e.addAtom('Ar',60.1,60.2,30.1,{thermal:false,v:[.02,.03,.04]});
  const K=e.kinetic();e._wallImpulse=0;e._reflectWalls();
  assert.deepEqual([...e.pos.slice(0,3)],[60,60,30]);near(e.voidHeat,K);near(e._wallImpulse,e.mass[0]*.09*1e4,1e-8);
});
test('Rejected integration restores the absorbed-energy ledger',()=>{
  const e=chamber({voidVelocity:true});e.addAtom('Ar',60.1,30,0,{thermal:false,v:[.02,0,0]});
  e.voidHeat=7;e._save();e._wallImpulse=0;e._reflectWalls();assert.ok(e.voidHeat>7);e._load();near(e.voidHeat,7);
});

test('Radiating heat scales the whole sample at once, so bonds are never strained', () => {
  // an N2 running hot: radiating must not pull one end of the bond harder than the other
  const e = chamber({ voidTemperature: true, thermostat: false, T: 50 });
  e.addAtom('N', 30.0, 30, 0, { thermal: false, v: [0.012, 0, 0] });
  e.addAtom('N', 31.1, 30, 0, { thermal: false, v: [0.012, 0, 0] });
  e.setBondOrder(0, 1, 3);
  const bond = () => Math.abs(e.pos[0] - e.pos[3]);
  for (let s = 0; s < 4000; s++) e.step();
  assert.ok(bond() > 0.9 && bond() < 1.4, `the bond survived, length ${bond().toFixed(3)} A`);
  assert.ok(e.voidHeat > 0, 'and the excess was radiated');
  assert.ok(e.temperature() <= 50 + 1e-6, `held at its setting, ${e.temperature().toFixed(2)} K`);
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
test('Pressure control preserves pinned fragments and cannot pass walls through them',()=>{
  const e=chamber({pressureControl:true,pressureTarget:200,pressureTau:100});
  e.addAtom('N',2,10,0,{thermal:false});e.addAtom('N',3.1,10,0,{thermal:false});e.setBondOrder(0,1,3);e.pinned[0]=1;
  e.computeForces();const positions=[...e.pos.slice(0,6)];for(let n=0;n<100;n++)e._barostat();
  assert.deepEqual([...e.pos.slice(0,6)],positions);assert.ok(e.box.x0<=2);assert.ok(e.box.y0<=10);
});
test('A void wall only acts on contact, never on atoms merely near it', () => {
  for (const ch of ['voidPressure', 'voidVelocity']) {
    // gliding parallel to a face it never touches
    const glide = chamber({ [ch]: true, thermostat: false });
    glide.addAtom('Ar', 59.5, 30, 0, { thermal: false, v: [0, 0.01, 0] });
    glide.step();
    assert.deepEqual([...glide.vel.slice(0, 3)], [0, 0.01, 0], `${ch} ignores an atom gliding past`);
    // approaching, still short of the face
    const near = chamber({ [ch]: true, thermostat: false });
    near.addAtom('Ar', 59.5, 30, 0, { thermal: false, v: [0.002, 0, 0] });
    near.step();
    assert.equal(near.vel[0], 0.002, `${ch} lets an atom reach the wall before acting`);
    assert.equal(near.voidHeat, 0);
  }
});

test('A molecule resting near a wall is left exactly as a reflecting chamber leaves it', () => {
  const settle = opts => {
    const e = new Engine({ width: 40, height: 30, depth: 20, T: 0, wallT: 0, seed: 3, thermostat: false, ...opts });
    e.addAtom('O', 20, 15, 0, { thermal: false });
    e.addAtom('H', 20.76, 15.59, 0, { thermal: false });
    e.addAtom('H', 19.24, 15.59, 0, { thermal: false });
    e.refresh();
    for (let s = 0; s < 4000; s++) { e.step(); e.vel.fill(0); }   // relax into its own geometry
    for (let i = 0; i < e.N; i++) e.pos[3 * i] += 18.74;          // 0.5 A clear of the +x face
    e.prev.set(e.pos); e.vel.fill(0); e.touch(); e.refresh();
    let peak = 0;
    for (let s = 0; s < 3000; s++) { e.step(); peak = Math.max(peak, e.kinetic()); }
    return peak;
  };
  const closed = settle({});
  for (const ch of ['voidPressure', 'voidVelocity']) {
    const open = settle({ [ch]: true });
    near(open, closed, 1e-18, );
    assert.ok(open < 1e-12, `${ch} sets nothing in motion, peak ${open}`);
  }
});

console.log(`${count} boundary checks passed.`);
