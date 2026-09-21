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
    const e=chamber({voidPressure:true});
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
test('Bar damping removes only outward normal energy at all six soft faces',()=>{
  for(let axis=0;axis<3;axis++) for(const sign of [-1,1]){
    const e=chamber({boundsMode:'forcefield',voidPressure:true,voidPressureTau:80});
    const xyz=[30,30,0],v=[.01,.02,.03];xyz[axis]=(axis===2?0:30)+sign*31;v[axis]=sign*.02;
    e.addAtom('Ar',...xyz,{thermal:false,v});const K=e.kinetic();
    e._voidBath(40);
    for(let d=0;d<3;d++) near(e.vel[d],d===axis?v[d]*Math.exp(-.25):v[d]);
    near(e.kinetic()+e.voidHeat,K);
    e.vel[axis]=-sign*.02;const before=[...e.vel.slice(0,3)];e._voidBath(40);
    assert.deepEqual([...e.vel.slice(0,3)],before,'inward and tangential velocities are untouched');
  }
});
test('Combined damping is monotone, leaves pinned atoms fixed and respects response/depth',()=>{
  const e=chamber({boundsMode:'forcefield',voidTemperature:true,voidPressure:true,voidTau:40,voidPressureTau:80,voidDepth:2});
  e.addAtom('Ar',61,61,0,{thermal:false,v:[.02,.03,.01]});
  e.addAtom('Ar',-1,40,0,{thermal:false,v:[-.02,0,0]});e.pinned[1]=1;
  const initial=e.kinetic();e._voidBath(20);
  near(e.vel[0],.02*Math.exp(-.375));near(e.vel[1],.03*Math.exp(-.375));near(e.vel[2],.01*Math.exp(-.25));
  near(e.vel[3],-.02);near(e.kinetic()+e.voidHeat,initial);
  const s=e.snapshot();e.voidPressureTau=1;e.restore(s);assert.equal(e.voidPressureTau,80);
  const json=e.toJSON();assert.equal(json.dampingVersion,2);assert.equal(json.voidPressureTau,80);assert.equal(json.voidTau,40);assert.equal(json.voidDepth,2);near(json.voidHeat,e.voidHeat);
});
test('Damping cannot hide soft-wall reaction force or affect solid-wall motion',()=>{
  for(const mode of ['solid','forcefield']){
    const e=chamber({boundsMode:mode,voidTemperature:true,voidPressure:true});e.addAtom('Ar',61,30,0,{thermal:false,v:[.02,0,0]});
    e.computeForces();const force=e.wallForce,K=e.kinetic();e._voidBath(20);e.computeForces();near(e.wallForce,force);
    if(mode==='solid'){near(e.kinetic(),K);near(e.voidHeat,0);}else assert.ok(force>0);
  }
});
console.log(`${count} boundary checks passed.`);
