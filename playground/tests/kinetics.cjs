const assert = require('node:assert/strict');
const { Engine } = require('../engine.js');
const K = require('../kinetics.js');
const { geo } = require('./reactions.cjs');
let count = 0;
function test(name, run) { run(); console.log('PASS', name); count++; }
function scene(parts, T = 298) {
  const e = new Engine({ width: 30, height: 30, depth: 30, T, thermostat: false });
  for (const [name, off, orders] of parts) {
    const b = e.N;
    for (const [s, x, y, z] of geo(name)) e.addAtom(s, 15 + off[0] + x, 15 + off[1] + y, off[2] + (z || 0), { thermal: false });
    e.touch(); e.refresh();
    for (const [a, c, o] of orders || []) e.setBondOrder(b + a, b + c, o);
  }
  e.touch(); e.refresh();
  return e;
}
const near = (x, want, tol, what) => assert.ok(Math.abs(x - want) <= tol, `${what}: ${x.toFixed(1)} kJ/mol, expected ${want} ± ${tol}`);

test('H + H₂ goes over a barrier close to the real one, 40 kJ/mol', () => {
  const e = scene([['H2', [0, 0, 0]], ['H', [geo('H2')[1][1] + 2.6, 0, 0]]], 0);
  const r = K.scanPair(e, 2, 1);
  assert.ok(r.ok);
  near(r.Ea, 40, 6, 'H + H2 barrier');
  assert.deepEqual(r.products, ['H2', 'H·']);
});

test('H + CH₄ abstracts a hydrogen rather than sticking as CH₅', () => {
  const e = scene([['CH4', [0, 0, 0]], ['H', [0.63 * 2.6, 0.63 * 2.6, 0.63 * 2.6]]], 0);
  const r = K.scanPair(e, 5, 1);
  assert.ok(r.ok); assert.deepEqual(r.products, ['CH3·', 'H2']);
  near(r.Ea, 47, 8, 'H + CH4 barrier');
});

test('A chlorine atom sitting in the plane of ethene is turned to attack the π bond from above', () => {
  const e = scene([['C2H4', [0, 0, 0], [[0, 1, 2]]]]);
  e.addAtom('Cl', 21, 15, 0, { thermal: false }); e.touch(); e.refresh();
  const r = K.scanPair(e, 6, 0);
  assert.ok(r.ok, r.reason); assert.deepEqual(r.products, ['C2H4Cl·']);
  assert.ok(r.Ea < 15, 'entrance barrier ' + r.Ea.toFixed(1));
  assert.ok(r.dE < -80, 'addition releases ' + (-r.dE).toFixed(0));
});

test('Breaking bonds costs their dissociation energies', () => {
  near(K.scanPair(scene([['CH4', [0, 0, 0]]], 0), 0, 1).Ea, 439, 6, 'CH3–H');
  const eth = scene([['C2H4', [0, 0, 0], [[0, 1, 2]]]], 0), r = K.scanPair(eth, 0, 1);
  near(r.Ea, 728, 20, 'H2C=CH2'); assert.deepEqual(r.products, ['CH2·', 'CH2·']);
});

test('Two hydrogens of one molecule leave together as H₂, both C–H bonds breaking as H–H forms', () => {
  const e = scene([['C2H6', [0, 0, 0]]], 0);
  const r = K.scanPair(e, 2, 6);
  assert.ok(r.ok, r.reason); assert.equal(r.channel.type, 'swap');
  assert.deepEqual(r.products, ['C2H4', 'H2']);
  near(r.Ea, 400, 30, 'ethane 1,2-elimination');
});

test('Asking about two atoms that cannot bond says so instead of inventing a time', () => {
  const e = new Engine({ width: 30, height: 30, depth: 30, T: 298, thermostat: false });
  e.addAtom('Ar', 10, 15, 0, { thermal: false }); e.addAtom('Ar', 14, 15, 0, { thermal: false }); e.touch(); e.refresh();
  const r = K.scanPair(e, 0, 1);
  assert.equal(r.ok, false); assert.match(r.reason, /do not form a covalent bond/);
});

test('A forecast is calculated on a copy: the live chamber does not move', () => {
  const e = scene([['C2H4', [0, 0, 0], [[0, 1, 2]]], ['CH3', [0, 0, 3.4]]]);
  const pos = e.pos.slice(0, 3 * e.N), pN = e.pN.slice(0, e.nPairs), E = e.computeForces();
  K.scanPair(e, 6, 0, { T: 298 });
  assert.deepEqual(Array.from(e.pos.subarray(0, 3 * e.N)), Array.from(pos));
  assert.deepEqual(Array.from(e.pN.subarray(0, e.nPairs)), Array.from(pN));
  assert.ok(Math.abs(e.computeForces() - E) < 1e-9);
});

test('Times follow transition-state theory: slower when colder, a barrierless pair limited by meeting', () => {
  const e = scene([['CH3', [0, 0, 0]], ['CH3', [3.6, 0, 0]]]);
  const r = K.scanPair(e, 0, 4);
  assert.ok(r.ok); assert.ok(r.Ea < 1);
  const f = K.forecast(e, r, { T: 298 });
  assert.ok(f.ok && f.bimolecular);
  assert.ok(f.halfLife > 1e-12 && f.halfLife < 1e-6, 'CH3 + CH3 in this chamber: ' + K.humanTime(f.halfLife));
  const c = scene([['CH4', [0, 0, 0]]]), b = K.scanPair(c, 0, 1);
  const t = [500, 1000, 2000].map(T => K.forecast(c, b, { T }).halfLife);
  assert.ok(t[0] > t[1] && t[1] > t[2]);
  assert.ok(K.forecast(c, b, { T: 298 }).halfLife > 3.15e16, 'methane does not fall apart at room temperature');
  assert.ok(f.fastest < f.halfLife && f.halfLife < f.slowest);
});

test('A saved scene gives the same forecast in a worker as on the page', () => {
  const e = scene([['C2H4', [0, 0, 0], [[0, 1, 2]]]]);
  e.addAtom('Cl', 15.6, 15, 3.4, { thermal: false }); e.touch(); e.refresh();
  const direct = K.scanPair(e, 6, 0, { T: e.T });
  const viaScene = K.study(JSON.parse(JSON.stringify(e.toJSON())), 6, 0);
  assert.ok(viaScene.ok); assert.ok(Math.abs(viaScene.Ea - direct.Ea) < 0.5, `${viaScene.Ea} vs ${direct.Ea}`);
  assert.deepEqual(viaScene.products, direct.products);
});

test('Durations read the way a person says them', () => {
  assert.equal(K.humanTime(2.5e-11), '25.0 ps');
  assert.equal(K.humanTime(90), '1.50 min');
  assert.equal(K.humanTime(7200), '2.00 h');
  assert.equal(K.humanTime(Infinity), 'never');
  assert.match(K.humanTime(1e30), /universe/);
});

function playEvent(e, i, j, steps) {
  const r = K.study(JSON.parse(JSON.stringify(e.toJSON())), i, j);
  assert.ok(r.ok, r.reason);
  K.applyEvent(e, r.event, r.event.atoms);
  for (let s = 0; s < steps; s++) e.step();
  return { r, got: e.fragments().list.map(f => e.formulaOf(f) + (e.isRadical(f) ? '·' : '')).sort() };
}

test('Skipping ahead puts the pair just past the barrier, and the live simulation finishes the reaction itself', () => {
  let made = 0;
  for (let seed = 1; seed <= 4; seed++) {
    const e = new Engine({ width: 30, height: 30, depth: 30, T: 298, seed, thermostatMode: 'kelvin' });
    for (const [s, x, y, z] of geo('C2H4')) e.addAtom(s, 12 + x, 15 + y, z || 0, { thermal: false });
    e.touch(); e.refresh(); e.setBondOrder(0, 1, 2);
    e.addAtom('Cl', 20, 15, 0, { thermal: false }); e.touch(); e.refresh();
    const { r, got } = playEvent(e, 6, 0, 3000);
    assert.deepEqual(r.products, ['C2H4Cl·']);
    if (got.join('+') === 'C2H4Cl·') made++;
  }
  assert.ok(made >= 3, made + ' of 4 additions completed');
});

test('An abstraction played from its moment ends as the forecast said: CH₃· and H₂', () => {
  let made = 0;
  for (let seed = 1; seed <= 4; seed++) {
    const e = new Engine({ width: 30, height: 30, depth: 30, T: 600, seed, thermostatMode: 'kelvin' });
    for (const [s, x, y, z] of geo('CH4')) e.addAtom(s, 15 + x, 15 + y, z || 0, { thermal: false });
    e.addAtom('H', 19, 15, 0, { thermal: false }); e.touch(); e.refresh();
    const { got } = playEvent(e, 5, 1, 2000);
    if (got.join('+') === 'CH3·+H2') made++;
  }
  assert.ok(made >= 3, made + ' of 4 abstractions completed');
});

test('Skipped time is part of the state: undo, rewind and saved scenes keep it', () => {
  const e = new Engine({ width: 20, height: 20, depth: 20, T: 298, thermostat: false });
  e.addAtom('H', 5, 5, 0, { thermal: false }); e.touch(); e.refresh();
  const before = e.snapshot();
  e.skipped = 3600;
  assert.equal(e.toJSON().skipped, 3600);
  e.restore(before); assert.equal(e.skipped, 0);
  e.skipped = 42; const s2 = e.snapshot(); e.skipped = 0; e.restore(s2); assert.equal(e.skipped, 42);
});

test('What happens next in Cl· + ethene + methane: the addition, by far, and never a reaction that changes nothing', () => {
  const e = new Engine({ width: 30, height: 30, depth: 12, T: 298, thermostat: false });
  for (const [s, x, y, z] of geo('C2H4')) e.addAtom(s, 12 + x, 15 + y, z || 0, { thermal: false });
  e.touch(); e.refresh(); e.setBondOrder(0, 1, 2);
  e.addAtom('Cl', 18, 15, 0, { thermal: false });
  for (const [s, x, y, z] of geo('CH4')) e.addAtom(s, 22 + x, 22 + y, z || 0, { thermal: false });
  e.touch(); e.refresh();
  const sv = K.survey(JSON.parse(JSON.stringify(e.toJSON())));
  assert.equal(sv.events[0].label, 'Cl· + C2H4 → C2H4Cl·');
  assert.ok(sv.events[0].share > 0.8, 'share ' + sv.events[0].share);
  for (const x of sv.events) assert.notEqual(x.reactants.map(f => f.replace('·', '')).sort().join(), x.products.map(f => f.replace('·', '')).sort().join());
  const abstraction = sv.events.find(x => x.label === 'Cl· + CH4 → CH3· + HCl');
  assert.ok(abstraction && abstraction.share < 0.05);
  const nx = K.pickNext(sv, 0.5, 0);
  assert.equal(nx.pick, sv.events[0]);
  assert.ok(Math.abs(nx.wait - Math.LN2 / sv.total) / nx.wait < 1e-9);
  const idx = sv.events[0].event.atoms;
  assert.ok(K.stillThere(e, sv.events[0].event, idx, sv.events[0].reactants));
});

test('A chamber of stable molecules at room temperature has nothing to skip to', () => {
  const e = new Engine({ width: 30, height: 30, depth: 12, T: 298, thermostat: false });
  for (const [dx, dy] of [[8, 8], [20, 20]]) for (const [s, x, y, z] of geo('CH4')) e.addAtom(s, dx + x, dy + y, z || 0, { thermal: false });
  e.touch(); e.refresh();
  const sv = K.survey(JSON.parse(JSON.stringify(e.toJSON())));
  assert.ok(Math.LN2 / sv.total > 3.15e16, 'methane at 298 K reacted in ' + K.humanTime(Math.LN2 / sv.total));
});

console.log(count + ' kinetics checks passed.');
