/* What a chamber actually does, run the way someone at a bench would run it:
 *   node playground/tests/lab.cjs
 * Every line is ten runs of the same experiment with different thermal seeds, because whether any
 * one collision ends in a bond is a coin toss here as it is in a flask. The hard assertions are
 * only the things that must never be otherwise — a settled mixture must not light itself, two
 * radicals must find each other. The rest is printed so the sandbox's range is visible, including
 * where it stops. */
const assert = require('node:assert/strict');
const { Engine, KB } = require('../engine.js');

const products = e => {
  const c = {};
  for (const f of e.fragments().list) { const s = e.formulaOf(f); c[s] = (c[s] || 0) + 1; }
  return Object.keys(c).sort().map(k => (c[k] > 1 ? c[k] + ' ' : '') + k).join(' + ');
};
/* The spark tool: a hard outward kick on the atoms under the pointer, and a short window in which
   the thermostat stands back and lets it do its work. */
function spark(e, list, T = 25000) {
  let cx = 0, cy = 0;
  for (const i of list) { cx += e.pos[3 * i]; cy += e.pos[3 * i + 1]; }
  cx /= list.length; cy /= list.length;
  for (const i of list) {
    const k = 3 * i;
    let ox = e.pos[k] - cx, oy = e.pos[k + 1] - cy, on = Math.hypot(ox, oy);
    if (on < 1e-6) { ox = e.gauss(); oy = e.gauss(); on = Math.hypot(ox, oy) || 1; }
    const sp = Math.sqrt(3 * KB * T / (e.mass[i] * 1e4));
    e.vel[k] += ox / on * sp; e.vel[k + 1] += oy / on * sp;
  }
  e.sparkHold = e.time + 2000; e.touch();
}
const put = {
  H:   (e, x, y) => { e.addAtom('H', x, y, 0); },
  Na:  (e, x, y) => { e.addAtom('Na', x, y, 0); },
  H2:  (e, x, y) => { e.addAtom('H', x, y, 0); e.addAtom('H', x + .74, y, 0); e.setBondOrder(e.N - 2, e.N - 1, 1); },
  O2:  (e, x, y) => { e.addAtom('O', x, y, 0); e.addAtom('O', x + 1.21, y, 0); e.setBondOrder(e.N - 2, e.N - 1, 2); },
  N2:  (e, x, y) => { e.addAtom('N', x, y, 0); e.addAtom('N', x + 1.1, y, 0); e.setBondOrder(e.N - 2, e.N - 1, 3); },
  Cl2: (e, x, y) => { e.addAtom('Cl', x, y, 0); e.addAtom('Cl', x + 1.99, y, 0); e.setBondOrder(e.N - 2, e.N - 1, 1); },
  Na2: (e, x, y) => { e.addAtom('Na', x, y, 0); e.addAtom('Na', x + 3.08, y, 0); e.setBondOrder(e.N - 2, e.N - 1, 1); },
  HCl: (e, x, y) => { e.addAtom('Cl', x, y, 0); e.addAtom('H', x + 1.27, y, 0); },
  NH3: (e, x, y) => { e.addAtom('N', x, y, 0); e.addAtom('H', x + 1, y + .3, 0); e.addAtom('H', x - .5, y + .9, 0); e.addAtom('H', x - .5, y - .9, .4); },
  CH4: (e, x, y) => { const d = .63; e.addAtom('C', x, y, 0); for (const [a, b, c] of [[d, d, d], [-d, -d, d], [-d, d, -d], [d, -d, -d]]) e.addAtom('H', x + a, y + b, c); },
  C2H4:(e, x, y) => { e.addAtom('C', x, y, 0); e.addAtom('C', x + 1.34, y, 0); e.setBondOrder(e.N - 2, e.N - 1, 2);
    e.addAtom('H', x - .55, y + .94, 0); e.addAtom('H', x - .55, y - .94, 0); e.addAtom('H', x + 1.89, y + .94, 0); e.addAtom('H', x + 1.89, y - .94, 0); },
};
const RUNS = 10;
function experiment(layout, { steps = 40000, ignite = null, ...opts } = {}) {
  const tally = new Map();
  for (let s = 0; s < RUNS; s++) {
    const e = new Engine({ width: 30, height: 24, depth: 14, T: 300, seed: 70 + s * 29, thermostatMode: 'kelvin', ...opts });
    layout(e); e.touch(); e.refresh();
    if (ignite) spark(e, ignite(e));
    for (let i = 0; i < steps; i++) e.step();
    const p = products(e);
    tally.set(p, (tally.get(p) || 0) + 1);
  }
  return [...tally.entries()].sort((a, b) => b[1] - a[1]);
}
const show = (name, rows) =>
  console.log('  ' + name.padEnd(30) + rows.slice(0, 3).map(([p, n]) => n + '/' + RUNS + '  ' + p).join('   ·   '));
const howOften = (rows, want) => (rows.find(([p]) => p === want) || [null, 0])[1];

console.log('\nWhat a closed chamber does in 40 ps, ' + RUNS + ' runs each. Kelvin stat at 300 K unless stated.\n');

console.log('Mixtures that must sit there, and do:');
{
  const cases = [
    ['H2 + O2', e => { put.H2(e, 8, 8); put.H2(e, 8, 16); put.O2(e, 20, 12); }, '2 H2 + O2'],
    ['N2 + 3 H2 (no catalyst)', e => { put.N2(e, 8, 12); put.H2(e, 18, 7); put.H2(e, 18, 13); put.H2(e, 18, 19); }, '3 H2 + N2'],
    ['C2H4 + H2 (no catalyst)', e => { put.C2H4(e, 8, 12); put.H2(e, 20, 12); }, 'C2H4 + H2'],
    ['CH4 + O2', e => { put.CH4(e, 8, 12); put.O2(e, 20, 12); }, 'CH4 + O2'],
  ];
  for (const [name, layout, intact] of cases) {
    const rows = experiment(layout);
    show(name, rows);
    assert.equal(howOften(rows, intact), RUNS, name + ' lit itself at room temperature');
  }
}
console.log('\nThings that go on their own:');
{
  const two = experiment(e => { put.H(e, 10, 12); put.H(e, 18, 13); });
  show('H· + H· → H2', two);
  assert.ok(howOften(two, 'H2') >= 8, 'two hydrogen radicals should pair up');
  const salt = experiment(e => { put.Na(e, 10, 12); e.addAtom('Cl', 16, 12.4, 0); });
  show('Na + Cl → NaCl', salt);
  assert.ok(howOften(salt, 'NaCl') >= 8, 'sodium and chlorine should make salt');
  show('Na + Cl2 → NaCl + Cl', experiment(e => { put.Na(e, 10, 12); put.Cl2(e, 18, 12.4); }, { steps: 60000 }));
}
console.log('\nThings that need lighting (adiabatic, as a closed vessel is):');
{
  const cold = experiment(e => { put.H2(e, 8, 8); put.H2(e, 8, 16); put.O2(e, 20, 12); }, { thermostat: false, steps: 60000 });
  const lit = experiment(e => { put.H2(e, 8, 8); put.H2(e, 8, 16); put.O2(e, 20, 12); }, { thermostat: false, steps: 60000, ignite: () => [0, 1] });
  show('2 H2 + O2, unlit', cold);
  show('2 H2 + O2, sparked', lit);
  assert.equal(howOften(cold, '2 H2 + O2'), RUNS);
  assert.ok(howOften(lit, '2 H2 + O2') <= RUNS - 7, 'a spark should start something');
  show('CH4 + 2 O2, sparked', experiment(e => { put.CH4(e, 8, 12); put.O2(e, 20, 8); put.O2(e, 20, 17); }, { thermostat: false, steps: 60000, ignite: () => [0, 1, 2, 3, 4] }));
  const salt = experiment(e => { put.Na2(e, 9, 12); put.Cl2(e, 20, 12); }, { thermostat: false, steps: 60000, ignite: () => [2, 3] });
  show('Na2 + Cl2, sparked', salt);
  assert.ok(howOften(salt, '2 NaCl') >= 8, 'a lit sodium/chlorine mixture should make salt');
}
console.log('\nWhere this model stops:');
{
  show('Na2 + Cl2, unlit', experiment(e => { put.Na2(e, 9, 12); put.Cl2(e, 20, 12); }, { steps: 60000 }));
  console.log('    ^ real sodium and chlorine need no spark. Nothing here transfers an electron between');
  console.log('      two closed shells, so every reaction has to start at a radical. The energies are');
  console.log('      right (2 NaCl is 506 kJ/mol below Na2 + Cl2, against about 518 measured) and a');
  console.log('      spark starts it, but it will not light itself.');
  show('HCl + NH3 → NH4Cl', experiment(e => { put.HCl(e, 9, 12); put.NH3(e, 18, 13); }, { steps: 60000 }));
  console.log('    ^ the white smoke of an acid meeting a base is a proton moving between two ions.');
  console.log('      This engine shares electrons between neutral atoms; it has no proton to move.');
}
console.log('');
