const assert = require('node:assert/strict');
const { Engine } = require('../engine.js');
let count = 0;
function test(name, run) { run(); console.log('PASS', name); count++; }
const q = { thermal: false };

function ethene(phi) {
  const e = new Engine({ width: 200, height: 200, depth: 200, T: 0, thermostat: false });
  e.box = { x0: -100, x1: 100, y0: -100, y1: 100, z0: -100, z1: 100 };
  const c = Math.cos(phi), s = Math.sin(phi);
  for (const a of [['C', -.67, 0, 0], ['C', .67, 0, 0], ['H', -1.24, .93, 0], ['H', -1.24, -.93, 0], ['H', 1.24, .93 * c, .93 * s], ['H', 1.24, -.93 * c, -.93 * s]]) e.addAtom(a[0], a[1], a[2], a[3], q);
  e.touch(); e.refresh(); e.setBondOrder(0, 1, 2);
  for (let k = 0; k < 400; k++) e.computeForces(1);
  return e;
}

function dihedral(e, a, b, c, d) {
  const p = i => [e.pos[3 * i], e.pos[3 * i + 1], e.pos[3 * i + 2]];
  const sub = (x, y) => x.map((v, k) => v - y[k]), cr = (x, y) => [x[1] * y[2] - x[2] * y[1], x[2] * y[0] - x[0] * y[2], x[0] * y[1] - x[1] * y[0]], dot = (x, y) => x[0] * y[0] + x[1] * y[1] + x[2] * y[2];
  const b1 = sub(p(b), p(a)), b2 = sub(p(c), p(b)), b3 = sub(p(d), p(c)), n1 = cr(b1, b2), n2 = cr(b2, b3);
  return Math.atan2(Math.hypot(...b2) * dot(b1, n2), dot(n1, n2)) * 180 / Math.PI;
}

function numericError(e) {
  e.computeForces();
  const F = Float64Array.from(e.frc.subarray(0, 3 * e.N));
  let worst = 0;
  for (let i = 0; i < 3 * e.N; i++) {
    const x = e.pos[i], h = 1e-5;
    e.pos[i] = x + h; const Ep = e.computeForces();
    e.pos[i] = x - h; const Em = e.computeForces();
    e.pos[i] = x;
    worst = Math.max(worst, Math.abs(-(Ep - Em) / (2 * h) - F[i]));
  }
  return worst;
}

test('Twisting a C=C bond costs the π-bond energy (ethene, exp. ≈ 272 kJ/mol at 90°)', () => {
  const E0 = ethene(0).Epot, E30 = ethene(Math.PI / 6).Epot, E90 = ethene(Math.PI / 2).Epot;
  assert.ok(E90 - E0 > 220 && E90 - E0 < 290, `90° twist ${(E90 - E0).toFixed(1)} kJ/mol`);
  assert.ok(Math.abs((E30 - E0) / (E90 - E0) - 0.25) < 0.03, `30° twist is not sin²-shaped: ${(E30 - E0).toFixed(1)}`);
});

test('Torsion forces are the exact gradient of the energy, including a radical mid-attack', () => {
  const twisted = ethene(0.6);
  twisted.addAtom('C', 2.9, 0.3, 0.5, q); twisted.addAtom('H', 3.2, 1.2, 0.9, q); twisted.touch(); twisted.refresh();
  assert.ok(numericError(twisted) < 1e-4, 'twisted ethene with a neighbour');
  const e = ethene(0.2);
  e.addAtom('Cl', -0.8, 0.1, 2.3, q); e.touch(); e.refresh();
  for (let k = 0; k < 300; k++) e.computeForces(1);
  assert.ok(numericError(e) < 1e-4, 'Cl approaching the π bond');
});

test('cis-2-butene stays cis at 600 K: a double bond does not rotate on a picosecond timescale', () => {
  const e = new Engine({ width: 30, height: 30, depth: 30, T: 600, seed: 5, thermostatMode: 'kelvin' });
  const A = [['C', -.67, 0, 0], ['C', .67, 0, 0], ['C', -1.45, 1.25, 0], ['C', 1.45, 1.25, 0], ['H', -1.2, -.95, 0], ['H', 1.2, -.95, 0],
    ['H', -2.5, 1.0, 0], ['H', -1.25, 1.85, .89], ['H', -1.25, 1.85, -.89], ['H', 2.5, 1.0, 0], ['H', 1.25, 1.85, .89], ['H', 1.25, 1.85, -.89]];
  for (const a of A) e.addAtom(a[0], a[1] + 15, a[2] + 15, a[3], q);
  e.touch(); e.refresh(); e.setBondOrder(0, 1, 2); e.minimize(500, 0.5); e.refresh(); e.thermalize(600);
  let worst = 0;
  for (let i = 0; i < 10000; i++) { e.step(); if (i % 10 === 0) worst = Math.max(worst, Math.abs(dihedral(e, 2, 0, 1, 3))); }
  assert.ok(worst < 70, `C–C=C–C dihedral reached ${worst.toFixed(0)}°`);
});

test('Single bonds still rotate freely: ethane has no π torsion', () => {
  const e = new Engine({ width: 200, height: 200, depth: 200, T: 0, thermostat: false });
  e.box = { x0: -100, x1: 100, y0: -100, y1: 100, z0: -100, z1: 100 };
  for (const a of [['C', 0, 0, 0], ['C', 1.54, 0, 0], ['H', -.36, 1.03, 0], ['H', -.36, -.51, .89], ['H', -.36, -.51, -.89], ['H', 1.9, -1.03, 0], ['H', 1.9, .51, .89], ['H', 1.9, .51, -.89]]) e.addAtom(a[0], a[1], a[2], a[3], q);
  e.touch(); e.refresh();
  assert.equal(e._torsions(), 0);
});

console.log(count + ' torsion checks passed.');
