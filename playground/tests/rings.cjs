const assert = require('node:assert/strict');
const { Engine } = require('../engine.js');
let count = 0;
function test(name, run) { run(); console.log('PASS', name); count++; }
const q = { thermal: false };
const box = () => { const e = new Engine({ width: 200, height: 200, depth: 200, T: 0, thermostat: false }); e.box = { x0: -100, x1: 100, y0: -100, y1: 100, z0: -100, z1: 100 }; return e; };
const atomE = s => { const e = box(); e.addAtom(s, 0, 0, 0, q); e.refresh(); return e.computeForces(); };
const EC = atomE('C'), EH = atomE('H');
function cycloalkane(n) {
  const R = 1.52 / (2 * Math.sin(Math.PI / n)), e = box(), C = [];
  for (let k = 0; k < n; k++) { const a = 2 * Math.PI * k / n; C.push([R * Math.cos(a), R * Math.sin(a)]); e.addAtom('C', C[k][0], C[k][1], 0, q); }
  for (let k = 0; k < n; k++) { const a = 2 * Math.PI * k / n, ux = Math.cos(a), uy = Math.sin(a); for (const z of [0.89, -0.89]) e.addAtom('H', C[k][0] + 0.63 * ux, C[k][1] + 0.63 * uy, z, q); }
  e.touch(); e.refresh();
  for (let k = 0; k < 4; k++) { e.minimize(4000, 0.02); e.refresh(); for (let m = 0; m < 100; m++) e.computeForces(1); }
  return { atomization: n * EC + 2 * n * EH - e.computeForces(), cc: Math.hypot(e.pos[0] - e.pos[3], e.pos[1] - e.pos[4], e.pos[2] - e.pos[5]), frags: e.fragments().list.length };
}

test('Cyclopropane carries its real ring strain, not the full cost of bending tetrahedral carbon to 60°', () => {
  const r = cycloalkane(3);
  assert.equal(r.frags, 1);
  assert.ok(Math.abs(r.atomization - 3401) < 70, `atomization ${r.atomization.toFixed(0)} kJ/mol (exp. 3401)`);
  assert.ok(Math.abs(r.cc - 1.51) < 0.05, `C–C ${r.cc.toFixed(3)} Å (exp. 1.510)`);
});

test('Larger rings, never fitted, keep their measured atomization energies', () => {
  for (const [n, exp] of [[4, 4600], [5, 5850], [6, 7055]]) {
    const r = cycloalkane(n);
    assert.ok(Math.abs(r.atomization - exp) / exp < 0.03, `C${n}H${2 * n}: ${r.atomization.toFixed(0)} kJ/mol (exp. ${exp})`);
  }
});

console.log(count + ' ring checks passed.');
