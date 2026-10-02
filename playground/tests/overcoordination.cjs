/* Diagnostic, not a validation pass: how far the model lets an atom take more partners than its
 * valence. Real CH5 and five-coordinate carbon do not exist; both scans should rise steeply.
 *   node playground/tests/overcoordination.cjs
 * Known issue: both currently fall into a bound well. Stretched bonds count as less used valence
 * (satF) than the Morse attraction they still provide, so the neighbours of a stretched centre look
 * open-shell and bond to one another. Hot hydrocarbons reach these wells within a picosecond. */
const { Engine } = require('../engine.js');
const { geo } = require('./reactions.cjs');
const q = { thermal: false };

function scan(label, mol, attacker, rs) {
  const M = geo(mol), A = geo(attacker), h0 = M[1];
  const L = Math.hypot(h0[1] - M[0][1], h0[2] - M[0][2], h0[3] - M[0][3]);
  const u = [(h0[1] - M[0][1]) / L, (h0[2] - M[0][2]) / L, (h0[3] - M[0][3]) / L];
  const rows = [];
  let Einf;
  for (const r of [6, ...rs]) {
    const e = new Engine({ width: 200, height: 200, depth: 200, T: 0, thermostat: false });
    e.box = { x0: -100, x1: 100, y0: -100, y1: 100, z0: -100, z1: 100 };
    for (const [s, x, y, z] of M) e.addAtom(s, x, y, z, q);
    const base = e.N, ax = M[0][1] - u[0] * r, ay = M[0][2] - u[1] * r, az = M[0][3] - u[2] * r;
    for (const [s, x, y, z] of A) e.addAtom(s, ax + (x - A[0][1]), ay + (y - A[0][2]), az + (z - A[0][3]), q);
    e.pinned[0] = e.pinned[base] = 1; e.touch(); e.refresh();
    for (let k = 0; k < 2; k++) e.minimize(3000, 0.05);
    e.refresh(); for (let k = 0; k < 200; k++) e.computeForces(1);
    const E = e.computeForces();
    if (Einf === undefined) { Einf = E; continue; }
    rows.push({ 'C···X (Å)': r, 'ΔE (kJ/mol)': +(E - Einf).toFixed(1), fragments: e.fragments().list.map(f => e.formulaOf(f)).join(' + ') });
  }
  console.log('\n' + label); console.table(rows);
  return Math.min(...rows.map(x => x['ΔE (kJ/mol)']));
}

const ch5 = scan('H pushed onto the back of CH4 (real: repulsive, CH5 does not exist)', 'CH4', 'H', [2.0, 1.8, 1.6, 1.45, 1.3, 1.2, 1.1]);
const c2 = scan('CH3 pushed onto the back of CH4 (real: repulsive)', 'CH4', 'CH3', [2.6, 2.2, 2.0, 1.8, 1.6]);
console.log(`\nDeepest well: CH5 ${ch5.toFixed(1)} kJ/mol, five-coordinate carbon ${c2.toFixed(1)} kJ/mol (both should be > 0).`);
