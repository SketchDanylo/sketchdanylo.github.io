/* Diagnostic, not a validation pass: report the modeled entrance path for
 * Cl· + ethene. No product, velocity kick or bond assignment is injected.
 * Reference PES: https://doi.org/10.1021/jp001221u */
const { Engine } = require('../engine.js');
const { geo, E } = require('./reactions.cjs');
const atoms = geo('C2H4'), reference = E('C2H4');
const rows = [];
for (let r = 4; r >= 1.7; r -= .1) {
  const e = new Engine({ T: 0, thermostat: false });
  e.box = { x0: -100, x1: 100, y0: -100, y1: 100, z0: -100, z1: 100 };
  for (const a of atoms) e.addAtom(...a, { thermal: false });
  e.addAtom('Cl', atoms[0][1], atoms[0][2], atoms[0][3] + r, { thermal: false });
  e.pinned[0] = e.pinned[6] = 1; e.refresh();
  for (let k = 0; k < 3; k++) e.minimize(2000, .02);
  const pair = (i, j) => e.pairMap.get(i * 1048576 + j);
  rows.push({ r: +r.toFixed(2), dE: +(e.Epot - reference).toFixed(2), cc: +e.pN[pair(0, 1)].toFixed(3), cCl: +e.bondStrength(pair(0, 6)).toFixed(3) });
}
console.table(rows);
console.log('Constrained entrance barrier (kJ/mol):', Math.max(...rows.filter(x => x.r >= 1.9).map(x => x.dE)));
