const assert = require('node:assert/strict');
const { Engine } = require('../engine.js');
const K = require('../kinetics.js');
const { geo, T } = require('./reactions.cjs');
T.CH3Cl = [['C', 0, 0, 0], ['Cl', 1.78, 0, 0], ['H', -0.36, 1.03, 0], ['H', -0.36, -0.51, 0.89], ['H', -0.36, -0.51, -0.89]];

function scene(parts) {
  const e = new Engine({ width: 30, height: 30, depth: 30, T: 0, thermostat: false });
  for (const [name, off, orders] of parts) {
    const b = e.N;
    for (const [s, x, y, z] of geo(name)) e.addAtom(s, 15 + off[0] + x, 15 + off[1] + y, off[2] + (z || 0), { thermal: false });
    e.touch(); e.refresh();
    for (const [a, c, o] of orders || []) e.setBondOrder(b + a, b + c, o);
  }
  e.touch(); e.refresh();
  return e;
}
const along = (name, k, d) => { const a = geo(name)[k]; const n = Math.hypot(a[1], a[2], a[3] || 0) || 1; return [a[1] / n * d + a[1], a[2] / n * d + a[2], (a[3] || 0) / n * d + (a[3] || 0)]; };
const out = (name, k, d) => { const g = geo(name), a = g[k], o = g[0], v = [a[1] - o[1], a[2] - o[2], (a[3] || 0) - (o[3] || 0)], n = Math.hypot(...v); return [a[1] + v[0] / n * d, a[2] + v[1] / n * d, (a[3] || 0) + v[2] / n * d]; };
const end = (name, d) => { const a = geo(name)[1]; return [a[1] + d, a[2], a[3] || 0]; };
const CC = [[0, 1, 2]];

const BARRIERS = [
  ['H + H₂ → H₂ + H', [['H2', [0, 0, 0]], ['H', [geo('H2')[1][1] + 2.6, 0, 0]]], 2, 1, ['H2', 'H·'], 40],
  ['H + CH₄ → H₂ + CH₃', [['CH4', [0, 0, 0]], ['H', along('CH4', 1, 1.6)]], 5, 1, ['CH3·', 'H2'], 50],
  ['Cl + CH₄ → HCl + CH₃', [['CH4', [0, 0, 0]], ['Cl', along('CH4', 1, 2.2)]], 5, 1, ['CH3·', 'HCl'], 11],
  ['Cl + H₂ → HCl + H', [['H2', [0, 0, 0]], ['Cl', [geo('H2')[1][1] + 3.0, 0, 0]]], 2, 1, ['H·', 'HCl'], 19],
  ['H + Cl₂ → HCl + Cl', [['Cl2', [0, 0, 0]], ['H', [geo('Cl2')[1][1] + 2.6, 0, 0]]], 2, 1, ['Cl·', 'HCl'], 5],
  ['CH₃ + Cl₂ → CH₃Cl + Cl', [['Cl2', [0, 0, 0]], ['CH3', [geo('Cl2')[1][1] + 3.2, 0, 0]]], 2, 1, ['CH3Cl', 'Cl·'], 2],
  ['OH + CH₄ → H₂O + CH₃', [['CH4', [0, 0, 0]], ['OH', along('CH4', 1, 2.0)]], 5, 1, ['CH3·', 'H2O'], 15],
  ['H + C₂H₄ → C₂H₅', [['C2H4', [0, 0, 0], CC], ['H', [0, 0, 3.0]]], 6, 0, ['C2H5·'], 9],
  ['Cl + C₂H₄ → C₂H₄Cl', [['C2H4', [0, 0, 0], CC], ['Cl', [0, 0, 3.2]]], 6, 0, ['C2H4Cl·'], 0],
  ['CH₃ + C₂H₄ → C₃H₇', [['C2H4', [0, 0, 0], CC], ['CH3', [0, 0, 3.4]]], 6, 0, ['C3H7·'], 31],
  ['CH₃ + CH₃ → C₂H₆', [['CH3', [0, 0, 0]], ['CH3', [3.6, 0, 0]]], 0, 4, ['C2H6'], 0],
  ['H + O₂ → OH + O', [['O2', [0, 0, 0], [[0, 1, 2]]], ['H', [geo('O2')[1][1] + 2.6, 0, 0]]], 2, 1, ['HO·', 'O·'], 70],
  ['H + C₂H₆ → H₂ + C₂H₅', [['C2H6', [0, 0, 0]], ['H', out('C2H6', 2, 1.6)]], 8, 2, ['C2H5·', 'H2'], 38],
  ['Cl + C₂H₆ → HCl + C₂H₅', [['C2H6', [0, 0, 0]], ['Cl', out('C2H6', 2, 2.2)]], 8, 2, ['C2H5·', 'HCl'], 1],
  ['OH + C₂H₆ → H₂O + C₂H₅', [['C2H6', [0, 0, 0]], ['OH', out('C2H6', 2, 2.0)]], 8, 2, ['C2H5·', 'H2O'], 8],
  ['OH + H₂ → H₂O + H', [['H2', [0, 0, 0]], ['OH', end('H2', 2.4)]], 2, 1, ['H2O', 'H·'], 15],
  ['H + H₂O → H₂ + OH', [['H2O', [0, 0, 0]], ['H', out('H2O', 1, 1.6)]], 3, 1, ['H2', 'HO·'], 76],
  ['H + HCl → H₂ + Cl', [['HCl', [0, 0, 0]], ['H', [-2.0, 0, 0]]], 2, 0, ['Cl·', 'H2'], 15],
  ['F + H₂ → HF + H', [['H2', [0, 0, 0]], ['F', end('H2', 2.4)]], 2, 1, ['H·', 'HF'], 4],
  ['F + CH₄ → HF + CH₃', [['CH4', [0, 0, 0]], ['F', out('CH4', 1, 2.0)]], 5, 1, ['CH3·', 'HF'], 2],
  ['H + F₂ → HF + F', [['F2', [0, 0, 0]], ['H', end('F2', 2.4)]], 2, 1, ['F·', 'HF'], 10],
  ['CH₃ + H₂ → CH₄ + H', [['H2', [0, 0, 0]], ['CH3', end('H2', 2.6)]], 2, 1, ['CH4', 'H·'], 45],
  ['OH + NH₃ → H₂O + NH₂', [['NH3', [0, 0, 0]], ['OH', out('NH3', 1, 2.0)]], 4, 1, ['H2O', 'NH2·'], 6],
];
const BONDS = [
  ['H–H', [['H2', [0, 0, 0]]], 0, 1, 436],
  ['CH₃–H', [['CH4', [0, 0, 0]]], 0, 1, 439],
  ['C₂H₅–H', [['C2H6', [0, 0, 0]]], 0, 2, 423],
  ['CH₃–CH₃', [['C2H6', [0, 0, 0]]], 0, 1, 377],
  ['C₂H₃–H (vinyl)', [['C2H4', [0, 0, 0], CC]], 0, 2, 465],
  ['C₆H₅–H', [['C6H6', [0, 0, 0], [[0, 1, 1.5], [1, 2, 1.5], [2, 3, 1.5], [3, 4, 1.5], [4, 5, 1.5], [5, 0, 1.5]]]], 0, 6, 473],
  ['HO–H', [['H2O', [0, 0, 0]]], 0, 1, 497],
  ['Cl–Cl', [['Cl2', [0, 0, 0]]], 0, 1, 243],
  ['H–Cl', [['HCl', [0, 0, 0]]], 0, 1, 432],
  ['O=O', [['O2', [0, 0, 0], [[0, 1, 2]]]], 0, 1, 498],
  ['H₂C=CH₂', [['C2H4', [0, 0, 0], CC]], 0, 1, 728],
  ['CH₃–Cl', [['CH3Cl', [0, 0, 0]]], 0, 1, 350],
  ['H–F', [['HF', [0, 0, 0]]], 0, 1, 570],
  ['F–F', [['F2', [0, 0, 0]]], 0, 1, 159],
  ['H₂N–H', [['NH3', [0, 0, 0]]], 0, 1, 450],
  ['HCO–H (formyl)', [['CH2O', [0, 0, 0], [[0, 1, 2]]]], 0, 2, 369],
];

const rows = [];
for (const [name, parts, i, j, products, real] of BARRIERS) {
  const r = K.scanPair(scene(parts), i, j);
  assert.ok(r.ok, name + ': ' + r.reason);
  const want = products.slice().sort().join(' + ');
  const hit = r.products.join(' + ') === want ? r : (r.alternatives || []).find(a => a.products.join(' + ') === want);
  assert.ok(hit, name + ' was not found; the scan gave ' + r.products.join(' + '));
  rows.push(['barrier', name + (hit === r ? '' : ' *'), hit.Ea, real]);
}
for (const [name, parts, i, j, real] of BONDS) {
  const r = K.scanPair(scene(parts), i, j);
  assert.ok(r.ok && r.kind === 'break', name + ': ' + r.reason);
  rows.push(['bond', name, r.Ea, real]);
}
const pad = (s, n) => String(s).padEnd(n), num = x => x.toFixed(0).padStart(5);
console.log(pad('', 8) + pad('reaction', 28) + 'model  real  error');
for (const [kind, name, model, real] of rows) console.log(pad(kind, 8) + pad(name, 28) + num(model) + ' ' + num(real) + ' ' + num(model - real) + (Math.abs(model - real) > 15 ? '  ✗' : ''));
const rms = list => Math.sqrt(list.reduce((s, r) => s + (r[2] - r[3]) ** 2, 0) / list.length);
console.log('rms error: barriers ' + rms(rows.filter(r => r[0] === 'barrier')).toFixed(1) + ' kJ/mol, bonds ' + rms(rows.filter(r => r[0] === 'bond')).toFixed(1) + ' kJ/mol');
console.log('* not the fastest route from these reactants; read from the slower routes the scan also found.');
console.log('benchmark: ' + rows.length + ' reactions scanned, products as expected.');
