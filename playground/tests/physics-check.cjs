/* Physics validation for the Chem Playground engine: node playground/tests/physics-check.cjs [--strict] */
const { Engine } = require('../engine.js');
const assert = require('node:assert/strict');

function eng(opts = {}) {
  const e = new Engine({ width: 200, height: 200, depth: 200, T: 0, thermostat: false, ...opts });
  e.box = { x0: -100, x1: 100, y0: -100, y1: 100, z0: -100, z1: 100 };
  e.recording = false;
  return e;
}
function build(atoms, opts) {
  const e = eng(opts);
  for (const a of atoms) e.addAtom(a[0], a[1], a[2], a[3] || 0, { thermal: false, charge: a[4] || 0 });
  e.refresh();
  return e;
}
const dist = (e, i, j) => Math.hypot(e.pos[3 * i] - e.pos[3 * j], e.pos[3 * i + 1] - e.pos[3 * j + 1], e.pos[3 * i + 2] - e.pos[3 * j + 2]);
function angle(e, j, i, k) {
  const u = [0, 1, 2].map(d => e.pos[3 * j + d] - e.pos[3 * i + d]), v = [0, 1, 2].map(d => e.pos[3 * k + d] - e.pos[3 * i + d]);
  const c = (u[0] * v[0] + u[1] * v[1] + u[2] * v[2]) / Math.hypot(...u) / Math.hypot(...v);
  return Math.acos(Math.max(-1, Math.min(1, c))) * 180 / Math.PI;
}
function energy(e) { e.needRebuild = true; e._checkRebuild(); const save = e.pN.slice(0, e.nPairs), c0 = e.cos0.slice(0, e.N); const E = e.computeForces(); e.pN.set(save); e.cos0.set(c0); return E; }
const jitter = (atoms, s = 0.15, seed = 1) => { let x = seed; const r = () => ((x = (x * 16807) % 2147483647) / 2147483647 - 0.5) * 2 * s; return atoms.map(a => [a[0], a[1] + r(), a[2] + r(), (a[3] || 0) + r(), a[4]]); };

const results = [];
const report = (name, got, want, unit, tol) => {
  const ok = tol == null ? true : Math.abs(got - want) <= tol;
  results.push({ name, got: +got.toFixed(3), want, unit, ok });
};

/* 1. analytic forces vs numerical gradient (frozen lagged state) */
function gradCheck(label, atoms) {
  const e = build(jitter(atoms, 0.25, 7));
  // relax lagged variables into a consistent state
  for (let k = 0; k < 200; k++) e.computeForces();
  const pN = e.pN.slice(0, e.nPairs), c0 = e.cos0.slice(0, e.N);
  const freeze = () => { e.pN.set(pN); e.cos0.set(c0); };
  freeze(); e.computeForces(); freeze();
  const F = e.frc.slice(0, 3 * e.N);
  // hold charges & polar-H radii (lagged) fixed as well by evaluating at nearby points only
  let maxErr = 0, maxF = 0; const h = 1e-5;
  for (let k = 0; k < 3 * e.N; k++) {
    const x0 = e.pos[k];
    e.pos[k] = x0 + h; freeze(); const Ep = e.computeForces();
    e.pos[k] = x0 - h; freeze(); const Em = e.computeForces();
    e.pos[k] = x0;
    const num = -(Ep - Em) / (2 * h);
    maxErr = Math.max(maxErr, Math.abs(num - F[k])); maxF = Math.max(maxF, Math.abs(F[k]));
  }
  freeze();
  results.push({ name: 'grad ' + label, got: +maxErr.toExponential(2), want: '< 1e-2·|F|max (' + maxF.toFixed(0) + ')', unit: 'kJ/mol/Å', ok: maxErr < 1e-2 * Math.max(1, maxF) });
}

const H2O = [['O', 0, 0, 0], ['H', 0.96, 0, 0], ['H', -0.24, 0.93, 0]];
const CH4 = [['C', 0, 0, 0], ['H', 0.63, 0.63, 0.63], ['H', -0.63, -0.63, 0.63], ['H', -0.63, 0.63, -0.63], ['H', 0.63, -0.63, -0.63]];
gradCheck('H2O', H2O);
gradCheck('CH4', CH4);
gradCheck('H+H2 (TS region)', [['H', 0, 0, 0], ['H', 0.95, 0, 0], ['H', 1.95, 0.05, 0]]);
gradCheck('ethene', [['C', 0, 0, 0], ['C', 1.34, 0, 0], ['H', -0.55, 0.93, 0], ['H', -0.55, -0.93, 0], ['H', 1.89, 0.93, 0], ['H', 1.89, -0.93, 0]]);
gradCheck('water dimer', [...H2O, ['O', 2.9, 0, 0.1], ['H', 3.3, 0.8, 0.2], ['H', 3.3, -0.8, 0.2]]);
gradCheck('Na/Cl triangle cluster', [['Na', 0, 0, 0], ['Na', 3.3, 0.2, 0], ['Cl', 1.6, 2.6, 0.3], ['Cl', 1.4, -2.7, -0.2], ['Na', 4.6, 2.9, 0.4]]);
gradCheck('H + H2O (attack)', [...H2O, ['H', 1.9, 0.6, 0.4]]);
gradCheck('ozone-like O3 ring strain', [['O', 0, 0, 0], ['O', 1.3, 0, 0], ['O', 0.65, 1.1, 0]]);

/* 2. structures after minimisation */
function relaxed(atoms) { const e = build(atoms); for (let k = 0; k < 4; k++) e.minimize(3000, 0.05); return e; }
{ const e = relaxed([['H', 0, 0, 0], ['H', 0.9, 0, 0]]); report('H–H length', dist(e, 0, 1), 0.741, 'Å', 0.01); report('H2 bond energy', -e.Epot, 436, 'kJ/mol', 5); }
{ const e = relaxed(H2O); report('O–H length', dist(e, 0, 1), 0.96, 'Å', 0.02); report('H–O–H angle', angle(e, 1, 0, 2), 104.5, '°', 3); }
{ const e = relaxed(CH4); report('C–H length', dist(e, 0, 1), 1.09, 'Å', 0.02); report('H–C–H angle', angle(e, 1, 0, 2), 109.47, '°', 2); report('CH4 atomisation', -e.Epot, 1652, 'kJ/mol', 60); }
{ const e = relaxed([['C', 0, 0, 0], ['C', 1.5, 0, 0], ['H', -0.4, 1, 0], ['H', -0.4, -0.5, 0.9], ['H', -0.4, -0.5, -0.9], ['H', 1.9, -1, 0], ['H', 1.9, 0.5, 0.9], ['H', 1.9, 0.5, -0.9]]); report('ethane C–C', dist(e, 0, 1), 1.54, 'Å', 0.03); }
{ const e = relaxed([['C', 0, 0, 0], ['C', 1.4, 0, 0], ['H', -0.55, 0.93, 0.1], ['H', -0.55, -0.93, 0], ['H', 1.95, 0.93, 0], ['H', 1.95, -0.93, -0.1]]); report('ethene C=C', dist(e, 0, 1), 1.34, 'Å', 0.03); report('ethene H–C=C', angle(e, 2, 0, 1), 121.5, '°', 5); }
{ const e = relaxed([['C', 0, 0, 0], ['C', 1.25, 0.05, 0], ['H', -1.05, 0.1, 0], ['H', 2.3, -0.1, 0.05]]); report('ethyne C≡C', dist(e, 0, 1), 1.20, 'Å', 0.03); report('ethyne H–C≡C', angle(e, 2, 0, 1), 180, '°', 3); }
{ const e = relaxed([['O', 0, 0, 0], ['O', 1.3, 0, 0]]); report('O=O length', dist(e, 0, 1), 1.21, 'Å', 0.02); report('O2 bond energy', -e.Epot, 498, 'kJ/mol', 10); }
{ const e = relaxed([['N', 0, 0, 0], ['N', 1.2, 0, 0]]); report('N≡N length', dist(e, 0, 1), 1.10, 'Å', 0.02); report('N2 bond energy', -e.Epot, 945, 'kJ/mol', 15); }
{ const e = relaxed([['C', 0, 0, 0], ['O', 1.25, 0, 0], ['O', -1.25, 0.1, 0]]); report('CO2 C=O', dist(e, 0, 1), 1.16, 'Å', 0.06); report('CO2 angle', angle(e, 1, 0, 2), 180, '°', 3); }
{
  const ring = []; for (let k = 0; k < 6; k++) { const a = k * Math.PI / 3; ring.push(['C', 1.4 * Math.cos(a), 1.4 * Math.sin(a), 0]); }
  for (let k = 0; k < 6; k++) { const a = k * Math.PI / 3; ring.push(['H', 2.48 * Math.cos(a), 2.48 * Math.sin(a), 0]); }
  const e = relaxed(ring); report('benzene C–C', dist(e, 0, 1), 1.39, 'Å', 0.03); report('benzene C–C–C', angle(e, 5, 0, 1), 120, '°', 2);
}
{ const e = relaxed([['N', 0, 0, 0], ['H', 1.0, 0, 0.3], ['H', -0.5, 0.87, 0.3], ['H', -0.5, -0.87, 0.3]]); report('NH3 H–N–H', angle(e, 1, 0, 2), 106.8, '°', 3); }

/* 3. intermolecular: water dimer, methane dimer */
function bindingEnergy(A, B, shift) {
  const eA = relaxed(A), eB = relaxed(B);
  const EA = eA.Epot, EB = eB.Epot;
  const both = [...A, ...B.map(b => [b[0], b[1] + shift[0], b[2] + shift[1], (b[3] || 0) + shift[2]])];
  const e = relaxed(both);
  return { dE: e.Epot - EA - EB, e };
}
{ const r = bindingEnergy(H2O, [['O', 0, 0, 0], ['H', 0.6, 0.75, 0], ['H', 0.6, -0.75, 0]], [2.95, 0.3, 0]); report('water dimer binding', r.dE, -21, 'kJ/mol', 8); let oo = dist(r.e, 0, 3); report('water dimer O···O', oo, 2.91, 'Å', 0.25); }
{ const r = bindingEnergy(CH4, CH4, [3.9, 0, 0]); report('methane dimer binding', r.dE, -2.2, 'kJ/mol', 2.5); }

/* 4. Reaction barriers live in playground/tests/reactions.cjs (relaxed minimum-energy paths).
      A symmetric-line scan underestimates them, so only the relaxed CH4 abstraction is checked here. */
/* 5. CH4 + H abstraction barrier: collinear H···H–C, relaxed scan (C, H_b, H_a pinned; CH3 umbrella relaxes) */
{
  const e0 = relaxed(CH4);
  const hx = e0.pos.slice(3, 6), cx = e0.pos.slice(0, 3); const u = [0, 1, 2].map(d => hx[d] - cx[d]); const L = Math.hypot(...u); const n = u.map(v => v / L);
  const atoms0 = []; for (let i = 0; i < 5; i++) atoms0.push([i ? 'H' : 'C', e0.pos[3 * i], e0.pos[3 * i + 1], e0.pos[3 * i + 2]]);
  const at = (rch, rhh) => {
    const pts = atoms0.map(a => a.slice());
    pts[1] = ['H', cx[0] + n[0] * rch, cx[1] + n[1] * rch, cx[2] + n[2] * rch];
    pts.push(['H', cx[0] + n[0] * (rch + rhh), cx[1] + n[1] * (rch + rhh), cx[2] + n[2] * (rch + rhh)]);
    const e = build(pts); e.pinned[0] = e.pinned[1] = e.pinned[5] = 1;
    e.minimize(600, 0.2); return energy(e);
  };
  // minimum-energy path over the coordinate s = r(C–H) − r(H···H), then its maximum
  let pathMax = -Infinity;
  for (let s = -1.0; s <= 1.0; s += 0.1) {
    let best = Infinity;
    for (let rch = 1.05; rch <= 2.3; rch += 0.05) { const rhh = rch - s; if (rhh < 0.7 || rhh > 2.3) continue; best = Math.min(best, at(rch, rhh)); }
    pathMax = Math.max(pathMax, best);
  }
  const far = build([...atoms0, ['H', 40, 40, 40]]);
  report('H + CH4 → H2 + CH3 barrier', pathMax - energy(far), 50, 'kJ/mol', 15);
}
/* 6. energy conservation (NVE, 2 ps, water + methane at 300 K) */
{
  const e = eng({ T: 300 });
  const mols = [H2O, CH4.map(a => [a[0], a[1], a[2] + 4, a[3]]), H2O.map(a => [a[0], a[1] + 3.2, a[2] + 1, a[3]]), CH4.map(a => [a[0], a[1] - 3.8, a[2] + 0.5, a[3] + 0.5])];
  for (const m of mols) for (const a of m) e.addAtom(a[0], a[1], a[2], a[3], { thermal: false });
  e.refresh(); e.minimize(2000, 0.1); e.thermalize(300); e.zeroMomentum();
  e.needForces = true;
  const E0 = e.Epot + e.kinetic(); let Emax = 0;
  for (let k = 0; k < 2000; k++) { e.step(); const Et = e.Epot + e.kinetic(); Emax = Math.max(Emax, Math.abs(Et - E0)); }
  report('NVE |ΔE| max over 2 ps', Emax, 0, 'kJ/mol', 3);
}
/* 7. deterministic step back */
{
  const e = eng({ T: 300, thermostat: true }); e.recording = true;
  for (const a of H2O) e.addAtom(a[0], a[1], a[2], a[3]);
  for (const a of CH4) e.addAtom(a[0], a[1] + 3.5, a[2], a[3]);
  e.refresh();
  for (let k = 0; k < 150; k++) e.step();
  const snap = e.pos.slice(0, 3 * e.N); e.step(); e.step();
  e.stepBack(); e.stepBack();
  let d = 0; for (let k = 0; k < snap.length; k++) d = Math.max(d, Math.abs(snap[k] - e.pos[k]));
  report('step back ×2 exactness', d, 0, 'Å', 1e-12);
}
/* 8. speed */
{
  const e = new Engine({ width: 50, height: 32, depth: 10, T: 300 }); e.recording = true;
  let s = 3; const r = () => ((s = (s * 16807) % 2147483647) / 2147483647);
  for (let m = 0; m < 30; m++) { const x = 5 + (m % 6) * 8, y = 5 + Math.floor(m / 6) * 5.5, M = m % 3 ? H2O : CH4; for (const a of M) e.addAtom(a[0], a[1] + x, a[2] + y, a[3]); }
  e.refresh(); e.minimize(500, 1);
  const t0 = performance.now(); let n = 0;
  while (performance.now() - t0 < 1500) { e.step(); n++; }
  results.push({ name: 'speed (' + e.N + ' atoms)', got: Math.round(n / 1.5), want: 'steps/s', unit: '', ok: true });
  results.push({ name: 'T after run', got: +e.temperature().toFixed(1), want: 300, unit: 'K', ok: true });
  results.push({ name: 'P (wall)', got: +e.pressureEMA.toFixed(1), want: '-', unit: 'bar', ok: true });
}

console.table(results.map(r => ({ check: r.name, got: r.got, want: r.want, unit: r.unit, ok: r.ok ? '✓' : '✗' })));
const bad = results.filter(r => !r.ok);
if (process.argv.includes('--strict')) assert.equal(bad.length, 0, bad.map(b => b.name).join(', '));
