/* Textbook-reaction checks for the Chem Playground engine: node playground/tests/reactions.cjs [--dyn]
 * A. elementary steps: barrier (relaxed collinear scan) and reaction energy vs literature
 * B. overall reaction energies from relaxed molecules vs standard gas-phase ΔH
 * C. (--dyn) molecular mixtures run at several temperatures; products are counted */
const { Engine } = require('../engine.js');

function eng(opts = {}) {
  const e = new Engine({ T: 0, thermostat: false, seed: 1, ...opts });
  e.box = { x0: -200, x1: 200, y0: -200, y1: 200, z0: -200, z1: 200 };
  e.recording = false;
  return e;
}
function build(atoms, pins = []) {
  const e = eng();
  for (const a of atoms) e.addAtom(a[0], a[1], a[2], a[3] || 0, { thermal: false });
  e.refresh();
  for (const i of pins) e.pinned[i] = 1;
  return e;
}
function energyOnly(e) { for (let k = 0; k < 400; k++) e.computeForces(1); return e.Epot; } // converge lagged bond orders at fixed geometry
function relax(e, it = 3) { for (let k = 0; k < it; k++) e.minimize(3000, 0.02); return e.Epot; }
const atomsOf = e => Array.from({ length: e.N }, (_, i) => [Engine.ELEMENTS[e.type[i]].sym, e.pos[3 * i], e.pos[3 * i + 1], e.pos[3 * i + 2]]);

/* molecule templates (rough geometry; relaxed before use) */
const T = {
  H: [['H', 0, 0, 0]], O: [['O', 0, 0, 0]], Cl: [['Cl', 0, 0, 0]], F: [['F', 0, 0, 0]], N: [['N', 0, 0, 0]], C: [['C', 0, 0, 0]],
  H2: [['H', 0, 0, 0], ['H', 0.74, 0, 0]], O2: [['O', 0, 0, 0], ['O', 1.21, 0, 0]], N2: [['N', 0, 0, 0], ['N', 1.1, 0, 0]],
  Cl2: [['Cl', 0, 0, 0], ['Cl', 1.99, 0, 0]], F2: [['F', 0, 0, 0], ['F', 1.42, 0, 0]],
  HCl: [['H', 0, 0, 0], ['Cl', 1.27, 0, 0]], HF: [['H', 0, 0, 0], ['F', 0.92, 0, 0]], OH: [['O', 0, 0, 0], ['H', 0.97, 0, 0]],
  H2O: [['O', 0, 0, 0], ['H', 0.96, 0, 0], ['H', -0.24, 0.93, 0]],
  HO2: [['O', 0, 0, 0], ['O', 1.33, 0, 0], ['H', -0.3, 0.93, 0]],
  NH3: [['N', 0, 0, 0], ['H', 1.0, 0, 0.3], ['H', -0.5, 0.87, 0.3], ['H', -0.5, -0.87, 0.3]],
  CH4: [['C', 0, 0, 0], ['H', 0.63, 0.63, 0.63], ['H', -0.63, -0.63, 0.63], ['H', -0.63, 0.63, -0.63], ['H', 0.63, -0.63, -0.63]],
  CH3: [['C', 0, 0, 0], ['H', 1.08, 0, 0], ['H', -0.54, 0.93, 0], ['H', -0.54, -0.93, 0.05]],
  CO2: [['C', 0, 0, 0], ['O', 1.16, 0, 0], ['O', -1.16, 0.05, 0]], CO: [['C', 0, 0, 0], ['O', 1.13, 0, 0]],
  C2H4: [['C', 0, 0, 0], ['C', 1.34, 0, 0], ['H', -0.55, 0.93, 0], ['H', -0.55, -0.93, 0], ['H', 1.89, 0.93, 0], ['H', 1.89, -0.93, 0]],
  C2H6: [['C', 0, 0, 0], ['C', 1.54, 0, 0], ['H', -0.4, 1, 0], ['H', -0.4, -0.5, 0.9], ['H', -0.4, -0.5, -0.9], ['H', 1.94, -1, 0], ['H', 1.94, 0.5, 0.9], ['H', 1.94, 0.5, -0.9]],
  NaCl: [['Na', 0, 0, 0], ['Cl', 2.36, 0, 0]], Na: [['Na', 0, 0, 0]]
};
const Emol = {};
function E(name) { if (Emol[name] === undefined) Emol[name] = relax(build(T[name])); return Emol[name]; }
const relaxedAtoms = {};
function geo(name) { if (!relaxedAtoms[name]) { const e = build(T[name]); relax(e); relaxedAtoms[name] = atomsOf(e); } return relaxedAtoms[name]; }

/* A. collinear atom transfer A + B–C → A–B + C with B the transferred atom */
const FAST = { on: false };
function transfer(attacker, aIdx, mol, bIdx, cIdx) {
  const M = geo(mol), X = geo(attacker);
  const B = M[bIdx], C = M[cIdx];
  const u = [B[1] - C[1], B[2] - C[2], B[3] - C[3]], L = Math.hypot(...u), n = u.map(v => v / L);
  const scanE = (rbc, rab) => {
    const pts = M.map(a => a.slice());
    // move B along the C→B axis to r_BC; drag the rest of B's side? B is terminal in all tests
    pts[bIdx] = [B[0], C[1] + n[0] * rbc, C[2] + n[1] * rbc, C[3] + n[2] * rbc];
    const ax = C[1] + n[0] * (rbc + rab), ay = C[2] + n[1] * (rbc + rab), az = C[3] + n[2] * (rbc + rab);
    // attacker: its atom aIdx at (ax,ay,az); other attacker atoms continue outward along n
    const A0 = X[aIdx];
    const extra = X.map((a, k) => { if (k === aIdx) return [a[0], ax, ay, az]; const d = Math.hypot(a[1] - A0[1], a[2] - A0[2], a[3] - A0[3]); return [a[0], ax + n[0] * d * 0.34, ay + n[1] * d * 0.94 + 0.3, az + n[2] * d * 0.34]; });
    const all = [...pts, ...extra], aGlobal = M.length + aIdx;
    const e = build(all, [cIdx, bIdx, aGlobal]);
    if (all.length > 3) { if (FAST.on) e.minimize(400, 0.2); else relax(e, 1); }
    else energyOnly(e);
    return e.Epot;
  };
  const Er = E(mol) + E(attacker);
  let pathMax = -Infinity;
  const r0 = Math.hypot(...u);
  const ds = FAST.on ? 0.15 : 0.1, dr = FAST.on ? 0.07 : 0.05;
  for (let s = -1.2; s <= 1.2001; s += ds) {
    let best = Infinity;
    for (let rbc = r0 * 0.95; rbc <= r0 + 1.6; rbc += dr) { const rab = rbc - s; if (rab < 0.7 || rab > 3.2) continue; best = Math.min(best, scanE(rbc, rab)); }
    pathMax = Math.max(pathMax, best); FAST.well = Math.min(FAST.well, best - Er);
  }
  return pathMax - Er;
}
function resetCache() { for (const k in Emol) delete Emol[k]; for (const k in relaxedAtoms) delete relaxedAtoms[k]; }
function elementary() {
  const out = {}; FAST.well = 0;
  out.HH2 = transfer('H', 0, 'H2', 1, 0);
  out.ClH2 = transfer('Cl', 0, 'H2', 1, 0);
  out.well = FAST.well; // deepest complex along the thermoneutral H + H2 and Cl + H2 paths (should be ≳ −15)
  out.HCH4 = transfer('H', 0, 'CH4', 1, 0);
  out.HCl2 = transfer('H', 0, 'Cl2', 1, 0);
  out.FH2 = transfer('F', 0, 'H2', 1, 0);
  out.OH2 = transfer('O', 0, 'H2', 1, 0);
  out.OHH2 = transfer('OH', 0, 'H2', 1, 0);
  out.HO2t = transfer('H', 0, 'O2', 1, 0);
  out.HO2a = addition();
  return out;
}
function addition() {
  const O2 = geo('O2'); let maxE = -Infinity;
  for (let r = 3.2; r >= 0.9; r -= 0.1) {
    const hx = O2[0][1] + r * Math.cos(1.9), hy = O2[0][2] + r * Math.sin(1.9);
    const e = build([...O2, ['H', hx, hy, 0]], [0, 2]); if (FAST.on) e.minimize(300, 0.3); else relax(e, 1);
    maxE = Math.max(maxE, e.Epot - E('O2') - E('H'));
  }
  return maxE;
}
function dimer(A, B, shift) {
  const both = [...geo(A), ...geo(B).map(b => [b[0], b[1] + shift[0], b[2] + shift[1], b[3] + shift[2]])];
  const e = build(both); relax(e, 2); return e.Epot - E(A) - E(B);
}
module.exports = { FAST, elementary, dimer, dE: (r, p) => p.reduce((s, x) => s + E(x), 0) - r.reduce((s, x) => s + E(x), 0), resetCache, E, geo, T };
if (require.main !== module) return;
const rows = [];
function step(label, barrier, dH, lit) { rows.push({ reaction: label, 'barrier (kJ/mol)': +barrier.toFixed(1), 'lit. Ea': lit[0], 'ΔE model': +dH.toFixed(1), 'lit. ΔH': lit[1] }); }
function dE(reac, prod) { return prod.reduce((s, x) => s + E(x), 0) - reac.reduce((s, x) => s + E(x), 0); }

step('H + H2 → H2 + H', transfer('H', 0, 'H2', 1, 0), 0, [40, 0]);
step('H + CH4 → H2 + CH3', transfer('H', 0, 'CH4', 1, 0), dE(['H', 'CH4'], ['H2', 'CH3']), [50, 3]);
step('Cl + H2 → HCl + H', transfer('Cl', 0, 'H2', 1, 0), dE(['Cl', 'H2'], ['HCl', 'H']), [23, 4]);
step('H + Cl2 → HCl + Cl', transfer('H', 0, 'Cl2', 1, 0), dE(['H', 'Cl2'], ['HCl', 'Cl']), [8, -189]);
step('F + H2 → HF + H', transfer('F', 0, 'H2', 1, 0), dE(['F', 'H2'], ['HF', 'H']), [4, -130]);
step('O + H2 → OH + H', transfer('O', 0, 'H2', 1, 0), dE(['O', 'H2'], ['OH', 'H']), [37, 8]);
step('OH + H2 → H2O + H', transfer('OH', 0, 'H2', 1, 0), dE(['OH', 'H2'], ['H2O', 'H']), [17, -62]);
step('H + O2 → OH + O', transfer('H', 0, 'O2', 1, 0), dE(['H', 'O2'], ['OH', 'O']), [70, 70]);
// radical addition H + O2 → HO2 (should be barrierless, ≈ −205)
{
  const O2 = geo('O2'); let maxE = -Infinity, minE = Infinity;
  for (let r = 3.2; r >= 0.9; r -= 0.1) {
    const hx = O2[0][1] + r * Math.cos(1.9), hy = O2[0][2] + r * Math.sin(1.9);
    const e = build([...O2, ['H', hx, hy, 0]], [0, 2]); relax(e, 1);
    const v = e.Epot - E('O2') - E('H'); maxE = Math.max(maxE, v); minE = Math.min(minE, v);
  }
  step('H + O2 → HO2 (addition)', maxE, dE(['H', 'O2'], ['HO2']), [0, -205]);
}
console.log('\nA. Elementary steps'); console.table(rows); rows.length = 0;

const overall = [
  ['2 H2 + O2 → 2 H2O', ['H2', 'H2', 'O2'], ['H2O', 'H2O'], -483.6],
  ['H2 + Cl2 → 2 HCl', ['H2', 'Cl2'], ['HCl', 'HCl'], -184.6],
  ['H2 + F2 → 2 HF', ['H2', 'F2'], ['HF', 'HF'], -546],
  ['N2 + 3 H2 → 2 NH3', ['N2', 'H2', 'H2', 'H2'], ['NH3', 'NH3'], -91.8],
  ['CH4 + 2 O2 → CO2 + 2 H2O', ['CH4', 'O2', 'O2'], ['CO2', 'H2O', 'H2O'], -802.3],
  ['2 CO + O2 → 2 CO2', ['CO', 'CO', 'O2'], ['CO2', 'CO2'], -566],
  ['C2H4 + H2 → C2H6', ['C2H4', 'H2'], ['C2H6'], -136.3],
  ['Na + ½ Cl2 → NaCl(g)', ['Na', 'Cl'], ['NaCl'], -412 + 0],
];
const t2 = overall.map(([l, r, p, lit]) => ({ reaction: l, 'ΔE model (kJ/mol)': +dE(r, p).toFixed(1), 'lit. ΔH (kJ/mol)': lit }));
console.log('\nB. Overall reaction energies (Na row: Na + Cl atoms → NaCl, lit. = bond energy)'); console.table(t2);

if (process.argv.includes('--dyn')) {
  function mix(recipe, Tk, ps, seed = 3) {
    const e = new Engine({ width: 40, height: 24, depth: 12, T: Tk, seed }); e.recording = false;
    let s = seed * 7919; const r = () => ((s = (s * 16807) % 2147483647) / 2147483647);
    for (const [name, n] of recipe) {
      const g = geo(name);
      for (let k = 0; k < n; k++) {
        let x, y, z, ok = false, tries = 0;
        while (!ok && tries++ < 500) {
          x = 3 + r() * 34; y = 3 + r() * 18; z = (r() - 0.5) * 6; ok = true;
          for (let i = 0; i < e.N; i++) if (Math.hypot(e.pos[3 * i] - x, e.pos[3 * i + 1] - y, e.pos[3 * i + 2] - z) < 3.6) { ok = false; break; }
        }
        const th = r() * 6.283, c = Math.cos(th), sn = Math.sin(th);
        for (const a of g) e.addAtom(a[0], x + a[1] * c - a[2] * sn, y + a[1] * sn + a[2] * c, z + a[3]);
      }
    }
    e.refresh();
    let err = null;
    for (let k = 0; k < ps * 1000; k++) { try { e.step(); } catch (x) { err = x.message; break; } }
    const fr = e.fragments(), cnt = {};
    for (const g of fr.list) { const f = e.formulaOf(g) + (e.isRadical(g) ? '·' : ''); cnt[f] = (cnt[f] || 0) + 1; }
    return (err ? '[STOPPED: ' + err + '] ' : '') + Object.entries(cnt).sort((a, b) => b[1] - a[1]).map(([f, n]) => n + ' ' + f).join(', ') + '  (T ' + e.temperature().toFixed(0) + ' K)';
  }
  const ps = +(process.argv.find(a => a.startsWith('--ps=')) || '--ps=30').slice(5);
  console.log('\nC. Mixtures, ' + ps + ' ps each');
  for (const [label, recipe, temps] of [
    ['H2 + O2', [['H2', 12], ['O2', 6]], [300, 1500, 3000]],
    ['H2 + Cl2', [['H2', 10], ['Cl2', 10]], [300, 1500, 3000]],
    ['CH4 + O2', [['CH4', 5], ['O2', 10]], [300, 2000, 3500]],
    ['H + Cl2 (radical)', [['H', 6], ['Cl2', 10]], [300]],
    ['Na + Cl2', [['Na', 10], ['Cl2', 5]], [300, 1000]],
  ]) for (const Tk of temps) console.log(label.padEnd(18), String(Tk).padStart(5) + ' K →', mix(recipe, Tk, ps));
}
