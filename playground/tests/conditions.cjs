/* Right conditions and wrong ones. Each mixture is run where real chemistry says it must not react
 * on this timescale (negative controls) and, where the real reaction is fast enough to see in
 * picoseconds, where it must. No atom may ever end over its valence.
 *   node playground/tests/conditions.cjs
 * Timescale matters: the simulator runs picoseconds, so a reaction whose real initiation takes
 * microseconds (H2 + Cl2 at 1200 K: Cl2 splits once per ~0.3 ms) correctly shows nothing. The
 * positive controls here are the fast steps real chemistry does in picoseconds. */
const assert = require('node:assert/strict');
const { Engine } = require('../engine.js');
let count = 0;
function test(name, run) { run(); console.log('PASS', name); count++; }
const TPL = {
  H2: [['H', -.37, 0, 0], ['H', .37, 0, 0]], Cl2: [['Cl', -1, 0, 0], ['Cl', 1, 0, 0]], F2: [['F', -.71, 0, 0], ['F', .71, 0, 0]],
  O2: [['O', -.6, 0, 0], ['O', .6, 0, 0]], N2: [['N', -.55, 0, 0], ['N', .55, 0, 0]],
  CH4: [['C', 0, 0, 0], ['H', .63, .63, .63], ['H', -.63, -.63, .63], ['H', -.63, .63, -.63], ['H', .63, -.63, -.63]],
  CH3: [['C', 0, 0, 0], ['H', 1.08, 0, 0], ['H', -.54, .93, 0], ['H', -.54, -.93, 0]],
  C2H4: [['C', -.67, 0, 0], ['C', .67, 0, 0], ['H', -1.22, .93, 0], ['H', -1.22, -.93, 0], ['H', 1.22, .93, 0], ['H', 1.22, -.93, 0]],
  Na: [['Na', 0, 0, 0]], Ar: [['Ar', 0, 0, 0]], He: [['He', 0, 0, 0]], H: [['H', 0, 0, 0]], F: [['F', 0, 0, 0]],
};
const ORD = { O2: [[0, 1, 2]], N2: [[0, 1, 3]], C2H4: [[0, 1, 2]] };
const VAL = { H: 1, C: 4, N: 3, O: 2, F: 1, Cl: 1, Na: 6, Ar: 0, He: 0 };   // Na bonds ionically: Na2Cl2 is a real ring
function run(mix, T, ps, seed, L = 40) {
  const e = new Engine({ width: L, height: L, depth: L, T, seed, thermostatMode: 'kelvin' });
  let s = seed * 7919 + 17; const rng = () => (s = (s * 16807) % 2147483647) / 2147483647;
  const centers = [];
  for (const [mol, n] of Object.entries(mix)) for (let k = 0; k < n; k++) {
    let c, tries = 0;
    do { c = [3 + rng() * (L - 6), 3 + rng() * (L - 6), e.box.z0 + 3 + rng() * (L - 6)]; tries++; } while (centers.some(d => Math.hypot(d[0] - c[0], d[1] - c[1], d[2] - c[2]) < 4.2) && tries < 500);
    centers.push(c);
    const u = [rng() - .5, rng() - .5, rng() - .5], un = Math.hypot(...u) || 1, base = e.N;
    // a linear template is turned to a random axis; anything else is placed as drawn (whole, never mixed)
    const linear = TPL[mol].every(([, , y, z]) => !y && !z);
    for (const [sym, x, y, z] of TPL[mol]) linear ? e.addAtom(sym, c[0] + x * u[0] / un, c[1] + x * u[1] / un, c[2] + x * u[2] / un, { thermal: false }) : e.addAtom(sym, c[0] + x, c[1] + y, c[2] + z, { thermal: false });
    for (const [i, j, o] of ORD[mol] || []) e.setBondOrder(base + i, base + j, o);
  }
  e.touch(); e.refresh(); e.minimize(300, 1); e.refresh(); e.thermalize(T);
  for (let i = 0; i < ps * 1000; i++) e.step();
  const n = new Float64Array(e.N); for (const b of e.bonds()) { n[b.i] += b.order; n[b.j] += b.order; }
  const over = []; for (let k = 0; k < e.N; k++) { const el = e.formulaOf([k]); if (n[k] > VAL[el] + 0.6) over.push(el + k); }
  const c = {}; for (const f of e.fragments().list) { const k = e.formulaOf(f); c[k] = (c[k] || 0) + 1; }
  return { c, over };
}
const same = (c, mix) => Object.keys(c).length === Object.keys(mix).length && Object.entries(mix).every(([m, n]) => c[m] === n);

test('Mixtures that real chemistry keeps unreacted stay unreacted, whatever the temperature allows', () => {
  const cases = [
    ['H2 + O2 at room temperature (needs a spark)', { H2: 10, O2: 5 }, 300],
    ['H2 + Cl2 in the dark at room temperature (needs light)', { H2: 8, Cl2: 8 }, 300],
    ['CH4 + O2 at room temperature', { CH4: 4, O2: 8 }, 300],
    ['N2 + H2 at Haber temperature with no catalyst', { N2: 5, H2: 15 }, 700],
    ['N2 at 3000 K: N≡N does not fuse', { N2: 16 }, 3000],
    ['noble gases with H2 at 3000 K', { He: 6, Ar: 6, H2: 6 }, 3000],
    ['ethylene + H2 at 600 K with no catalyst', { C2H4: 6, H2: 6 }, 600],
    ['methane at 1000 K', { CH4: 8 }, 1000],
  ];
  const bad = [];
  for (const [label, mix, T] of cases) for (const seed of [3, 7]) {
    const r = run(mix, T, 60, seed);
    if (!same(r.c, mix) || r.over.length) bad.push(`${label} (seed ${seed}): ${JSON.stringify(r.c)}${r.over.length ? ' over valence ' + r.over.join(',') : ''}`);
  }
  assert.equal(bad.length, 0, bad.join('\n'));
});

test('Fast reactions that real chemistry runs in picoseconds happen here too', () => {
  // sodium burns in chlorine; methyl radicals recombine; an H atom strips Cl2 (k ≈ 2e-11 cm3/s at 300 K)
  const na = run({ Na: 8, Cl2: 4 }, 300, 60, 3);
  assert.ok((na.c.NaCl || 0) >= 5, 'Na + Cl2 → NaCl: ' + JSON.stringify(na.c));
  const me = run({ CH3: 8 }, 300, 60, 3);
  assert.ok((me.c.C2H6 || 0) >= 3, 'CH3 + CH3 → C2H6: ' + JSON.stringify(me.c));
  const h = run({ H: 3, Cl2: 30 }, 300, 60, 3);
  assert.ok((h.c.HCl || 0) >= 2, 'H + Cl2 → HCl + Cl: ' + JSON.stringify(h.c));
  const f = run({ F: 3, H2: 30 }, 1000, 60, 3);
  assert.ok((f.c.HF || 0) >= 2, 'F + H2 → HF + H at 1000 K: ' + JSON.stringify(f.c));
  for (const r of [na, me, h, f]) assert.equal(r.over.length, 0, 'an atom over its valence: ' + r.over.join(','));
});

console.log(count + ' condition checks passed.');
