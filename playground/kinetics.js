(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory(require('./engine.js'));
  else root.ChemKinetics = factory(root.ChemEngine);
})(typeof self !== 'undefined' ? self : this, function (CE) {
'use strict';
const { Engine, ELEMENTS, PAIR } = CE;
const NT = ELEMENTS.length;
const R = 0.008314462618, KB_SI = 1.380649e-23, H_SI = 6.62607015e-34, AMU = 1.66053906660e-27, LN2 = Math.log(2);
const RESTRAINT_K = 4000, BREAK_REACH = 3.2, START_GAP = 4.2;

function blank(T) {
  const e = new Engine({ width: 200, height: 200, depth: 200, T: T || 0, thermostat: false });
  e.box = { x0: -100, x1: 100, y0: -100, y1: 100, z0: -100, z1: 100 };
  return e;
}

function isolate(src, atoms, place, T) {
  const e = blank(T), map = new Map();
  const at = i => place && place.has(i) ? place.get(i) : [src.pos[3 * i], src.pos[3 * i + 1], src.pos[3 * i + 2]];
  let cx = 0, cy = 0, cz = 0;
  for (const i of atoms) { const p = at(i); cx += p[0]; cy += p[1]; cz += p[2]; }
  cx /= atoms.length; cy /= atoms.length; cz /= atoms.length;
  for (const i of atoms) {
    const p = at(i);
    map.set(i, e.addAtom(ELEMENTS[src.type[i]].sym, p[0] - cx, p[1] - cy, p[2] - cz, { thermal: false, charge: src.formal[i], V: src.val[i] }));
  }
  e.touch(); e.refresh();
  for (let p = 0; p < src.nPairs; p++) {
    const a = map.get(src.pI[p]), b = map.get(src.pJ[p]);
    if (a !== undefined && b !== undefined && src.pN[p] > 1) e.setBondOrder(a, b, src.pN[p]);
  }
  e.refresh();
  return { e, map };
}

function settle(e, iters) {
  e.minimize(iters, 0.05);
  for (let k = 0; k < 120; k++) e.computeForces(1);
  e.minimize(Math.max(200, iters >> 1), 0.05);
  for (let k = 0; k < 60; k++) e.computeForces(1);
  return e.computeForces() - (e.Erestraint || 0);
}

function dist(e, i, j) { return Math.hypot(e.pos[3 * j] - e.pos[3 * i], e.pos[3 * j + 1] - e.pos[3 * i + 1], e.pos[3 * j + 2] - e.pos[3 * i + 2]); }

function fragmentFormulas(e) {
  const f = e.fragments().list.map(g => e.formulaOf(g) + (e.isRadical(g) ? '·' : ''));
  f.sort();
  return f;
}

function relaxedEnergyOf(src, atoms, T) {
  const { e } = isolate(src, atoms, null, T);
  return { E: settle(e, 3000), formula: fragmentFormulas(e) };
}

function pairParams(src, i, j) { return PAIR[src.type[i] * NT + src.type[j]]; }

function bonded(src, i, j) {
  const a = Math.min(i, j), b = Math.max(i, j), p = src.pairMap.get(a * 1048576 + b);
  return p !== undefined && src.bondStrength(p) > 0.25;
}

function runChannel(src, atoms, i, j, E0, ch, opts) {
  const { e, map } = isolate(src, atoms, ch.place, opts.T);
  const a = map.get(i), b = map.get(j), k = ch.k === undefined ? undefined : map.get(ch.k);
  const pp = pairParams(src, i, j);
  const coord = () => k === undefined ? dist(e, a, b) : dist(e, a, b) - dist(e, b, k);
  const c0 = coord();
  let c1;
  if (ch.type === 'break') c1 = c0 + BREAK_REACH;
  else if (ch.type === 'form') c1 = pp.re[1];
  else c1 = pp.re[1] - (pairParams(src, j, ch.k).re[1] + BREAK_REACH);
  const steps = opts.steps || Math.max(14, Math.ceil(Math.abs(c1 - c0) / 0.08));
  const iters = opts.iters || 800, path = [];
  e.restraints = [k === undefined ? { i: a, j: b, r: c0, k: RESTRAINT_K } : { i: a, j: b, m: b, n: k, r: c0, k: RESTRAINT_K }];
  for (let s = 0; s <= steps; s++) {
    const target = c0 + (c1 - c0) * s / steps;
    e.restraints[0].r = target;
    const E = settle(e, iters);
    path.push({ x: coord(), r: dist(e, a, b), E: E - E0, pos: e.pos.slice(0, 3 * e.N) });
  }
  let top = 0;
  for (let s = 1; s < path.length; s++) if (path[s].E > path[top].E) top = s;
  let Eend;
  if (ch.type === 'break') Eend = path[path.length - 1].E;
  else { e.restraints = null; Eend = settle(e, 3000) - E0; }
  const nowBonded = bondedIn(e, a, b);
  const done = ch.type === 'break' ? !nowBonded : nowBonded && (k === undefined || !bondedIn(e, b, k));
  return {
    channel: ch, done, Ea: Math.max(0, path[top].E), dE: Eend, barrierAtEnd: top === path.length - 1,
    start: path[0].pos, end: e.pos.slice(0, 3 * e.N), E0,
    products: fragmentFormulas(e), path: path.map(p => ({ x: p.x, r: p.r, E: p.E })),
    ts: path[top].pos ? { pos: path[top].pos, map } : null, tsDistance: path[top].r,
    leaving: ch.k
  };
}

const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]], dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const unit = a => { const n = Math.hypot(a[0], a[1], a[2]); return n > 1e-9 ? [a[0] / n, a[1] / n, a[2] / n] : null; };
const P3 = (src, i) => [src.pos[3 * i], src.pos[3 * i + 1], src.pos[3 * i + 2]];

function neighbours(src, a) {
  const out = [];
  for (let p = 0; p < src.nPairs; p++) {
    if (src.bondStrength(p) <= 0.25) continue;
    if (src.pI[p] === a) out.push(src.pJ[p]); else if (src.pJ[p] === a) out.push(src.pI[p]);
  }
  return out;
}

function faces(src, a) {
  const nb = neighbours(src, a), here = P3(src, a), us = nb.map(b => unit(sub(P3(src, b), here))).filter(Boolean);
  const out = [];
  if (!us.length) return out;
  const sum = us.reduce((s, u) => [s[0] + u[0], s[1] + u[1], s[2] + u[2]], [0, 0, 0]);
  if (Math.hypot(...sum) > 0.3) out.push(unit(sum.map(x => -x)));
  if (us.length >= 2) { const n = unit(cross(us[0], us[1])); if (n) out.push(n, n.map(x => -x)); }
  return out;
}

function rotateOnto(a, b) {
  const v = cross(a, b), c = dot(a, b);
  if (c < -0.999999) { const t = unit(cross(a, Math.abs(a[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0])); return p => { const d = dot(t, p); return [2 * d * t[0] - p[0], 2 * d * t[1] - p[1], 2 * d * t[2] - p[2]]; }; }
  const k = 1 / (1 + c);
  return p => [
    (c + k * v[0] * v[0]) * p[0] + (k * v[0] * v[1] - v[2]) * p[1] + (k * v[0] * v[2] + v[1]) * p[2],
    (k * v[1] * v[0] + v[2]) * p[0] + (c + k * v[1] * v[1]) * p[1] + (k * v[1] * v[2] - v[0]) * p[2],
    (k * v[2] * v[0] - v[1]) * p[0] + (k * v[2] * v[1] + v[0]) * p[1] + (c + k * v[2] * v[2]) * p[2]];
}

function pose(src, frag, x, y, dir) {
  const target = P3(src, y).map((v, d) => v + dir[d] * START_GAP), px = P3(src, x);
  const own = faces(src, x)[0], R = own ? rotateOnto(own, dir.map(v => -v)) : (p => p);
  const place = new Map();
  for (const a of frag) { const r = R(sub(P3(src, a), px)); place.set(a, [target[0] + r[0], target[1] + r[1], target[2] + r[2]]); }
  return place;
}

function directions(list) {
  const out = [];
  for (const d of list) if (d && !out.some(o => dot(o, d) > 0.9)) out.push(d);
  return out;
}

function bondedIn(e, a, b) {
  const p = e.pairMap.get(Math.min(a, b) * 1048576 + Math.max(a, b));
  return p !== undefined && e.bondStrength(p) > 0.25;
}

function distinct(list, best) {
  const seen = new Set([best.products.join('+')]), out = [];
  for (const r of list) { const key = r.products.join('+'); if (!seen.has(key)) { seen.add(key); out.push(r); } }
  return out;
}

function scanPair(src, i, j, opts = {}) {
  if (src.needRebuild || src.needForces) { src._checkRebuild(); src.computeForces(); }
  const frag = src.fragments(), ci = frag.comp[i], cj = frag.comp[j];
  const atoms = ci === cj ? frag.list[ci].slice() : frag.list[ci].concat(frag.list[cj]);
  const pp = pairParams(src, i, j);
  if (i === j) return { ok: false, reason: 'pick two different atoms' };
  if (!pp || !pp.bond) return { ok: false, reason: ELEMENTS[src.type[i]].sym + ' and ' + ELEMENTS[src.type[j]].sym + ' do not form a covalent bond in this model' };
  const kind = bonded(src, i, j) ? 'break' : ci === cj ? 'close' : 'form';
  const reactants = ci === cj ? [relaxedEnergyOf(src, frag.list[ci], opts.T)] : [relaxedEnergyOf(src, frag.list[ci], opts.T), relaxedEnergyOf(src, frag.list[cj], opts.T)];
  const E0 = reactants.reduce((s, x) => s + x.E, 0);
  const channels = [];
  const attacker = x => frag.list[frag.comp[x]];
  if (kind === 'break' || kind === 'close') channels.push({ type: kind === 'break' ? 'break' : 'form' });
  else {
    for (const [x, y] of [[i, j], [j, i]]) {
      const user = unit(sub(P3(src, x), P3(src, y)));
      for (const d of directions([user, ...faces(src, y)])) channels.push({ type: 'form', from: x, at: y, place: pose(src, attacker(x), x, y, d) });
      for (const k of neighbours(src, y)) {
        if (k === x) continue;
        for (const d of directions([unit(sub(P3(src, y), P3(src, k))), user])) channels.push({ type: 'transfer', at: y, k, from: x, place: pose(src, attacker(x), x, y, d) });
      }
    }
  }
  const run = (ch, o) => { const [x, y] = ch.from !== undefined ? [ch.from, ch.at] : [i, j]; return runChannel(src, atoms, x, y, E0, ch, o); };
  const coarse = channels.length > 1 ? { steps: 10, iters: 300, T: opts.T } : opts;
  let results = [];
  for (const ch of channels) {
    results.push(run(ch, coarse));
    if (opts.onStep) opts.onStep(results.length / (channels.length + 2));
  }
  if (channels.length > 1) {
    const good = results.filter(r => r.done).sort((p, q) => p.Ea - q.Ea);
    const keep = good.filter(r => r.Ea <= good[0].Ea + 25).slice(0, 3);
    const refined = keep.map(r => run(r.channel, opts));
    results = refined.concat(results.filter(r => !keep.includes(r)));
  }
  const ok = results.filter(r => r.done).sort((p, q) => p.Ea - q.Ea);
  if (!ok.length) return { ok: false, reason: kind === 'break' ? 'the bond would not come apart' : 'no path forms this bond: the atoms are saturated and nothing they hold lets go', tried: results.length };
  const best = ok[0];
  return {
    ok: true, kind, i, j, atoms, Ea: best.Ea, dE: best.dE, barrierAtEnd: best.barrierAtEnd,
    reactants: reactants.flatMap(x => x.formula), products: best.products, path: best.path, ts: best.ts, tsDistance: best.tsDistance,
    channel: { type: best.channel.type, k: best.channel.k }, bimolecular: ci !== cj, T: opts.T,
    alternatives: distinct(ok.slice(1), best).map(r => ({ channel: { type: r.channel.type, k: r.channel.k }, Ea: r.Ea, dE: r.dE, products: r.products }))
  };
}

function collisionRate(src, fragA, fragB, T) {
  const mass = list => list.reduce((s, i) => s + src.mass[i], 0);
  const radius = list => {
    let cx = 0, cy = 0, cz = 0; for (const i of list) { cx += src.pos[3 * i]; cy += src.pos[3 * i + 1]; cz += src.pos[3 * i + 2]; }
    cx /= list.length; cy /= list.length; cz /= list.length;
    let r = 0; for (const i of list) r = Math.max(r, Math.hypot(src.pos[3 * i] - cx, src.pos[3 * i + 1] - cy, src.pos[3 * i + 2] - cz) + ELEMENTS[src.type[i]].rvdw * 0.5);
    return r;
  };
  const mA = mass(fragA), mB = mass(fragB), mu = mA * mB / (mA + mB) * AMU;
  const d = (radius(fragA) + radius(fragB)) * 1e-10;
  const v = Math.sqrt(8 * KB_SI * T / (Math.PI * mu));
  return Math.PI * d * d * v * 1e6;
}

function partnerDensity(src, fragB) {
  const f = src.fragments(), want = src.formulaOf(fragB);
  let n = 0; for (const g of f.list) if (src.formulaOf(g) === want) n++;
  const b = src.box, V = (b.x1 - b.x0) * (b.y1 - b.y0) * (b.z1 - b.z0) * 1e-24;
  return { n, perCm3: n / V };
}

const UNCERTAINTY = { Ea: 10, A: 10 };

function forecast(src, scan, opts = {}) {
  const T = opts.T ?? src.T, RT = R * T;
  if (!scan.ok || T <= 0) return { ok: false, reason: scan.reason || 'temperature is zero, so nothing activated ever happens' };
  let A, unit, per, note;
  if (!scan.bimolecular) {
    A = scan.kind === 'break' && scan.barrierAtEnd ? 10 ** 15.5 : KB_SI * T / H_SI;
    per = 1; unit = 's⁻¹';
    note = scan.kind === 'break' && scan.barrierAtEnd ? 'simple bond fission: loose transition state, A ≈ 10^15.5 s⁻¹' : 'unimolecular, tight transition state: A ≈ kT/h';
  } else {
    const f = src.fragments(), fa = f.list[f.comp[scan.i]], fb = f.list[f.comp[scan.j]];
    const steric = opts.steric ?? 0.1;
    A = steric * collisionRate(src, fa, fb, T);
    const dens = opts.partnerPerCm3 != null ? { n: null, perCm3: opts.partnerPerCm3 } : partnerDensity(src, fb);
    per = dens.perCm3; unit = 'cm³ s⁻¹';
    note = 'bimolecular: collision rate × steric factor ' + steric + ', partner density ' + dens.perCm3.toExponential(2) + ' cm⁻³ (' + (dens.n ?? '?') + ' in the chamber)';
  }
  const k = A * Math.exp(-scan.Ea / RT);
  const kLo = A / UNCERTAINTY.A * Math.exp(-(scan.Ea + UNCERTAINTY.Ea) / RT);
  const kHi = A * UNCERTAINTY.A * Math.exp(-Math.max(0, scan.Ea - UNCERTAINTY.Ea) / RT);
  const life = x => LN2 / (x * per);
  return { ok: true, T, Ea: scan.Ea, dE: scan.dE, A, unit, k, kLo, kHi, halfLife: life(k), fastest: life(kHi), slowest: life(kLo), note, bimolecular: scan.bimolecular };
}

function humanTime(s) {
  if (!isFinite(s)) return 'never';
  const units = [[1e-12, 'fs', 1e-15], [1e-9, 'ps', 1e-12], [1e-6, 'ns', 1e-9], [1e-3, 'µs', 1e-6], [1, 'ms', 1e-3], [60, 's', 1], [3600, 'min', 60], [86400, 'h', 3600], [86400 * 365.25, 'days', 86400]];
  for (const [lim, name, div] of units) if (s < lim) return fmt(s / div) + ' ' + name;
  const y = s / (86400 * 365.25);
  if (y < 1e4) return fmt(y) + ' years';
  if (y < 1.4e10) return '10^' + Math.log10(y).toFixed(0) + ' years';
  return 'longer than the age of the universe';
}
function fmt(x) { return x >= 100 ? x.toFixed(0) : x >= 10 ? x.toFixed(1) : x.toFixed(2); }

function fromScene(scene) {
  const e = new Engine({ T: scene.T || 0, thermostat: false });
  if (scene.box) e.box = { ...scene.box };
  for (const a of scene.atoms) e.addAtom(a[0], a[1], a[2], a[3], { thermal: false, charge: a[7] || 0, V: a[8] });
  e.touch(); e.refresh();
  for (const [i, j, o] of scene.bonds || []) e.setBondOrder(i, j, o);
  e.refresh();
  return e;
}

function study(scene, i, j, temps) {
  const e = fromScene(scene), T = scene.T || 298;
  const scan = scanPair(e, i, j, { T });
  if (!scan.ok) return scan;
  const at = t => { const f = forecast(e, scan, { T: t }); return f.ok ? { T: t, halfLife: f.halfLife, fastest: f.fastest, slowest: f.slowest } : { T: t, halfLife: Infinity }; };
  const now = forecast(e, scan, { T });
  return {
    ok: true, kind: scan.kind, channel: scan.channel, Ea: scan.Ea, dE: scan.dE, reactants: scan.reactants, products: scan.products,
    path: scan.path.map(p => p.E), bimolecular: scan.bimolecular, alternatives: scan.alternatives,
    now, table: (temps || [200, 298, 500, 1000, 1500, 2500]).map(at), T
  };
}

return { scanPair, forecast, humanTime, isolate, settle, fromScene, study, UNCERTAINTY };
});
