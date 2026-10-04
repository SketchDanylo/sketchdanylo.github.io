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
  for (let n = 0; n < 3 * e.N; n++) e.pos[n] += 0.03 * Math.sin(12.9898 * n + 4.1414);
  e.touch(); e.refresh();
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
  const ka = ch.ka === undefined ? undefined : map.get(ch.ka);
  const pp = pairParams(src, i, j);
  const terms = ch.type === 'swap' ? [[a, b, 1], [a, ka, -0.5], [b, k, -0.5]] : k === undefined ? [[a, b, 1]] : [[a, b, 1], [b, k, -1]];
  const coord = () => terms.reduce((s, [x, y, w]) => s + w * dist(e, x, y), 0);
  const c0 = coord();
  let c1;
  if (ch.type === 'break') c1 = c0 + BREAK_REACH;
  else if (ch.type === 'form') c1 = pp.re[1];
  else if (ch.type === 'swap') c1 = pp.re[1] - 0.5 * (pairParams(src, i, ch.ka).re[1] + pairParams(src, j, ch.k).re[1] + 2 * BREAK_REACH);
  else c1 = pp.re[1] - (pairParams(src, j, ch.k).re[1] + BREAK_REACH);
  const steps = opts.steps || Math.max(14, Math.ceil(Math.abs(c1 - c0) / 0.08));
  const iters = opts.iters || 800, path = [];
  const restrain = r => { e.restraints = [{ terms, r, k: RESTRAINT_K }]; };
  const point = s => {
    restrain(c0 + (c1 - c0) * s / steps);
    const E = settle(e, iters);
    const bo = []; for (let p = 0; p < e.nPairs; p++) if (e.pN[p] > 1.001) bo.push([e.pI[p], e.pJ[p], e.pN[p]]);
    return { x: coord(), r: dist(e, a, b), E: E - E0, pos: e.pos.slice(0, 3 * e.N), bo, bonds: bondsOf(e) };
  };
  for (let s = 0; s <= steps; s++) path.push(point(s));
  let Eend, done, products, endPos;
  if (ch.type === 'break') {
    Eend = path[path.length - 1].E;
    done = !bondedIn(e, a, b) && freeValence(e).every(x => x > -0.5);
    products = fragmentFormulas(e); endPos = e.pos.slice(0, 3 * e.N);
  } else {
    e.restraints = null; Eend = settle(e, 3000) - E0;
    const sane = freeValence(e).every(x => x > -0.5);
    done = sane && bondedIn(e, a, b) && (k === undefined || !bondedIn(e, b, k)) && (ka === undefined || !bondedIn(e, a, ka));
    products = fragmentFormulas(e); endPos = e.pos.slice(0, 3 * e.N);
    if (done) {
      const parts = e.fragments().list;
      if (parts.length > 1) {
        const apart = parts.reduce((sum, g) => sum + relaxedEnergyOf(e, g, opts.T).E, 0) - E0;
        if (Eend < apart - 60) done = false;
        else Eend = apart;
      }
    }
    if (done && !opts.oneWay) {
      const known = new Set([...path[0].bonds, ...bondsOf(e)]);
      for (let s = steps; s >= 0; s--) { const q = point(s); if (q.E < path[s].E && [...q.bonds].every(x => known.has(x))) path[s] = q; }
    }
  }
  e.restraints = null;
  let top = 0;
  for (let s = 1; s < path.length; s++) if (path[s].E > path[top].E) top = s;
  return {
    channel: ch, done, Ea: Math.max(0, path[top].E), dE: Eend, barrierAtEnd: top === path.length - 1,
    event: committed(path, top, j),
    start: path[0].pos, end: endPos, E0,
    products, path: path.map(p => ({ x: p.x, r: p.r, E: p.E })),
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

function committed(path, top, at) {
  const last = path.length - 1, drop = Math.max(8, 0.3 * (path[top].E - path[last].E));
  let c = top + 1;
  while (c < last && path[c].E > path[top].E - drop) c++;
  c = Math.min(c, last);
  const from = path[Math.max(0, c - 1)].pos, to = path[c].pos;
  return { pos: Array.from(to), push: Array.from(to, (v, q) => v - from[q]), bo: path[c].bo, at };
}

function bondedIn(e, a, b) {
  const p = e.pairMap.get(Math.min(a, b) * 1048576 + Math.max(a, b));
  return p !== undefined && e.bondStrength(p) > 0.25;
}

function bondsOf(e) {
  const out = new Set();
  for (let p = 0; p < e.nPairs; p++) if (e.bondStrength(p) > 0.25) out.add(Math.min(e.pI[p], e.pJ[p]) * 1048576 + Math.max(e.pI[p], e.pJ[p]));
  return out;
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
  const relaxed = list => {
    if (!opts.cache) return relaxedEnergyOf(src, list, opts.T);
    const key = list.join(',');
    if (!opts.cache.has(key)) opts.cache.set(key, relaxedEnergyOf(src, list, opts.T));
    return opts.cache.get(key);
  };
  const reactants = ci === cj ? [relaxed(frag.list[ci])] : [relaxed(frag.list[ci]), relaxed(frag.list[cj])];
  const free = freeValence(src);
  const E0 = reactants.reduce((s, x) => s + x.E, 0);
  const channels = [];
  const attacker = x => frag.list[frag.comp[x]];
  if (kind === 'break') channels.push({ type: 'break' });
  else if (kind === 'close') {
    channels.push({ type: 'form' });
    if (free[i] <= 0.5 && free[j] <= 0.5) for (const ka of neighbours(src, i)) for (const kb of neighbours(src, j)) if (ka !== j && kb !== i) channels.push({ type: 'swap', from: i, at: j, ka, k: kb });
  } else {
    for (const [x, y] of [[i, j], [j, i]]) {
      const user = unit(sub(P3(src, x), P3(src, y)));
      for (const d of directions([user, ...faces(src, y)])) channels.push({ type: 'form', from: x, at: y, place: pose(src, attacker(x), x, y, d) });
      if (free[y] > 0.5) continue;
      for (const k of neighbours(src, y)) {
        if (k === x) continue;
        for (const d of directions([unit(sub(P3(src, y), P3(src, k))), user])) channels.push({ type: 'transfer', at: y, k, from: x, place: pose(src, attacker(x), x, y, d) });
      }
    }
  }
  const run = (ch, o) => { const [x, y] = ch.from !== undefined ? [ch.from, ch.at] : [i, j]; return runChannel(src, atoms, x, y, E0, ch, o); };
  const coarse = channels.length > 1 || opts.quick ? { steps: channels.length > 1 ? 10 : 16, iters: 300, T: opts.T } : opts;
  let results = [];
  for (const ch of channels) {
    results.push(run(ch, coarse));
    if (opts.onStep) opts.onStep(results.length / (channels.length + 2));
  }
  if (channels.length > 1 && !opts.quick) {
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
    event: { atoms, pos: best.event.pos, push: best.event.push, bo: best.event.bo, anchor: ci !== cj ? frag.list[frag.comp[best.event.at]] : atoms },
    alternatives: distinct(ok.slice(1), best).map(r => ({ channel: { type: r.channel.type, k: r.channel.k }, Ea: r.Ea, dE: r.dE, products: r.products, barrierAtEnd: r.barrierAtEnd, event: { atoms, pos: r.event.pos, push: r.event.push, bo: r.event.bo, anchor: ci !== cj ? frag.list[frag.comp[r.event.at]] : atoms } }))
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
    path: scan.path.map(p => p.E), bimolecular: scan.bimolecular, alternatives: scan.alternatives, event: scan.event,
    now, table: (temps || [200, 298, 500, 1000, 1500, 2500]).map(at), T
  };
}

function siteKey(src, frag, a) {
  const nb = neighbours(src, a).map(b => ELEMENTS[src.type[b]].sym).sort().join('');
  return src.formulaOf(frag) + ':' + ELEMENTS[src.type[a]].sym + '(' + nb + ')';
}

function freeValence(src) {
  const used = new Float64Array(src.N);
  for (let p = 0; p < src.nPairs; p++) if (src.bondStrength(p) > 0.25) { used[src.pI[p]] += src.pN[p]; used[src.pJ[p]] += src.pN[p]; }
  return Array.from(used, (u, i) => src.val[i] - u);
}

function candidates(src, cap = 24) {
  const frag = src.fragments(), free = freeValence(src), out = [], seen = new Map(), pairSeen = new Set();
  const pi = new Set();
  for (let p = 0; p < src.nPairs; p++) if (src.bondStrength(p) > 0.25 && src.pN[p] > 1.2) { pi.add(src.pI[p]); pi.add(src.pJ[p]); }
  const add = (i, j, kind, key, prio) => {
    const d = Math.hypot(src.pos[3 * j] - src.pos[3 * i], src.pos[3 * j + 1] - src.pos[3 * i + 1], src.pos[3 * j + 2] - src.pos[3 * i + 2]);
    const have = seen.get(key);
    if (have) { have.mult++; if (d < have.d) { have.i = i; have.j = j; have.d = d; } return; }
    const c = { i, j, kind, key, prio, mult: 1, d };
    seen.set(key, c); out.push(c);
  };
  const L = frag.list, open = L.map(f => f.some(x => free[x] > 0.5));
  for (const first of [true, false]) for (let A = 0; A < L.length; A++) for (let B = 0; B < L.length; B++) {
    if (A === B || (!first && open[B])) continue;
    for (const r of L[A]) {
      const radical = free[r] > 0.5, pib = pi.has(r);
      if (radical !== first || (!radical && !pib)) continue;
      for (const t of L[B]) {
        const pp = pairParams(src, r, t);
        if (!pp || !pp.bond) continue;
        if (!radical && ELEMENTS[src.type[t]].sym !== 'H') continue;
        const pk = Math.min(r, t) + ',' + Math.max(r, t);
        if (pairSeen.has(pk)) continue;
        pairSeen.add(pk);
        const ka = siteKey(src, L[A], r), kb = siteKey(src, L[B], t);
        add(r, t, 'meet', ka < kb ? ka + '|' + kb : kb + '|' + ka, radical ? 0 : 1);
      }
    }
  }
  for (let p = 0; p < src.nPairs; p++) {
    if (src.bondStrength(p) <= 0.25) continue;
    const i = src.pI[p], j = src.pJ[p], f = L[frag.comp[i]];
    const ki = siteKey(src, f, i), kj = siteKey(src, f, j);
    add(i, j, 'break', 'break ' + (ki < kj ? ki + '–' + kj : kj + '–' + ki) + ' ' + src.pN[p].toFixed(1), 2);
  }
  out.sort((x, y) => x.prio - y.prio || x.d - y.d);
  return { list: out.slice(0, cap), dropped: Math.max(0, out.length - cap) };
}

function chamberVolume(e) { const b = e.box; return (b.x1 - b.x0) * (b.y1 - b.y0) * (b.z1 - b.z0) * 1e-24; }

function scanCandidate(e, c, T, quick, cache) {
  const scan = scanPair(e, c.i, c.j, { T, quick, cache });
  if (!scan.ok) return null;
  const f = forecast(e, scan, { T, partnerPerCm3: 1 });
  if (!f.ok) return null;
  const V = chamberVolume(e), one = (s, kRate) => ({ label: scan.reactants.join(' + ') + ' → ' + s.products.join(' + '), i: c.i, j: c.j, Ea: s.Ea, dE: s.dE, rate: scan.bimolecular ? kRate * c.mult / V : kRate * c.mult, kRate, bimolecular: scan.bimolecular, refined: !quick, reactants: scan.reactants, products: s.products, event: s.event, kind: scan.kind });
  const main = one(scan, f.k);
  main.also = (scan.alternatives || []).filter(a => !(bare(a.products) === bare(scan.reactants))).map(a => {
    const g = forecast(e, { ok: true, Ea: a.Ea, kind: a.channel.type === 'break' ? 'break' : scan.kind, barrierAtEnd: a.barrierAtEnd, bimolecular: scan.bimolecular, i: scan.i, j: scan.j }, { T, partnerPerCm3: 1 });
    return g.ok ? one(a, g.k) : null;
  }).filter(Boolean);
  return main;
}

const bare = list => list.map(f => f.replace(/·/g, '')).sort().join(' + ');

function combine(results, T) {
  const events = [], seen = new Map(), flat = [];
  for (const r of results) if (r) { flat.push(r); for (const a of r.also || []) flat.push({ ...a, cand: r.cand, alt: true }); }
  for (const r of flat) {
    if (!r) continue;
    if (bare(r.reactants) === bare(r.products)) continue;
    const key = bare(r.reactants) + ' → ' + bare(r.products);
    const prev = seen.get(key);
    if (prev) { prev.rate += r.rate; if (r.refined && !prev.refined || r.refined === prev.refined && r.Ea < prev.Ea) Object.assign(prev, { Ea: r.Ea, event: r.event, refined: r.refined, cand: r.cand }); continue; }
    const x = { ...r }; seen.set(key, x); events.push(x);
  }
  const total = events.reduce((s, x) => s + x.rate, 0);
  for (const x of events) { x.share = total > 0 ? x.rate / total : 0; x.halfLife = LN2 / x.rate; }
  events.sort((x, y) => y.rate - x.rate);
  return { ok: true, T, events, total };
}

const REFINE = 3, COOL = [80, 160, 300];

function reuse(cached, c, V) {
  const one = x => ({ ...x, i: c.i, j: c.j, event: null, cand: c, rate: x.bimolecular ? x.kRate * c.mult / V : x.kRate * c.mult });
  return { ...one(cached), also: (cached.also || []).map(one) };
}

function survey(scene, opts = {}) {
  const e = fromScene(scene), T = scene.T || 298, cache = new Map();
  const { list, dropped } = candidates(e, opts.cap);
  const quick = list.map((c, n) => { if (opts.onProgress) opts.onProgress(n, list.length + REFINE); return { c, r: scanCandidate(e, c, T, true, cache) }; });
  const best = quick.filter(q => q.r).sort((x, y) => x.r.Ea - y.r.Ea).slice(0, REFINE);
  for (const q of best) q.r = scanCandidate(e, q.c, T, false, cache);
  const out = combine(quick.map(q => q.r), T);
  out.dropped = dropped; out.scanned = list.length;
  return out;
}

function pickNext(sv, u1 = Math.random(), u2 = Math.random()) {
  if (!sv.total || !isFinite(sv.total)) return null;
  const wait = -Math.log(1 - u1) / sv.total;
  let acc = 0, pick = sv.events[sv.events.length - 1];
  for (const x of sv.events) { acc += x.share; if (u2 < acc) { pick = x; break; } }
  return { wait, pick };
}

function stillThere(eng, ev, idx, reactants) {
  const f = eng.fragments(), comps = new Set(idx.map(i => f.comp[i]));
  return bare([...comps].map(c => eng.formulaOf(f.list[c]))) === bare(reactants);
}

function formed(eng, idx, products) {
  const f = eng.fragments(), comps = new Set(idx.map(i => f.comp[i]));
  return bare([...comps].map(c => eng.formulaOf(f.list[c]))) === bare(products) ? [...comps].flatMap(c => f.list[c]) : null;
}

function verifyEvent(scene, ev, products, T, seeds = 3, fs = 400) {
  let ok = 0;
  for (let s = 1; s <= seeds; s++) {
    const e = fromScene(scene);
    e.T = T; e.thermostat = true; e.thermostatMode = 'kelvin'; e.rngState = s * 7919 + 13;
    applyEvent(e, ev, ev.atoms);
    const cool = COOL.slice();
    for (let k = 1; k <= fs; k++) {
      e.step();
      if (cool.length && k >= cool[0]) { cool.shift(); const f0 = e.fragments(), cs = new Set(ev.atoms.map(i => f0.comp[i])); e.thermalize(T, [...cs].flatMap(c => f0.list[c])); }
    }
    if (formed(e, ev.atoms, products)) ok++;
  }
  return ok;
}

function applyEvent(eng, ev, idx) {
  const anchor = new Set(ev.anchor);
  let cx = 0, cy = 0, cz = 0, ix = 0, iy = 0, iz = 0, n = 0;
  ev.atoms.forEach((k, m) => {
    if (!anchor.has(k)) return;
    const i = idx[m];
    cx += eng.pos[3 * i]; cy += eng.pos[3 * i + 1]; cz += eng.pos[3 * i + 2];
    ix += ev.pos[3 * m]; iy += ev.pos[3 * m + 1]; iz += ev.pos[3 * m + 2]; n++;
  });
  const dx = (cx - ix) / n, dy = (cy - iy) / n, dz = (cz - iz) / n;
  ev.atoms.forEach((k, m) => {
    const i = idx[m];
    eng.pos[3 * i] = ev.pos[3 * m] + dx; eng.pos[3 * i + 1] = ev.pos[3 * m + 1] + dy; eng.pos[3 * i + 2] = ev.pos[3 * m + 2] + dz;
    eng.born[i] = eng.time;
  });
  eng.touch(); eng.refresh();
  for (const [a, b, o] of ev.bo) eng.setBondOrder(idx[a], idx[b], o);
  eng.thermalize(eng.T, idx);
  let mD = 0;
  idx.forEach((i, m) => { mD += eng.mass[i] * (ev.push[3 * m] ** 2 + ev.push[3 * m + 1] ** 2 + ev.push[3 * m + 2] ** 2); });
  if (mD > 0) {
    const sc = Math.sqrt(2 * CE.KB * Math.max(eng.T, 50) / (CE.KEU * mD));
    idx.forEach((i, m) => { for (let d = 0; d < 3; d++) eng.vel[3 * i + d] += sc * ev.push[3 * m + d]; });
  }
  eng.refresh();
}

return { scanPair, forecast, humanTime, isolate, settle, fromScene, study, applyEvent, verifyEvent, stillThere, candidates, scanCandidate, combine, survey, pickNext, reuse, REFINE, COOL, UNCERTAINTY };
});
