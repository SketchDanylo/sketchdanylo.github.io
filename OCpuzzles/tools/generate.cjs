/* Chem Puzzles generator.
 *   node OCpuzzles/tools/generate.cjs
 * Builds every puzzle offline and writes them into OCpuzzles/puzzles.js.
 * Names come from the site's own naming engine (native-namer.js) and are
 * translated to German; each name is parsed back by an independent German parser
 * and must give the same RDKit canonical structure. Drawings, R/S and E/Z come
 * from RDKit (vendored in the naming-engine branch; override with RDKIT_JS).
 * Task text is German, theme names are English (UI).
 */
'use strict';
const fs = require('fs'), path = require('path');
const C = require('./chem.cjs');
const HAND = require('./handmade.cjs');

/* ---------- deterministic randomness ---------- */
function mulberry(a) { return () => { a |= 0; a = a + 0x6D2B79F5 | 0; let t = Math.imul(a ^ a >>> 15, 1 | a); t = t + Math.imul(t ^ t >>> 7, 61 | t) ^ t; return ((t ^ t >>> 14) >>> 0) / 4294967296; }; }
const rnd = mulberry(20261126);
const ri = (a, b) => a + Math.floor(rnd() * (b - a + 1));
const pick = a => a[Math.floor(rnd() * a.length)];
const chance = p => rnd() < p;
const shuffle = a => { a = a.slice(); for (let i = a.length - 1; i > 0; i--) { const j = Math.floor(rnd() * (i + 1)); [a[i], a[j]] = [a[j], a[i]]; } return a; };
const clampC = c => Math.max(100, Math.min(3000, Math.round(c / 10) * 10));
const jit = (c, d = 60) => c + ri(-d, d);

/* ---------- molecule table and puzzle list ---------- */
const M = [], MI = new Map();
function mol(sm, opts) {
  const key = sm + (opts ? '|H' : '');
  if (MI.has(key)) return MI.get(key);
  const d = opts ? depictH(sm) : C.depict(sm);
  if (!d) { MI.set(key, -1); return -1; }
  M.push([d.a.map(a => [a[0], Math.round(a[1] * 100), Math.round(a[2] * 100)].concat(a[3] ? [a[3]] : [])), d.b]);
  MI.set(key, M.length - 1); return M.length - 1;
}
function depictH(sm) { // keep explicit [H]
  const RD = C.RD, m = RD.get_mol(sm, JSON.stringify({ removeHs: false }));
  m.set_new_coords(true); const block = m.get_molblock(); m.delete();
  const L = block.split('\n'); const na = parseInt(L[3].slice(0, 3)), nb = parseInt(L[3].slice(3, 6));
  const a = [], b = [];
  for (let i = 0; i < na; i++) { const s = L[4 + i]; a.push([s.slice(31, 34).trim(), +s.slice(0, 10), -(+s.slice(10, 20))]); }
  for (let i = 0; i < nb; i++) { const s = L[4 + na + i]; b.push([+s.slice(0, 3) - 1, +s.slice(3, 6) - 1, +s.slice(6, 9), +s.slice(9, 12)]); }
  let sum = 0; b.forEach(x => sum += Math.hypot(a[x[0]][1] - a[x[1]][1], a[x[0]][2] - a[x[1]][2]));
  const k = b.length / sum, xs = a.map(x => x[1]), ys = a.map(x => x[2]);
  const cx = (Math.min(...xs) + Math.max(...xs)) / 2, cy = (Math.min(...ys) + Math.max(...ys)) / 2;
  a.forEach(x => { x[1] = Math.round((x[1] - cx) * k * 100) / 100; x[2] = Math.round((x[2] - cy) * k * 100) / 100; });
  for (let i = 0; i < a.length; i++) for (let j = i + 1; j < a.length; j++) if (Math.hypot(a[i][1] - a[j][1], a[i][2] - a[j][2]) < 0.5) return null;
  return { a, b: b.map(x => x[3] ? x : x.slice(0, 3)) };
}

const P = [], KEYS = new Set(), COUNT = {};
function put(fam, p) {
  const key = fam + '|' + p.q + '|' + JSON.stringify([p.m, p.o, p.rg, p.m2, p.items, p.tm]);
  if (KEYS.has(key)) return false;
  KEYS.add(key); COUNT[fam] = (COUNT[fam] || 0) + 1;
  p.i = fam + '-' + COUNT[fam]; p.c = clampC(p.c); P.push(p); return true;
}
function mc(fam, p, correct, wrongs) {
  const ws = [...new Set(wrongs.filter(w => w != null && w !== correct))].slice(0, 3);
  if (ws.length < 2) return false;
  const o = shuffle([correct, ...ws]);
  return put(fam, Object.assign(p, { k: 'mc', o, a: o.indexOf(correct) }));
}

/* ---------- graph helpers ---------- */
const clone = g => g.map(a => ({ el: a.el, ch: a.ch || 0, nb: a.nb.map(e => ({ n: e.n, o: e.o })) }));
function removeAtom(g, i) {
  const out = [];
  g.forEach((a, k) => { if (k !== i) out.push({ el: a.el, ch: a.ch || 0, nb: a.nb.filter(e => e.n !== i).map(e => ({ n: e.n > i ? e.n - 1 : e.n, o: e.o })) }); });
  return out;
}
const valOk = g => g.every((a, i) => a.nb.reduce((s, e) => s + e.o, 0) <= (C.VAL[a.el] + (a.ch > 0 && a.el === 'N' ? 1 : 0) - (a.ch < 0 && a.el === 'O' ? 1 : 0)));
const heavyC = (g, i) => g[i].nb.filter(e => g[e.n].el === 'C').length;
const FR = { Me: 'C', Et: 'CC', Pr: 'CCC', iPr: 'C(C)C', tBu: 'C(C)(C)C', F: 'F', Cl: 'Cl', Br: 'Br', I: 'I', OH: 'O', NH2: 'N', OMe: 'OC', NO2: '[N+](=O)[O-]', COOH: 'C(=O)O', CHO: 'C=O', CN: 'C#N', Ac: 'C(C)=O' };

function build(spec) {
  const g = C.G(); let A;
  if (spec.ring) {
    A = []; for (let i = 0; i < spec.ring; i++) A.push(C.add(g, 'C'));
    for (let i = 0; i < spec.ring; i++) C.bond(g, A[i], A[(i + 1) % spec.ring], spec.arom && i % 2 === 0 ? 2 : 1);
  } else A = C.chain(g, spec.n);
  for (const [p, o] of spec.uns || []) { if (p < 1 || p >= A.length) return null; C.setOrder(g, A[p - 1], A[p], o); }
  for (const [p, t] of spec.fg || []) {
    const c = A[p - 1]; if (c === undefined) return null;
    if (t === 'ol') C.bond(g, c, C.add(g, 'O'));
    else if (t === 'nh2') C.bond(g, c, C.add(g, 'N'));
    else if (t === 'on' || t === 'al') C.bond(g, c, C.add(g, 'O'), 2);
    else if (t === 'cooh') { C.bond(g, c, C.add(g, 'O'), 2); C.bond(g, c, C.add(g, 'O')); }
    else if (t === 'cn') C.bond(g, c, C.add(g, 'N'), 3);
  }
  for (const [p, f] of spec.subs || []) { const c = A[p - 1]; if (c === undefined) return null; C.attach(g, c, FR[f]); }
  return valOk(g) ? { g, A } : null;
}

/* ---------- naming with independent check ---------- */
const NAME = new Map(), FAILS = [];
function nm(g) {
  let sm; try { sm = C.smiles(g); } catch (e) { return null; }
  if (NAME.has(sm)) return NAME.get(sm);
  let n = null;
  try { n = C.nameDE(g); } catch (e) { n = null; }
  if (n) { const err = C.roundTrip(n, sm); if (err) { FAILS.push(err); n = null; } }
  NAME.set(sm, n); return n;
}
const can = sm => { try { return C.canon(sm); } catch (e) { return null; } };
const smOf = g => C.smiles(g);
const fmt = f => f.replace(/(\d+)/g, d => d.split('').map(x => '₀₁₂₃₄₅₆₇₈₉'[x]).join(''));

/* name tokenizer for distractors */
const PRE_RE = /^(?:(\d+(?:,\d+)*)-)?(di|tri|tetra)?(\(propan-2-yl\)|propan-2-yl|methyl|ethyl|propyl|butyl|fluor|chlor|brom|iod|hydroxy|amino|oxo|methoxy|ethoxy|cyano|nitro|phenyl)/i;
function splitName(n) {
  let s = n, pre = [], m;
  while ((m = s.match(PRE_RE))) { pre.push({ locs: m[1] ? m[1].split(',').map(Number) : null, mult: m[2] || '', sub: m[3].toLowerCase() }); s = s.slice(m[0].length).replace(/^-(?=\d)/, ''); if (!s) return null; }
  return { pre, parent: s };
}
const cap = s => s.replace(/^((?:[\d,]+-)?)([a-zäöü])/, (m, p, c) => p + c.toUpperCase());
function joinName(pre, parent) {
  let out = '';
  pre.forEach((p, i) => { const t = (p.locs ? p.locs.join(',') + '-' : '') + p.mult + p.sub; out += (i && /^\d/.test(t) ? '-' : '') + t; });
  out += (pre.length && /^\d/.test(parent) ? '-' : '') + (pre.length ? parent.toLowerCase() : parent);
  return cap(out);
}
const STEMLEN = { meth: 1, eth: 2, prop: 3, but: 4, pent: 5, hex: 6, hept: 7, oct: 8, non: 9, dec: 10 };
function reversedNumbering(n) {
  const sp = splitName(n); if (!sp) return null;
  const par = sp.parent.toLowerCase();
  if (/benzol|phenol|anilin|benzoe|benzal|cyclo|oat$|säure$|al$|nitril$|amid$/.test(par)) return null;
  const st = Object.keys(STEMLEN).filter(k => par.startsWith(k)).sort((a, b) => b.length - a.length)[0]; if (!st) return null;
  const L = STEMLEN[st];
  const rev = x => L + 1 - x;
  const pre = sp.pre.map(p => ({ ...p, locs: p.locs ? p.locs.map(rev).sort((a, b) => a - b) : null }));
  let parent = sp.parent.replace(/-(\d+(?:,\d+)*)-((?:di|tri)?(?:en|in))/g, (m, l, t) => '-' + l.split(',').map(x => L - +x).sort((a, b) => a - b).join(',') + '-' + t)
    .replace(/-(\d+(?:,\d+)*)-((?:di|tri)?(?:ol|on|amin))$/, (m, l, t) => '-' + l.split(',').map(x => rev(+x)).sort((a, b) => a - b).join(',') + '-' + t);
  const out = joinName(pre, parent);
  return out === n ? null : out;
}
function nonAlphabetical(n) {
  const sp = splitName(n); if (!sp || sp.pre.length < 2) return null;
  const out = joinName(sp.pre.slice().reverse(), sp.parent);
  return out === n ? null : out;
}

/* ---------- random structures ---------- */
function randSpec(L) {
  if (L >= 2 && chance(0.22)) return randRing(L);
  const n = ri(3, Math.min(8, 4 + L));
  const s = { n, uns: [], fg: [], subs: [] };
  const f = pick(L === 0 ? ['', '', 'ol', 'en'] : ['', 'ol', 'ol', 'on', 'al', 'cooh', 'nh2', 'cn', 'en', 'en', 'in']);
  if (f === 'ol' || f === 'nh2') s.fg.push([ri(1, n), f]);
  else if (f === 'on') s.fg.push([ri(2, n - 1), 'on']);
  else if (['al', 'cooh', 'cn'].includes(f)) s.fg.push([1, f]);
  else if (f === 'en') s.uns.push([ri(1, n - 1), 2]);
  else if (f === 'in') s.uns.push([ri(1, n - 1), 3]);
  if (L >= 3 && chance(0.35)) { if (f !== 'en' && f !== 'in') s.uns.push([ri(1, n - 1), 2]); else s.fg.push([ri(1, n), 'ol']); }
  const ns = L === 0 ? ri(0, 1) : ri(L >= 2 ? 1 : 0, Math.min(3, L + 1));
  const pool = L <= 1 ? ['Me', 'Me', 'Et', 'Cl', 'Br'] : ['Me', 'Me', 'Me', 'Et', 'Et', 'Pr', 'iPr', 'F', 'Cl', 'Br', 'I'];
  for (let k = 0; k < ns; k++) { const t = pick(pool); s.subs.push([/^[A-Z][a-z]?$/.test(t) && t.length <= 2 && !['Me', 'Et', 'Pr'].includes(t) ? ri(1, n) : ri(2, Math.max(2, n - 1)), t]); }
  return s;
}
function randRing(L) {
  if (chance(0.55)) { // benzene
    const k = ri(1, L >= 3 ? 3 : 2), s = { ring: 6, arom: true, subs: [] }, used = new Set();
    for (let i = 0; i < k; i++) { let p = i === 0 ? 1 : ri(2, 6); while (used.has(p)) p = ri(2, 6); used.add(p); s.subs.push([p, pick(['Me', 'Et', 'Cl', 'Br', 'NO2', 'OH', 'NH2', 'COOH', 'CHO', 'OMe', 'CN', 'Ac'])]); }
    return s;
  }
  const r = pick([5, 6, 6]), s = { ring: r, subs: [], fg: [] };
  if (chance(0.5)) s.fg.push([1, pick(['ol', 'on', 'nh2'])]);
  const k = ri(1, 2), used = new Set([1]);
  for (let i = 0; i < k; i++) { let p = s.fg.length ? ri(2, r) : (i === 0 ? 1 : ri(2, r)); if (used.has(p) && s.fg.length) continue; used.add(p); s.subs.push([p, pick(['Me', 'Me', 'Et', 'Cl', 'Br', 'iPr'])]); }
  return s;
}
function perturb(spec) {
  const s = JSON.parse(JSON.stringify(spec));
  const len = s.ring || s.n;
  const all = [...(s.subs || []).map(x => ['subs', x]), ...(s.fg || []).filter(x => !['al', 'cooh', 'cn'].includes(x[1])).map(x => ['fg', x]), ...(s.uns || []).map(x => ['uns', x])];
  if (!all.length) { if (!s.ring && chance(0.5)) { s.n = Math.max(2, s.n + pick([-1, 1])); return s; } return null; }
  const [kind, x] = pick(all);
  if (kind === 'uns') x[0] = ri(1, len - 1); else x[0] = ri(s.ring ? 1 : (kind === 'subs' && ['Me', 'Et', 'Pr', 'iPr'].includes(x[1]) ? 1 : 1), len);
  if (chance(0.25) && kind === 'subs' && ['Me', 'Et'].includes(x[1])) x[1] = x[1] === 'Me' ? 'Et' : 'Me';
  return s;
}
function scoreName(n, spec, g) {
  const locs = (n.match(/\d+/g) || []).length, sp = splitName(n);
  const kinds = sp ? new Set(sp.pre.map(p => p.sub)).size : 0;
  let s = 0.25 * g.length + 0.55 * locs + 0.9 * Math.max(0, kinds - 1) + (/(di|tri)/.test(n) ? 0.6 : 0);
  if (/-\d+-(?:di|tri)?(en|in)/.test(n) && /(ol|on|al|säure|nitril|amin)$/.test(n)) s += 1.2;
  if (sp && spec && !spec.ring) { const par = sp.parent.toLowerCase(); const st = Object.keys(STEMLEN).filter(k => par.startsWith(k)).sort((a, b) => b.length - a.length)[0]; if (st && STEMLEN[st] !== spec.n) s += 1.4; }
  if (spec && spec.ring && spec.arom) s += 0.4;
  if (spec && spec.ring && !spec.arom) s += 0.8;
  return s;
}

/* ======================================================================
   FAMILIES
   ====================================================================== */
const T = { el: 'Elements', nom: 'Nomenclature', str: 'Structure', fg: 'Functional groups', bond: 'Bonding', st: 'Stereochemistry', ab: 'Acids & bases', phys: 'Physical properties', mech: 'Mechanisms', rx: 'Reactions', syn: 'Synthesis' };

/* ---------- 1. elements ---------- */
const EL = [
  [1, 'H', 'Wasserstoff', 1, 1], [2, 'He', 'Helium', 8, 1], [3, 'Li', 'Lithium', 1, 1], [4, 'Be', 'Beryllium', 2, 2], [5, 'B', 'Bor', 3, 2], [6, 'C', 'Kohlenstoff', 4, 1], [7, 'N', 'Stickstoff', 5, 1], [8, 'O', 'Sauerstoff', 6, 1], [9, 'F', 'Fluor', 7, 1], [10, 'Ne', 'Neon', 8, 1],
  [11, 'Na', 'Natrium', 1, 1], [12, 'Mg', 'Magnesium', 2, 1], [13, 'Al', 'Aluminium', 3, 1], [14, 'Si', 'Silicium', 4, 2], [15, 'P', 'Phosphor', 5, 1], [16, 'S', 'Schwefel', 6, 1], [17, 'Cl', 'Chlor', 7, 1], [18, 'Ar', 'Argon', 8, 1], [19, 'K', 'Kalium', 1, 1], [20, 'Ca', 'Calcium', 2, 1],
  [21, 'Sc', 'Scandium', 0, 3], [22, 'Ti', 'Titan', 0, 2], [23, 'V', 'Vanadium', 0, 3], [24, 'Cr', 'Chrom', 0, 2], [25, 'Mn', 'Mangan', 0, 2], [26, 'Fe', 'Eisen', 0, 1], [27, 'Co', 'Cobalt', 0, 2], [28, 'Ni', 'Nickel', 0, 2], [29, 'Cu', 'Kupfer', 0, 1], [30, 'Zn', 'Zink', 0, 1],
  [31, 'Ga', 'Gallium', 3, 3], [32, 'Ge', 'Germanium', 4, 3], [33, 'As', 'Arsen', 5, 2], [34, 'Se', 'Selen', 6, 3], [35, 'Br', 'Brom', 7, 1], [36, 'Kr', 'Krypton', 8, 2], [37, 'Rb', 'Rubidium', 1, 3], [38, 'Sr', 'Strontium', 2, 3], [39, 'Y', 'Yttrium', 0, 3], [40, 'Zr', 'Zirconium', 0, 3],
  [41, 'Nb', 'Niob', 0, 3], [42, 'Mo', 'Molybdän', 0, 3], [44, 'Ru', 'Ruthenium', 0, 3], [45, 'Rh', 'Rhodium', 0, 3], [46, 'Pd', 'Palladium', 0, 2], [47, 'Ag', 'Silber', 0, 1], [48, 'Cd', 'Cadmium', 0, 3], [49, 'In', 'Indium', 3, 3], [50, 'Sn', 'Zinn', 4, 2], [51, 'Sb', 'Antimon', 5, 3],
  [52, 'Te', 'Tellur', 6, 3], [53, 'I', 'Iod', 7, 1], [54, 'Xe', 'Xenon', 8, 2], [55, 'Cs', 'Caesium', 1, 3], [56, 'Ba', 'Barium', 2, 2], [74, 'W', 'Wolfram', 0, 2], [78, 'Pt', 'Platin', 0, 2], [79, 'Au', 'Gold', 0, 1], [80, 'Hg', 'Quecksilber', 0, 2], [82, 'Pb', 'Blei', 4, 1],
  [83, 'Bi', 'Bismut', 5, 3], [86, 'Rn', 'Radon', 8, 3], [92, 'U', 'Uran', 0, 2],
];
const ROMAN = ['', 'I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII'];
function famElements() {
  const base = { 1: 100, 2: 230, 3: 420 };
  const near = (e, key) => shuffle(EL.filter(x => x !== e)).sort((a, b) => (b[1][0] === e[1][0]) - (a[1][0] === e[1][0]) || Math.abs(a[0] - e[0]) - Math.abs(b[0] - e[0])).slice(0, 6).map(x => x[key]);
  for (const e of EL) {
    mc('el', { c: jit(base[e[4]], 40), t: T.el, q: `Welches Element hat das Symbol <b>${e[1]}</b>?` }, e[2], shuffle(near(e, 2)).slice(0, 3));
    mc('el', { c: jit(base[e[4]] + 60, 40), t: T.el, q: `Welches Symbol hat das Element <b>${e[2]}</b>?` }, e[1], shuffle(near(e, 1)).slice(0, 3));
    if (e[0] <= 36 && e[3]) mc('el', { c: jit(e[0] <= 20 ? 330 : 470, 40), t: T.el, q: `In welcher Hauptgruppe steht <b>${e[2]}</b>?` }, ROMAN[e[3]], shuffle(ROMAN.slice(1).filter(r => r !== ROMAN[e[3]])).slice(0, 3));
    if (e[0] <= 20) mc('el', { c: jit(220 + e[0] * 6, 30), t: T.el, q: `Welches Element hat die Ordnungszahl <b>${e[0]}</b>?` }, e[2], EL.filter(x => x[0] !== e[0] && Math.abs(x[0] - e[0]) <= 3).map(x => x[2]));
    if (e[0] <= 20 && e[3]) mc('el', { c: jit(e[0] <= 10 ? 380 : 450, 40), t: T.el, q: `Wie viele Valenzelektronen hat ein <b>${e[2]}</b>-Atom?` }, String(e[3]), shuffle(['1', '2', '3', '4', '5', '6', '7', '8'].filter(x => x !== String(e[3]))).slice(0, 3));
  }
}

/* ---------- 2. homologous series basics ---------- */
const STEM = ['', 'Meth', 'Eth', 'Prop', 'But', 'Pent', 'Hex', 'Hept', 'Oct', 'Non', 'Dec'];
function famBasics() {
  for (let n = 1; n <= 10; n++) {
    const others = shuffle([1, 2, 3, 4, 5, 6, 7, 8, 9, 10].filter(k => k !== n && Math.abs(k - n) <= 3)).slice(0, 3);
    mc('ba', { c: jit(110 + n * 18, 20), t: T.nom, q: `Wie heißt das unverzweigte Alkan mit <b>${n}</b> C-Atom${n > 1 ? 'en' : ''}?` }, STEM[n] + 'an', others.map(k => STEM[k] + 'an'));
    mc('ba', { c: jit(160 + n * 15, 20), t: T.str, q: `Welche Summenformel hat <b>${STEM[n]}an</b>?`, f: 1 }, fmt(`C${n}H${2 * n + 2}`), [fmt(`C${n}H${2 * n}`), fmt(`C${n}H${2 * n + 1}`), fmt(`C${n + 1}H${2 * n + 2}`), n > 1 ? fmt(`C${n}H${2 * n - 2}`) : null]);
    const g = C.chain(C.G(), 0); // placeholder
    const sm = 'C'.repeat(n); const mi = mol(sm);
    if (n >= 2 && mi >= 0) mc('ba', { c: jit(120 + n * 20, 20), t: T.nom, q: 'Wie heißt diese Verbindung?', m: mi }, STEM[n] + 'an', others.map(k => STEM[k] + 'an'));
    if (n >= 2 && n <= 6) {
      const sme = 'C=C' + 'C'.repeat(n - 2), mie = mol(sme), ne = nm(C.parseSmiles(sme));
      if (mie >= 0 && ne) mc('ba', { c: jit(260 + n * 20, 20), t: T.nom, q: 'Wie heißt diese Verbindung?', m: mie }, ne, [STEM[n] + 'an', n > 2 ? STEM[n] + '-1-in' : null, STEM[Math.min(10, n + 1)] + '-1-en', n > 2 ? STEM[n - 1] + '-1-en' : 'Propen'].map(x => x && (x.endsWith('-1-in') ? x.replace('-1-in', '-1-in') : x)));
      const smo = 'C'.repeat(n) + 'O', mio = mol(smo), no = nm(C.parseSmiles(smo));
      if (mio >= 0 && no) mc('ba', { c: jit(220 + n * 20, 20), t: T.nom, q: 'Wie heißt diese Verbindung?', m: mio }, no, [STEM[n] + 'an', STEM[Math.min(10, n + 1)] + (n + 1 > 2 ? 'an-1-ol' : 'anol'), STEM[n] + 'al', n > 2 ? STEM[n - 1] + 'an-1-ol' : 'Methanol'].map(x => x.replace(/^(Meth|Eth)an-1-ol$/, '$1anol')));
    }
  }
}

/* ---------- 3. naming pool ---------- */
const POOL = []; // {spec, g, sm, name, score}
function fillPool(target) {
  const seen = new Set(); let tries = 0;
  while (POOL.length < target && tries < target * 30) {
    tries++;
    const L = pick([0, 1, 1, 2, 2, 2, 3, 3, 4]);
    const spec = randSpec(L); const b = build(spec); if (!b) continue;
    if (b.g.length > 13 || b.g.length < 3) continue;
    const sm = smOf(b.g), c = can(sm); if (!c || seen.has(c)) continue;
    const n = nm(b.g); if (!n) continue;
    const mi = mol(sm); if (mi < 0) continue;
    seen.add(c); POOL.push({ spec, g: C.parseSmiles(sm), sm, can: c, name: n, mi, score: scoreName(n, spec, b.g) });
  }
}
function isomerNames(e, want = 6) {
  const out = new Map();
  for (let t = 0; t < 40 && out.size < want; t++) {
    const s = perturb(e.spec); if (!s) continue; const b = build(s); if (!b) continue;
    const sm = smOf(b.g), c = can(sm); if (!c || c === e.can) continue;
    const n = nm(b.g); if (!n || n === e.name) continue;
    out.set(c, { name: n, sm, g: b.g, can: c });
  }
  return [...out.values()];
}
const celoName = s => 200 + 150 * s;
function famNaming() {
  let mcN = 0, mmN = 0, txN = 0;
  for (const e of POOL) {
    const isos = isomerNames(e);
    const base = celoName(e.score);
    const r = rnd();
    if (r < 0.22) { // choose the name
      const w = [reversedNumbering(e.name), nonAlphabetical(e.name), ...shuffle(isos.map(x => x.name))].filter(Boolean);
      // a distractor must not name the same molecule correctly
      if (mc('nm', { c: jit(base), t: T.nom, q: 'Wie heißt diese Verbindung nach IUPAC?', m: e.mi }, e.name, w)) mcN++;
    } else if (r < 0.42) { // choose the structure
      const opts = isos.filter(x => mol(x.sm) >= 0).slice(0, 3);
      if (opts.length < 3) continue;
      const o = shuffle([e.mi, ...opts.map(x => mol(x.sm))]);
      if (put('sm', { c: jit(base - 40), t: T.str, k: 'mm', q: `Welche Struktur ist <b>${e.name}</b>?`, o, a: o.indexOf(e.mi) })) mmN++;
    } else { // type the name
      const acc = [e.name, C.schoolForm(e.name)].filter(Boolean);
      if (/Toluol|Anilin|Phenol/.test(e.name)) {}
      if (put('tn', { c: jit(base + 170), t: T.nom, k: 'txt', q: 'Benenne die Verbindung nach IUPAC.', m: e.mi, mt: e.mi, acc, ans: e.name })) txN++;
    }
  }
  return { mcN, mmN, txN };
}

/* ---------- 4. formula, class, hybridisation, oxidation number ---------- */
function famFormula(n) {
  for (const e of shuffle(POOL).slice(0, n)) {
    const f = C.formula(e.g);
    const m = f.match(/C(\d*)H(\d*)(.*)/); if (!m) continue;
    const c = +(m[1] || 1), h = +(m[2] || 1), rest = m[3];
    const v = [[c, h + 2], [c, h - 2], [c + 1, h + 2], [c - 1, h], [c, h + 1]].filter(([a, b]) => a > 0 && b >= 0).map(([a, b]) => fmt('C' + (a > 1 ? a : '') + 'H' + (b > 1 ? b : b === 1 ? '' : '0') + rest).replace('H₀', ''));
    if (rest) v.push(fmt('C' + (c > 1 ? c : '') + 'H' + (h > 1 ? h : '') + rest.replace(/O(\d*)/, (x, d) => 'O' + ((+(d || 1)) + 1))));
    if (chance(0.5)) { put('ft', { c: jit(430 + 28 * e.g.length + (e.spec.ring ? 120 : 0)), t: T.str, k: 'txt', q: 'Gib die Summenformel der Verbindung an.', m: e.mi, ff: f, ans: fmt(f), ph: 'z. B. C5H12O' }); continue; }
    mc('fo', { c: jit(300 + 28 * e.g.length + (e.spec.ring ? 120 : 0)), t: T.str, q: 'Welche Summenformel hat die Verbindung?', m: e.mi, f: 1 }, fmt(f), shuffle(v));
  }
}
const CLASS = {
  Alkan: ['CCCC', 'CC(C)CC', 'CCCCCC', 'CC(C)(C)C', 'CCC(C)CC'], Alken: ['C=CCC', 'CC=CC', 'CC(C)=CC', 'C=C(C)CC'], Alkin: ['C#CCC', 'CC#CC', 'C#CC(C)C'],
  Alkohol: ['CCCO', 'CC(O)C', 'CC(C)(C)O', 'OCCCC', 'CC(O)CC'], Aldehyd: ['CCC=O', 'CC(C)C=O', 'CCCCC=O'], Keton: ['CC(=O)C', 'CCC(=O)C', 'CCC(=O)CC', 'CC(C)C(=O)C'],
  Carbonsäure: ['CC(=O)O', 'CCC(=O)O', 'CC(C)C(=O)O', 'CCCCC(=O)O'], Ester: ['CC(=O)OC', 'CC(=O)OCC', 'CCC(=O)OC', 'CC(C)C(=O)OCC', 'CCCC(=O)OC'],
  Ether: ['COC', 'CCOCC', 'CCOC', 'CC(C)OC', 'CCCOC'], Amin: ['CCN', 'CC(C)N', 'CCCN', 'CCNCC', 'CN(C)C'], Amid: ['CC(N)=O', 'CCC(N)=O', 'CC(=O)NC', 'CC(C)C(N)=O'],
  Nitril: ['CC#N', 'CCC#N', 'CC(C)C#N'], Halogenalkan: ['CCCBr', 'CC(Cl)C', 'CC(C)(C)Br', 'CCCCI', 'ClCCC'],
};
const CONFUSE = { Aldehyd: ['Keton', 'Carbonsäure', 'Alkohol'], Keton: ['Aldehyd', 'Ester', 'Ether'], Ester: ['Ether', 'Carbonsäure', 'Keton'], Ether: ['Ester', 'Alkohol', 'Keton'], Amin: ['Amid', 'Nitril', 'Alkohol'], Amid: ['Amin', 'Ester', 'Carbonsäure'], Nitril: ['Amin', 'Alkin', 'Amid'], Carbonsäure: ['Ester', 'Aldehyd', 'Alkohol'], Alkohol: ['Ether', 'Aldehyd', 'Carbonsäure'], Alken: ['Alkan', 'Alkin', 'Aromat'], Alkin: ['Alken', 'Nitril', 'Alkan'], Alkan: ['Alken', 'Halogenalkan', 'Alkin'], Halogenalkan: ['Alkan', 'Alkohol', 'Ether'] };
const CLEVEL = { Alkan: 300, Alken: 330, Alkin: 360, Alkohol: 340, Halogenalkan: 400, Aldehyd: 480, Keton: 480, Carbonsäure: 450, Ether: 580, Ester: 640, Amin: 560, Amid: 720, Nitril: 680 };
function famClass() {
  for (const cls in CLASS) for (const sm of CLASS[cls]) {
    const mi = mol(sm); if (mi < 0) continue;
    mc('cl', { c: jit(CLEVEL[cls], 50), t: T.fg, q: 'Zu welcher Stoffklasse gehört die Verbindung?', m: mi }, cls, CONFUSE[cls]);
    const n = nm(C.parseSmiles(sm));
    if (n) mc('cl', { c: jit(CLEVEL[cls] + 80, 50), t: T.fg, q: `Zu welcher Stoffklasse gehört <b>${n}</b>?` }, cls, CONFUSE[cls]);
  }
  // functional group spotting in larger molecules
  const FGQ = [['Hydroxygruppe', a => a.el === 'O' && a.nb.length === 1 && a.nb[0].o === 1], ['Carbonylgruppe', a => a.el === 'O' && a.nb.length === 1 && a.nb[0].o === 2]];
  for (const e of shuffle(POOL).slice(0, 220)) {
    const g = e.g, idx = [];
    g.forEach((a, i) => { if (a.el === 'O' && a.nb.length === 1 && a.nb[0].o === 1 && !g[a.nb[0].n].nb.some(x => g[x.n].el === 'O' && x.o === 2)) idx.push([i, 'OH']); });
    g.forEach((a, i) => { if (a.el === 'N' && a.nb.length === 1 && a.nb[0].o === 1) idx.push([i, 'NH2']); });
    if (idx.length !== 1 || e.spec.ring && e.spec.arom) continue;
    const [i, kind] = idx[0];
    put('fp', { c: jit(kind === 'OH' ? 420 : 520, 50), t: T.fg, k: 'pick', q: kind === 'OH' ? 'Klicke auf das Atom der Hydroxygruppe.' : 'Klicke auf das Stickstoffatom der Aminogruppe.', m: e.mi, a: [i] });
  }
}
function hyb(g, i) { const o = g[i].nb.map(e => e.o); if (o.includes(3) || o.filter(x => x === 2).length >= 2) return 'sp'; if (o.includes(2)) return 'sp²'; return 'sp³'; }
function famHybrid(n) {
  let k = 0;
  for (const e of shuffle(POOL)) {
    if (k >= n) break;
    const cs = e.g.map((a, i) => i).filter(i => e.g[i].el === 'C');
    const want = pick(['sp', 'sp²', 'sp²', 'sp³']);
    const cand = cs.filter(i => hyb(e.g, i) === want); if (!cand.length) continue;
    const i = pick(cand), h = hyb(e.g, i);
    const special = e.g[i].nb.some(x => e.g[x.n].el === 'N' && x.o === 3) ? 200 : e.g[i].nb.some(x => e.g[x.n].el === 'O' && x.o === 2) ? 120 : 0;
    if (mc('hy', { c: jit(720 + (h === 'sp' ? 120 : h === 'sp²' ? 60 : 0) + special + 8 * e.g.length, 50), t: T.bond, q: 'Welche Hybridisierung hat das markierte C-Atom?', m: e.mi, hl: [i] }, h, ['sp', 'sp²', 'sp³'])) k++;
  }
}
function oxNum(g, i) {
  let ox = -C.hcount(g, i);
  for (const e of g[i].nb) { const el = g[e.n].el; if (['O', 'N', 'F', 'Cl', 'Br', 'I', 'S'].includes(el)) ox += e.o; }
  return ox;
}
const sgn = v => v > 0 ? '+' + v : v < 0 ? '−' + (-v) : '0';
function famOx(n) {
  let k = 0;
  for (const e of shuffle(POOL)) {
    if (k >= n) break;
    const cs = e.g.map((a, i) => i).filter(i => e.g[i].el === 'C');
    const het = cs.filter(i => e.g[i].nb.some(x => e.g[x.n].el !== 'C'));
    const i = het.length && chance(0.75) ? pick(het) : pick(cs);
    const v = oxNum(e.g, i);
    const w = [v - 1, v + 1, v - 2, v + 2, -v].filter(x => x !== v && x >= -4 && x <= 4).map(sgn);
    const hard = e.g[i].nb.some(x => e.g[x.n].el !== 'C') ? 140 : 0;
    if (mc('ox', { c: jit(1000 + hard + 15 * e.g.length + Math.abs(v) * 25, 50), t: T.bond, q: 'Welche Oxidationszahl hat das markierte C-Atom?', m: e.mi, hl: [i] }, sgn(v), shuffle(w))) k++;
  }
}

/* ---------- 5. stereochemistry ---------- */
function centersOf(sm) { try { return C.stereoTags(sm).CIP_atoms.map(x => x[0]); } catch (e) { return []; } }
function famStereoCenters(n) {
  let k = 0, z = 0;
  for (const e of shuffle(POOL)) {
    if (k >= n) break;
    const cs = centersOf(e.sm).filter(i => e.g[i].el === 'C');
    if (cs.length >= 1 && cs.length <= 3) {
      if (chance(0.7)) { if (put('sc', { c: jit(1120 + 140 * (cs.length - 1) + 10 * e.g.length + (e.spec.ring ? 120 : 0), 50), t: T.st, k: 'pick', q: cs.length === 1 ? 'Klicke auf das Stereozentrum.' : `Markiere alle ${cs.length === 2 ? 'zwei' : 'drei'} Stereozentren.`, m: e.mi, a: cs })) k++; }
      else if (mc('sc', { c: jit(1180 + 120 * cs.length + (e.spec.ring ? 120 : 0), 50), t: T.st, q: 'Wie viele Stereozentren besitzt die Verbindung?', m: e.mi }, String(cs.length), ['0', '1', '2', '3', '4'])) k++;
    } else if (cs.length === 0 && z < n * 0.15 && e.g.length >= 6) {
      if (mc('sc', { c: jit(1250, 50), t: T.st, q: 'Wie viele Stereozentren besitzt die Verbindung?', m: e.mi }, '0', ['1', '2', '3'])) { k++; z++; }
    }
  }
}
// substituent fragments with an increasing CIP rank (verified by RDKit, not by this order)
const GRP = { H: '', Me: 'C', Et: 'CC', Pr: 'CCC', iPr: 'C(C)C', vin: 'C=C', tBu: 'C(C)(C)C', CH2OH: 'CO', CHO: 'C=O', COOH: 'C(=O)O', NH2: 'N', OH: 'O', OMe: 'OC', F: 'F', Cl: 'Cl', Br: 'Br', I: 'I', CN: 'C#N', CH2Cl: 'CCl' };
const firstAtom = k => k === 'H' ? 'H' : GRP[k].match(/^(Cl|Br|[A-Z])/)[1];
function famEZ(n) {
  const keys = Object.keys(GRP).filter(k => k !== 'COOH');
  let k = 0, tries = 0;
  while (k < n && tries < n * 40) {
    tries++;
    const L = [pick(keys), pick(keys)], R = [pick(keys), pick(keys)];
    if (L[0] === 'H' && L[1] === 'H' || R[0] === 'H' && R[1] === 'H') continue;
    const none = L[0] === L[1] || R[0] === R[1];
    if (none && !chance(0.12)) continue;
    // p (on left, "/") must be non-H; r (on right "/") non-H
    let [p, q] = L[0] === 'H' ? [L[1], L[0]] : L; let [r, s] = R[0] === 'H' ? [R[1], R[0]] : R;
    const sm = `${GRP[p]}/C${q === 'H' ? '' : '(' + GRP[q] + ')'}=C${s === 'H' ? '' : '(' + GRP[s] + ')'}/${GRP[r]}`;
    let tags; try { tags = C.stereoTags(sm); } catch (e) { continue; }
    const lab = tags.CIP_bonds.length ? tags.CIP_bonds[0][2] : '';
    let ans = lab === '(E)' ? 'E' : lab === '(Z)' ? 'Z' : (none ? 'keine E/Z-Isomerie' : null);
    if (!ans) continue;
    // geometry check: p and r must be drawn trans
    const nP = C.parseSmiles(GRP[p]).length, nQ = q === 'H' ? 0 : C.parseSmiles(GRP[q]).length, nS = s === 'H' ? 0 : C.parseSmiles(GRP[s]).length;
    const iP = 0, iCL = nP, iCR = nP + 1 + nQ, iR = iCR + 1 + nS;
    const mi = mol(sm); if (mi < 0) continue;
    const A = M[mi][0]; const side = (a) => { const [x1, y1] = [A[iCL][1], A[iCL][2]], [x2, y2] = [A[iCR][1], A[iCR][2]]; return Math.sign((x2 - x1) * (A[a][2] - y1) - (y2 - y1) * (A[a][1] - x1)); };
    if (!none && side(iP) === side(iR)) { FAILS.push('EZ geometry ' + sm); continue; }
    const simple = (q === 'H' && s === 'H');
    const sameFirst = (x, y) => x !== 'H' && y !== 'H' && firstAtom(x) === firstAtom(y);
    const deep = sameFirst(L[0], L[1]) || sameFirst(R[0], R[1]);
    const c = none ? 1250 : simple ? 1180 : deep ? 1560 : 1360;
    if (mc('ez', { c: jit(c + (/CHO|COOH|CH2OH|vin|iPr|tBu|CH2Cl/.test([p, q, r, s].join()) ? 60 : 0), 50), t: T.st, q: 'Wie ist die Doppelbindung konfiguriert?', m: mi }, ans, ['E', 'Z', 'keine E/Z-Isomerie'])) k++;
  }
}
function famRS(n) {
  const keys = Object.keys(GRP).filter(k => k !== 'H');
  let k = 0, tries = 0;
  while (k < n && tries < n * 40) {
    tries++;
    const withH = chance(0.75);
    const gs = withH ? ['H', ...shuffle(keys).slice(0, 3)] : shuffle(keys).slice(0, 4);
    if (new Set(gs).size < 4) continue;
    const tag = chance(0.5) ? '@' : '@@';
    const [a, b, c2, d] = gs.slice(1).concat([]); // a,b,c2 + maybe d
    const explicitH = withH && chance(0.6);
    let sm;
    if (withH) sm = explicitH ? `${GRP[gs[1]]}[C${tag}]([H])(${GRP[gs[2]]})${GRP[gs[3]]}` : `${GRP[gs[1]]}[C${tag}H](${GRP[gs[2]]})${GRP[gs[3]]}`;
    else sm = `${GRP[gs[0]]}[C${tag}](${GRP[gs[1]]})(${GRP[gs[2]]})${GRP[gs[3]]}`;
    let tags; try { tags = C.stereoTags(sm); } catch (e) { continue; }
    const t = tags.CIP_atoms.find(x => x[1] === '(R)' || x[1] === '(S)'); if (!t) continue;
    if (tags.CIP_atoms.length !== 1) continue;
    const mi = explicitH ? mol(sm, true) : mol(sm); if (mi < 0) continue;
    if (!M[mi][1].some(b => b[3])) continue; // need a wedge
    const firsts = gs.filter(x => x !== 'H').map(firstAtom);
    const deep = firsts.length !== new Set(firsts).size;
    const deeper = ['iPr', 'vin', 'CH2OH', 'CHO', 'CH2Cl'].filter(x => gs.includes(x)).length >= 2;
    const c = 1420 + (deep ? 180 : 0) + (deeper ? 120 : 0) + (withH && !explicitH ? 130 : 0) + (!withH ? 200 : 0);
    if (mc('rs', { c: jit(c, 50), t: T.st, q: 'Ist das Stereozentrum R- oder S-konfiguriert?', m: mi }, t[1][1], ['R', 'S'].filter(x => x !== t[1][1]).concat(['Kein Stereozentrum']))) k++;
  }
}

/* ---------- 6. acids, bases, boiling points ---------- */
const ACIDS = [['Ethan', 50], ['Ammoniak', 38], ['Ethin', 25], ['Aceton (α-H)', 19.3], ['Ethanol', 16], ['Wasser', 15.7], ['Phenol', 10], ['Ammonium-Ion', 9.25], ['4-Nitrophenol', 7.15], ['Kohlensäure', 6.35], ['Essigsäure', 4.76], ['Benzoesäure', 4.2], ['Ameisensäure', 3.75], ['Fluorwasserstoff', 3.17], ['Chloressigsäure', 2.87], ['Dichloressigsäure', 1.35], ['Trichloressigsäure', 0.66], ['Chlorwasserstoff', -7], ['Schwefelsäure', -3], ['4-Chlorbutansäure', 4.52], ['2-Chlorbutansäure', 2.86], ['Pentan-2,4-dion', 9], ['Malonsäurediethylester', 13], ['Cyclopentadien', 16], ['Pikrinsäure', 0.4], ['Blausäure (HCN)', 9.2], ['Schwefelwasserstoff', 7]];
const BASES = [['Triethylamin', 10.75], ['Methylamin', 10.64], ['Ammoniak', 9.25], ['Pyridin', 5.25], ['Anilin', 4.6], ['4-Nitroanilin', 1.0], ['Acetamid', -0.5], ['Pyrrol', -3.8], ['Imidazol', 6.95], ['Guanidin', 13.6], ['Wasser', -1.7], ['Hydroxid-Ion', 15.7], ['Acetat-Ion', 4.76], ['Amid-Ion (NH₂⁻)', 38], ['Ethanolat-Ion', 16]];
function spaced(list, n, gap) {
  for (let t = 0; t < 200; t++) { const s = shuffle(list).slice(0, n); const v = s.map(x => x[1]).sort((a, b) => a - b); if (v.every((x, i) => !i || x - v[i - 1] >= gap)) return s; }
  return null;
}
function famAcidBase(n) {
  let k = 0, t = 0;
  while (k < n && t < n * 20) {
    t++;
    const base = chance(0.4), list = base ? BASES : ACIDS;
    const s = spaced(list, pick([3, 4, 4]), 1.4); if (!s) continue;
    const sorted = s.slice().sort((a, b) => a[1] - b[1]);
    const most = chance(0.75);
    const ans = base ? (most ? sorted[sorted.length - 1] : sorted[0]) : (most ? sorted[0] : sorted[sorted.length - 1]);
    const runner = base ? (most ? sorted[sorted.length - 2] : sorted[1]) : (most ? sorted[1] : sorted[sorted.length - 2]);
    const gap = Math.abs(ans[1] - runner[1]);
    const c = 1000 + (gap < 2.5 ? 260 : gap < 5 ? 140 : 0) + (s.length - 3) * 60 + (most ? 0 : 80) + (base ? 120 : 0);
    const q = base ? (most ? 'Welche Verbindung ist die stärkste Base?' : 'Welche Verbindung ist die schwächste Base?') : (most ? 'Welche Verbindung ist die stärkste Säure?' : 'Welche Verbindung ist die schwächste Säure?');
    if (put('ab', { c: jit(c, 40), t: T.ab, k: 'mc', q, o: s.map(x => x[0]), a: s.indexOf(ans), e: s.slice().sort((a, b) => a[1] - b[1]).map(x => `${x[0]}: pKs ${base ? '(korr. Säure) ' : ''}${String(x[1]).replace('.', ',')}`).join(' · ') })) k++;
    if (chance(0.3)) {
      const s2 = spaced(list, 4, 1.4); if (!s2) continue;
      const items = s2.slice().sort((a, b) => base ? b[1] - a[1] : a[1] - b[1]).map(x => x[0]);
      put('ao', { c: jit(1300 + (base ? 150 : 0), 60), t: T.ab, k: 'ord', q: base ? 'Ordne nach abnehmender Basizität.' : 'Ordne nach abnehmender Acidität.', items, ends: base ? ['stärkste Base', 'schwächste Base'] : ['stärkste Säure', 'schwächste Säure'] });
    }
  }
}
const BP = [['Methan', -162], ['Ethan', -89], ['Propan', -42], ['Butan', -1], ['Pentan', 36], ['Hexan', 69], ['Heptan', 98], ['Octan', 126], ['2-Methylbutan', 28], ['2,2-Dimethylpropan', 10], ['2-Methylpentan', 60], ['Methanol', 65], ['Ethanol', 78], ['Propan-1-ol', 97], ['Butan-1-ol', 118], ['Pentan-1-ol', 138], ['2-Methylpropan-2-ol', 82], ['Diethylether', 35], ['Dimethylether', -24], ['Aceton', 56], ['Propanal', 48], ['Ethanal', 20], ['Essigsäure', 118], ['Wasser', 100], ['Ethan-1,2-diol', 197], ['Propan-1,2,3-triol', 290], ['Chlormethan', -24], ['Dichlormethan', 40], ['Trichlormethan', 61], ['Tetrachlormethan', 77], ['Ethylacetat', 77], ['Methylamin', -6], ['Ethylamin', 17], ['Trimethylamin', 3], ['Benzol', 80], ['Toluol', 111], ['Phenol', 182], ['Butanal', 75]];
function famBP(n) {
  let k = 0, t = 0;
  while (k < n && t < n * 30) {
    t++;
    const s = spaced(BP, pick([2, 2, 3, 3, 4]), 14); if (!s) continue;
    const sorted = s.slice().sort((a, b) => b[1] - a[1]);
    const hi = chance(0.7), ans = hi ? sorted[0] : sorted[sorted.length - 1], run = hi ? sorted[1] : sorted[sorted.length - 2];
    const gap = Math.abs(ans[1] - run[1]);
    const c = 430 + (s.length - 2) * 90 + (gap < 25 ? 240 : gap < 50 ? 120 : 0) + (hi ? 0 : 60);
    if (put('bp', { c: jit(c, 40), t: T.phys, k: 'mc', q: hi ? 'Welche Verbindung hat die höchste Siedetemperatur?' : 'Welche Verbindung hat die niedrigste Siedetemperatur?', o: s.map(x => x[0]), a: s.indexOf(ans), e: sorted.map(x => `${x[0]}: ${x[1]} °C`).join(' · ') })) k++;
  }
}

/* ---------- 7. substitution / elimination mechanisms ---------- */
function famMech() {
  const SUB = { methyl: ['Brommethan', 'Iodmethan', 'Chlormethan'], prim: ['1-Brombutan', '1-Brompropan', '1-Chlorpentan', '1-Iodhexan', '1-Brom-2-methylpropan'], sec: ['2-Brombutan', '2-Brompentan', '3-Brompentan', 'Bromcyclohexan', '2-Iodpropan'], tert: ['2-Brom-2-methylpropan', '2-Chlor-2-methylbutan', '3-Brom-3-methylpentan', '1-Brom-1-methylcyclohexan'], talc: ['2-Methylbutan-2-ol', '2-Methylpropan-2-ol', '1-Methylcyclohexan-1-ol'] };
  const RULES = [
    ['methyl', 'NaOH (aq)', 'SN2', 1280], ['methyl', 'NaCN in DMSO', 'SN2', 1300], ['prim', 'NaI in Aceton', 'SN2', 1320], ['prim', 'NaCN in DMSO', 'SN2', 1320], ['prim', 'NaN₃ in DMF', 'SN2', 1350],
    ['prim', 'KOtBu, Δ', 'E2', 1560], ['sec', 'NaCN in DMSO', 'SN2', 1500], ['sec', 'NaN₃ in DMF', 'SN2', 1500], ['sec', 'NaOEt in EtOH, Δ', 'E2', 1480], ['sec', 'KOtBu, Δ', 'E2', 1450],
    ['tert', 'NaOEt in EtOH, Δ', 'E2', 1420], ['tert', 'H₂O, 25 °C', 'SN1', 1380], ['tert', 'Ethanol, 25 °C', 'SN1', 1450], ['tert', 'KOH in Ethanol, Δ', 'E2', 1430], ['talc', 'konz. H₂SO₄, Δ', 'E1', 1550], ['talc', 'konz. H₃PO₄, Δ', 'E1', 1580],
  ];
  const WHY = { SN2: 'Gutes Nucleophil, wenig gehindertes C-Atom, polar-aprotisches oder wässriges Medium → SN2.', E2: 'Starke Base und Wärme (oder sperrige Base) → E2.', SN1: 'Tertiäres Substrat, schwaches Nucleophil, polar-protisches Lösungsmittel → SN1 (Solvolyse).', E1: 'Tertiärer Alkohol + starke Säure + Wärme → Carbokation → E1.' };
  for (const [cls, cond, ans, c] of RULES) for (const s of SUB[cls]) put('me', { c: jit(c, 50), t: T.mech, k: 'mc', q: `<b>${s}</b> + ${cond}: Welcher Mechanismus dominiert?`, o: ['SN1', 'SN2', 'E1', 'E2'], a: ['SN1', 'SN2', 'E1', 'E2'].indexOf(ans), e: WHY[ans] });
}

/* ---------- 8. reactions on graphs ---------- */
function alkeneAtoms(g) { for (let i = 0; i < g.length; i++) for (const e of g[i].nb) if (e.o === 2 && g[i].el === 'C' && g[e.n].el === 'C' && i < e.n) return [i, e.n]; return null; }
function addAcross(g, a, b, Xa, Xb) { const h = clone(g); C.setOrder(h, a, b, 1); if (Xa) C.attach(h, a, Xa); if (Xb) C.attach(h, b, Xb); return valOk(h) ? h : null; }
const subst = (g, i) => g[i].nb.filter(e => g[e.n].el === 'C').length;
function markov(g) { const ab = alkeneAtoms(g); if (!ab) return null; const [a, b] = ab; const ha = C.hcount(g, a), hb = C.hcount(g, b); if (ha === hb) return null; return ha < hb ? { more: a, less: b } : { more: b, less: a }; }
function hetero(g, el, test) { return g.map((a, i) => i).filter(i => g[i].el === el && test(i)); }
function ohCarbon(g) { const o = hetero(g, 'O', i => g[i].nb.length === 1 && g[i].nb[0].o === 1); if (o.length !== 1) return null; const c = g[o[0]].nb[0].n; if (g[c].nb.some(e => g[e.n].el === 'O' && e.o === 2)) return null; return { o: o[0], c }; }
function carbonyl(g) { const o = hetero(g, 'O', i => g[i].nb.length === 1 && g[i].nb[0].o === 2); if (o.length !== 1) return null; const c = g[o[0]].nb[0].n; const others = g[c].nb.filter(e => e.n !== o[0]); if (others.some(e => g[e.n].el !== 'C')) return null; return { o: o[0], c, ald: C.hcount(g, c) > 0 }; }
function acid(g) { for (let c = 0; c < g.length; c++) { const os = g[c].nb.filter(e => g[e.n].el === 'O'); if (os.length === 2 && os.some(e => e.o === 2) && os.some(e => e.o === 1 && g[e.n].nb.length === 1)) return { c, oh: os.find(e => e.o === 1).n, oo: os.find(e => e.o === 2).n }; } return null; }
function halide(g) { const x = hetero(g, 'Br', () => true).concat(hetero(g, 'Cl', () => true), hetero(g, 'I', () => true)); if (x.length !== 1) return null; return { x: x[0], c: g[x[0]].nb[0].n }; }
function nitrile(g) { const n = hetero(g, 'N', i => g[i].nb.length === 1 && g[i].nb[0].o === 3); return n.length === 1 ? { n: n[0], c: g[n[0]].nb[0].n } : null; }
const hasStereoBond = sm => { try { return C.stereoTags(sm).CIP_bonds.length > 0; } catch (e) { return true; } };
function eliminate(g, c, leave, mode) { // remove leaving atom index `leave` on carbon c; mode zaitsev|hofmann
  const betas = g[c].nb.filter(e => g[e.n].el === 'C' && C.hcount(g, e.n) > 0).map(e => e.n);
  if (!betas.length) return null;
  const prods = betas.map(b => { let h = clone(g); C.setOrder(h, c, b, 2); h = removeAtom(h, leave); const bi = b > leave ? b - 1 : b, ci = c > leave ? c - 1 : c; return { h, deg: subst(h, ci) + subst(h, bi) }; });
  const uniq = new Map(); prods.forEach(p => { const cc = can(smOf(p.h)); if (cc && !uniq.has(cc)) uniq.set(cc, p); });
  const list = [...uniq.values()]; if (!list.length) return null;
  list.sort((x, y) => y.deg - x.deg);
  const target = mode === 'zaitsev' ? list[0] : list[list.length - 1];
  if (list.length > 1 && list.filter(p => p.deg === target.deg).length > 1) return null;
  if (hasStereoBond(smOf(target.h))) return null;
  return { main: target.h, others: list.filter(p => p !== target).map(p => p.h) };
}
const R = { // reagent -> function(g) => product graph | 'none' | null(not applicable)
  hbr: g => { const m = markov(g); return m ? addAcross(g, m.more, m.less, 'Br', null) : null; },
  hcl: g => { const m = markov(g); return m ? addAcross(g, m.more, m.less, 'Cl', null) : null; },
  hbrp: g => { const m = markov(g); return m ? addAcross(g, m.less, m.more, 'Br', null) : null; },
  hyd: g => { const m = markov(g); return m ? addAcross(g, m.more, m.less, 'O', null) : null; },
  bh3: g => { const m = markov(g); return m ? addAcross(g, m.less, m.more, 'O', null) : null; },
  br2: g => { const ab = alkeneAtoms(g); return ab ? addAcross(g, ab[0], ab[1], 'Br', 'Br') : null; },
  h2: g => { const ab = alkeneAtoms(g); return ab ? addAcross(g, ab[0], ab[1], null, null) : null; },
  kmno4c: g => { const ab = alkeneAtoms(g); return ab ? addAcross(g, ab[0], ab[1], 'O', 'O') : null; },
  pcc: g => { const x = ohCarbon(g); if (!x) return null; const h = C.hcount(g, x.c); if (h === 0) return 'none'; const k = clone(g); C.setOrder(k, x.c, x.o, 2); return k; },
  cr: g => { const x = ohCarbon(g); if (!x) return null; const h = C.hcount(g, x.c); if (h === 0) return 'none'; const k = clone(g); C.setOrder(k, x.c, x.o, 2); if (h >= 2) C.bond(k, x.c, C.add(k, 'O')); return k; },
  pbr3: g => { const x = ohCarbon(g); if (!x || C.hcount(g, x.c) === 0) return null; const k = clone(g); k[x.o].el = 'Br'; return k; },
  dehyd: g => { const x = ohCarbon(g); if (!x) return null; const r = eliminate(g, x.c, x.o, 'zaitsev'); return r ? r.main : null; },
  nabh4: g => { const x = carbonyl(g); if (!x) return null; const k = clone(g); C.setOrder(k, x.c, x.o, 1); return k; },
  lahacid: g => { const x = acid(g); if (!x) return null; return removeAtom(g, x.oo); },
  naoh: g => { const x = halide(g); if (!x || subst(g, x.c) > 1) return null; const k = clone(g); k[x.x].el = 'O'; return k; },
  nacn: g => { const x = halide(g); if (!x || subst(g, x.c) > 2) return null; let k = clone(g); k = removeAtom(k, x.x); const c = x.c > x.x ? x.c - 1 : x.c; const n = C.add(k, 'C'); C.bond(k, c, n); C.bond(k, n, C.add(k, 'N'), 3); return k; },
  nh3: g => { const x = halide(g); if (!x || subst(g, x.c) > 1) return null; const k = clone(g); k[x.x].el = 'N'; return k; },
  naoet: g => { const x = halide(g); if (!x || subst(g, x.c) < 2) return null; const r = eliminate(g, x.c, x.x, 'zaitsev'); return r ? r.main : null; },
  kotbu: g => { const x = halide(g); if (!x || subst(g, x.c) < 2) return null; const r = eliminate(g, x.c, x.x, 'hofmann'); return r ? r.main : null; },
  mgco2: g => { const x = halide(g); if (!x) return null; let k = removeAtom(g, x.x); const c = x.c > x.x ? x.c - 1 : x.c; const a = C.add(k, 'C'); C.bond(k, c, a); C.bond(k, a, C.add(k, 'O'), 2); C.bond(k, a, C.add(k, 'O')); return k; },
  lahcn: g => { const x = nitrile(g); if (!x) return null; const k = clone(g); C.setOrder(k, x.c, x.n, 1); return k; },
  hydcn: g => { const x = nitrile(g); if (!x) return null; const k = clone(g); k[x.n].el = 'O'; C.setOrder(k, x.c, x.n, 2); C.bond(k, x.c, C.add(k, 'O')); return k; },
  meoh: g => esterify(g, 'C'), etoh: g => esterify(g, 'CC'),
  memgbr: g => grig(g, 'C'), etmgbr: g => grig(g, 'CC'), phmgbr: null,
};
function esterify(g, alk) { const x = acid(g); if (!x) return null; const k = clone(g); C.attach(k, x.oh, alk); return k; }
function grig(g, alk) { const x = carbonyl(g); if (!x) return null; const k = clone(g); C.setOrder(k, x.c, x.o, 1); C.attach(k, x.c, alk); return valOk(k) ? k : null; }
const RG = {
  hbr: 'HBr', hcl: 'HCl', hbrp: 'HBr, ROOR (Peroxide)', hyd: 'H₂O, H₂SO₄ (kat.)', bh3: 'BH₃·THF; dann H₂O₂, NaOH', br2: 'Br₂ in CH₂Cl₂', h2: 'H₂, Pd/C', kmno4c: 'KMnO₄ (kalt, verdünnt)',
  pcc: 'PCC in CH₂Cl₂', cr: 'K₂Cr₂O₇, H₂SO₄', pbr3: 'PBr₃', dehyd: 'konz. H₂SO₄, Δ', nabh4: 'NaBH₄ in MeOH', lahacid: 'LiAlH₄; dann H₃O⁺', naoh: 'NaOH (aq)', nacn: 'NaCN in DMSO', nh3: 'NH₃ (Überschuss)',
  naoet: 'NaOEt in EtOH, Δ', kotbu: 'KOtBu, Δ', mgco2: 'Mg, Ether; dann CO₂; dann H₃O⁺', lahcn: 'LiAlH₄; dann H₂O', hydcn: 'H₃O⁺, Δ', meoh: 'CH₃OH, H₂SO₄ (kat.)', etoh: 'C₂H₅OH, H₂SO₄ (kat.)', memgbr: 'CH₃MgBr; dann H₃O⁺', etmgbr: 'C₂H₅MgBr; dann H₃O⁺',
};
const RC = { h2: 900, br2: 950, hbr: 1100, hcl: 1100, hyd: 1180, kmno4c: 1300, bh3: 1380, hbrp: 1450, pcc: 1300, cr: 1260, pbr3: 1250, dehyd: 1400, nabh4: 1150, lahacid: 1380, naoh: 1150, nacn: 1260, nh3: 1300, naoet: 1480, kotbu: 1650, mgco2: 1720, lahcn: 1460, hydcn: 1430, meoh: 1240, etoh: 1260, memgbr: 1500, etmgbr: 1520 };
const RWHY = { hbr: 'Markownikow-Addition über das stabilere Carbokation.', hcl: 'Markownikow-Addition über das stabilere Carbokation.', hbrp: 'Radikalische Addition: Br geht an das C-Atom mit mehr H (anti-Markownikow).', hyd: 'Säurekatalysierte Hydratisierung nach Markownikow.', bh3: 'Hydroborierung/Oxidation: OH am weniger substituierten C (anti-Markownikow).', br2: 'Elektrophile Addition: je ein Br an beide C-Atome.', h2: 'Katalytische Hydrierung der Doppelbindung.', kmno4c: 'syn-Dihydroxylierung: zwei OH-Gruppen.', pcc: 'PCC oxidiert primäre Alkohole nur bis zum Aldehyd, sekundäre zum Keton; tertiäre gar nicht.', cr: 'Starkes Oxidationsmittel: primär → Carbonsäure, sekundär → Keton, tertiär → keine Reaktion.', pbr3: 'OH wird durch Br ersetzt.', dehyd: 'Dehydratisierung nach Saizew: die höher substituierte Doppelbindung.', nabh4: 'Reduktion der C=O-Gruppe zum Alkohol.', lahacid: 'LiAlH₄ reduziert Carbonsäuren zu primären Alkoholen.', naoh: 'SN2 am primären C-Atom.', nacn: 'SN2: Cyanid ersetzt das Halogen (Kettenverlängerung um ein C).', nh3: 'SN2 mit Ammoniak (Überschuss) → primäres Amin.', naoet: 'E2 mit kleiner starker Base → Saizew-Produkt.', kotbu: 'E2 mit sperriger Base → Hofmann-Produkt (weniger substituiert).', mgco2: 'Grignard-Reagenz + CO₂ → Carbonsäure mit einem C-Atom mehr.', lahcn: 'Nitril → primäres Amin.', hydcn: 'Saure Hydrolyse: Nitril → Carbonsäure.', meoh: 'Fischer-Veresterung.', etoh: 'Fischer-Veresterung.', memgbr: 'Grignard-Addition an die C=O-Gruppe.', etmgbr: 'Grignard-Addition an die C=O-Gruppe.' };
const APPLIES = { alkene: ['hbr', 'hcl', 'hbrp', 'hyd', 'bh3', 'br2', 'h2', 'kmno4c'], alcohol: ['pcc', 'cr', 'pbr3', 'dehyd'], carbonyl: ['nabh4', 'memgbr', 'etmgbr'], acid: ['lahacid', 'meoh', 'etoh'], halide: ['naoh', 'nacn', 'nh3', 'naoet', 'kotbu', 'mgco2'], nitrile: ['lahcn', 'hydcn'] };
function substrate(kind) {
  const n = ri(3, 6), s = { n, uns: [], fg: [], subs: [] };
  if (kind === 'alkene') { s.uns.push([ri(1, n - 1), 2]); if (chance(0.7)) s.subs.push([ri(1, n), 'Me']); if (chance(0.25)) s.subs.push([ri(1, n), 'Me']); }
  if (kind === 'alcohol') { const p = ri(1, n); s.fg.push([p, 'ol']); if (chance(0.45)) s.subs.push([p, 'Me']); if (chance(0.3)) s.subs.push([ri(2, n), 'Me']); }
  if (kind === 'carbonyl') { if (chance(0.5)) s.fg.push([1, 'al']); else s.fg.push([ri(2, n - 1), 'on']); if (chance(0.4)) s.subs.push([ri(2, n - 1), 'Me']); }
  if (kind === 'acid') { s.fg.push([1, 'cooh']); if (chance(0.4)) s.subs.push([ri(2, n), 'Me']); }
  if (kind === 'halide') { const p = ri(1, n); s.subs.push([p, pick(['Br', 'Br', 'Cl', 'I'])]); if (chance(0.4)) s.subs.push([p === 1 ? 2 : p, 'Me']); if (chance(0.25)) s.subs.push([ri(2, n), 'Me']); }
  if (kind === 'nitrile') { s.fg.push([1, 'cn']); if (chance(0.5)) s.subs.push([ri(2, n), 'Me']); }
  return build(s);
}
const RXN = []; // {sub, subSm, r, prod, prodName, kind}
function famReactions(nProd, nReag) {
  let made = 0, t = 0; const kinds = Object.keys(APPLIES);
  const pending = [];
  while (pending.length < nProd + nReag && t < (nProd + nReag) * 40) {
    t++;
    const kind = pick(kinds), b = substrate(kind); if (!b) continue;
    const g = b.g, sm = smOf(g), sn = nm(g); if (!sn) continue;
    const outs = {};
    for (const r of APPLIES[kind]) { const p = R[r](g); if (p === null || p === undefined) continue; if (p === 'none') { outs[r] = 'none'; continue; } const ps = smOf(p), pc = can(ps), pn = nm(p); if (!pc || !pn) continue; outs[r] = { g: p, sm: ps, can: pc, name: pn }; }
    const rs = Object.keys(outs); if (rs.length < 3) continue;
    const r = pick(rs); pending.push({ g, sm, sn, kind, r, outs });
  }
  for (const x of pending) {
    const out = x.outs[x.r];
    const mi = mol(x.sm); if (mi < 0) continue;
    const ansName = out === 'none' ? 'Keine Reaktion' : out.name;
    const others = Object.entries(x.outs).filter(([r, o]) => r !== x.r && o !== 'none' && (out === 'none' || o.can !== out.can)).map(([r, o]) => o.name);
    const isoNames = out !== 'none' && others.length < 3 ? [] : [];
    const asReagent = made >= nProd || chance(0.38);
    if (!asReagent && out !== 'none' && chance(0.38)) {
      const pm = mol(out.sm);
      if (pm >= 0 && put('rt', { c: jit(RC[x.r] + 240, 60), t: T.rx, k: 'txt', q: 'Benenne das Hauptprodukt.', m: mi, rg: RG[x.r], mt: pm, ans: out.name, e: RWHY[x.r] })) made++;
      RXN.push(x); continue;
    }
    if (!asReagent) {
      const ws = shuffle(others); if (out !== 'none' && x.kind === 'alcohol' && chance(0.3)) ws.push('Keine Reaktion');
      if (chance(0.5) && out !== 'none') {
        const opts = [out.sm, ...Object.entries(x.outs).filter(([r, o]) => r !== x.r && o !== 'none' && o.can !== out.can).map(([r, o]) => o.sm)];
        const uniq = [...new Map(opts.map(s => [can(s), s])).values()].slice(0, 4);
        if (uniq.length >= 3) { const o = shuffle(uniq.map(s => mol(s))); if (!o.includes(-1)) { if (put('rp', { c: jit(RC[x.r] - 30, 60), t: T.rx, k: 'mm', q: 'Welches ist das Hauptprodukt?', m: mi, rg: RG[x.r], o, a: o.indexOf(mol(out.sm)), e: RWHY[x.r] })) made++; continue; } }
      }
      if (mc('rx', { c: jit(RC[x.r] + (out === 'none' ? 120 : 0), 60), t: T.rx, q: 'Welches ist das Hauptprodukt?', m: mi, rg: RG[x.r], e: RWHY[x.r] }, ansName, ws)) made++;
    } else {
      if (out === 'none') continue;
      const pm = mol(out.sm); if (pm < 0) continue;
      const wrong = Object.entries(x.outs).filter(([r, o]) => r !== x.r && (o === 'none' || o.can !== out.can)).map(([r]) => RG[r]);
      mc('rg', { c: jit(RC[x.r] + 160, 60), t: T.rx, q: 'Mit welchem Reagenz gelingt diese Umsetzung?', m: mi, m2: pm, e: RWHY[x.r] }, RG[x.r], shuffle(wrong));
    }
    RXN.push(x);
  }
}

/* ---------- 9. two steps, Grignard disconnections, reagent sequences ---------- */
const CHAINS = [['alkene', 'hbr', 'nacn'], ['alkene', 'hbr', 'naoet'], ['alkene', 'bh3', 'pcc'], ['alkene', 'bh3', 'cr'], ['alkene', 'hyd', 'cr'], ['alkene', 'hbrp', 'nacn'], ['alkene', 'hbrp', 'mgco2'], ['alcohol', 'pbr3', 'nacn'], ['alcohol', 'pbr3', 'mgco2'], ['alcohol', 'pcc', 'memgbr'], ['alcohol', 'pcc', 'etmgbr'], ['carbonyl', 'nabh4', 'pbr3'], ['carbonyl', 'memgbr', 'dehyd'], ['acid', 'lahacid', 'pcc'], ['halide', 'nacn', 'lahcn'], ['halide', 'nacn', 'hydcn'], ['halide', 'naoet', 'bh3'], ['halide', 'kotbu', 'hbr'], ['nitrile', 'hydcn', 'lahacid'], ['acid', 'meoh', null]];
function run(g, seq) { let cur = g; for (const r of seq) { const p = R[r](cur); if (!p || p === 'none') return null; cur = p; } return cur; }
function famTwoStep(nProd, nSeq) {
  let a = 0, b = 0, t = 0;
  while ((a < nProd || b < nSeq) && t < (nProd + nSeq) * 60) {
    t++;
    const ch = pick(CHAINS.filter(c => c[2])); const s = substrate(ch[0]); if (!s) continue;
    const g = s.g, end = run(g, [ch[1], ch[2]]); if (!end) continue;
    const esm = smOf(end), ec = can(esm), en = nm(end); if (!ec || !en) continue;
    const ssm = smOf(g), mi = mol(ssm); if (mi < 0 || !nm(g)) continue;
    // alternatives: swap regiochem / other reagents in either step
    const alts = new Map();
    const swaps = { hbr: ['hbrp'], hbrp: ['hbr'], bh3: ['hyd'], hyd: ['bh3'], pcc: ['cr'], cr: ['pcc'], naoet: ['kotbu', 'naoh'], kotbu: ['naoet'], nacn: ['naoh', 'nh3'], mgco2: ['nacn'], memgbr: ['etmgbr'], etmgbr: ['memgbr'], lahcn: ['hydcn'], hydcn: ['lahcn'], nabh4: ['memgbr'], pbr3: ['dehyd'], dehyd: ['pbr3'], lahacid: ['meoh'] };
    const seqs = [];
    for (const x of swaps[ch[1]] || []) seqs.push([x, ch[2]]);
    for (const y of swaps[ch[2]] || []) seqs.push([ch[1], y]);
    seqs.push([ch[2], ch[1]]);
    for (const sq of seqs) { const p = run(g, sq); const key = sq.join('>'); if (p) { const c = can(smOf(p)); if (c && c !== ec) alts.set(key, { sq, p, c, name: nm(p) }); } else alts.set(key, { sq, p: null }); }
    const altList = [...alts.values()];
    if (a < nProd && chance(0.4)) {
      const pm = mol(esm);
      if (pm >= 0 && put('tt', { c: jit(2200 + (ch.includes('mgco2') || ch.includes('memgbr') || ch.includes('kotbu') ? 150 : 0), 80), t: T.syn, k: 'txt', q: 'Benenne das Produkt nach beiden Schritten.', m: mi, rg: `1) ${RG[ch[1]]}  →  2) ${RG[ch[2]]}`, mt: pm, ans: en })) a++;
      continue;
    }
    if (a < nProd && chance(0.55)) {
      const ws = altList.filter(x => x.p && x.name).map(x => x.name);
      const mid = run(g, [ch[1]]); if (mid) { const mn = nm(mid); if (mn && can(smOf(mid)) !== ec) ws.push(mn); }
      if (mc('tp', { c: jit(2050 + (ch.includes('mgco2') || ch.includes('memgbr') || ch.includes('kotbu') ? 150 : 0), 80), t: T.syn, q: 'Welches Produkt entsteht nach beiden Schritten?', m: mi, rg: `1) ${RG[ch[1]]}  →  2) ${RG[ch[2]]}` }, en, shuffle(ws))) a++;
    } else if (b < nSeq) {
      const pm = mol(esm); if (pm < 0) continue;
      const fmtSeq = sq => `1) ${RG[sq[0]]}  →  2) ${RG[sq[1]]}`;
      const ws = altList.map(x => fmtSeq(x.sq));
      if (mc('sq', { c: jit(2350 + (ch.includes('mgco2') || ch.includes('kotbu') ? 150 : 0) + (ch[0] === 'alkene' ? 60 : 0), 100), t: T.syn, q: 'Mit welcher Reagenzienfolge gelingt die Umwandlung?', m: mi, m2: pm, e: `${RWHY[ch[1]]} ${RWHY[ch[2]]}` }, fmtSeq([ch[1], ch[2]]), shuffle(ws))) b++;
    }
  }
}
const GRIG = { C: 'Methylmagnesiumbromid', CC: 'Ethylmagnesiumbromid', CCC: 'Propylmagnesiumbromid', 'C(C)C': 'Isopropylmagnesiumbromid' };
function famGrignardRetro(n) {
  let k = 0, t = 0;
  while (k < n && t < n * 50) {
    t++;
    const b = substrate('carbonyl'); if (!b) continue; const cn = nm(b.g); if (!cn) continue;
    const R0 = pick(Object.keys(GRIG)); const p = (() => { const x = carbonyl(b.g); if (!x) return null; const h = clone(b.g); C.setOrder(h, x.c, x.o, 1); C.attach(h, x.c, R0); return valOk(h) ? h : null; })(); if (!p) continue;
    const pc = can(smOf(p)), pn = nm(p); if (!pc || !pn) continue;
    // wrong pairs: other Grignard with same carbonyl, same Grignard with isomeric carbonyls
    const wrong = [];
    for (const R1 of Object.keys(GRIG)) { if (R1 === R0) continue; const x = carbonyl(b.g); const h = clone(b.g); C.setOrder(h, x.c, x.o, 1); C.attach(h, x.c, R1); const c1 = can(smOf(h)); if (c1 && c1 !== pc) wrong.push(`${GRIG[R1]} + ${cn}`); }
    for (let tt = 0; tt < 8; tt++) { const b2 = substrate('carbonyl'); if (!b2) continue; const n2 = nm(b2.g); if (!n2 || n2 === cn) continue; const x = carbonyl(b2.g); const h = clone(b2.g); C.setOrder(h, x.c, x.o, 1); C.attach(h, x.c, R0); const c1 = can(smOf(h)); if (c1 && c1 !== pc) wrong.push(`${GRIG[R0]} + ${n2}`); }
    const pm = mol(smOf(p)); if (pm < 0) continue;
    if (mc('gr', { c: jit(2150 + (C.hcount(p, carbonyl(b.g).c) === 0 ? 120 : 0), 80), t: T.syn, q: `Welche Kombination liefert nach wässriger Aufarbeitung <b>${pn}</b>?`, m: pm, e: 'Das C-Atom, das die OH-Gruppe trägt, war das Carbonyl-C; einer seiner Reste kam aus dem Grignard-Reagenz.' }, `${GRIG[R0]} + ${cn}`, shuffle(wrong))) k++;
  }
}
function famSNStereo(n) {
  const NU = [['C#N', 'NaCN in DMSO'], ['[N-]=[N+]=[N-]', 'NaN₃ in DMF'], ['I', 'NaI in Aceton'], ['O', 'NaOH (aq, verd.)'], ['SC', 'NaSCH₃ in DMF']];
  let k = 0, t = 0;
  while (k < n && t < n * 40) {
    t++;
    const r1 = pick(['C', 'CC', 'CCC']), r2 = pick(['CC', 'CCC', 'C(C)C', 'CCCC']); if (r1 === r2) continue;
    const X = pick(['Br', 'Cl', 'I']); const [nu, cond] = pick(NU); if (nu === X) continue;
    if (nu.startsWith('[N-]')) continue; // azide names are outside the parser
    const tag = pick(['@', '@@']); const inv = tag === '@' ? '@@' : '@';
    const sub = `${r1}[C${tag}H](${X})${r2}`, prod = `${r1}[C${inv}H](${nu})${r2}`;
    let ts, tp; try { ts = C.stereoTags(sub).CIP_atoms; tp = C.stereoTags(prod).CIP_atoms; } catch (e) { continue; }
    if (ts.length !== 1 || tp.length !== 1) continue;
    const ls = ts[0][1].slice(1, 2), lp = tp[0][1].slice(1, 2);
    const flat = s => s.replace(/@/g, '').replace('[CH]', 'C');
    const sn = nm(C.parseSmiles(flat(sub))), pn = nm(C.parseSmiles(flat(prod))); if (!sn || !pn) continue;
    const other = lp === 'R' ? 'S' : 'R';
    const q = `(${ls})-${sn} reagiert mit ${cond}. Welches Produkt entsteht?`;
    if (mc('ss', { c: jit(2200 + (ls === lp ? 0 : 150), 80), t: T.st, q, e: `SN2 → Walden-Umkehr am Stereozentrum. Der Deskriptor ergibt sich danach neu aus den CIP-Prioritäten: hier (${lp}).` }, `(${lp})-${pn}`, [`(${other})-${pn}`, `racemisches ${pn}`, `(${ls === lp ? other : ls})-${sn}`])) k++;
  }
  // SN1 at tertiary centres -> racemate
  for (const [s, p] of [['3-Brom-3-methylhexan', '3-Methylhexan-3-ol'], ['3-Chlor-3-methylheptan', '3-Methylheptan-3-ol'], ['3-Brom-3-methylheptan', '3-Methylheptan-3-ol'], ['4-Brom-4-methyloctan', '4-Methyloctan-4-ol']])
    for (const L of ['R', 'S']) mc('ss', { c: jit(2350, 60), t: T.st, q: `(${L})-${s} wird in Wasser/Aceton erwärmt (Solvolyse). Welches Substitutionsprodukt entsteht?`, e: 'SN1: Das planare Carbokation wird von beiden Seiten angegriffen → (nahezu) racemisches Produkt.' }, `racemisches ${p}`, [`(R)-${p}`, `(S)-${p}`, `(${L})-${p} (Retention)`]);
}

/* ---------- 9b. three-step routes and "which starting material?" ---------- */
const KIND_OF = g => { if (alkeneAtoms(g)) return 'alkene'; if (acid(g)) return 'acid'; if (ohCarbon(g)) return 'alcohol'; if (carbonyl(g)) return 'carbonyl'; if (halide(g)) return 'halide'; if (nitrile(g)) return 'nitrile'; return null; };
function walk(g, steps) {
  const seq = []; let cur = g;
  for (let k = 0; k < steps; k++) {
    const kind = KIND_OF(cur); if (!kind) return null;
    const opts = shuffle(APPLIES[kind].filter(r => !seq.includes(r)));
    let ok = false;
    for (const r of opts) { const p = R[r](cur); if (!p || p === 'none' || !valOk(p) || p.length > 14) continue; const n = nm(p); if (!n) continue; seq.push(r); cur = p; ok = true; break; }
    if (!ok) return null;
  }
  return { seq, end: cur };
}
function famThreeStep(nProd, nSeq, steps = 3) {
  let a = 0, b = 0, t = 0;
  const swaps = { hbr: ['hbrp'], hbrp: ['hbr'], bh3: ['hyd'], hyd: ['bh3'], pcc: ['cr'], cr: ['pcc'], naoet: ['kotbu'], kotbu: ['naoet'], nacn: ['naoh'], naoh: ['nacn'], mgco2: ['nacn'], memgbr: ['etmgbr'], etmgbr: ['memgbr'], lahcn: ['hydcn'], hydcn: ['lahcn'], pbr3: ['dehyd'], dehyd: ['pbr3'], h2: ['br2'], br2: ['h2'], kmno4c: ['hyd'], nh3: ['naoh'], meoh: ['etoh'], etoh: ['meoh'], lahacid: ['meoh'], nabh4: ['memgbr'], hcl: ['hbr'] };
  while ((a < nProd || b < nSeq) && t < (nProd + nSeq) * 80) {
    t++;
    const s = substrate(pick(Object.keys(APPLIES))); if (!s) continue;
    const w = walk(s.g, steps); if (!w) continue;
    const ec = can(smOf(w.end)), en = nm(w.end); if (!ec || !en) continue;
    const ssm = smOf(s.g), mi = mol(ssm); if (mi < 0 || !nm(s.g)) continue;
    const alts = [];
    for (let k = 0; k < steps; k++) for (const x of swaps[w.seq[k]] || []) { const sq = w.seq.slice(); sq[k] = x; alts.push(sq); }
    for (let k = 0; k + 1 < steps; k++) { const sq = w.seq.slice(); [sq[k], sq[k + 1]] = [sq[k + 1], sq[k]]; alts.push(sq); }
    const res = alts.map(sq => { const p = run(s.g, sq); return { sq, p, c: p ? can(smOf(p)) : null, name: p ? nm(p) : null }; }).filter(x => x.c !== ec);
    const hard = w.seq.filter(r => ['kotbu', 'bh3', 'hbrp', 'mgco2', 'memgbr', 'etmgbr', 'lahcn'].includes(r)).length;
    const fmtSeq = sq => sq.map((r, i) => `${i + 1}) ${RG[r]}`).join('  →  ');
    if (a < nProd && chance(0.5)) {
      const ws = res.filter(x => x.name).map(x => x.name);
      if (mc('t' + steps, { c: jit((steps === 3 ? 2420 : 2700) + 70 * hard, 90), t: T.syn, q: `Welches Produkt entsteht nach allen ${steps === 3 ? 'drei' : 'vier'} Schritten?`, m: mi, rg: fmtSeq(w.seq) }, en, shuffle(ws))) a++;
    } else if (b < nSeq) {
      const pm = mol(smOf(w.end)); if (pm < 0) continue;
      if (mc('s' + steps, { c: jit((steps === 3 ? 2620 : 2870) + 70 * hard, 110), t: T.syn, q: 'Mit welcher Reagenzienfolge gelingt die Umwandlung?', m: mi, m2: pm, e: w.seq.map(r => RWHY[r]).join(' ') }, fmtSeq(w.seq), shuffle(res.map(x => fmtSeq(x.sq))))) b++;
    }
  }
}
function famReverse(n) {
  let k = 0, t = 0;
  while (k < n && t < n * 60) {
    t++;
    const kind = pick(Object.keys(APPLIES)), s = substrate(kind); if (!s) continue;
    const r = pick(APPLIES[kind]); const p = R[r](s.g); if (!p || p === 'none') continue;
    const pc = can(smOf(p)), pm = mol(smOf(p)); if (!pc || pm < 0 || !nm(p) || !nm(s.g)) continue;
    const sc = can(smOf(s.g));
    // distractor starting materials: isomers of the substrate that give a different product with the same reagent
    const ds = new Map();
    for (let tt = 0; tt < 30 && ds.size < 3; tt++) {
      const s2 = substrate(kind); if (!s2) continue; const c2 = can(smOf(s2.g)); if (!c2 || c2 === sc || ds.has(c2)) continue;
      if (C.formula(s2.g) !== C.formula(s.g) && chance(0.7)) continue;
      const p2 = R[r](s2.g); if (p2 && p2 !== 'none' && can(smOf(p2)) === pc) continue;
      if (mol(smOf(s2.g)) < 0) continue; ds.set(c2, smOf(s2.g));
    }
    if (ds.size < 3) continue;
    const o = shuffle([mol(smOf(s.g)), ...[...ds.values()].map(x => mol(x))]);
    if (put('re', { c: jit(RC[r] + 420, 70), t: T.rx, k: 'mm', q: 'Welches Edukt liefert mit diesem Reagenz das gezeigte Produkt?', m: pm, rg: RG[r], rev: 1, o, a: o.indexOf(mol(smOf(s.g))), e: RWHY[r] })) k++;
  }
}

/* ---------- 10. curly-arrow templates ---------- */
function famArrows() {
  // proton transfer
  const BASE = [['HO', '−', 3], ['H₃N', '', 1], ['H₂O', '', 2], ['CH₃O', '−', 3], ['H₂N', '−', 2], ['EtO', '−', 3]];
  const ACID = [['Cl'], ['Br'], ['I']];
  for (const [bt, bq, nl] of BASE) for (const [x] of ACID) {
    const lps = nl === 1 ? [0] : nl === 2 ? [-60, 60] : [-90, 90, 0];
    put('ar', { c: jit(bq ? 1020 : 1080, 40), t: T.mech, k: 'arr', q: `Zeichne die Elektronenpfeile für den Protonentransfer: ${bt}${bq ? '⁻' : ''} + H–${x}.`,
      tm: { atoms: [{ id: 'B', t: bt, x: 60, y: 80, q: bq, lp: lps }, { id: 'H', t: 'H', x: 160, y: 80 }, { id: 'X', t: x, x: 220, y: 80, lp: [-90, 90, 0] }], bonds: [{ id: 'hx', a: 'H', b: 'X' }] },
      sol: [['lp:B>a:H', 'b:hx>a:X']], e: 'Das freie Elektronenpaar der Base greift das H an; die Bindungselektronen H–X wandern zum X.' });
  }
  // SN2
  const NU = [['HO', '−'], ['N≡C', '−'], ['I', '−'], ['HS', '−'], ['CH₃O', '−'], ['CH₃S', '−'], ['H₃N', '']];
  const RS = [['H', 'H', 'H'], ['H', 'H', 'CH₃'], ['H', 'H', 'C₂H₅'], ['H', 'H', 'C₆H₅']];
  const LG = ['Cl', 'Br', 'I', 'OTs'];
  for (const [nt, nq] of NU) for (const rs of RS) for (const lg of LG) {
    if (nt === 'I' && lg === 'I') continue;
    if (!chance(0.45)) continue;
    put('ar', { c: jit(1320 + (rs[2] !== 'H' ? 40 : 0) + (lg === 'OTs' ? 60 : 0), 40), t: T.mech, k: 'arr', q: `SN2: ${nt}${nq ? '⁻' : ''} reagiert mit dem gezeigten Substrat. Zeichne die Elektronenpfeile.`,
      tm: { atoms: [{ id: 'N', t: nt, x: 50, y: 90, q: nq, lp: nq ? [-90, 90, 0] : [0] }, { id: 'C', t: 'C', x: 150, y: 90 }, { id: 'a', t: rs[0], x: 150, y: 44 }, { id: 'b', t: rs[1], x: 122, y: 128 }, { id: 'c', t: rs[2], x: 180, y: 128 }, { id: 'L', t: lg, x: 222, y: 90, lp: lg === 'OTs' ? [] : [-90, 90, 0] }],
        bonds: [{ id: 'ca', a: 'C', b: 'a' }, { id: 'cb', a: 'C', b: 'b' }, { id: 'cc', a: 'C', b: 'c' }, { id: 'cl', a: 'C', b: 'L' }] },
      sol: [['lp:N>a:C', 'b:cl>a:L']], e: 'Rückseitenangriff: freies Elektronenpaar des Nucleophils → C; gleichzeitig wandert das Bindungselektronenpaar C–X zur Abgangsgruppe.' });
  }
  // nucleophilic addition to C=O
  const CN = [['H', '−', 'Hydrid (aus NaBH₄)', 1], ['N≡C', '−', 'Cyanid', 3], ['HO', '−', 'Hydroxid', 3], ['CH₃O', '−', 'Methanolat', 3]];
  const CO = [['H', 'H'], ['H', 'CH₃'], ['CH₃', 'CH₃'], ['H', 'C₆H₅'], ['CH₃', 'C₂H₅'], ['H', 'C₂H₅']];
  for (const [nt, nq, nn, nl] of CN) for (const [r1, r2] of CO) {
    put('ar', { c: jit(1420 + (nt === 'H' ? 60 : 0), 40), t: T.mech, k: 'arr', q: `Nucleophile Addition: ${nn} greift die Carbonylgruppe an. Zeichne die Elektronenpfeile.`,
      tm: { atoms: [{ id: 'N', t: nt, x: 50, y: 110, q: nq, lp: nl === 1 ? [0] : [-90, 90, 0] }, { id: 'C', t: 'C', x: 150, y: 110 }, { id: 'O', t: 'O', x: 150, y: 50, lp: [-150, -30] }, { id: 'a', t: r1, x: 190, y: 146 }, { id: 'b', t: r2, x: 112, y: 152 }],
        bonds: [{ id: 'co', a: 'C', b: 'O', o: 2 }, { id: 'ca', a: 'C', b: 'a' }, { id: 'cb', a: 'C', b: 'b' }] },
      sol: [['lp:N>a:C', 'b:co>a:O']], e: 'Das Nucleophil greift das elektrophile Carbonyl-C an; die π-Elektronen der C=O-Bindung klappen zum Sauerstoff (Alkoxid).' });
  }
  // Grignard addition (bond as source)
  for (const [r1, r2] of CO) put('ar', { c: jit(1560, 40), t: T.mech, k: 'arr', q: 'Grignard-Addition: CH₃MgBr greift die Carbonylgruppe an. Zeichne die Elektronenpfeile.',
    tm: { atoms: [{ id: 'Mg', t: 'MgBr', x: 30, y: 60 }, { id: 'Me', t: 'H₃C', x: 70, y: 110 }, { id: 'C', t: 'C', x: 160, y: 110 }, { id: 'O', t: 'O', x: 160, y: 50, lp: [-150, -30] }, { id: 'a', t: r1, x: 200, y: 146 }, { id: 'b', t: r2, x: 124, y: 152 }],
      bonds: [{ id: 'mg', a: 'Mg', b: 'Me' }, { id: 'co', a: 'C', b: 'O', o: 2 }, { id: 'ca', a: 'C', b: 'a' }, { id: 'cb', a: 'C', b: 'b' }] },
    sol: [['b:mg>a:C', 'b:co>a:O']], e: 'Die polare C–Mg-Bindung liefert das nucleophile C: Bindungselektronen C–Mg → Carbonyl-C; π-Elektronen C=O → O.' });
  // E2
  const BE = [['HO', 'Hydroxid'], ['EtO', 'Ethanolat'], ['tBuO', 'tert-Butanolat']];
  const LE = ['Br', 'Cl', 'I', 'OTs'];
  const SE = ['CH₃', 'H', 'C₂H₅'];
  for (const [bt, bn] of BE) for (const lg of LE) for (const se of SE) {
    put('ar', { c: jit(1640 + (lg === 'OTs' ? 40 : 0), 40), t: T.mech, k: 'arr', q: `E2: ${bn} eliminiert H–${lg === 'OTs' ? 'OTs' : lg}. Zeichne alle Elektronenpfeile.`,
      tm: { atoms: [{ id: 'O', t: bt, x: 64, y: 40, q: '−', lp: [-90, 90, 8] }, { id: 'Hb', t: 'H', x: 150, y: 62 }, { id: 'Cb', t: 'C', x: 150, y: 110 }, { id: 'h2', t: 'H', x: 106, y: 110 }, { id: 'h3', t: 'H', x: 150, y: 158 }, { id: 'Ca', t: 'C', x: 210, y: 110 }, { id: 'Ha', t: 'H', x: 210, y: 62 }, { id: 'Me', t: se, x: 268, y: 110 }, { id: 'L', t: lg, x: 210, y: 162, lp: lg === 'OTs' ? [] : [180, 0, 90] }],
      bonds: [{ id: 'chb', a: 'Cb', b: 'Hb' }, { id: 'c2', a: 'Cb', b: 'h2' }, { id: 'c3', a: 'Cb', b: 'h3' }, { id: 'cc', a: 'Cb', b: 'Ca' }, { id: 'cha', a: 'Ca', b: 'Ha' }, { id: 'cme', a: 'Ca', b: 'Me' }, { id: 'cl', a: 'Ca', b: 'L' }] },
      sol: [['lp:O>a:Hb', 'b:chb>b:cc', 'b:cl>a:L']], e: 'Base → β-H; die C–H-Bindungselektronen bilden die π-Bindung; die C–X-Bindung bricht heterolytisch.' });
  }
  // alkene + HX, step 1
  const AL = [['H', 'H', 'H', 'H'], ['CH₃', 'H', 'H', 'H'], ['CH₃', 'CH₃', 'H', 'H'], ['CH₃', 'H', 'CH₃', 'H'], ['C₂H₅', 'H', 'H', 'H'], ['CH₃', 'CH₃', 'CH₃', 'H']];
  for (const al of AL) for (const x of ['Br', 'Cl', 'I']) put('ar', { c: jit(1240, 40), t: T.mech, k: 'arr', q: `Elektrophile Addition von H–${x}, 1. Schritt: Zeichne die Elektronenpfeile für die Protonierung.`,
    tm: { atoms: [{ id: 'C1', t: 'C', x: 70, y: 100 }, { id: 'C2', t: 'C', x: 130, y: 100 }, { id: 'a', t: al[0], x: 40, y: 64 }, { id: 'b', t: al[1], x: 40, y: 136 }, { id: 'c', t: al[2], x: 160, y: 136 }, { id: 'd', t: al[3], x: 160, y: 64 }, { id: 'H', t: 'H', x: 205, y: 40 }, { id: 'X', t: x, x: 262, y: 40, lp: [-90, 90, 0] }],
      bonds: [{ id: 'cc', a: 'C1', b: 'C2', o: 2 }, { id: 'ca', a: 'C1', b: 'a' }, { id: 'cb', a: 'C1', b: 'b' }, { id: 'c2c', a: 'C2', b: 'c' }, { id: 'cd', a: 'C2', b: 'd' }, { id: 'hx', a: 'H', b: 'X' }] },
    sol: [['b:cc>a:H', 'b:hx>a:X']], e: 'Die π-Elektronen greifen das Proton an, die H–X-Bindung bricht: Carbokation + Halogenid.' });
  // SN1 step 2 / carbocation capture
  const NUC = [['H₂O', '', [-60, 60]], ['Br', '−', [-90, 90, 0]], ['CH₃OH', '', [-60, 60]], ['Cl', '−', [-90, 90, 0]]];
  const CAT = [['CH₃', 'CH₃', 'CH₃'], ['CH₃', 'CH₃', 'C₂H₅'], ['CH₃', 'C₂H₅', 'C₆H₅'], ['CH₃', 'H', 'C₆H₅']];
  for (const [nt, nq, lp] of NUC) for (const ct of CAT) put('ar', { c: jit(1140, 40), t: T.mech, k: 'arr', q: `SN1, 2. Schritt: ${nt}${nq ? '⁻' : ''} fängt das Carbokation ab. Zeichne den Elektronenpfeil.`,
    tm: { atoms: [{ id: 'N', t: nt, x: 50, y: 100, q: nq, lp: lp.map(a => a + 0) }, { id: 'C', t: 'C', x: 160, y: 100, q: '+' }, { id: 'a', t: ct[0], x: 160, y: 50 }, { id: 'b', t: ct[1], x: 205, y: 128 }, { id: 'c', t: ct[2], x: 120, y: 140 }],
      bonds: [{ id: 'ca', a: 'C', b: 'a' }, { id: 'cb', a: 'C', b: 'b' }, { id: 'cc', a: 'C', b: 'c' }] },
    sol: [['lp:N>a:C']], e: 'Ein freies Elektronenpaar des Nucleophils bildet die neue Bindung zum positiv geladenen C-Atom.' });
}

/* ---------- 11. hand-written puzzles ---------- */
function famHand() {
  const map = { text: 'txt', choice: 'mc', order: 'ord', match: 'mat', click: 'pick', arrows: 'arr', synth: 'syn' };
  for (const h of HAND) {
    const p = { c: h.celo, t: { syn1: T.syn, syn2: T.syn, syn3: T.syn, syn4: T.syn, syn5: T.syn, syn6: T.syn, syn7: T.syn, syn8: T.syn }[h.id] || ({ text: T.nom, choice: T.rx, order: T.phys, match: T.rx, click: T.st, arrows: T.mech }[h.type]), k: map[h.type], q: h.task, e: h.explain, hint: h.hint };
    if (h.formula) p.f2 = h.formula;
    if (h.mol) p.tm = h.mol;
    if (h.type === 'text') { p.acc = h.accept; p.ans = h.explain.split(/[:(–]/)[0].trim(); const sm = { propan: 'CCC', '2mb': 'CC(C)CC', buten: 'CC=CC', dibr: 'BrCCBr', grig: 'CC(C)(C)O' }[h.id]; if (sm) { p.mt = mol(sm); p.ans = { propan: 'Propan', '2mb': '2-Methylbutan', buten: 'But-2-en', dibr: '1,2-Dibromethan', grig: '2-Methylpropan-2-ol' }[h.id]; } }
    if (h.type === 'choice') { p.o = h.options; p.a = h.correct[0]; }
    if (h.type === 'order') { p.items = h.items; p.ends = h.ends; }
    if (h.type === 'match') p.pairs = h.pairs;
    if (h.type === 'click') p.a = h.answer;
    if (h.type === 'arrows') p.sol = h.solutions;
    if (h.type === 'synth') p.syn = { start: h.start, target: h.target, nodes: h.nodes, reagents: h.reagents, edges: h.edges, maxSteps: h.maxSteps, minYield: h.minYield };
    // hand puzzles keep their ids
    KEYS.add('hand|' + h.id); p.i = h.id; p.c = clampC(p.c); P.push(p);
  }
  const t = { id: 'iupac' }; void t;
}
// fix theme for hand puzzles by content
function themeFix() {
  const th = { br: T.el, propan: T.nom, ord6: T.el, alkmatch: T.nom, fg: T.fg, sdp: T.phys, '2mb': T.nom, buten: T.nom, iso: T.str, rtyp: T.rx, dibr: T.rx, hcl: T.mech, acid: T.ab, markov: T.rx, stereo: T.st, sn1: T.mech, cip: T.st, sn2: T.mech, e2: T.mech, da: T.mech, nitro: T.rx, grig: T.rx, claisen: T.mech, sn2stereo: T.st };
  for (const p of P) if (th[p.i]) p.t = th[p.i];
}

/* ======================================================================
   RUN
   ====================================================================== */
(async () => {
  await C.initRDKit();
  famHand(); themeFix();
  famElements(); famBasics();
  fillPool(1500);
  const nr = famNaming();
  famFormula(300); famClass(); famHybrid(170); famOx(230);
  famStereoCenters(320); famEZ(260); famRS(300);
  famAcidBase(230); famBP(130); famMech();
  famReactions(700, 280); famTwoStep(300, 170); famThreeStep(150, 150); famThreeStep(70, 90, 4); famReverse(220); famGrignardRetro(130); famSNStereo(110);
  famArrows();

  // stats
  const bins = {}; for (const p of P) { const b = Math.floor(p.c / 100) * 100; bins[b] = (bins[b] || 0) + 1; }
  const mid = P.filter(p => p.c >= 900 && p.c <= 1900).length;
  console.log('puzzles', P.length, 'molecules', M.length, '900-1900:', mid, (100 * mid / P.length).toFixed(1) + '%');
  console.log(Object.entries(bins).map(([k, v]) => k + ':' + v).join('  '));
  console.log(Object.entries(COUNT).map(([k, v]) => k + ':' + v).join('  '), nr);
  if (FAILS.length) { console.log('round-trip/geometry rejections:', FAILS.length); console.log(FAILS.slice(0, 15).join('\n')); }

  const out = '/* Generated by tools/generate.cjs — do not edit by hand. */\nwindow.CP_DATA=' + JSON.stringify({ v: 1, M, P }) + ';\n';
  fs.writeFileSync(path.join(__dirname, '..', 'puzzles.js'), out);
  // stamp the data version into the page and the service worker so caches stay consistent
  const ver = require('crypto').createHash('sha1').update(out).digest('hex').slice(0, 10);
  const ver2 = require('crypto').createHash('sha1').update(out + fs.readFileSync(path.join(__dirname, '..', 'namer-de.js'), 'utf8') + fs.readFileSync(path.join(__dirname, '..', 'index.html'), 'utf8').replace(/\?v=[\w-]+/g, '')).digest('hex').slice(0, 10);
  for (const [f, re, rep] of [['index.html', /(puzzles|namer-de)\.js\?v=[\w-]+/g, '$1.js?v=' + ver2], ['sw.js', /(puzzles|namer-de)\.js\?v=[\w-]+/g, '$1.js?v=' + ver2], ['sw.js', /VERSION = '[^']+'/, "VERSION = 'cp-" + ver2 + "'"]]) {
    const fp = path.join(__dirname, '..', f); if (fs.existsSync(fp)) fs.writeFileSync(fp, fs.readFileSync(fp, 'utf8').replace(re, rep));
  }
  console.log('data version', ver);
  console.log('wrote puzzles.js', (out.length / 1024).toFixed(0) + ' KB');
})();
