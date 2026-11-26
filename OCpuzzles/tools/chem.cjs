/* Chem Puzzles build helpers: molecular graphs, SMILES I/O, naming (German),
 * an independent German name parser for round-trip checks, and RDKit depictions.
 * Dev-time only. The page never loads this file.
 */
'use strict';
const path = require('path');
const { name: nameEN } = require('./native-namer.js');

const RDKIT_JS = process.env.RDKIT_JS ||
  path.resolve(__dirname, '../../../naming-engine/src/vendor/rdkit/dist/RDKit_minimal.js');
let RD = null;
async function initRDKit() { if (!RD) RD = await require(RDKIT_JS)(); return RD; }

/* ---------------- graphs ---------------- */
// atom: {el, nb:[{n,o}], ch?}
const VAL = { C: 4, N: 3, O: 2, F: 1, Cl: 1, Br: 1, I: 1, S: 2 };
function G() { return []; }
function add(g, el, ch = 0) { g.push({ el, nb: [], ch }); return g.length - 1; }
function bond(g, a, b, o = 1) { g[a].nb.push({ n: b, o }); g[b].nb.push({ n: a, o }); }
function hcount(g, i) {
  const a = g[i], s = a.nb.reduce((t, e) => t + e.o, 0);
  let v = VAL[a.el]; if (a.el === 'N' && a.ch === 1) v = 4; if (a.el === 'O' && a.ch === -1) v = 1; if (a.el === 'O' && a.ch === 1) v = 3;
  return Math.max(0, v - s);
}
function chain(g, n) { const ids = []; for (let i = 0; i < n; i++) { ids.push(add(g, 'C')); if (i) bond(g, ids[i - 1], ids[i]); } return ids; }
function setOrder(g, a, b, o) { g[a].nb.find(e => e.n === b).o = o; g[b].nb.find(e => e.n === a).o = o; }

/* graph -> SMILES (input order is kept: atom 0 first, DFS) */
function smiles(g) {
  const seen = new Array(g.length).fill(false), ring = new Map(); let label = 1;
  // find ring closures with a first DFS
  const order = [], parent = new Array(g.length).fill(-1), closures = [];
  (function dfs(i) { seen[i] = true; order.push(i); for (const e of g[i].nb) { if (e.n === parent[i]) continue; if (seen[e.n]) { if (order.indexOf(e.n) < order.indexOf(i)) closures.push([e.n, i, e.o]); continue; } parent[e.n] = i; dfs(e.n); } })(0);
  if (order.length !== g.length) throw new Error('disconnected graph');
  const openAt = new Map(), closeAt = new Map();
  for (const [a, b, o] of closures) { const d = label++; (openAt.get(a) || openAt.set(a, []).get(a)).push([d, o]); (closeAt.get(b) || closeAt.set(b, []).get(b)).push([d, o]); }
  const sym = i => { const a = g[i]; const org = ['C', 'N', 'O', 'F', 'Cl', 'Br', 'I', 'S'].includes(a.el); if (!a.ch && org) return a.el; const h = hcount(g, i); return '[' + a.el + (h ? 'H' + (h > 1 ? h : '') : '') + (a.ch > 0 ? '+' : a.ch < 0 ? '-' : '') + ']'; };
  const bs = o => o === 2 ? '=' : o === 3 ? '#' : '';
  const done = new Array(g.length).fill(false);
  function write(i, from) {
    done[i] = true; let s = sym(i);
    for (const [d, o] of openAt.get(i) || []) s += bs(o) + (d < 10 ? d : '%' + d);
    for (const [d] of closeAt.get(i) || []) s += d < 10 ? d : '%' + d;
    const kids = g[i].nb.filter(e => e.n !== from && parent[e.n] === i);
    kids.forEach((e, k) => { const t = bs(e.o) + write(e.n, i); s += k < kids.length - 1 ? '(' + t + ')' : t; });
    return s;
  }
  return write(0, -1);
}
function parseSmiles(sm) {
  const g = [], stack = [], rings = new Map(); let at = -1, order = 1;
  const toks = sm.match(/\[[^\]]+\]|Cl|Br|[A-Z][a-z]?|[()=#]|%\d\d|\d/g);
  if (toks.join('') !== sm) throw new Error('unsupported SMILES ' + sm);
  for (const t of toks) {
    if (t === '(') { stack.push(at); continue; } if (t === ')') { at = stack.pop(); continue; }
    if (t === '=') { order = 2; continue; } if (t === '#') { order = 3; continue; }
    if (/^%?\d/.test(t)) { if (rings.has(t)) { const [o2, oo] = rings.get(t); bond(g, at, o2, Math.max(order, oo)); rings.delete(t); } else rings.set(t, [at, order]); order = 1; continue; }
    const el = t.match(/[A-Z][a-z]?/)[0], ch = t.includes('+') ? 1 : t.includes('-') ? -1 : 0;
    const i = add(g, el, ch); if (at >= 0) bond(g, at, i, order); at = i; order = 1;
  }
  return g;
}
const forNamer = g => g.map((a, i) => ({ el: a.el, h: hcount(g, i), charge: a.ch || 0, nb: a.nb.map(e => ({ n: e.n, o: e.o })) }));

/* ---------------- RDKit wrappers ---------------- */
function canon(sm) { const m = RD.get_mol(sm); if (!m || !m.is_valid()) throw new Error('RDKit rejects ' + sm); const s = m.get_smiles(); m.delete(); return s; }
function stereoTags(sm) { const m = RD.get_mol(sm); const t = JSON.parse(m.get_stereo_tags()); m.delete(); return t; }
function formula(g) {
  const c = {}; g.forEach((a, i) => { c[a.el] = (c[a.el] || 0) + 1; const h = hcount(g, i); if (h) c.H = (c.H || 0) + h; });
  const keys = Object.keys(c).filter(k => k !== 'C' && k !== 'H').sort();
  return ['C', 'H', ...keys].filter(k => c[k]).map(k => k + (c[k] > 1 ? c[k] : '')).join('');
}
/* 2D depiction from RDKit: atoms [el,x,y,ch,hLabel] bonds [a,b,order,stereo(1 wedge,6 hash)] */
function depict(sm) {
  const m = RD.get_mol(sm);
  if (!m || !m.is_valid()) throw new Error('RDKit rejects ' + sm);
  const mb = m.get_new_coords(true) || m.get_molblock();
  m.set_new_coords(true);
  try { m.straighten_depiction(); } catch (e) {}
  const block = m.get_molblock(); m.delete();
  const L = block.split('\n'); const [na, nbd] = [parseInt(L[3].slice(0, 3)), parseInt(L[3].slice(3, 6))];
  const atoms = [], bonds = [];
  for (let i = 0; i < na; i++) { const s = L[4 + i]; atoms.push([s.slice(31, 34).trim(), +s.slice(0, 10), -(+s.slice(10, 20))]); }
  for (let i = 0; i < nbd; i++) { const s = L[4 + na + i]; bonds.push([+s.slice(0, 3) - 1, +s.slice(3, 6) - 1, +s.slice(6, 9), +s.slice(9, 12)]); }
  for (const l of L) if (l.startsWith('M  CHG')) { const p = l.slice(6).trim().split(/\s+/).map(Number); for (let k = 0; k < p[0]; k++) atoms[p[1 + 2 * k] - 1][3] = p[2 + 2 * k]; }
  // normalise: bond length 1, centre 0, 2 decimals
  let sum = 0; bonds.forEach(b => sum += Math.hypot(atoms[b[0]][1] - atoms[b[1]][1], atoms[b[0]][2] - atoms[b[1]][2]));
  const k = bonds.length ? bonds.length / sum : 1;
  const xs = atoms.map(a => a[1]), ys = atoms.map(a => a[2]);
  const cx = (Math.min(...xs) + Math.max(...xs)) / 2, cy = (Math.min(...ys) + Math.max(...ys)) / 2;
  atoms.forEach(a => { a[1] = Math.round((a[1] - cx) * k * 100) / 100; a[2] = Math.round((a[2] - cy) * k * 100) / 100; });
  // overlap check
  for (let i = 0; i < atoms.length; i++) for (let j = i + 1; j < atoms.length; j++)
    if (Math.hypot(atoms[i][1] - atoms[j][1], atoms[i][2] - atoms[j][2]) < 0.55) return null;
  return { a: atoms.map(a => a[3] ? a : a.slice(0, 3)), b: bonds.map(b => b[3] ? b : b.slice(0, 3)) };
}

/* ---------------- names ---------------- */
function toGerman(en) {
  let s = en.replace(/benzen-1-amine$/, 'aniline').replace(/benzene-1-carboxylic acid$/, 'benzoic acid').replace(/benzen-1-ol$/, 'phenol').replace(/benzene-1-carbaldehyde$/, 'benzaldehyde').replace(/benzene-1-carbonitrile$/, 'benzonitrile');
  const special = { 'benzoic acid': 'benzoesäure', 'benzaldehyde': 'benzaldehyd', 'benzene': 'benzol', 'phenol': 'phenol', 'benzonitrile': 'benzonitril', 'aniline': 'anilin' };
  for (const k in special) s = s.replace(new RegExp(k + '$'), special[k]);
  s = s.replace(/fluoro/g, 'fluor').replace(/chloro/g, 'chlor').replace(/bromo/g, 'brom').replace(/iodo/g, 'iod');
  s = s.replace(/oic acid$/, 'säure').replace(/enitrile$/, 'nitril').replace(/nitrile$/, 'nitril')
    .replace(/amine$/, 'amin').replace(/amide$/, 'amid').replace(/oate$/, 'oat');
  s = s.replace(/idene/g, 'iden').replace(/-yn-/g, '-in-').replace(/-en-(\d)/g, '-en-$1');
  s = s.replace(/benzene/g, 'benzol').replace(/ane\b/g, 'an').replace(/ene\b/g, 'en').replace(/yne\b/g, 'in').replace(/one\b/g, 'on');
  s = s.replace(/^(\w+) (\S+oat)$/, (m, a, b) => a + (/^\d/.test(b) ? '-' : '') + b);            // ester: "ethyl ethanoate" -> "ethylethanoat"
  if (/ /.test(s)) return null;
  return s.replace(/^((?:\(\d?[EZRS](?:,\d?[EZRS])*\)-)?(?:[\d,]+-)?(?:N(?:,N)?-)?\(?)([a-z])/, (m, p, c) => p + c.toUpperCase());
}
function nameDE(g) {
  const r = nameEN(forNamer(g));
  if (r.status === 'unsupported' || !r.name) return null;
  if (r.coveredAtoms !== undefined && r.coveredAtoms !== g.length) return null;
  return toGerman(r.name);
}
/* Traditional German school form: "Butan-2-ol" -> "2-Butanol", "2-Methylpent-1-en" -> "2-Methyl-1-penten" */
function schoolForm(n) {
  const m = n.match(/^(.*?)((?:[Cc]yclo)?(?:[Mm]eth|[Ee]th|[Pp]rop|[Bb]ut|[Pp]ent|[Hh]ex|[Hh]ept|[Oo]ct|[Nn]on|[Dd]ec)(?:an|a)?)-([\d,]+)-((?:di|tri)?(?:ol|on|en|in|amin))$/);
  if (!m) return null;
  const [, pre, stem, locs, suf] = m;
  const word = (stem + suf).toLowerCase();
  return pre ? pre + '-' + locs + '-' + word : locs + '-' + word.charAt(0).toUpperCase() + word.slice(1);
}

/* ---------------- independent German name parser ---------------- */
const STEMS = ['', 'Meth', 'Eth', 'Prop', 'But', 'Pent', 'Hex', 'Hept', 'Oct', 'Non', 'Dec'];
const MULT = { di: 2, tri: 3, tetra: 4 };
const PREFIX = { methyl: 'C', ethyl: 'CC', propyl: 'CCC', butyl: 'CCCC', 'propan-2-yl': 'C(C)C', fluor: 'F', chlor: 'Cl', brom: 'Br', iod: 'I', hydroxy: 'O', amino: 'N', oxo: '=O', methoxy: 'OC', ethoxy: 'OCC', cyano: 'C#N', nitro: 'N(=O)=O', phenyl: 'c' };
function parseGerman(nm) {
  // returns graph or throws
  let s = nm.toLowerCase().replace(/\((propan-2-yl)\)/g, '$1');
  const g = G();
  const pre = [];
  // ester: alkyl + alkanoat
  let esterAlkyl = null;
  const em = s.match(/^(methyl|ethyl|propyl|butyl)(.+oat)$/); if (em) { esterAlkyl = em[1]; s = em[2].replace(/^-/, ''); }
  // prefixes
  const prefRe = /^(?:(\d+(?:,\d+)*)-)?(di|tri|tetra)?(propan-2-yl|methyl|ethyl|propyl|butyl|fluor|chlor|brom|iod|hydroxy|amino|oxo|methoxy|ethoxy|cyano|nitro|phenyl)/;
  let m;
  while ((m = s.match(prefRe))) {
    const locs = m[1] ? m[1].split(',').map(Number) : [null];
    const k = m[2] ? MULT[m[2]] : 1;
    if (locs.length !== k && !(locs[0] === null && k === 1)) throw new Error('multiplier mismatch ' + nm);
    locs.forEach(l => pre.push([l, m[3]]));
    s = s.slice(m[0].length).replace(/^-(?=\d)/, '');
    if (!s) throw new Error('no parent ' + nm);
  }
  let atoms = null, ring = false;
  // parents
  if (/^benzol/.test(s) || /^phenol$/.test(s) || /^benzoesäure$/.test(s) || /^benzaldehyd$/.test(s) || /^benzonitril$/.test(s) || /^anilin$/.test(s)) {
    ring = true; atoms = [];
    for (let i = 0; i < 6; i++) atoms.push(add(g, 'C'));
    for (let i = 0; i < 6; i++) bond(g, atoms[i], atoms[(i + 1) % 6], i % 2 ? 1 : 2);
    if (s === 'phenol') { bond(g, atoms[0], add(g, 'O')); }
    if (s === 'anilin') { bond(g, atoms[0], add(g, 'N')); }
    if (s === 'benzoesäure') { const c = add(g, 'C'); bond(g, atoms[0], c); bond(g, c, add(g, 'O'), 2); bond(g, c, add(g, 'O')); }
    if (s === 'benzaldehyd') { const c = add(g, 'C'); bond(g, atoms[0], c); bond(g, c, add(g, 'O'), 2); }
    if (s === 'benzonitril') { const c = add(g, 'C'); bond(g, atoms[0], c); bond(g, c, add(g, 'N'), 3); }
    s = '';
  } else {
    const cyc = s.startsWith('cyclo'); if (cyc) s = s.slice(5);
    const st = STEMS.map(x => x.toLowerCase()).map((x, i) => [x, i]).filter(([x]) => x && s.startsWith(x)).sort((a, b) => b[0].length - a[0].length)[0];
    if (!st) throw new Error('no stem ' + nm);
    const n = st[1]; s = s.slice(st[0].length);
    atoms = chain(g, n); if (cyc) { bond(g, atoms[n - 1], atoms[0]); ring = true; }
    // unsaturation: "an" | "a-1,3-dien" | "-2-en" | "-1-in" | "-1-en-3-in"
    let un = s.match(/^a?-(\d+(?:,\d+)*)-(di|tri)?(en|in)(?:-(\d+(?:,\d+)*)-(di|tri)?(in))?/);
    if (un) {
      const apply = (locs, t) => locs.split(',').map(Number).forEach(l => setOrder(g, atoms[l - 1], atoms[l % n], t === 'en' ? 2 : 3));
      apply(un[1], un[3]); if (un[4]) apply(un[4], un[6]);
      s = s.slice(un[0].length);
    } else if (s.startsWith('an')) s = s.slice(2);
    else if (s.startsWith('en') && n === 2) { setOrder(g, atoms[0], atoms[1], 2); s = s.slice(2); }
    else if (s.startsWith('in') && n === 2) { setOrder(g, atoms[0], atoms[1], 3); s = s.slice(2); }
    else throw new Error('no saturation ' + nm);
    // suffix
    const sf = s.match(/^(?:-(\d+(?:,\d+)*)-)?(di|tri)?(ol|on|al|säure|nitril|amin|amid|oat)?$/);
    if (!sf) throw new Error('bad suffix ' + nm + ' / ' + s);
    const kind = sf[3], k = sf[2] ? MULT[sf[2]] : 1;
    let locs = sf[1] ? sf[1].split(',').map(Number) : null;
    if (kind && !locs) { if (['al', 'säure', 'nitril', 'amid', 'oat'].includes(kind)) locs = k === 2 ? [1, n] : [1]; else if (n <= 2 || ring) locs = [1]; else throw new Error('missing locant ' + nm); }
    if (locs && locs.length !== k) throw new Error('suffix multiplier mismatch ' + nm);
    for (const l of locs || []) {
      const c = atoms[l - 1];
      if (kind === 'ol') bond(g, c, add(g, 'O'));
      else if (kind === 'on') bond(g, c, add(g, 'O'), 2);
      else if (kind === 'al') bond(g, c, add(g, 'O'), 2);
      else if (kind === 'säure') { bond(g, c, add(g, 'O'), 2); bond(g, c, add(g, 'O')); }
      else if (kind === 'nitril') bond(g, c, add(g, 'N'), 3);
      else if (kind === 'amin') bond(g, c, add(g, 'N'));
      else if (kind === 'amid') { bond(g, c, add(g, 'O'), 2); bond(g, c, add(g, 'N')); }
      else if (kind === 'oat') { bond(g, c, add(g, 'O'), 2); const o = add(g, 'O'); bond(g, c, o); attach(g, o, { methyl: 'C', ethyl: 'CC', propyl: 'CCC', butyl: 'CCCC' }[esterAlkyl]); }
    }
  }
  for (const [l, p] of pre) {
    const c = l === null ? atoms[0] : atoms[l - 1];
    if (c === undefined) throw new Error('bad locant ' + nm);
    if (p === 'oxo') { bond(g, c, add(g, 'O'), 2); continue; }
    if (p === 'phenyl') { const r = []; for (let i = 0; i < 6; i++) r.push(add(g, 'C')); for (let i = 0; i < 6; i++) bond(g, r[i], r[(i + 1) % 6], i % 2 ? 1 : 2); bond(g, c, r[0]); continue; }
    attach(g, c, PREFIX[p]);
  }
  if (s && !ring && s.length) { /* handled */ }
  return g;
}
function attach(g, at, frag) {
  if (!frag) throw new Error('unknown fragment');
  const sub = parseSmiles(frag.startsWith('=') ? frag.slice(1) : frag);
  const off = g.length; sub.forEach(a => g.push({ el: a.el, ch: a.ch, nb: a.nb.map(e => ({ n: e.n + off, o: e.o })) }));
  bond(g, at, off, frag.startsWith('=') ? 2 : 1);
}
const NameDE = require('../namer-de.js');
function chargedNitro(g) { for (let i = 0; i < g.length; i++) if (g[i].el === 'N') { const dbl = g[i].nb.filter(e => g[e.n].el === 'O' && e.o === 2); if (dbl.length === 2) { g[i].ch = 1; const o = dbl[1].n; dbl[1].o = 1; g[o].nb.find(e => e.n === i).o = 1; g[o].ch = -1; } } return g; }
function readGraph(nm) { return chargedNitro(NameDE.parse(nm).g); }
function roundTrip(nm, sm) {
  let g; try { g = readGraph(nm); } catch (e) { return 'parse: ' + e.message + ' :: ' + nm; }
  let a, b; try { a = canon(smiles(g)); b = canon(sm); } catch (e) { return 'rdkit: ' + e.message; }
  return a === b ? null : 'mismatch ' + nm + ' -> ' + a + ' vs ' + b;
}

module.exports = { readGraph, chargedNitro, attach, VAL, initRDKit, G, add, bond, chain, setOrder, hcount, smiles, parseSmiles, canon, stereoTags, formula, depict, nameDE, toGerman, schoolForm, parseGerman, roundTrip, get RD() { return RD; } };
