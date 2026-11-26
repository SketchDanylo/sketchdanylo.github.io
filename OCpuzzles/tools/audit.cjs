/* Audit of the shipped puzzles.js, using only the shipped data:
 * drawings are turned back into molfiles for RDKit, names are parsed by the
 * independent German parser. Nothing here calls the generator or the namer.
 *   node OCpuzzles/tools/audit.cjs
 */
'use strict';
const fs = require('fs'), path = require('path');
const C = require('./chem.cjs');
const src = fs.readFileSync(path.join(__dirname, '..', 'puzzles.js'), 'utf8');
const window = {}; eval(src); const { M, P } = window.CP_DATA;

function molblock(m) {
  const [A, B] = m; const L = ['', '  audit', ''];
  L.push(String(A.length).padStart(3) + String(B.length).padStart(3) + '  0  0  0  0  0  0  0  0999 V2000');
  for (const a of A) L.push((a[1] / 100).toFixed(4).padStart(10) + (-a[2] / 100).toFixed(4).padStart(10) + '0.0000'.padStart(10) + ' ' + a[0].padEnd(3) + ' 0  0  0  0  0  0  0  0  0  0  0  0');
  for (const b of B) L.push(String(b[0] + 1).padStart(3) + String(b[1] + 1).padStart(3) + String(b[2]).padStart(3) + String(b[3] || 0).padStart(3));
  const ch = A.map((a, i) => [i + 1, a[3]]).filter(x => x[1]);
  if (ch.length) L.push('M  CHG' + String(ch.length).padStart(3) + ch.map(([i, c]) => String(i).padStart(4) + String(c).padStart(4)).join(''));
  L.push('M  END'); return L.join('\n');
}
const cache = new Map();
function rd(mi) {
  if (cache.has(mi)) return cache.get(mi);
  const m = C.RD.get_mol(molblock(M[mi]), JSON.stringify({ removeHs: true }));
  if (!m || !m.is_valid()) throw new Error('molfile invalid for molecule ' + mi);
  const r = { can: m.get_smiles(), tags: JSON.parse(m.get_stereo_tags()) };
  // a version without stereo, for constitution comparisons
  const flat = C.RD.get_mol(r.can.replace(/[@\/\\]/g, '').replace(/\[C(H?)\]/g, (x, h) => 'C'));
  r.flat = flat ? flat.get_smiles() : r.can; if (flat) flat.delete();
  m.delete(); cache.set(mi, r); return r;
}
const NameDE = require('../namer-de.js');
const canName = n => { try { return C.canon(C.smiles(C.readGraph(n))); } catch (e) { return null; } };
const plain = s => s.replace(/<[^>]+>/g, '');
// graph from depiction (for hybridisation / oxidation number / formula)
function graphOf(mi) { const [A, B] = M[mi]; const g = A.map(a => ({ el: a[0], ch: a[3] || 0, nb: [] })); for (const b of B) { g[b[0]].nb.push({ n: b[1], o: b[2] }); g[b[1]].nb.push({ n: b[0], o: b[2] }); } return g; }

const bad = [], stats = {};
const fail = (p, msg) => bad.push(p.i + ' [' + p.c + ']: ' + msg);
const ok = fam => stats[fam] = (stats[fam] || 0) + 1;

(async () => {
  await C.initRDKit();
  const ids = new Set();
  for (const p of P) {
    const fam = p.i.includes('-') ? p.i.split('-')[0] : p.i;
    try {
      if (ids.has(p.i)) fail(p, 'duplicate id'); ids.add(p.i);
      if (!(p.c >= 100 && p.c <= 3000)) fail(p, 'celo out of range');
      if (p.k === 'mc') { if (!Array.isArray(p.o) || p.o.length < 2 || p.a < 0 || p.a >= p.o.length) fail(p, 'bad options'); if (new Set(p.o).size !== p.o.length) fail(p, 'duplicate options'); }
      if (p.k === 'mm') { if (new Set(p.o).size !== p.o.length || p.o.some(x => x < 0 || x >= M.length)) fail(p, 'bad structure options'); const cs = p.o.map(x => rd(x).can); if (new Set(cs).size !== cs.length) fail(p, 'two options are the same molecule'); }
      const ans = p.k === 'mc' ? p.o[p.a] : null;

      if (fam === 'nm') { const c = canName(ans); if (c !== rd(p.m).can) fail(p, 'answer name does not match drawing: ' + ans);
        for (const w of p.o) if (w !== ans) { const cw = canName(w); if (cw === c) stats.misnumbered = (stats.misnumbered || 0) + 1; } }
      if (fam === 'tn') { if (canName(p.ans) !== rd(p.m).can) fail(p, 'typed answer does not match drawing: ' + p.ans); }
      if (p.k === 'txt' && p.mt !== undefined) { // the page's own check must accept the official answer and the school form
        const tg = NameDE.fromDepiction(M[p.mt]);
        for (const a of [p.ans, ...(p.acc || []).filter(x => /[a-z]/i.test(x) && x.length > 3 && !/^[a-z0-9]+$/.test(x))]) { const r = NameDE.read(a); if (!r.ok || !NameDE.same(r.g, tg)) fail(p, 'page reader rejects accepted answer ' + a); }
        if (canName(p.ans) !== rd(p.mt).flat && canName(p.ans) !== rd(p.mt).can) fail(p, 'typed answer vs target (RDKit): ' + p.ans);
      }
      if (fam === 'ft') { const f = C.formula(graphOf(p.m)); if (f !== p.ff) fail(p, 'typed formula ' + p.ff + ' vs ' + f); }
      if (fam === 'sm') { const n = plain(p.q).match(/ist (.+)\?$/)[1]; if (canName(n) !== rd(p.o[p.a]).can) fail(p, 'structure answer does not match ' + n); for (let k = 0; k < p.o.length; k++) if (k !== p.a && rd(p.o[k]).can === rd(p.o[p.a]).can) fail(p, 'distractor equals answer'); }
      if (['rx', 'tp', 't3', 't4'].includes(fam)) { const cs = p.o.filter(x => x !== 'Keine Reaktion').map(canName); if (cs.some(x => !x)) fail(p, 'unparseable option ' + p.o.join(' | ')); if (new Set(cs).size !== cs.length) fail(p, 'two options name the same molecule'); }
      if (fam === 'gr') { const n = plain(p.q).match(/liefert nach wässriger Aufarbeitung (.+)\?$/)[1]; if (canName(n) !== rd(p.m).can) fail(p, 'target name does not match drawing'); }
      if (fam === 'fo') { const g = graphOf(p.m); const f = C.formula(g).replace(/(\d+)/g, d => d.split('').map(x => '₀₁₂₃₄₅₆₇₈₉'[x]).join('')); if (f !== ans) fail(p, 'formula ' + f + ' vs ' + ans); }
      if (fam === 'hy') { const g = graphOf(p.m), i = p.hl[0]; const o = g[i].nb.map(e => e.o); const h = o.includes(3) || o.filter(x => x === 2).length >= 2 ? 'sp' : o.includes(2) ? 'sp²' : 'sp³'; if (h !== ans) fail(p, 'hybridisation'); }
      if (fam === 'ox') { const g = graphOf(p.m), i = p.hl[0]; let v = -C.hcount(g, i); for (const e of g[i].nb) if (g[e.n].el !== 'C' && g[e.n].el !== 'H') v += e.o; const s = v > 0 ? '+' + v : v < 0 ? '−' + (-v) : '0'; if (s !== ans) fail(p, 'oxidation number ' + s + ' vs ' + ans); }
      if (fam === 'rs') { const t = rd(p.m).tags.CIP_atoms.find(x => /[RS]/.test(x[1])); if (!t || t[1][1] !== ans) fail(p, 'R/S from drawing: ' + (t && t[1]) + ' vs ' + ans); }
      if (fam === 'ez') { const t = rd(p.m).tags.CIP_bonds; const lab = t.length ? t[0][2][1] : null; const want = ans === 'keine E/Z-Isomerie' ? null : ans; if (lab !== want) fail(p, 'E/Z from drawing: ' + lab + ' vs ' + ans); }
      if (fam === 'sc' && p.k === 'pick') { const m = C.RD.get_mol(rd(p.m).flat); const tags = JSON.parse(m.get_stereo_tags()); m.delete();
        const n = tags.CIP_atoms.length; if (n !== p.a.length) fail(p, 'stereocentre count ' + n + ' vs ' + p.a.length);
        const g = graphOf(p.m); for (const i of p.a) if (g[i].el !== 'C' || C.hcount(g, i) > 1) fail(p, 'picked atom cannot be a stereocentre'); }
      if (fam === 'sc' && p.k === 'mc') { const m = C.RD.get_mol(rd(p.m).flat); const n = JSON.parse(m.get_stereo_tags()).CIP_atoms.length; m.delete(); if (String(n) !== ans) fail(p, 'stereocentre count ' + n + ' vs ' + ans); }
      if (p.k === 'arr') { const ids2 = new Set(p.tm.atoms.map(a => a.id)), bids = new Set(p.tm.bonds.map(b => b.id)); for (const s of p.sol.flat()) { const [a, b] = s.split('>'); for (const x of [a, b]) { const [k, id] = x.split(':'); if (k === 'b' ? !bids.has(id) : !ids2.has(id)) fail(p, 'arrow refers to missing ' + x); } } }
      if (fam === 'fp') { const g = graphOf(p.m), i = p.a[0]; const want = /Hydroxy/.test(p.q) ? 'O' : 'N'; if (g[i].el !== want || g[i].nb.length !== 1 || g[i].nb[0].o !== 1) fail(p, 'picked atom is not the asked group'); }
      if (p.k === 'pick' && p.m !== undefined) for (const i of p.a) if (i >= M[p.m][0].length) fail(p, 'pick index out of range');
      ok(fam);
    } catch (e) { fail(p, 'exception ' + e.message); }
  }
  console.log('checked', P.length, 'puzzles;', Object.entries(stats).map(([k, v]) => k + ':' + v).join(' '));
  console.log(bad.length ? bad.length + ' problems:\n' + bad.slice(0, 40).join('\n') : 'no problems');
  process.exitCode = bad.length ? 1 : 0;
})();
