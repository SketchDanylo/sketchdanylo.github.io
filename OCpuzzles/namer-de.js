/* Chem Puzzles name reader: German (and English) IUPAC-style names -> molecular graph.
 * Used by the page to understand typed answers and by tools/ to check generated names.
 *   NameDE.read(text)        -> { ok, g, ast, text, fixed } | { ok:false, error }
 *   NameDE.same(g1, g2)      -> true when both graphs are the same constitution
 * Graph: [{ el, ch, nb:[{n,o}] }], hydrogens implicit.
 */
(function (root, factory) { const api = factory(); if (typeof module === 'object' && module.exports) module.exports = api; else root.NameDE = api; })(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';
  const VAL = { C: 4, N: 3, O: 2, F: 1, Cl: 1, Br: 1, I: 1, S: 2 };
  const add = (g, el, ch = 0) => { g.push({ el, ch, nb: [] }); return g.length - 1; };
  const bond = (g, a, b, o = 1) => { g[a].nb.push({ n: b, o }); g[b].nb.push({ n: a, o }); };
  const setOrder = (g, a, b, o) => { const e1 = g[a].nb.find(e => e.n === b), e2 = g[b].nb.find(e => e.n === a); if (!e1) throw new E('bond'); e1.o = o; e2.o = o; };
  function hcount(g, i) { const a = g[i], s = a.nb.reduce((t, e) => t + e.o, 0); let v = VAL[a.el]; if (a.el === 'N' && a.ch === 1) v = 4; if (a.el === 'O' && a.ch === -1) v = 1; if (a.el === 'O' && a.ch === 1) v = 3; return Math.max(0, v - s); }
  class E extends Error {}
  const fail = m => { throw new E(m); };

  /* tiny SMILES reader for fragments (attachment = atom 0) */
  function frag(sm) {
    const g = [], st = [], rings = {}; let at = -1, o = 1;
    for (const t of sm.match(/Cl|Br|[A-Z]|[()=#]|\d/g)) {
      if (t === '(') { st.push(at); continue; } if (t === ')') { at = st.pop(); continue; }
      if (t === '=') { o = 2; continue; } if (t === '#') { o = 3; continue; }
      if (/\d/.test(t)) { if (rings[t] !== undefined) { bond(g, at, rings[t][0], Math.max(o, rings[t][1])); delete rings[t]; } else rings[t] = [at, o]; o = 1; continue; }
      const i = add(g, t); if (at >= 0) bond(g, at, i, o); at = i; o = 1;
    }
    return g;
  }
  function graft(g, at, sub, subAt = 0, order = 1) { const off = g.length; sub.forEach(a => g.push({ el: a.el, ch: a.ch || 0, nb: a.nb.map(e => ({ n: e.n + off, o: e.o })) })); bond(g, at, off + subAt, order); return off; }

  const STEMS = [['undec', 11], ['dodec', 12], ['meth', 1], ['eth', 2], ['prop', 3], ['but', 4], ['pent', 5], ['hex', 6], ['hept', 7], ['oct', 8], ['non', 9], ['dec', 10]];
  const MULT = { di: 2, tri: 3, tetra: 4, penta: 5, hexa: 6, bis: 2, tris: 3, tetrakis: 4 };
  // simple substituents: fragment SMILES (atom 0 attaches); '=' prefix = double bond to the parent
  const SIMPLE = {
    fluor: 'F', chlor: 'Cl', brom: 'Br', iod: 'I', hydroxy: 'O', amino: 'N', oxo: '=O', nitro: 'N(=O)=O', cyano: 'C#N',
    methoxy: 'OC', ethoxy: 'OCC', propoxy: 'OCCC', butoxy: 'OCCCC', phenoxy: 'OC1=CC=CC=C1', formyl: 'C=O', carboxy: 'C(=O)O', acetyl: 'C(C)=O',
    vinyl: 'C=C', ethenyl: 'C=C', ethinyl: 'C#C', methyliden: '=C', methylen: '=C', mercapto: 'S', sulfanyl: 'S',
  };
  const SIMPLE_KEYS = Object.keys(SIMPLE).sort((a, b) => b.length - a.length);

  /* ---------- text clean-up: English, school and trivial forms ---------- */
  function prep(raw) {
    let s = raw.toLowerCase().normalize('NFC').trim();
    s = s.replace(/[–—‐‑−]/g, '-').replace(/[\[{]/g, '(').replace(/[\]}]/g, ')').replace(/saeure/g, 'säure').replace(/\s*-\s*/g, '-');
    s = s.replace(/^(?:\((?:\d*[rsez],?\s*)+\)-|(?:cis|trans|rac|\(±\)|\(\+\)|\(-\))-)/, '');
    // English
    const EN = [[/benzoic acid$/, 'benzoesäure'], [/oic acid$/, 'säure'], [/ic acid$/, 'säure'], [/fluoro/g, 'fluor'], [/chloro/g, 'chlor'], [/bromo/g, 'brom'], [/iodo/g, 'iod'], [/nitrile$/, 'nitril'], [/amine$/, 'amin'], [/amide$/, 'amid'], [/oate$/, 'oat'], [/aldehyde$/, 'aldehyd'], [/benzene/g, 'benzol'], [/aniline/g, 'anilin'], [/toluene/g, 'toluol'], [/ylidene/g, 'yliden'], [/ethynyl/g, 'ethinyl'],
      [/ane\b/g, 'an'], [/ene\b/g, 'en'], [/yne\b/g, 'in'], [/-yn-/g, '-in-'], [/one\b/g, 'on'], [/anenitril/g, 'annitril'], [/enenitril/g, 'ennitril']];
    for (const [a, b] of EN) s = s.replace(a, b);
    // trivial / older German names
    const TRIV = [[/^isopropanol$/, 'propan-2-ol'], [/^aceton$/, 'propanon'], [/^essigsäure$/, 'ethansäure'], [/^ameisensäure$/, 'methansäure'], [/^propionsäure$/, 'propansäure'], [/^buttersäure$/, 'butansäure'], [/^formaldehyd$/, 'methanal'], [/^acetaldehyd$/, 'ethanal'], [/^ethylen$/, 'ethen'], [/^acetylen$/, 'ethin'], [/^propylen$/, 'propen'], [/^glycerin$/, 'propan-1,2,3-triol'], [/^(?:ethylen)?glykol$/, 'ethan-1,2-diol'], [/^(?:tert|t)-?butanol$/, '2-methylpropan-2-ol'], [/^isobutan$/, '2-methylpropan'], [/^neopentan$/, '2,2-dimethylpropan'], [/^isopentan$/, '2-methylbutan']];
    for (const [a, b] of TRIV) s = s.replace(a, b);
    s = s.replace(/toluol/g, 'methylbenzol').replace(/acetat$/, 'ethanoat').replace(/formiat$/, 'methanoat').replace(/propionat$/, 'propanoat').replace(/butyrat$/, 'butanoat');
    s = s.replace(/(?:\(|\b)?(?:iso-?propyl)\)?/g, '(propan-2-yl)').replace(/\(?sec-?butyl\)?/g, '(butan-2-yl)').replace(/\(?(?:tert|t)-?butyl\)?/g, '(2-methylpropan-2-yl)').replace(/\(?iso-?butyl\)?/g, '(2-methylpropyl)');
    s = s.replace(/\s+/g, '');
    s = s.replace(/^(.+?)säure(methyl|ethyl|propyl|butyl)ester$/, '$2$1oat');
    s = s.replace(/(\d)(?=[a-zäöü(])/g, '$1-').replace(/([a-zäöü)])(?=\d)/g, '$1-').replace(/--+/g, '-').replace(/^-|-$/g, '');
    return s;
  }

  /* ---------- parser ---------- */
  function locList(s, i) { const m = s.slice(i).match(/^((?:\d+|n)(?:,(?:\d+|n))*)-/); return m ? { locs: m[1].split(',').map(x => x === 'n' ? 'N' : +x), len: m[0].length } : null; }
  function stemAt(s, i) { for (const [t, n] of STEMS.slice().sort((a, b) => b[0].length - a[0].length)) if (s.startsWith(t, i)) return { n, len: t.length }; return null; }
  // a substituent starting at s[i]; returns { sub:{g,A,at,order}, len, key }
  function subAt(s, i) {
    if (s[i] === '(') { // bracketed compound substituent
      let d = 0, j = i; for (; j < s.length; j++) { if (s[j] === '(') d++; else if (s[j] === ')') { d--; if (!d) break; } }
      if (d) fail('Klammer nicht geschlossen');
      return { sub: compound(s.slice(i + 1, j)), len: j - i + 1 };
    }
    for (const k of SIMPLE_KEYS) if (s.startsWith(k, i)) { const sm = SIMPLE[k]; const g = frag(sm.replace(/^=/, '')); return { sub: { g, A: [0], at: 0, order: sm.startsWith('=') ? 2 : 1, name: k }, len: k.length }; }
    if (s.startsWith('phenyl', i)) return { sub: ring6('phenyl'), len: 6 };
    const cy = s.startsWith('cyclo', i) ? 5 : 0;
    const st = stemAt(s, i + cy);
    if (st) {
      const r = s.slice(i + cy + st.len);
      let m = r.match(/^yl/);
      if (m) return { sub: alkyl(st.n, 1, !!cy), len: cy + st.len + 2 };
      m = r.match(/^an-(\d+)-yl/);
      if (m && !cy) return { sub: alkyl(st.n, +m[1], false), len: cy + st.len + m[0].length };
    }
    return null;
  }
  function alkyl(n, attach, cyc) {
    const g = []; const A = []; for (let k = 0; k < n; k++) { A.push(add(g, 'C')); if (k) bond(g, A[k - 1], A[k]); }
    if (cyc) { if (n < 3) fail('Ring zu klein'); bond(g, A[n - 1], A[0]); }
    if (attach < 1 || attach > n) fail('Lokant passt nicht');
    return { g, A, at: A[attach - 1], order: 1, name: 'alkyl' + n + (cyc ? 'c' : '') + '@' + attach };
  }
  function ring6(name) { const g = frag('C1=CC=CC=C1'); return { g, A: [0, 1, 2, 3, 4, 5], at: 0, order: 1, name }; }
  function compound(inner) { // e.g. "1-methylethyl", "2-bromethyl", "propan-2-yl", "4-nitrophenyl"
    const { pre, i } = prefixes(inner, 0);
    let sub;
    if (i === inner.length && pre.length && pre[pre.length - 1].locs.length === 1 && pre[pre.length - 1].locs[0] === null) sub = pre.pop().sub; // last item is the substituent's own parent
    else { const base = subAt(inner, i); if (!base || i + base.len !== inner.length) fail('Substituent nicht lesbar: ' + inner); sub = base.sub; }
    applyPrefixes(sub.g, sub.A, pre, false);
    return sub;
  }
  function prefixes(s, i) {
    const pre = [];
    for (;;) {
      let j = i, locs = null, mult = 1;
      const L = locList(s, j); if (L) { locs = L.locs; j += L.len; }
      const mm = s.slice(j).match(/^(tetrakis|tris|bis|tetra|penta|hexa|tri|di)(?=[a-z(])/);
      let sub = null;
      if (mm) { const t = subAt(s, j + mm[1].length); if (t) { mult = MULT[mm[1]]; sub = t; j += mm[1].length; } }
      if (!sub) sub = subAt(s, j);
      if (!sub) return { pre, i, pending: L && !pre.some(p => p.pendingOnly) ? null : null };
      // a plain locant group in front of the parent ("2-butanol") is not a prefix: subAt failed in that case
      j += sub.len;
      if (locs && locs.length !== mult) fail('Anzahl der Lokanten passt nicht zu ' + (mult > 1 ? 'di/tri' : 'einem Substituenten'));
      if (!locs && mult > 1) locs = new Array(mult).fill(null);
      pre.push({ locs: locs || [null], sub: sub.sub, text: s.slice(i, j) });
      if (s[j] === '-') j++;
      i = j;
    }
  }
  function applyPrefixes(g, A, pre, isParent, info = {}) {
    for (const p of pre) for (const l of p.locs) {
      let at;
      if (l === 'N') { if (info.nAtom === undefined) fail('N-Lokant ohne Stickstoff'); at = info.nAtom; }
      else if (l === null) { if (A.length <= 2 || info.ring || info.single) at = A[info.default || 0]; else fail('Lokant fehlt'); }
      else { if (l < 1 || l > A.length) fail('Lokant ' + l + ' gibt es nicht'); at = A[l - 1]; }
      graft(g, at, p.sub.g, p.sub.at, p.sub.order);
    }
  }
  const RETAINED = { phenol: ['O', 1], anilin: ['N', 1], benzoesäure: ['C(=O)O', 1], benzaldehyd: ['C=O', 1], benzonitril: ['C#N', 1], benzol: null };

  function parse(s) {
    const { pre, i } = prefixes(s, 0);
    let rest = s.slice(i), pending = null;
    const L = locList(rest, 0); if (L) { pending = L.locs; rest = rest.slice(L.len); }
    const g = [], ast = { pre: [], unsat: [], suffix: null, parent: null };
    let A, n, ring = false, nAtom, ester = null;
    // retained benzene parents
    const ret = Object.keys(RETAINED).find(k => rest === k);
    if (ret) {
      const r = ring6('benzol'); r.g.forEach(a => g.push(a)); A = r.A; ring = true; n = 6; ast.parent = 'benzol';
      if (RETAINED[ret]) { const f = frag(RETAINED[ret][0]); graft(g, A[0], f); if (ret === 'anilin') nAtom = g.length - 1; ast.suffix = { kind: ret, locs: [1] }; }
      if (pending) fail('Lokant vor dem Stammnamen passt hier nicht');
      rest = '';
    } else {
      const cyc = rest.startsWith('cyclo'); if (cyc) rest = rest.slice(5);
      const st = stemAt(rest, 0); if (!st) fail('Stammname nicht erkannt');
      n = st.n; rest = rest.slice(st.len);
      for (let k = 0; k < n; k++) { const a = add(g, 'C'); if (k) bond(g, a - 1, a); } A = g.map((_, k) => k);
      if (cyc) { if (n < 3) fail('Ring zu klein'); bond(g, A[n - 1], A[0]); ring = true; }
      ast.parent = (cyc ? 'cyclo' : 'chain') + n;
      // saturation / unsaturation
      let m;
      const unsat = [];
      if ((m = rest.match(/^a?-((?:\d+,)*\d+)-(di|tri)?(en|in)/))) { unsat.push([m[1].split(',').map(Number), m[3], m[2]]); rest = rest.slice(m[0].length); if ((m = rest.match(/^-((?:\d+,)*\d+)-(di|tri)?(in)/))) { unsat.push([m[1].split(',').map(Number), m[3], m[2]]); rest = rest.slice(m[0].length); } }
      else if ((m = rest.match(/^a?(di|tri)?(en|in)(?![a-zäöü]*yl)/)) && !rest.startsWith('an')) { unsat.push([null, m[2], m[1]]); rest = rest.slice(m[0].length); }
      else if (rest.startsWith('an')) rest = rest.slice(2);
      else fail('Endung -an/-en/-in fehlt');
      for (const u of unsat) {
        let locs = u[0]; const k = u[2] ? MULT[u[2]] : 1;
        if (!locs) { if (pending && pending.length === k) { locs = pending; pending = null; } else if (k === 1 && (n <= 3 || ring)) locs = [1]; else fail('Lokant der Mehrfachbindung fehlt'); }
        if (locs.length !== k) fail('Anzahl der Lokanten passt nicht');
        for (const l of locs) { if (l < 1 || l > n - (ring ? 0 : 1)) fail('Lokant ' + l + ' gibt es nicht'); setOrder(g, A[l - 1], A[l % n], u[1] === 'en' ? 2 : 3); }
        ast.unsat.push([locs.slice().sort((a, b) => a - b), u[1]]);
      }
      // suffix
      m = rest.match(/^(?:-((?:\d+,)*\d+)-)?(di|tri)?(ol|on|al|säure|nitril|amin|amid|oat)$/);
      if (rest && !m) fail('Endung nicht erkannt: ' + rest);
      if (m) {
        const kind = m[3], k = m[2] ? MULT[m[2]] : 1; let locs = m[1] ? m[1].split(',').map(Number) : null;
        if (!locs) {
          if (['al', 'säure', 'nitril', 'amid', 'oat'].includes(kind)) locs = k === 2 ? [1, n] : [1];
          else if (pending && pending.length === k) { locs = pending; pending = null; }
          else if (k === 1 && (n <= 2 || ring)) locs = [1];
          else if (kind === 'on' && k === 1 && (n === 3 || n === 4)) locs = [2];
          else fail('Lokant der Endung fehlt');
        }
        if (locs.length !== k) fail('Anzahl der Lokanten passt nicht');
        for (const l of locs) {
          if (l < 1 || l > n) fail('Lokant ' + l + ' gibt es nicht');
          const c = A[l - 1];
          if (kind === 'ol') bond(g, c, add(g, 'O'));
          else if (kind === 'on' || kind === 'al') bond(g, c, add(g, 'O'), 2);
          else if (kind === 'säure') { bond(g, c, add(g, 'O'), 2); bond(g, c, add(g, 'O')); }
          else if (kind === 'nitril') bond(g, c, add(g, 'N'), 3);
          else if (kind === 'amin') { nAtom = add(g, 'N'); bond(g, c, nAtom); }
          else if (kind === 'amid') { bond(g, c, add(g, 'O'), 2); nAtom = add(g, 'N'); bond(g, c, nAtom); }
          else if (kind === 'oat') { bond(g, c, add(g, 'O'), 2); ester = add(g, 'O'); bond(g, c, ester); }
        }
        if (kind === 'on') ast.ketoneC = locs.map(l => A[l - 1]);
        ast.suffix = { kind, locs: locs.slice().sort((a, b) => a - b) };
      }
      if (pending) fail('Lokant vor dem Stammnamen passt hier nicht');
    }
    // ester alkyl = first locant-free alkyl prefix
    let preList = pre;
    if (ester !== null) {
      const idx = pre.findIndex(p => p.locs.length === 1 && p.locs[0] === null && /^alkyl\d+@1$/.test(p.sub.name || ''));
      if (idx < 0) fail('Ester: Alkylgruppe fehlt');
      graft(g, ester, pre[idx].sub.g, pre[idx].sub.at); ast.ester = pre[idx].sub.name;
      preList = pre.filter((_, k) => k !== idx);
    }
    applyPrefixes(g, A, preList, true, { nAtom, ring, single: n <= 2 });
    for (const c of ast.ketoneC || []) if (g[c].nb.filter(e => g[e.n].el === 'C').length < 2) fail('Keton am Kettenende wäre ein Aldehyd');
    delete ast.ketoneC;
    for (let k = 0; k < g.length; k++) { const s2 = g[k].nb.reduce((t, e) => t + e.o, 0); if (s2 > VAL[g[k].el] + (g[k].el === 'N' ? 2 : 0)) fail('Zu viele Bindungen an einem Atom'); }
    ast.pre = preList.map(p => [subKey(p.sub), p.locs.map(l => l === null ? 0 : l).sort((a, b) => a - b).join(',')]).sort();
    ast.order = preList.map(p => p.text.replace(/^[\d,n-]+/, '').replace(/^(?:di|tri|tetra|bis|tris)/, '').replace(/^\(/, ''));
    return { g, ast };
  }
  function subKey(sub) { // canonical string of a substituent, rooted at its attachment atom
    const g = sub.g, seen = new Set();
    const walk = (i, from) => { seen.add(i); const kids = g[i].nb.filter(e => e.n !== from && !seen.has(e.n)).map(e => e.o + walk(e.n, i)).sort(); return g[i].el + (kids.length ? '(' + kids.join(',') + ')' : ''); };
    return sub.order + walk(sub.at, -1);
  }

  /* ---------- typo tolerance ---------- */
  const VOCAB = ['meth', 'eth', 'prop', 'but', 'pent', 'hex', 'hept', 'oct', 'non', 'dec', 'yl', 'an', 'en', 'in', 'a', 'ol', 'on', 'al', 'säure', 'nitril', 'amin', 'amid', 'oat', 'di', 'tri', 'tetra', 'bis', 'tris', 'cyclo', 'benzol', 'phenol', 'anilin', 'benzoesäure', 'benzaldehyd', 'benzonitril', 'phenyl', ...Object.keys(SIMPLE)];
  function lev(a, b) { const m = a.length, n = b.length; if (!m) return n; if (!n) return m; let prev = Array.from({ length: n + 1 }, (_, j) => j); for (let i = 1; i <= m; i++) { const cur = [i]; for (let j = 1; j <= n; j++) cur[j] = Math.min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + (a[i - 1] === b[j - 1] ? 0 : 1)); prev = cur; } return prev[n]; }
  function fixRun(w) { // best segmentation of a letter run into vocabulary words, with edit cost
    const N = w.length, best = new Array(N + 1).fill(null); best[0] = { cost: 0, out: '' };
    for (let i = 0; i < N; i++) { if (!best[i]) continue;
      for (const t of VOCAB) for (let len = Math.max(1, t.length - 2); len <= t.length + 2 && i + len <= N; len++) {
        const piece = w.slice(i, i + len); const c = lev(piece, t); if (c > (t.length <= 3 ? 0 : t.length <= 5 ? 1 : 2)) continue;
        const cost = best[i].cost + c; if (!best[i + len] || cost < best[i + len].cost) best[i + len] = { cost, out: best[i].out + t };
      } }
    return best[N];
  }
  function fuzzy(s) {
    let total = 0, ok = true;
    const out = s.replace(/[a-zäöü]+/g, w => { const r = fixRun(w); if (!r) { ok = false; return w; } total += r.cost; return r.out; });
    return ok && total > 0 && total <= 3 ? out : null;
  }

  /* ---------- public ---------- */
  function read(text) {
    const s = prep(text || '');
    if (!s) return { ok: false, error: 'leer' };
    try { const r = parse(s); return { ok: true, g: r.g, ast: r.ast, text: s, fixed: false }; }
    catch (e) {
      if (!(e instanceof E)) throw e;
      const f = fuzzy(s);
      if (f) { try { const r = parse(f); return { ok: true, g: r.g, ast: r.ast, text: f, fixed: true }; } catch (e2) { if (!(e2 instanceof E)) throw e2; } }
      return { ok: false, error: e.message, text: s };
    }
  }
  function normNitro(g) { // charged nitro -> N(=O)=O
    for (let i = 0; i < g.length; i++) if (g[i].el === 'N' && g[i].ch === 1) { const om = g[i].nb.find(e => g[e.n].el === 'O' && g[e.n].ch === -1 && e.o === 1); if (om) { g[i].ch = 0; g[om.n].ch = 0; setOrder(g, i, om.n, 2); } }
    return g;
  }
  function same(g1, g2) {
    g1 = normNitro(g1.map(a => ({ el: a.el, ch: a.ch || 0, nb: a.nb.map(e => ({ ...e })) }))).filter(a => a.el !== 'H');
    g2 = normNitro(g2.map(a => ({ el: a.el, ch: a.ch || 0, nb: a.nb.map(e => ({ ...e })) }))).filter(a => a.el !== 'H');
    if (g1.length !== g2.length) return false;
    const h1 = g1.map((_, i) => hcount(g1, i)), h2 = g2.map((_, i) => hcount(g2, i));
    aromatize(g1); aromatize(g2);
    const col = g => { const hh = g === g1 ? h1 : h2; let c = g.map((a, i) => a.el + a.ch + ':' + a.nb.length + ':' + hh[i] + ':' + a.nb.map(e => e.o).sort().join('')); for (let r = 0; r < 4; r++) c = g.map((a, i) => c[i] + '|' + a.nb.map(e => e.o + c[e.n]).sort().join(',')); return c; };
    const c1 = col(g1), c2 = col(g2);
    if (c1.slice().sort().join() !== c2.slice().sort().join()) return false;
    const map = new Array(g1.length).fill(-1), used = new Array(g2.length).fill(false);
    const order = g1.map((_, i) => i).sort((a, b) => g1[b].nb.length - g1[a].nb.length);
    const fits = (i, j) => { if (c1[i] !== c2[j]) return false; for (const e of g1[i].nb) { const k = map[e.n]; if (k >= 0) { const f = g2[j].nb.find(x => x.n === k); if (!f || f.o !== e.o) return false; } } return true; };
    const go = d => { if (d === order.length) return true; const i = order[d]; for (let j = 0; j < g2.length; j++) if (!used[j] && fits(i, j)) { map[i] = j; used[j] = true; if (go(d + 1)) return true; map[i] = -1; used[j] = false; } return false; };
    return go(0);
  }
  function aromatize(g) { // alternating 6-rings get bond order 4, so both Kekulé forms compare equal
    const n = g.length;
    for (let s0 = 0; s0 < n; s0++) {
      const walk = (path) => { const i = path[path.length - 1]; if (path.length === 6) { const f = g[i].nb.find(e => e.n === path[0]); if (!f) return false; const os = path.map((a, k) => g[a].nb.find(e => e.n === path[(k + 1) % 6]).o); if (os.every((o, k) => o === 4 || o === (os[0] === 4 ? 4 : (k % 2 ? 3 - os[0] : os[0]))) && os.includes(2) || os.every(o => o === 4)) { path.forEach((a, k) => { const b = path[(k + 1) % 6]; g[a].nb.find(e => e.n === b).o = 4; g[b].nb.find(e => e.n === a).o = 4; }); return true; } return false; }
        for (const e of g[i].nb) if (!path.includes(e.n) && e.n > path[0] && (e.o === 1 || e.o === 2 || e.o === 4)) { path.push(e.n); if (walk(path)) return true; path.pop(); } return false; };
      if (g[s0].el === 'C') walk([s0]);
    }
    return g;
  }
  const fromDepiction = m => { const [A, B] = m; const idx = []; const g = []; A.forEach((a, i) => { if (a[0] !== 'H') { idx[i] = g.length; g.push({ el: a[0], ch: a[3] || 0, nb: [] }); } }); for (const b of B) if (idx[b[0]] !== undefined && idx[b[1]] !== undefined) bond(g, idx[b[0]], idx[b[1]], b[2]); return g; };
  const formula = g => { const c = {}; g.forEach((a, i) => { c[a.el] = (c[a.el] || 0) + 1; const h = hcount(g, i); if (h) c.H = (c.H || 0) + h; }); return c; };
  return { read, same, prep, parse: s => parse(prep(s)), fromDepiction, formula, hcount, normNitro };
});
