/* Offline, connectivity-based organic nomenclature.
 * References and explicit boundaries: docs/naming-engine.md.
 * Pure graph input; no DOM, network, formula matching, or coordinate-based stereo.
 */
(function (root, factory) {
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.NativeNamer = api;
})(typeof globalThis !== 'undefined' ? globalThis : this, function () {
  'use strict';
  class Unsupported extends Error {}
  const fail = message => { throw new Unsupported(message); };
  const SMALL = ['', 'meth', 'eth', 'prop', 'but', 'pent', 'hex', 'hept', 'oct', 'non',
    'dec', 'undec', 'dodec', 'tridec', 'tetradec', 'pentadec', 'hexadec', 'heptadec', 'octadec', 'nonadec',
    'icos', 'henicos', 'docos', 'tricos', 'tetracos', 'pentacos', 'hexacos', 'heptacos', 'octacos', 'nonacos'];
  const UNITS = ['', 'hen', 'do', 'tri', 'tetra', 'penta', 'hexa', 'hepta', 'octa', 'nona'];
  const TENS = ['', '', 'cos', 'triacont', 'tetracont', 'pentacont', 'hexacont', 'heptacont', 'octacont', 'nonacont'];
  function stem(n) {
    if (SMALL[n]) return SMALL[n];
    if (n >= 30 && n < 100) return UNITS[n % 10] + TENS[Math.floor(n / 10)];
    if (n === 100) return 'hect';
    fail('Carbon parents longer than 100 atoms are not yet supported.');
  }
  function multi(n, complex = false) {
    if (n === 1) return '';
    if (complex) return n === 2 ? 'bis' : n === 3 ? 'tris' : multi(n) + 'kis';
    return ['', '', 'di', 'tri', 'tetra'][n] || stem(n) + 'a';
  }
  const HALO = { F: 'fluoro', Cl: 'chloro', Br: 'bromo', I: 'iodo' };
  const SENIORITY = { acid: 100, carboxylate: 99, ester: 98, acylhalide: 97, sulfo: 95,
    sulfonate: 94, amide: 90, nitrile: 80, aldehyde: 70, ketone: 60, olate: 51,
    ol: 50, thiolate: 41, thiol: 40, amine: 30 };
  const TERMINAL = new Set(['acid', 'carboxylate', 'ester', 'acylhalide', 'amide', 'nitrile', 'aldehyde']);
  const SUFFIX = { acid: 'oic acid', carboxylate: 'oate', ester: 'oate', amide: 'amide',
    nitrile: 'nitrile', aldehyde: 'al', ketone: 'one', ol: 'ol', olate: 'olate',
    thiol: 'thiol', thiolate: 'thiolate', amine: 'amine', sulfo: 'sulfonic acid', sulfonate: 'sulfonate' };
  const EXO = { acid: 'carboxylic acid', carboxylate: 'carboxylate', ester: 'carboxylate',
    amide: 'carboxamide', nitrile: 'carbonitrile', aldehyde: 'carbaldehyde' };
  const cmp = (a, b) => {
    for (let i = 0; i < Math.max(a.length, b.length); i++) {
      if ((a[i] ?? Infinity) !== (b[i] ?? Infinity)) return (a[i] ?? Infinity) - (b[i] ?? Infinity);
    }
    return 0;
  };
  const alphabetical = s => s.replace(/\((?:\d+[RSrsEZ],?)+\)-?/g,'').replace(/\d+[a-z]?[,-]?/g, '').replace(/[^a-z]/gi, '').toLowerCase();
  const compound = s => /[\d()[\]{}, -]/.test(s)
    || /(?:amino|hydroxy|oxo|sulfanyl|oxy|carboxy|fluoro|chloro|bromo|iodo).*yl$|ylamino$|ylsulfanyl$/.test(s);
  const wrap = s => s.includes('{') ? '(' + s + ')' : s.includes('[') ? '{' + s + '}' : s.includes('(') ? '[' + s + ']' : '(' + s + ')';
  const groupName = s => compound(s) ? wrap(s) : s;
  const join = (prefix, base) => prefix + (/^\d/.test(base) && prefix ? '-' : '') + base;
  function prefixes(items, omitLocants = false) {
    const groups = new Map();
    items.forEach(p => {
      if (!groups.has(p.name)) groups.set(p.name, []);
      groups.get(p.name).push(p.loc);
    });
    return [...groups].sort((a, b) => alphabetical(a[0]).localeCompare(alphabetical(b[0])))
      .map(([name, locs]) => {
        locs.sort((a, b) => typeof a === 'number' && typeof b === 'number' ? a - b : String(a).localeCompare(String(b), 'en', { numeric: true }));
        const complex = compound(name) || (locs.length > 1 && /^(di|tri|tetra|dec)/.test(name));
        return (omitLocants ? '' : locs.join(',') + '-') + multi(locs.length, complex)
          + (complex ? wrap(name) : name);
      }).join('-');
  }
  function hydro(n, doubles, triples, cyclic) {
    let out = (cyclic ? 'cyclo' : '') + stem(n);
    if (!doubles.length && !triples.length) return out + 'ane';
    for (const [locs, ending] of [[doubles, 'ene'], [triples, 'yne']]) {
      if (!locs.length) continue;
      out = out.replace(/e$/, '');
      if (locs.length > 1 && !/[aeiou]$/.test(out)) out += 'a';
      out += (n === 2 ? '' : '-' + locs.join(',') + '-') + multi(locs.length) + ending;
    }
    return out;
  }
  // Atom and bond patterns in fixed IUPAC ring numbering. Both directions and all
  // rotations are matched; indicated H is preserved by the explicit Kekule pattern.
  const RINGS = [
    ['pyridine', 'NCCCCC', '212121'], ['pyridazine', 'NNCCCC', '121212'],
    ['pyrimidine', 'NCNCCC', '212121'], ['pyrazine', 'NCCNCC', '212121'],
    ['1H-pyrrole', 'NCCCC', '12121'], ['furan', 'OCCCC', '12121'], ['thiophene', 'SCCCC', '12121'],
    ['1H-imidazole', 'NCNCC', '12121'], ['1H-pyrazole', 'NNCCC', '12121'],
    ['1,3-oxazole', 'OCNCC', '12121'], ['1,2-oxazole', 'ONCCC', '12121'],
    ['1,3-thiazole', 'SCNCC', '12121'], ['1,2-thiazole', 'SNCCC', '12121'],
    ['aziridine', 'NCC', '111'], ['azetidine', 'NCCC', '1111'],
    ['pyrrolidine', 'NCCCC', '11111'], ['piperidine', 'NCCCCC', '111111'],
    ['azepane', 'NCCCCCC', '1111111'], ['oxirane', 'OCC', '111'],
    ['oxetane', 'OCCC', '1111'], ['oxolane', 'OCCCC', '11111'], ['oxane', 'OCCCCC', '111111'],
    ['thiirane', 'SCC', '111'], ['thietane', 'SCCC', '1111'],
    ['thiolane', 'SCCCC', '11111'], ['thiane', 'SCCCCC', '111111'],
    ['morpholine', 'OCCNCC', '111111'], ['piperazine', 'NCCNCC', '111111'],
    ['1,3-dioxolane', 'OCOCC', '11111'], ['1,4-dioxane', 'OCCOCC', '111111']
  ];

  class Engine {
    constructor(input, options) {
      if (input.some(a => a.isotope || a.stereo || a.chirality || (a.nb || []).some(e => e.stereo)))
        fail('Explicit isotopes or stereochemical annotations require another naming method.');
      this.g = input.map(a => ({ el: a.el, h: a.h ?? 0, charge: a.charge || 0,
        nb: (a.nb || []).map(e => ({ n: e.n, o: e.o })) }));
      this.stereo=options.stereo||{atoms:[],bonds:[]};
      if(this.stereo.atoms.some(x=>!this.g[x.atom]||this.g[x.atom].el!=='C'||!/^[RSrs]$/.test(x.label)))
        fail('The assigned heteroatom stereocenter needs a naming method not yet supported.');
      if(this.stereo.bonds.some(x=>!/^[EZ]$/.test(x.label)||!this.g[x.a]?.nb.some(e=>e.n===x.b&&e.o===2)))
        fail('Invalid stereochemical bond assignment.');
      this.work = 0;
      this.limit = options.maxWork ?? 250000;
      this.depth = 0;
      this.memo = new Map();
      this.validate();
      this.rings = this.findRings();
      this.ringAtoms = new Set(this.rings.flatMap(r => r.atoms));
    }
    tick(n = 1) { if ((this.work += n) > this.limit) fail('The structure exceeds the local search budget; no partial name was returned.'); }
    validate() {
      const g = this.g;
      if (!g.length) fail('Draw a molecule first.');
      if (g.length > 1024) fail('The offline engine currently accepts up to 1024 heavy atoms.');
      for (let i = 0; i < g.length; i++) {
        const a = g[i];
        if (!['C', 'N', 'O', 'S', 'P', 'F', 'Cl', 'Br', 'I'].includes(a.el)) fail('Unsupported element environment: ' + a.el + '.');
        if (!Number.isInteger(a.h) || a.h < 0 || !Number.isInteger(a.charge)) fail('Invalid hydrogen count or formal charge.');
        const neighbors = new Set();
        for (const e of a.nb) {
          if (!g[e.n] || e.n === i || neighbors.has(e.n) || ![1, 2, 3].includes(e.o)) fail('Invalid molecular graph or unsupported bond type.');
          if (!g[e.n].nb.some(r => r.n === i && r.o === e.o)) fail('The molecular graph has an unmatched bond.');
          neighbors.add(e.n);
        }
        const val = a.h + a.nb.reduce((s, e) => s + e.o, 0);
        const allowed = a.el === 'C' ? [4] : a.el === 'N' ? [a.charge === 1 ? 4 : 3]
          : a.el === 'O' ? [a.charge === -1 ? 1 : 2] : a.el === 'S' ? [2, 4, 6] : a.el === 'P' ? [5] : [1];
        if (!allowed.includes(val)) fail('Invalid or unsupported valence at atom ' + (i + 1) + ' (' + a.el + ').');
        if (a.charge && !((a.el === 'O' && a.charge === -1) || (a.el === 'N' && a.charge === 1)))
          fail('This charged atom environment is not supported locally.');
        if (a.el === 'N' && a.charge === 1 && !(a.nb.length === 3
          && a.nb.some(e => g[e.n].el === 'O' && e.o === 2 && !g[e.n].charge)
          && a.nb.some(e => g[e.n].el === 'O' && e.o === 1 && g[e.n].charge === -1)
          && a.nb.some(e => (g[e.n].el === 'C' || (g.length === 4 && g.every(x => x.el === 'N' || x.el === 'O'))) && e.o === 1)))
          fail('Charged nitrogen environments other than nitro groups require another naming method.');
      }
      if (this.component(0, new Set(g.map((_, i) => i))).size !== g.length) fail('Separate fragments and salts require a separate naming method.');
    }
    bond(a, b) { return this.g[a].nb.find(e => e.n === b)?.o || 0; }
    component(start, allowed, blocked = new Set()) {
      const seen = new Set(), stack = [start];
      while (stack.length) {
        this.tick();
        const i = stack.pop();
        if (seen.has(i) || blocked.has(i) || !allowed.has(i)) continue;
        seen.add(i);
        this.g[i].nb.forEach(e => { if (!seen.has(e.n)) stack.push(e.n); });
      }
      return seen;
    }
    findRings() {
      // Bridges partition a cactus graph into intact rings and acyclic atoms.
      const g = this.g, tin = [], low = [], bridges = new Set();
      let time = 0;
      const edge = (a, b) => Math.min(a, b) + ':' + Math.max(a, b);
      const visit = (i, parent) => {
        this.tick(); tin[i] = low[i] = ++time;
        for (const e of g[i].nb) {
          if (e.n === parent) continue;
          if (tin[e.n]) low[i] = Math.min(low[i], tin[e.n]);
          else {
            visit(e.n, i); low[i] = Math.min(low[i], low[e.n]);
            if (low[e.n] > tin[i]) bridges.add(edge(i, e.n));
          }
        }
      };
      visit(0, -1);
      const ringNB = g.map((a, i) => a.nb.filter(e => !bridges.has(edge(i, e.n))).map(e => e.n));
      const seen = new Set(), rings = [];
      for (let i = 0; i < g.length; i++) {
        if (!ringNB[i].length || seen.has(i)) continue;
        const members = new Set(), stack = [i];
        while (stack.length) {
          const a = stack.pop(); if (members.has(a)) continue;
          members.add(a); ringNB[a].forEach(n => stack.push(n));
        }
        if ([...members].some(a => ringNB[a].length !== 2)) {
          const candidates = this.polycycleCandidates(members, ringNB);
          members.forEach(a => seen.add(a));
          rings.push({ atoms: [...members], candidates }); continue;
        }
        const r = []; let p = -1, c = i;
        do { r.push(c); seen.add(c); const next = ringNB[c].find(n => n !== p); p = c; c = next; } while (c !== i);
        rings.push({ atoms: r });
      }
      return rings;
    }
    polycycleCandidates(members, ringNB) {
      const g = this.g, candidates = [];
      if ([...members].some(i => g[i].el !== 'C' || ringNB[i].some(j => this.bond(i, j) !== 1)))
        fail('Unsaturated or heteroatom-containing fused/bridged/spiro systems require another naming method.');
      const junctions = [...members].filter(i => ringNB[i].length !== 2);
      const walk = (start, next) => {
        const path = [start]; let prev = start, curr = next;
        while (true) {
          this.tick(); path.push(curr);
          if (junctions.includes(curr)) return path;
          const after = ringNB[curr].find(i => i !== prev); prev = curr; curr = after;
          if (path.length > members.size + 1) fail('Unsupported ring numbering.');
        }
      };
      if (junctions.length === 2 && junctions.every(i => ringNB[i].length === 3)) {
        for (const head of junctions) {
          const paths = ringNB[head].map(n => walk(head, n));
          if (paths.some(p => p.at(-1) === head)) fail('Unsupported ring assembly.');
          const lengths = paths.map(p => p.length - 2).sort((a, b) => b - a);
          for (let a = 0; a < 3; a++) for (let b = 0; b < 3; b++) {
            if (a === b) continue;
            const d = 3 - a - b, ordered = [paths[a], paths[b], paths[d]];
            if (cmp(ordered.map(p => p.length - 2), lengths)) continue;
            const arr = [...ordered[0], ...ordered[1].slice(1, -1).reverse(), ...ordered[2].slice(1, -1)];
            if (new Set(arr).size !== members.size) fail('Unsupported polycyclic overlap.');
            candidates.push({ arr, ring: true, poly: true, base: 'bicyclo[' + lengths.join('.') + ']' + stem(arr.length) + 'ane' });
          }
        }
      } else if (junctions.length === 1 && ringNB[junctions[0]].length === 4) {
        const center = junctions[0], paths = ringNB[center].map(n => walk(center, n));
        const lengths = [...new Set(paths.map(p => p.length - 2))].sort((a, b) => a - b);
        const min = lengths[0], max = lengths.at(-1);
        for (const small of paths) for (const large of paths) {
          if (small.length - 2 !== min || large.length - 2 !== max || small.slice(1, -1).some(i => large.slice(1, -1).includes(i))) continue;
          const arr = [...small.slice(1, -1), center, ...large.slice(1, -1)];
          candidates.push({ arr, ring: true, poly: true, base: 'spiro[' + [min, max].join('.') + ']' + stem(arr.length) + 'ane' });
        }
      }
      if (!candidates.length) fail('This polycyclic ring system is outside the local numbering rules.');
      return candidates;
    }
    groups(set) {
      const g = this.g, groups = [], used = new Set();
      const add = (type, c, atoms, extra = {}) => { groups.push({ type, anchors: [c], atoms, ...extra }); atoms.forEach(x => used.add(x)); };
      for (const c of set) {
        if (g[c].el !== 'C') continue;
        const nb = g[c].nb.filter(e => set.has(e.n));
        const oxo = nb.find(e => g[e.n].el === 'O' && e.o === 2 && g[e.n].nb.length === 1);
        if (oxo) {
          const acidO = nb.find(e => g[e.n].el === 'O' && e.o === 1 && g[e.n].nb.length === 1);
          const esterO = nb.find(e => g[e.n].el === 'O' && e.o === 1 && g[e.n].nb.length === 2 && !this.ringAtoms.has(e.n));
          const amideN = nb.find(e => g[e.n].el === 'N' && e.o === 1 && !this.ringAtoms.has(e.n) && !g[e.n].charge && g[e.n].nb.every(x => x.o === 1));
          const halo = nb.find(e => HALO[g[e.n].el]);
          if (acidO) add(g[acidO.n].charge === -1 ? 'carboxylate' : 'acid', c, [oxo.n, acidO.n]);
          else if (esterO) add('ester', c, [oxo.n, esterO.n], { link: esterO.n });
          else if (amideN) add('amide', c, [oxo.n, amideN.n], { link: amideN.n });
          else if (halo) add('acylhalide', c, [oxo.n, halo.n], { halogen: g[halo.n].el });
          else add(g[c].h ? 'aldehyde' : 'ketone', c, [oxo.n]);
        } else {
          const nitrile = nb.find(e => g[e.n].el === 'N' && e.o === 3 && g[e.n].nb.length === 1);
          if (nitrile) add('nitrile', c, [nitrile.n]);
        }
      }
      for (const c of set) {
        if (g[c].el !== 'C') continue;
        for (const e of g[c].nb) {
          if (!set.has(e.n) || used.has(e.n) || e.o !== 1) continue;
          const a = g[e.n];
          if ((a.el === 'O' || a.el === 'S') && a.nb.length === 1)
            add(a.el === 'O' ? (a.charge ? 'olate' : 'ol') : (a.charge ? 'thiolate' : 'thiol'), c, [e.n]);
          if (a.el === 'S' && a.nb.length === 4) {
            const oxy = a.nb.filter(x => g[x.n].el === 'O' && x.o === 2 && g[x.n].nb.length === 1);
            const oh = a.nb.find(x => g[x.n].el === 'O' && x.o === 1 && g[x.n].nb.length === 1);
            if (oxy.length === 2 && oh) add(g[oh.n].charge ? 'sulfonate' : 'sulfo', c, [e.n, oh.n, ...oxy.map(x => x.n)]);
          }
        }
      }
      for (const n of set) {
        const a = g[n];
        if (a.el !== 'N' || used.has(n) || this.ringAtoms.has(n) || a.charge) continue;
        if (a.nb.every(e => g[e.n].el === 'C' && e.o === 1)) {
          const anchors = a.nb.map(e => e.n).filter(c => set.has(c));
          if (anchors.length) add('amine', anchors[0], [n], { anchors, link: n });
        }
      }
      return groups;
    }
    ringCandidates(ring) {
      const out = [], n = ring.length, g = this.g;
      for (let start = 0; start < n; start++) for (const dir of [1, -1]) {
        const arr = ring.map((_, k) => ring[(start + dir * k + n * 2) % n]);
        const els = arr.map(i => g[i].el).join('');
        const orders = arr.map((i, k) => this.bond(i, arr[(k + 1) % n])).join('');
        if (/^C+$/.test(els)) {
          const benzene = n === 6 && /^(121212|212121)$/.test(orders);
          out.push({ arr, ring: true, benzene, base: benzene ? 'benzene' : null });
        } else {
          const alternating = s => /^(121212|212121)$/.test(s);
          const match = RINGS.find(([, es, os]) => es === els && (os === orders || (alternating(os) && alternating(orders))));
          if (match) out.push({ arr, ring: true, base: match[0], hetero: true });
        }
      }
      return out;
    }
    candidates(set, root) {
      const out = [], g = this.g;
      for (const ring of this.rings) {
        if (ring.atoms.every(i => set.has(i)) && (root === undefined || ring.atoms.includes(root)))
          out.push(...(ring.candidates || this.ringCandidates(ring.atoms)));
      }
      const carbons = [...set].filter(i => g[i].el === 'C' && !this.ringAtoms.has(i));
      const allowed = new Set(carbons);
      const neighbors = i => g[i].nb.filter(e => allowed.has(e.n));
      // Enumerate endpoint-to-endpoint tree paths, including the attachment as an
      // internal point for radicals (e.g. hexan-2-yl, not 1-methylpentyl).
      const ends = carbons.filter(i => neighbors(i).length <= 1);
      for (const start of ends) {
        const stack = [[start, -1, [start]]];
        while (stack.length) {
          this.tick();
          const [i, prev, arr] = stack.pop();
          const next = neighbors(i).filter(e => e.n !== prev);
          if (!next.length && (root === undefined || arr.includes(root))) out.push({ arr });
          for (const e of next) stack.push([e.n, i, [...arr, e.n]]);
        }
      }
      return out;
    }
    locate(cand, group) {
      for (const a of group.anchors) if (cand.arr.includes(a)) return { group, anchor: a, at: a, exo: false };
      if (cand.ring && TERMINAL.has(group.type)) {
        const c = group.anchors[0];
        if (this.ringAtoms.has(c)) return null;
        const edge = this.g[c].nb.find(e => cand.arr.includes(e.n) && e.o === 1);
        if (edge && this.g[c].nb.filter(e => this.g[e.n].el === 'C').length === 1 && EXO[group.type])
          return { group, anchor: c, at: edge.n, exo: true };
      }
      return null;
    }
    decorate(cand, top, root, mode) {
      const pos = new Map(cand.arr.map((i, k) => [i, k + 1]));
      const located = top.map(g => this.locate(cand, g)).filter(Boolean);
      if (located.some(x => !x.exo && TERMINAL.has(x.group.type) && !cand.ring
        && x.at !== cand.arr[0] && x.at !== cand.arr.at(-1))) return null;
      if (mode === 'acyl' && cand.arr[0] !== root) return null;
      const doubles = [], triples = [];
      if (!cand.base) for (let k = 0; k < cand.arr.length - (cand.ring ? 0 : 1); k++) {
        const order = this.bond(cand.arr[k], cand.arr[(k + 1) % cand.arr.length]);
        if (order === 2) doubles.push(k + 1);
        if (order === 3) triples.push(k + 1);
      }
      const suffixLocs = located.map(x => pos.get(x.at)).sort((a, b) => a - b);
      const used = new Set(cand.arr);
      located.forEach(x => { x.group.atoms.forEach(i => used.add(i)); if (x.exo) used.add(x.anchor); });
      const branchLocs = [];
      for (const c of cand.arr) for (const e of this.g[c].nb) if (!used.has(e.n)) branchLocs.push(pos.get(c));
      branchLocs.sort((a, b) => a - b);
      return { ...cand, pos, located, doubles, triples, suffixLocs, branchLocs,
        rank: [-located.length, cand.ring ? -1 : 0, -cand.arr.length, -(doubles.length + triples.length), -doubles.length],
        locRank: [mode === 'radical' ? [pos.get(root)] : suffixLocs,
          [...doubles, ...triples].sort((a, b) => a - b), doubles, branchLocs] };
    }
    compare(a, b) {
      let d = cmp(a.rank, b.rank);
      for (let i = 0; !d && i < a.locRank.length; i++) d = cmp(a.locRank[i], b.locRank[i]);
      return d;
    }
    describe(set, root, mode = 'parent') {
      this.tick();
      if (++this.depth > 64) fail('More than 64 nested substituent levels require another naming method.');
      try {
        const memoKey = mode + ':' + root + ':' + [...set].sort((a, b) => a - b).join(',');
        if (this.memo.has(memoKey)) return this.memo.get(memoKey);
        const groups = this.groups(set);
        const highest = mode === 'parent' ? Math.max(0, ...groups.map(g => SENIORITY[g.type])) : 0;
        const top = mode === 'parent' ? groups.filter(g => SENIORITY[g.type] === highest) : [];
        const cands = this.candidates(set, root).map(c => this.decorate(c, top, root, mode)).filter(Boolean);
        if (!cands.length) fail('No supported parent hydride for this ring or substituent.');
        cands.sort((a, b) => this.compare(a, b));
        const best = cands[0];
        if (top.length && !best.located.length) fail('The principal functional group has no supported parent.');
        // Fully resolve tied paths. In particular, numbering must not depend on
        // the order the user drew atoms or on neighbor insertion order.
        const tied = cands.filter(c => this.compare(c, best) === 0);
        const successes = [], errors = [];
        for (const cand of tied) {
          try { successes.push(this.assemble(set, cand, mode, root)); }
          catch (err) { if (!(err instanceof Unsupported)) throw err; errors.push(err); }
        }
        if (!successes.length) throw errors[0];
        successes.sort((a, b) => a.alpha.localeCompare(b.alpha, 'en', { numeric: true }) || a.name.localeCompare(b.name));
        this.memo.set(memoKey, successes[0]);
        return successes[0];
      } finally { this.depth--; }
    }
    takeBranch(start, parent, set, occupied) {
      const branch = this.component(start, set, occupied);
      if (!branch.size) fail('A substituent overlaps an already named group.');
      const exits = [];
      for (const i of branch) for (const e of this.g[i].nb) if (set.has(e.n) && !branch.has(e.n)) exits.push([i, e.n]);
      if (exits.length !== 1 || exits[0][0] !== start || exits[0][1] !== parent)
        fail('A bridging substituent cannot be expressed by this local naming method.');
      return branch;
    }
    radical(set, start) {
      const a = this.g[start];
      if (a.el === 'C') {
        const co = a.nb.find(e => set.has(e.n) && this.g[e.n].el === 'O' && e.o === 2 && this.g[e.n].nb.length === 1);
        const acidO = a.nb.find(e => set.has(e.n) && this.g[e.n].el === 'O' && e.o === 1 && this.g[e.n].nb.length === 1);
        if (co && acidO && set.size === 3 && !this.g[acidO.n].charge) return 'carboxy';
        const nit = a.nb.find(e => this.g[e.n].el === 'N' && e.o === 3);
        if (nit && set.size === 2) return 'cyano';
        return this.describe(set, start, co ? 'acyl' : 'radical').name;
      }
      if (this.ringAtoms.has(start)) return this.describe(set, start, 'radical').name;
      if (HALO[a.el] && set.size === 1) return HALO[a.el];
      if (a.el === 'N') {
        const oxy = a.nb.filter(e => set.has(e.n) && this.g[e.n].el === 'O');
        if (set.size === 3 && oxy.length === 2 && a.charge === 1
          && oxy.some(e => e.o === 2 && !this.g[e.n].charge)
          && oxy.some(e => e.o === 1 && this.g[e.n].charge === -1)) return 'nitro';
        if (a.charge || a.nb.some(e => e.o !== 1)) fail('Unsupported nitrogen substituent or charge.');
        const names = this.linkedParts(set, start);
        if (!names.length) return 'amino';
        return prefixes(names.map(name => ({ name, loc: 0 })), true) + 'amino';
      }
      if (a.el === 'O' || a.el === 'S') {
        if (a.nb.some(e => e.o !== 1)) fail('Unsupported oxygen or sulfur substituent.');
        const parts = this.linkedParts(set, start);
        if (!parts.length) return a.el === 'O' ? (a.charge ? 'oxido' : 'hydroxy') : 'sulfanyl';
        if (parts.length !== 1 || a.charge) fail('Unsupported heteroatom linkage.');
        const name = parts[0];
        if (a.el === 'S') return groupName(name) + 'sulfanyl';
        if (/^(meth|eth|prop|but|pent|hex|hept|oct|non|dec)yl$/.test(name) || name === 'phenyl') return name.replace(/yl$/, 'oxy');
        return groupName(name) + 'oxy';
      }
      fail('Unsupported substituent environment.');
    }
    linkedParts(set, link) {
      const used = new Set([link]), names = [];
      for (const e of this.g[link].nb) {
        if (!set.has(e.n) || used.has(e.n)) continue;
        if (e.o !== 1) fail('Unsupported multiple bond through a linking heteroatom.');
        const branch = this.takeBranch(e.n, link, set, used);
        names.push(this.radical(branch, e.n)); branch.forEach(i => used.add(i));
      }
      if (used.size !== set.size) fail('The substituent was not completely named.');
      return names;
    }
    assemble(set, c, mode, root) {
      this.tick();
      const g = this.g, used = new Set(c.arr), items = [], alkyls = [];
      const suffixType = c.located[0]?.group.type;
      if (c.located.some(x => x.exo !== c.located[0].exo)) fail('Mixed endocyclic and exocyclic suffix groups require another naming method.');
      const consume = i => { if (used.has(i)) fail('A functional group overlaps the parent.'); used.add(i); };
      for (const x of c.located) {
        const group = x.group;
        group.atoms.forEach(consume);
        if (x.exo) consume(x.anchor);
      }
      if (mode === 'acyl') {
        const oxo = g[root].nb.find(e => set.has(e.n) && g[e.n].el === 'O' && e.o === 2);
        if (!oxo) fail('Acyl carbonyl missing.');
        consume(oxo.n);
      }
      // N/O substituents of suffix groups are named before remaining branches.
      for (const x of c.located) {
        const { group } = x;
        if (group.link === undefined) continue;
        for (const e of g[group.link].nb) {
          if (!set.has(e.n) || used.has(e.n)) continue;
          const branch = this.takeBranch(e.n, group.link, set, used);
          const name = this.radical(branch, e.n);
          branch.forEach(i => used.add(i));
          if (group.type === 'ester') alkyls.push(name);
          else items.push({ name, loc: c.located.length > 1 ? 'N^' + c.pos.get(x.at) : 'N' });
        }
      }
      for (const a of c.arr) for (const e of g[a].nb) {
        if (!set.has(e.n) || used.has(e.n)) continue;
        const branch = this.takeBranch(e.n, a, set, used);
        let name;
        if (e.o === 2 && branch.size === 1 && g[e.n].el === 'O' && !g[e.n].charge) name = 'oxo';
        else if (e.o === 1) name = this.radical(branch, e.n);
        else if (e.o === 2 && branch.size === 1 && g[e.n].el === 'C' && g[e.n].h === 2) name = 'methylidene';
        else fail('This multiply bonded substituent is not supported locally.');
        branch.forEach(i => used.add(i)); items.push({ name, loc: c.pos.get(a) });
      }
      if (used.size !== set.size || [...set].some(i => !used.has(i))) fail('The complete structure could not be represented in the name.');
      let base = c.base || hydro(c.arr.length, c.doubles, c.triples, c.ring);
      const locs = c.suffixLocs, count = locs.length, exo = c.located[0]?.exo;
      if (mode === 'radical') {
        const loc = c.pos.get(root);
        if (c.benzene) base = 'phenyl';
        else if (!c.base && !c.doubles.length && !c.triples.length && loc === 1) base = (c.ring ? 'cyclo' : '') + stem(c.arr.length) + 'yl';
        else if (!c.ring && c.arr.length === 2) base = base.replace(/e$/, 'yl');
        else base = base.replace(/e$/, '') + '-' + loc + '-yl';
      } else if (mode === 'acyl') {
        base = base.replace(/e$/, '') + 'oyl';
      } else if (suffixType) {
        if (c.hetero && c.located.some(x => !x.exo && TERMINAL.has(suffixType))) fail('This heterocyclic suffix requires specialized nomenclature.');
        let suffix = exo ? EXO[suffixType] : SUFFIX[suffixType];
        if (suffixType === 'acylhalide') suffix = 'oyl ' + ({ F: 'fluoride', Cl: 'chloride', Br: 'bromide', I: 'iodide' }[c.located[0].group.halogen]);
        if (!suffix) fail('Unsupported suffix combination.');
        const multiplier = /^[aeiou]/.test(suffix) ? multi(count).replace(/a$/, '') : multi(count);
        suffix = multiplier + suffix;
        if (!exo && /^[aeiou]/.test(suffix)) base = base.replace(/e$/, '');
        const locanted = exo ? (count > 1 || items.length > 0 || c.hetero || c.poly) : !TERMINAL.has(suffixType) && (c.arr.length > 2 || count > 1 || c.hetero);
        base += (locanted ? '-' + locs.join(',') + '-' : '') + suffix;
      }
      const omit = items.every(x => typeof x.loc === 'number') && (c.arr.length === 1
        || (c.ring && !c.hetero && !c.poly && !suffixType && mode === 'parent' && items.length === 1)
        || (!c.ring && c.arr.length === 2 && !suffixType && mode === 'parent' && items.length === 1));
      let name = join(prefixes(items, omit), base);
      const retained = { 'benzen-1-ol': 'phenol', 'benzen-1-amine': 'aniline',
        'benzenecarboxylic acid': 'benzoic acid', 'benzenecarboxylate': 'benzoate',
        'benzenecarbaldehyde': 'benzaldehyde', 'benzenecarboxamide': 'benzamide', 'benzenecarbonitrile': 'benzonitrile' };
      name = retained[name] || name;
      const descriptors=this.stereo.atoms.filter(x=>c.pos.has(x.atom)).map(x=>({loc:c.pos.get(x.atom),label:x.label}));
      for(const b of this.stereo.bonds){
        if(c.pos.has(b.a)&&c.pos.has(b.b)) descriptors.push({loc:Math.min(c.pos.get(b.a),c.pos.get(b.b)),label:b.label});
      }
      descriptors.sort((a,b)=>a.loc-b.loc||a.label.localeCompare(b.label));
      if(descriptors.length) name='('+descriptors.map(x=>x.loc+x.label).join(',')+')-'+name;
      if (alkyls.length) {
        if (alkyls.length !== count || new Set(alkyls).size !== 1) fail('Mixed ester substituents require an additional naming method.');
        name = multi(alkyls.length, compound(alkyls[0])) + groupName(alkyls[0]) + ' ' + name;
      }
      const alpha = [...items].sort((a, b) => alphabetical(a.name).localeCompare(alphabetical(b.name)) || String(a.loc).localeCompare(String(b.loc), 'en', { numeric: true }))
        .map(x => alphabetical(x.name) + ':' + String(x.loc).padStart(5, '0')).join(';');
      return { name, alpha, coveredAtoms: used.size };
    }
    run() {
      const set = new Set(this.g.map((_, i) => i));
      const result = this.specialParent(set) || this.describe(set);
      // Conservative guard: never return a connectivity-only name when a
      // stereodescriptor falls outside the selected parent/substituent paths.
      if((result.name.match(/\d+[RSrsEZ](?=[,)])/g)||[]).length!==this.stereo.atoms.length+this.stereo.bonds.length)
        fail('Not every specified stereogenic unit can be represented by the current naming rules.');
      return { status: 'ok', name: result.name, coveredAtoms: result.coveredAtoms,
        atomCount: this.g.length, nomenclature: this.stereo.atoms.length||this.stereo.bonds.length?'systematic-stereochemical':'systematic-connectivity',
        notes: ['Offline systematic name. Only explicitly specified stereochemistry is assigned; preferred IUPAC names may differ.'] };
    }
    specialParent(set) {
      const g = this.g;
      if (g.length === 1 && !g[0].charge) {
        const hydride = { 'N3': 'ammonia', 'O2': 'water', 'S2': 'hydrogen sulfide' }[g[0].el + g[0].h];
        if (hydride) return { name: hydride, coveredAtoms: 1 };
        const halide = { F: 'fluoride', Cl: 'chloride', Br: 'bromide', I: 'iodide' }[g[0].el];
        if (halide && g[0].h === 1) return { name: 'hydrogen ' + halide, coveredAtoms: 1 };
      }
      if (g.length === 2 && g[0].el === g[1].el && g.every(a => !a.charge)) {
        const simple = { O1: 'hydrogen peroxide', O2: 'dioxygen', N3: 'dinitrogen', N1: 'hydrazine' }[g[0].el + this.bond(0, 1)];
        if (simple) return { name: simple, coveredAtoms: 2 };
      }
      const carbonDioxide = g.length === 3 && g.some(a => a.el === 'C' && a.nb.length === 2
        && a.nb.every(e => g[e.n].el === 'O' && e.o === 2));
      if (carbonDioxide) return { name: 'carbon dioxide', coveredAtoms: 3 };
      const sulfuric = g.length === 5 && g.some(a => a.el === 'S' && a.nb.length === 4
        && a.nb.every(e => g[e.n].el === 'O' && g[e.n].nb.length === 1)
        && a.nb.filter(e => e.o === 2).length === 2);
      if (sulfuric) {
        const h = g.reduce((sum, a) => sum + a.h, 0);
        return { name: h === 2 ? 'sulfuric acid' : h === 1 ? 'hydrogen sulfate' : 'sulfate', coveredAtoms: 5 };
      }
      const nitric = g.length === 4 && g.some(a => a.el === 'N' && a.charge === 1 && a.nb.length === 3
        && a.nb.every(e => g[e.n].el === 'O' && g[e.n].nb.length === 1)
        && a.nb.filter(e => e.o === 2).length === 1);
      if (nitric) {
        const h = g.reduce((sum, a) => sum + a.h, 0);
        return { name: h === 1 ? 'nitric acid' : 'nitrate', coveredAtoms: 4 };
      }
      const phosphorus = [...set].filter(i => g[i].el === 'P');
      if (phosphorus.length) {
        if (phosphorus.length !== 1) fail('Polyphosphates are not supported locally.');
        const p = phosphorus[0], nb = g[p].nb;
        const oxo = nb.filter(e => g[e.n].el === 'O' && e.o === 2 && g[e.n].nb.length === 1);
        const single = nb.filter(e => g[e.n].el === 'O' && e.o === 1 && !g[e.n].charge);
        if (nb.length !== 4 || oxo.length !== 1 || single.length !== 3) fail('Only neutral orthophosphate acids and esters are supported locally.');
        const used = new Set([p, ...nb.map(e => e.n)]), alkyls = []; let h = 0;
        for (const e of single) {
          const tail = g[e.n].nb.find(x => x.n !== p);
          if (!tail) { h++; continue; }
          const branch = this.takeBranch(tail.n, e.n, set, used);
          alkyls.push(this.radical(branch, tail.n)); branch.forEach(i => used.add(i));
        }
        if (used.size !== set.size) fail('Incomplete phosphate structure.');
        const parts = new Map(); alkyls.sort().forEach(n => parts.set(n, (parts.get(n) || 0) + 1));
        const words = [...parts].map(([n, count]) => multi(count, compound(n)) + groupName(n));
        if (h) words.push((h > 1 ? multi(h) : '') + 'hydrogen');
        return { name: h === 3 ? 'phosphoric acid' : words.join(' ') + ' phosphate', coveredAtoms: used.size };
      }
      const anhydrides = [...set].filter(i => g[i].el === 'O' && g[i].nb.length === 2
        && g[i].nb.every(e => e.o === 1 && g[e.n].el === 'C' && g[e.n].nb.some(x => g[x.n].el === 'O' && x.o === 2)));
      if (anhydrides.length) {
        if (anhydrides.length !== 1 || this.ringAtoms.has(anhydrides[0])) fail('Cyclic or repeated anhydrides require another naming method.');
        const o = anhydrides[0], used = new Set([o]), acids = [];
        for (const e of g[o].nb) {
          const branch = this.takeBranch(e.n, o, set, used);
          acids.push(this.describe(branch, e.n, 'acyl').name.replace(/oyl$/, 'oic'));
          branch.forEach(i => used.add(i));
        }
        if (used.size !== set.size) fail('Incomplete anhydride structure.');
        return { name: [...new Set(acids)].sort().join(' ') + ' anhydride', coveredAtoms: used.size };
      }
      return null;
    }
  }
  function name(graph, options = {}) {
    try { return new Engine(graph, options).run(); }
    catch (error) {
      if (error instanceof Unsupported) return { status: 'unsupported', reason: error.message, atomCount: graph.length };
      throw error;
    }
  }
  return Object.freeze({ name, stem, version: '2.1.0' });
});
