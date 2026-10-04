(function (root) {
'use strict';
const K = () => root.ChemKinetics;
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

class ForecastCard {
  constructor(actions) {
    this.actions = actions; this.worker = null; this.ids = []; this.req = 0; this.anchor = null;
    this.panel = document.createElement('aside');
    this.panel.className = 'atom-inspector forecast-card'; this.panel.inert = true;
    this.panel.setAttribute('role', 'dialog'); this.panel.setAttribute('aria-label', 'Reaction forecast');
    this.panel.innerHTML = `<div class="inspector-head"><div><p class="eyebrow">FORECAST</p><h2 id="fcTitle"></h2></div>
      <div class="inspector-acts"><button class="iconbtn" id="fcClose" aria-label="Close forecast" title="Close"><svg viewBox="0 0 16 16"><path d="M4.6 4.6l6.8 6.8M11.4 4.6l-6.8 6.8"/></svg></button></div></div>
      <div id="fcBody"></div>`;
    document.body.append(this.panel);
    this.el = id => this.panel.querySelector('#' + id);
    this.el('fcClose').onclick = () => this.close();
    this.panel.addEventListener('keydown', e => { if (e.key === 'Escape') { e.stopPropagation(); this.close(); } });
  }
  get isOpen() { return this.panel.classList.contains('open'); }
  picked() {
    const eng = this.actions.engine;
    return this.ids.map(id => eng.indexOfId(id)).filter(i => i >= 0);
  }
  choose(i) {
    const eng = this.actions.engine, id = eng.ids[i];
    if (this.ids.length >= 2) this.ids = [];
    if (this.ids[0] === id) return;
    this.ids.push(id);
    if (this.ids.length === 1) { this.req++; this.panel.classList.remove('open'); this.panel.inert = true; return; }
    this.run();
  }
  run() {
    const eng = this.actions.engine, [i, j] = this.picked();
    if (i === undefined || j === undefined) { this.close(); return; }
    const scene = eng.toJSON(), id = ++this.req;
    this.sceneIds = Array.from(eng.ids.subarray(0, eng.N)); this.result = null;
    this.show('Searching for the path…', '<div class="fc-busy"><i></i></div><p class="fc-note">Pulling the two atoms along the reaction in a private copy, relaxing everything else at each step. The simulation keeps running.</p>');
    if (!this.worker) {
      this.worker = new Worker('kinetics-worker.js?v=' + (this.actions.version || ''));
      this.worker.onmessage = ev => { if (ev.data.id === this.req) { this.result = ev.data.result; this.render(ev.data.result); } };
      this.worker.onerror = () => this.show('Forecast unavailable', '<p class="fc-note">The calculation could not start in this browser.</p>');
    }
    this.worker.postMessage({ id, scene, i, j });
  }
  show(title, body) {
    this.el('fcTitle').textContent = title;
    this.el('fcBody').innerHTML = body;
    this.panel.classList.add('open'); this.panel.inert = false;
  }
  render(r) {
    const pretty = this.actions.pretty || (x => x), side = list => list.map(pretty).join(' + ');
    if (!r.ok) { this.show('No reaction here', '<p class="fc-note">' + esc(r.reason) + '.</p>'); return; }
    const H = K().humanTime, now = r.now, T = Math.round(r.T);
    const what = r.kind === 'break' ? 'bond breaks' : r.channel.type === 'transfer' ? 'atom transfer' : r.bimolecular ? 'the two join' : 'ring closes';
    const verdict = !isFinite(now.halfLife) || now.halfLife > 3.15e16 ? 'Not on any human timescale' : now.halfLife < 1e-9 ? 'Practically at once' : 'Half-life ' + H(now.halfLife);
    const prof = this.profile(r.path, r.Ea);
    const rows = r.table.map(x => '<tr' + (Math.abs(x.T - T) < 1 ? ' class="now"' : '') + '><td>' + x.T + ' K</td><td>' + H(x.halfLife) + '</td></tr>').join('');
    const alt = (r.alternatives || []).slice(0, 2).map(a => '<li>' + esc(side(a.products)) + ' <span class="m">barrier ' + a.Ea.toFixed(0) + ' kJ/mol</span></li>').join('');
    const where = r.bimolecular ? 'at ' + T + ' K, for this pair at the density in this chamber' : 'at ' + T + ' K';
    const body = `<p class="fc-what">${esc(what)}</p>
      <div class="fc-big"><b>${esc(verdict)}</b><span>${esc(where)}</span></div>
      ${isFinite(now.halfLife) && now.halfLife < 3.15e16 ? `<p class="fc-range">likely between <b>${esc(H(now.fastest))}</b> and <b>${esc(H(now.slowest))}</b></p>` : ''}
      <dl class="fc-stats"><div><dt>Barrier</dt><dd>${r.Ea.toFixed(0)} kJ/mol</dd></div><div><dt>${r.dE <= 0 ? 'Releases' : 'Costs'}</dt><dd>${Math.abs(r.dE).toFixed(0)} kJ/mol</dd></div></dl>
      ${prof}
      ${r.event && isFinite(now.halfLife) && now.halfLife < 1e30 ? `<button class="fc-skip" id="fcSkip">Skip the wait · watch it happen</button><p class="fc-skipnote">Jumps the clock by a waiting time drawn at random from this rate, then plays the reaction itself in the live simulation from just past the top of the barrier.</p>` : ''}
      <table class="fc-table"><thead><tr><th>Temperature</th><th>Half-life</th></tr></thead><tbody>${rows}</tbody></table>
      ${alt ? '<p class="fc-sub">Also possible, more slowly</p><ul class="fc-alt">' + alt + '</ul>' : ''}
      <details class="quantum-notes"><summary>How this is estimated</summary>
      <p>The barrier is measured on this simulation's own energy surface: the two atoms are pulled through the reaction in a private copy while every other atom relaxes. Transfers are pulled along the difference between the bond forming and the bond breaking, and the lowest of the paths tried is kept.</p>
      <p>The time comes from transition-state theory: rate = A·e<sup>−barrier/RT</sup>. A is about 10<sup>15.5</sup> s⁻¹ for a bond simply coming apart, kT/h otherwise, and for two molecules the collision rate × 0.1. The range allows ±10 kJ/mol on the barrier and ×/÷10 on A.</p>
      <p>It is only as good as the model: where the model's barrier is wrong, so is the time. Known errors: C–H bonds next to another carbon are 30–40 kJ/mol too weak, hydrogen taken by Cl· costs about 20 kJ/mol too much, and CH₃· adds to C=C 20 kJ/mol too easily. Above about 2000 K the live simulation also decomposes hydrocarbons much faster than these numbers say; the forecast is the more trustworthy of the two there.</p></details>`;
    this.show(side(r.reactants) + ' → ' + side(r.products), body);
    const skip = this.el('fcSkip');
    if (skip) skip.onclick = () => {
      const res = this.result, ids = this.sceneIds; this.close();
      const side = list => list.map(this.actions.pretty || (x => x)).join(' + ');
      this.actions.skip?.(res.event, -Math.log(1 - Math.random()) * res.now.halfLife / Math.LN2, ids, side(res.reactants) + ' → ' + side(res.products), res.reactants, true);
    };
  }
  pool() {
    if (this.workers) return this.workers;
    const n = Math.max(1, Math.min(6, (navigator.hardwareConcurrency || 2) - 1));
    this.pending = new Map(); this.seq = 0;
    this.workers = Array.from({ length: n }, () => {
      const w = new Worker('kinetics-worker.js?v=' + (this.actions.version || ''));
      w.onmessage = ev => { const p = this.pending.get(ev.data.id); if (p) { this.pending.delete(ev.data.id); p(ev.data.result); } };
      w.busy = 0; return w;
    });
    return this.workers;
  }
  dispatch(msg) {
    const ws = this.pool(), w = ws.reduce((a, b) => (b.busy < a.busy ? b : a)), id = 'j' + (++this.seq);
    w.busy++;
    return new Promise(res => { this.pending.set(id, r => { w.busy--; res(r); }); w.postMessage({ ...msg, id }); });
  }
  async next() {
    const eng = this.actions.engine, Kn = K();
    if (!eng.N) return;
    if (eng.needRebuild || eng.needForces) eng.refresh();
    const key = ++this.req; this.ids = [];
    if (!this.auto) this.log = [];
    this.actions.pause?.();
    const scene = eng.toJSON(), ids = Array.from(eng.ids.subarray(0, eng.N)), T = scene.T || 298;
    const { list, dropped } = Kn.candidates(eng);
    if (!list.length) { this.show('Nothing here can react', '<p class="fc-note">No radical, no π bond and no bond to break was found.</p>'); return; }
    const b = eng.box, V = (b.x1 - b.x0) * (b.y1 - b.y0) * (b.z1 - b.z0) * 1e-24;
    const memo = this.memo || (this.memo = new Map()), mk = c => Math.round(T) + '|' + c.key;
    const fresh = list.filter(c => !(memo.get(mk(c)) || {}).refined);
    const total = Math.max(1, fresh.length + Math.min(Kn.REFINE, fresh.length));
    let done = 0;
    const bar = () => '<div class="fc-bar"><i style="width:' + (100 * done / total).toFixed(0) + '%"></i></div><p class="fc-note">Forecasting ' + list.length + ' possible reaction' + (list.length > 1 ? 's' : '') + ' in this chamber at ' + Math.round(T) + ' K' + (fresh.length < list.length ? ', ' + (list.length - fresh.length) + ' already known' : '') + ', ' + this.pool().length + ' at a time. The simulation is paused until you choose.</p>';
    this.show(this.auto ? 'Keeping going…' : 'What happens next?', bar() + this.logHtml());
    this.panel.style.left = Math.max(12, innerWidth - (this.panel.offsetWidth || 330) - 24) + 'px'; this.panel.style.top = '76px';
    this.lastScene = { scene, key, T };
    const job = (c, quick) => this.dispatch({ type: 'scan', key: 'sv' + key, scene, c, T, quick }).then(r => { done++; if (key === this.req) this.el('fcBody').innerHTML = bar() + this.logHtml(); if (r && r.ok !== false) { r.cand = c; const m = memo.get(mk(c)); if (!m || !m.refined || r.refined) memo.set(mk(c), { ...r, event: null, cand: null, also: (r.also || []).map(a => ({ ...a, event: null })) }); return r; } return null; });
    const quick = await Promise.all(list.map(c => { const m = memo.get(mk(c)); return m ? Promise.resolve(Kn.reuse(m, c, V)) : job(c, true); }));
    if (key !== this.req) return;
    const order = quick.map((r, n) => ({ r, n })).filter(x => x.r && !x.r.refined).sort((a, c) => a.r.Ea - c.r.Ea).slice(0, Kn.REFINE);
    const refined = await Promise.all(order.map(x => job(list[x.n], false)));
    if (key !== this.req) return;
    order.forEach((x, m) => { if (refined[m]) quick[x.n] = refined[m]; });
    this.renderSurvey(Kn.combine(quick, T), ids, dropped);
  }
  logHtml() {
    if (!this.log || !this.log.length) return '';
    const H = K().humanTime;
    return '<p class="fc-sub">So far</p><ol class="fc-log">' + this.log.slice(-8).map(x => '<li><span class="m">+' + esc(H(x.wait)) + '</span> ' + esc(x.label) + '</li>').join('') + '</ol>';
  }
  async go(sv, ids) {
    const side = list => list.map(this.actions.pretty || (x => x)).join(' + ');
    const ls = this.lastScene, key = this.req;
    const same = (x, want) => x.products.map(f => f.replace(/·/g, '')).sort().join() === want.products.map(f => f.replace(/·/g, '')).sort().join();
    let events = sv.events.slice(), tries = 0;
    while (events.length && tries++ < 4) {
      const total = events.reduce((a, x) => a + x.rate, 0);
      const nx = K().pickNext({ ...sv, events: events.map(x => ({ ...x, share: x.rate / total })), total });
      if (!nx) return false;
      let pick = nx.pick;
      this.show('Setting up the reaction…', '<div class="fc-busy"><i></i></div>' + this.logHtml());
      if (pick.cand && ls && !(pick.event && pick.refined)) {
        const r = await this.dispatch({ type: 'scan', key: 'sv' + ls.key, scene: ls.scene, c: pick.cand, T: ls.T, quick: false });
        if (key !== this.req) return false;
        pick = r && r.ok !== false ? [r, ...(r.also || [])].find(x => same(x, pick)) || null : null;
      }
      let good = !!(pick && pick.event);
      if (good && ls) {
        const n = await this.dispatch({ type: 'verify', scene: ls.scene, event: pick.event, products: pick.products, T: ls.T });
        if (key !== this.req) return false;
        good = n >= 2;
      }
      if (!good) { events = events.filter(x => x !== nx.pick); continue; }
      const label = side(pick.reactants) + ' → ' + side(pick.products);
      const ok = this.actions.skip?.(pick.event, nx.wait, ids, label, pick.reactants, !this.auto);
      if (ok === false) { this.failNote = 'The molecules changed before the reaction could be played. Search again.'; return false; }
      (this.log = this.log || []).push({ wait: nx.wait, label });
      return true;
    }
    this.failNote = 'The likeliest reactions did not end as forecast when tried on a copy, so nothing was skipped.';
    return false;
  }
  showWatching() {
    const n = this.log.length;
    this.show('Keeping going · ' + n + ' reaction' + (n > 1 ? 's' : ''), '<p class="fc-note" id="fcWatch">Watching the reaction</p><label class="fc-auto"><input type="checkbox" id="fcAuto" checked> Keep going on its own</label>' + this.logHtml());
    this.el('fcAuto').onchange = ev => { this.auto = ev.target.checked; if (!this.auto) clearInterval(this.timer); };
    this.watchThenNext();
  }
  watchThenNext() {
    const eng = this.actions.engine, t0 = eng.time, key = this.req, span = 3000;
    clearInterval(this.timer);
    this.timer = setInterval(() => {
      if (!this.auto || key !== this.req) { clearInterval(this.timer); return; }
      const left = Math.max(0, t0 + span - eng.time);
      const el = this.el('fcWatch'); if (el) el.textContent = (left / 1000).toFixed(1) + ' ps of the reaction left to watch before the next search';
      if (left <= 0) { clearInterval(this.timer); this.next(); }
    }, 200);
  }
  renderSurvey(sv, ids, dropped) {
    const H = K().humanTime, pretty = this.actions.pretty || (x => x), side = list => list.map(pretty).join(' + ');
    const wait = sv.total > 0 ? Math.LN2 / sv.total : Infinity, never = !isFinite(wait) || wait > 3.15e16;
    if (this.auto && !never && sv.events.length) {
      const key = this.req;
      this.go(sv, ids).then(ok => { if (ok) this.showWatching(); else if (key === this.req) { this.auto = false; this.renderSurvey(sv, ids, dropped); } });
      return;
    }
    if (this.auto && never) this.auto = false;
    const rows = sv.events.slice(0, 6).map(x => '<tr><td>' + esc(side(x.reactants)) + ' → ' + esc(side(x.products)) + '</td><td>' + (x.share >= 0.001 ? (100 * x.share).toFixed(x.share > 0.1 ? 0 : 1) + '%' : '<0.1%') + '</td></tr>').join('');
    const note = this.failNote ? '<p class="fc-note">' + esc(this.failNote) + '</p>' : '';
    this.failNote = null;
    const body = `${note}<div class="fc-big"><b>${never ? 'Nothing on any human timescale' : 'Next reaction in about ' + esc(H(wait))}</b><span>${never ? 'the fastest step found would take ' + esc(H(sv.events[0] ? sv.events[0].halfLife : Infinity)) : 'half-life of the whole chamber at ' + Math.round(sv.T) + ' K'}</span></div>
      ${never || !sv.events.length ? '' : '<button class="fc-skip" id="fcNext">Skip to it · watch it happen</button><label class="fc-auto"><input type="checkbox" id="fcAuto"> Keep going on its own</label><p class="fc-skipnote">Which reaction comes first, and when, is drawn at random from these rates, the way it would be in a real flask. Each reaction you skip to plays out in the live simulation; kept going, it watches each for 3 ps and then looks for the next.</p>'}
      ${this.logHtml()}
      <table class="fc-table fc-list"><thead><tr><th>Possible reaction</th><th>Chance next</th></tr></thead><tbody>${rows}</tbody></table>
      ${dropped ? '<p class="fc-note">' + dropped + ' more candidates were not checked.</p>' : ''}
      <details class="quantum-notes"><summary>What is searched</summary><p>A radical meeting any atom of another molecule, a π bond meeting the hydrogens of another molecule, and every distinct bond coming apart. Each is forecast as with the hourglass tool, counted as often as it occurs in the chamber, and its share is its rate over the total.</p><p>Not yet searched: two closed-shell molecules reacting with each other through their π systems (Diels–Alder, ene reactions), rearrangements inside one molecule, and anything ionic.</p></details>`;
    this.show('What happens next at ' + Math.round(sv.T) + ' K', body);
    const b = this.el('fcNext'), auto = this.el('fcAuto');
    if (b) b.onclick = () => {
      const key = this.req, fail = () => { if (key !== this.req) return; this.auto = false; this.renderSurvey(sv, ids, dropped); };
      if (auto && auto.checked) { this.auto = true; this.log = []; this.go(sv, ids).then(ok => ok ? this.showWatching() : fail()); return; }
      this.go(sv, ids).then(ok => ok ? this.close() : fail());
    };
  }
  profile(path, Ea) {
    if (!path || path.length < 2) return '';
    const lo = Math.min(0, ...path), hi = Math.max(1, ...path), W = 280, Hh = 70, pad = 6;
    const pts = path.map((e, k) => [pad + (W - 2 * pad) * k / (path.length - 1), pad + (Hh - 2 * pad) * (1 - (e - lo) / (hi - lo))]);
    const top = path.indexOf(Math.max(...path)), zero = pad + (Hh - 2 * pad) * (1 - (0 - lo) / (hi - lo));
    return `<svg class="fc-prof" viewBox="0 0 ${W} ${Hh}" aria-label="Energy along the reaction"><line x1="${pad}" x2="${W - pad}" y1="${zero}" y2="${zero}"/><polyline points="${pts.map(p => p.map(v => v.toFixed(1)).join(',')).join(' ')}"/>${Ea > 0.5 ? `<circle cx="${pts[top][0].toFixed(1)}" cy="${pts[top][1].toFixed(1)}" r="3"/>` : ''}</svg><p class="fc-axis"><span>reactants</span><span>energy along the path</span><span>products</span></p>`;
  }
  place(sx, sy, reach = 0) {
    const p = this.panel, w = p.offsetWidth || 320, h = p.offsetHeight || 380, gap = 30 + reach;
    let x = sx + gap, y = sy - h / 2;
    if (x + w > innerWidth - 12) x = sx - gap - w;
    x = Math.max(12, Math.min(innerWidth - w - 12, x));
    y = Math.max(12, Math.min(innerHeight - h - 12, y));
    p.style.left = x + 'px'; p.style.top = y + 'px';
  }
  close() {
    this.req++; this.ids = []; this.auto = false; clearInterval(this.timer);
    this.panel.classList.remove('open'); this.panel.inert = true;
    this.actions.focus?.();
  }
}
root.ForecastCard = ForecastCard;
})(typeof self !== 'undefined' ? self : this);
