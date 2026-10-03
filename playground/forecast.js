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
    this.show('Searching for the path…', '<div class="fc-busy"><i></i></div><p class="fc-note">Pulling the two atoms along the reaction in a private copy, relaxing everything else at each step. The simulation keeps running.</p>');
    if (!this.worker) {
      this.worker = new Worker('kinetics-worker.js?v=' + (this.actions.version || ''));
      this.worker.onmessage = ev => { if (ev.data.id === this.req) this.render(ev.data.result); };
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
      <table class="fc-table"><thead><tr><th>Temperature</th><th>Half-life</th></tr></thead><tbody>${rows}</tbody></table>
      ${alt ? '<p class="fc-sub">Also possible, more slowly</p><ul class="fc-alt">' + alt + '</ul>' : ''}
      <details class="quantum-notes"><summary>How this is estimated</summary>
      <p>The barrier is measured on this simulation's own energy surface: the two atoms are pulled through the reaction in a private copy while every other atom relaxes. Transfers are pulled along the difference between the bond forming and the bond breaking, and the lowest of the paths tried is kept.</p>
      <p>The time comes from transition-state theory: rate = A·e<sup>−barrier/RT</sup>. A is about 10<sup>15.5</sup> s⁻¹ for a bond simply coming apart, kT/h otherwise, and for two molecules the collision rate × 0.1. The range allows ±10 kJ/mol on the barrier and ×/÷10 on A.</p>
      <p>It is only as good as the model: where the model's barrier is wrong, so is the time. Known errors: aromatic and vinylic C–H bonds are 50–120 kJ/mol too weak, Cl + CH₄ comes out 20 kJ/mol too slow. Above about 2000 K the live simulation also decomposes hydrocarbons much faster than these numbers say; the forecast is the more trustworthy of the two there.</p></details>`;
    this.show(side(r.reactants) + ' → ' + side(r.products), body);
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
    this.req++; this.ids = [];
    this.panel.classList.remove('open'); this.panel.inert = true;
    this.actions.focus?.();
  }
}
root.ForecastCard = ForecastCard;
})(typeof self !== 'undefined' ? self : this);
