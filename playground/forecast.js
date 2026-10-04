(function (root) {
'use strict';
const K = () => root.ChemKinetics;

class ForecastCard {
  constructor(actions) {
    this.actions = actions; this.worker = null; this.ids = []; this.req = 0; this.fx = null;
  }
  get isOpen() { return false; }
  picked() {
    const eng = this.actions.engine;
    return this.ids.map(id => eng.indexOfId(id)).filter(i => i >= 0);
  }
  choose(i) {
    const eng = this.actions.engine, id = eng.ids[i];
    if (this.ids.length === 2 && this.result && this.ids.includes(id)) { this.skipProbe(); return; }
    if (this.ids.length >= 2) { this.ids = []; this.result = null; if (this.fx) this.fx.probe = null; }
    if (this.ids[0] === id) return;
    this.ids.push(id);
    if (this.ids.length === 1) { this.req++; this.actions.pause?.(); return; }
    this.run();
  }
  run() {
    const eng = this.actions.engine, [i, j] = this.picked();
    if (i === undefined || j === undefined) { this.close(); return; }
    const scene = eng.toJSON(), id = ++this.req, pair = this.ids.slice();
    this.sceneIds = Array.from(eng.ids.subarray(0, eng.N)); this.result = null;
    this.fx = { ...(this.fx || {}), scan: null, probe: null, focus: pair };
    if (!this.worker) {
      this.worker = new Worker('kinetics-worker.js?v=' + (this.actions.version || ''));
      this.worker.onmessage = ev => { if (ev.data.id === this.req) this.render(ev.data.result, pair); };
      this.worker.onerror = () => this.say('the forecast could not start in this browser');
    }
    this.worker.postMessage({ id, scene, i, j });
  }
  render(r, pair) {
    const H = K().humanTime, pretty = this.actions.pretty || (x => x);
    this.fx = { ...(this.fx || {}), focus: null };
    if (!r.ok) { this.ids = []; this.fx.fail = { atoms: pair, t0: performance.now() }; this.say(r.reason); return; }
    this.result = r;
    const h = r.now.halfLife, never = !isFinite(h) || h > 3.15e16;
    this.fx.probe = {
      atoms: pair, t0: performance.now(), path: r.path, Ea: r.Ea,
      text: never ? 'never' : '½ ' + H(h),
      sub: r.products.map(pretty).join(' + '),
      skip: !!r.event && !never
    };
  }
  skipProbe() {
    const res = this.result, ids = this.sceneIds, pair = this.ids.slice();
    if (!res || !res.event || !(res.now.halfLife < 3.15e16)) return;
    this.close();
    const side = list => list.map(this.actions.pretty || (x => x)).join(' + ');
    const wait = -Math.log(1 - Math.random()) * res.now.halfLife / Math.LN2;
    if (this.actions.skip?.(res.event, wait, ids, side(res.reactants) + ' → ' + side(res.products), res.reactants, true) !== false) this.jumped(pair, wait);
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
    this.actions.pause?.();
    const scene = eng.toJSON(), ids = Array.from(eng.ids.subarray(0, eng.N)), T = scene.T || 298;
    const { list } = Kn.candidates(eng);
    const fx = this.fx = { scan: list.map(c => ({ a: ids[c.i], b: ids[c.j], done: false, Ea: null })), jump: this.fx && this.fx.jump, t0: performance.now() };
    if (!list.length) { this.say('nothing here can react'); this.finish(); return; }
    const b = eng.box, V = (b.x1 - b.x0) * (b.y1 - b.y0) * (b.z1 - b.z0) * 1e-24;
    const memo = this.memo || (this.memo = new Map()), mk = c => Math.round(T) + '|' + c.key;
    this.lastScene = { scene, key, T };
    const mark = (n, r) => { if (n === undefined || !fx.scan) return; fx.scan[n].done = true; fx.scan[n].Ea = r && r.ok !== false ? r.Ea : null; };
    const job = (c, quick, n) => this.dispatch({ type: 'scan', key: 'sv' + key, scene, c, T, quick }).then(r => { mark(n, r); if (r && r.ok !== false) { r.cand = c; const m = memo.get(mk(c)); if (!m || !m.refined || r.refined) memo.set(mk(c), { ...r, event: null, cand: null, also: (r.also || []).map(a => ({ ...a, event: null })) }); return r; } return null; });
    const quick = await Promise.all(list.map((c, n) => { const m = memo.get(mk(c)); if (m) { mark(n, m); return Promise.resolve(Kn.reuse(m, c, V)); } return job(c, true, n); }));
    if (key !== this.req) return;
    const order = quick.map((r, n) => ({ r, n })).filter(x => x.r && !x.r.refined).sort((a, c) => a.r.Ea - c.r.Ea).slice(0, Kn.REFINE);
    const refined = await Promise.all(order.map(x => job(list[x.n], false)));
    if (key !== this.req) return;
    order.forEach((x, m) => { if (refined[m]) quick[x.n] = refined[m]; });
    this.decide(Kn.combine(quick, T), ids);
  }
  say(text) { this.fx = { ...(this.fx || {}), scan: null, focus: null, msg: { text, t0: performance.now() } }; }
  finish() { if (this.auto) { this.auto = false; clearInterval(this.timer); this.actions.autoChanged?.(false); } }
  stop() { this.req++; this.finish(); if (this.fx) { this.fx.scan = null; this.fx.focus = null; } }
  toggleAuto() { if (this.auto) { this.stop(); return; } this.auto = true; this.actions.autoChanged?.(true); this.next(); }
  async decide(sv, ids) {
    const H = K().humanTime, wait = sv.total > 0 ? Math.LN2 / sv.total : Infinity;
    if (!isFinite(wait) || wait > 3.15e16 || !sv.events.length) {
      this.say((sv.events[0] ? 'stable · nothing for ' + H(sv.events[0].halfLife) : 'nothing here can react') + (this.actions.absorbs?.() ? ' · L for light' : ''));
      this.finish(); return;
    }
    const key = this.req, ok = await this.go(sv, ids);
    if (key !== this.req) return;
    if (this.fx) { this.fx.scan = null; this.fx.focus = null; }
    if (ok && this.auto) this.watchThenNext(); else this.finish();
  }
  async go(sv, ids) {
    const ls = this.lastScene, key = this.req;
    const same = (x, want) => x.products.map(f => f.replace(/·/g, '')).sort().join() === want.products.map(f => f.replace(/·/g, '')).sort().join();
    const failed = new Set();
    let pair = null;
    for (let tries = 0; tries < 4; tries++) {
      const nx = K().pickNext(sv);
      if (!nx) return false;
      if (failed.has(nx.pick)) break;
      let pick = nx.pick;
      pair = [ids[pick.i], ids[pick.j]];
      if (this.fx) { this.fx.scan = null; this.fx.focus = pair; }
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
      if (!good) { failed.add(nx.pick); continue; }
      const side = list => list.map(this.actions.pretty || (x => x)).join(' + ');
      const ok = this.actions.skip?.(pick.event, nx.wait, ids, side(pick.reactants) + ' → ' + side(pick.products), pick.reactants, !this.auto);
      if (ok === false) break;
      this.jumped(pair, nx.wait);
      return true;
    }
    if (pair) this.fx = { ...(this.fx || {}), scan: null, focus: null, fail: { atoms: pair, t0: performance.now() } };
    return false;
  }
  jumped(pair, wait) {
    this.fx = { ...(this.fx || {}), scan: null, focus: null, jump: { atoms: pair, text: '+' + K().humanTime(wait), t0: performance.now() } };
  }
  watchThenNext() {
    const eng = this.actions.engine, t0 = eng.time, key = this.req, span = 3000;
    clearInterval(this.timer);
    this.timer = setInterval(() => {
      if (!this.auto || key !== this.req) { clearInterval(this.timer); return; }
      if (eng.time >= t0 + span) { clearInterval(this.timer); this.next(); }
    }, 200);
  }
  close() {
    this.req++; this.ids = []; this.result = null; this.finish();
    if (this.fx) { this.fx.probe = null; this.fx.focus = null; }
    this.actions.focus?.();
  }
}
root.ForecastCard = ForecastCard;
})(typeof self !== 'undefined' ? self : this);
