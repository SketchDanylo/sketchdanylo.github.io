/* Chem Playground — application: time loop, tools, conditions, molecules, keybinds. */
(function () {
'use strict';
const { Engine, ELEMENTS, BY_SYM, KB } = ChemEngine;
const $ = id => document.getElementById(id);
const store = {
  get(k, d) { try { const v = localStorage.getItem('cp.' + k); return v == null ? d : JSON.parse(v); } catch (e) { return d; } },
  set(k, v) { try { localStorage.setItem('cp.' + k, JSON.stringify(v)); } catch (e) { } }
};
const clamp = (x, a, b) => Math.max(a, Math.min(b, x));
const SUB = { 0: '₀', 1: '₁', 2: '₂', 3: '₃', 4: '₄', 5: '₅', 6: '₆', 7: '₇', 8: '₈', 9: '₉' };
const pretty = f => f.replace(/(\d+)(?=[+−-]$)/, m => m).replace(/([A-Za-z\)])(\d+)/g, (m, a, d) => a + d.split('').map(c => SUB[c]).join(''));
const esc = s => String(s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const STP_T = 298.15, ATM = 1.01325;
const SPS1 = 20000;                 // 1.0× = 20 000 steps (20 ps) of 1 fs per real second
const SPAN_MIN = 0.5, SPAN_MAX = 50; // visible width 0.05 – 5 nm

/* ======================= engine ======================= */
const savedBox = store.get('box', { w: 44, h: 26, d: 12 });
const eng = new Engine({ width: savedBox.w, height: savedBox.h, depth: savedBox.d, T: store.get('T', STP_T) });
eng.thermostat = store.get('thermostat', true);
eng.recording = true;

const canvas = $('field');
const R = new FieldRenderer(canvas);
R.mode = store.get('mode', 'density');
let rp = new Float64Array(0);        // interpolated render positions

/* ======================= time ======================= */
const time = { playing: false, speed: store.get('speed', 1), acc: 0, last: performance.now(), rateEMA: 0, stepsWindow: 0, windowStart: performance.now(), limited: false, manual: null };
function stepsPerSecond(s) {
  if (s <= 0) return 0;
  if (s <= 1) return Math.pow(10, Math.log10(SPS1) * (s - 0.1) / 0.9); // 0.1× → 1 step/s, 1× → 20 000 steps/s
  return SPS1 * s;
}
function fmtTime(fs) {
  if (fs < 1000) return fs.toFixed(0) + ' fs';
  if (fs < 1e6) return (fs / 1000).toFixed(3) + ' ps';
  return (fs / 1e6).toFixed(4) + ' ns';
}
function fmtRate(sps) {
  const fsps = sps; // fs per real second
  if (fsps < 1000) return fsps.toFixed(fsps < 10 ? 1 : 0) + ' fs/s';
  if (fsps < 1e6) return (fsps / 1000).toFixed(fsps < 1e4 ? 2 : 1) + ' ps/s';
  return (fsps / 1e6).toFixed(2) + ' ns/s';
}

/* ======================= undo ======================= */
const undoStack = [], redoStack = [];
function pushUndo() {
  undoStack.push({ s: eng.snapshot(), box: { ...eng.box } });
  if (undoStack.length > 60) undoStack.shift();
  redoStack.length = 0;
}
function applySnap(o) {
  eng.box = { ...o.box };
  eng.restore(o.s);
  eng.checkpoints.length = 0;
  edited();
}
function undo() { if (!undoStack.length) return toast('Nothing to undo'); redoStack.push({ s: eng.snapshot(), box: { ...eng.box } }); applySnap(undoStack.pop()); }
function redo() { if (!redoStack.length) return; undoStack.push({ s: eng.snapshot(), box: { ...eng.box } }); applySnap(redoStack.pop()); }
let editTick = 0;
function edited() { // topology changed by the user: no reaction/bond events for this change
  eng.touch(); eng.refresh();
  editTick++; bondTrack.reset = true; species.reset = true;
  selection.forEach(i => { if (i >= eng.N) selection.delete(i); });
  updateEmpty(); scheduleSave();
}

/* ======================= toasts ======================= */
function toast(msg, cls) {
  const t = document.createElement('div'); t.className = 'toast ' + (cls || ''); t.textContent = msg;
  $('toasts').appendChild(t);
  setTimeout(() => { t.classList.add('out'); setTimeout(() => t.remove(), 320); }, 2600);
}

/* ======================= camera ======================= */
function fitBox(animate) {
  const b = eng.box, W = R.W, H = R.H;
  const span = clamp(Math.max((b.x1 - b.x0) * 1.12, (b.y1 - b.y0) * 1.18 * W / H), SPAN_MIN, SPAN_MAX);
  camTo((b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2 + 1.2 * span / W * 20, span, animate);
}
let camAnim = null;
function camTo(cx, cy, span, animate) {
  if (!animate) { R.cam.cx = cx; R.cam.cy = cy; R.cam.span = span; camAnim = null; return; }
  camAnim = { from: { ...R.cam }, to: { cx, cy, span }, t0: performance.now(), dur: 420 };
}
function zoomAt(sx, sy, factor) {
  const [wx, wy] = R.toWorld(sx, sy);
  const span = clamp(R.cam.span * factor, SPAN_MIN, SPAN_MAX);
  const k = span / R.cam.span;
  R.cam.cx = wx - (wx - R.cam.cx) * k; R.cam.cy = wy - (wy - R.cam.cy) * k; R.cam.span = span;
  camAnim = null; clampCam();
}
function clampCam() {
  const b = eng.box, m = 6;
  R.cam.cx = clamp(R.cam.cx, b.x0 - m, b.x1 + m); R.cam.cy = clamp(R.cam.cy, b.y0 - m, b.y1 + m);
}
function resize() { R.resize(window.innerWidth, window.innerHeight); if (condR) condR.resize(condCanvasSize()[0], condCanvasSize()[1]); }
window.addEventListener('resize', resize);

/* ======================= scrub fields =======================
 A number you can drag (horizontal), scroll, nudge with arrows, or click to type. The underline
 shows the position inside a soft range; typed values may go beyond it. */
function scrub(el, o) {
  const toP = v => {
    if (o.map === 'log') return (Math.log(Math.max(o.min, v) + o.off) - Math.log(o.min + o.off)) / (Math.log(o.max + o.off) - Math.log(o.min + o.off));
    return (v - o.min) / (o.max - o.min);
  };
  const fromP = p => {
    if (o.map === 'log') return Math.exp(Math.log(o.min + o.off) + p * (Math.log(o.max + o.off) - Math.log(o.min + o.off))) - o.off;
    return o.min + p * (o.max - o.min);
  };
  el.tabIndex = 0;
  const api = {
    render() {
      if (el.classList.contains('editing')) return;
      const v = o.get(), p = toP(v);
      el.innerHTML = esc(o.fmt(v)) + (o.unit ? '<span class="u">' + o.unit + '</span>' : '');
      el.style.setProperty('--p', clamp(p, 0, 1));
      el.classList.toggle('over', p > 1.0001 || p < -0.0001);
    },
    set(v, commit) { v = o.snap ? o.snap(v) : v; if (o.hardMin != null) v = Math.max(o.hardMin, v); if (o.hardMax != null) v = Math.min(o.hardMax, v); o.set(v, commit); api.render(); }
  };
  let drag = null;
  el.addEventListener('pointerdown', e => {
    if (el.classList.contains('editing') || e.button !== 0) return;
    e.preventDefault(); e.stopPropagation();
    drag = { x: e.clientX, p: clamp(toP(o.get()), 0, 1), moved: false, v0: o.get() };
    el.setPointerCapture(e.pointerId);
  });
  el.addEventListener('pointermove', e => {
    if (!drag) return;
    const dx = e.clientX - drag.x;
    if (!drag.moved && Math.abs(dx) < 3) return;
    drag.moved = true; el.classList.add('dragging');
    const p = clamp(drag.p + dx / (e.shiftKey ? 900 : 180), 0, 1);
    api.set(fromP(p));
  });
  el.addEventListener('pointerup', e => {
    if (!drag) return;
    el.classList.remove('dragging');
    if (!drag.moved) startEdit(); else { api.set(o.get(), true); bump(); }
    drag = null;
  });
  el.addEventListener('wheel', e => {
    e.preventDefault(); e.stopPropagation();
    const p = toP(o.get()) - Math.sign(e.deltaY) * (e.shiftKey ? 0.005 : 0.025);
    api.set(fromP(clamp(p, 0, 1)), true);
  }, { passive: false });
  el.addEventListener('keydown', e => {
    if (el.classList.contains('editing')) return;
    if (e.key === 'Enter') { e.preventDefault(); startEdit(); }
    else if (e.key === 'ArrowUp' || e.key === 'ArrowRight' || e.key === 'ArrowDown' || e.key === 'ArrowLeft') {
      e.preventDefault(); e.stopPropagation();
      const d = (e.key === 'ArrowUp' || e.key === 'ArrowRight') ? 1 : -1;
      api.set(fromP(clamp(toP(o.get()) + d * (e.shiftKey ? 0.005 : 0.03), 0, 1)), true);
    }
  });
  function bump() { el.classList.remove('bump'); void el.offsetWidth; el.classList.add('bump'); }
  function startEdit() {
    el.classList.add('editing');
    el.contentEditable = 'true';
    el.textContent = o.edit ? o.edit(o.get()) : String(+o.get().toPrecision(6));
    const r = document.createRange(); r.selectNodeContents(el); const sel = getSelection(); sel.removeAllRanges(); sel.addRange(r);
    el.focus();
    const done = ok => {
      el.removeEventListener('keydown', onKey, true); el.removeEventListener('blur', onBlur);
      const txt = el.textContent.trim();
      el.contentEditable = 'false'; el.classList.remove('editing');
      if (ok) {
        const v = o.parse ? o.parse(txt) : parseFloat(txt.replace(',', '.'));
        if (Number.isFinite(v)) { api.set(v, true); bump(); } else toast('Enter a number, e.g. ' + o.fmt(o.get()) + (o.unit ? ' ' + o.unit : ''));
      }
      api.render();
    };
    const onKey = e => {
      e.stopPropagation();
      if (e.key === 'Enter') { e.preventDefault(); done(true); el.blur(); }
      else if (e.key === 'Escape') { e.preventDefault(); done(false); el.blur(); }
    };
    const onBlur = () => { if (el.classList.contains('editing')) done(true); };
    el.addEventListener('keydown', onKey, true); el.addEventListener('blur', onBlur);
  }
  api.startEdit = startEdit;
  api.render();
  return api;
}

/* ======================= temperature ======================= */
function blackbody(T) { // approximate sRGB of a blackbody (Tanner Helland fit), dim for cold
  if (T <= 0.5) return 'rgb(90,150,255)';
  const t = Math.max(1000, Math.min(40000, T)) / 100;
  let r, g, b;
  if (t <= 66) { r = 255; g = 99.47 * Math.log(t) - 161.12; b = t <= 19 ? 0 : 138.52 * Math.log(t - 10) - 305.04; }
  else { r = 329.7 * Math.pow(t - 60, -0.1332); g = 288.12 * Math.pow(t - 60, -0.0755); b = 255; }
  let f = 1;
  if (T < 1000) { // below incandescence: blend from cold blue to dull red
    const k = T / 1000; r = 90 + (r - 90) * k; g = 150 + (g - 150) * k; b = 255 + (b - 255) * k; f = 0.75 + 0.25 * k;
  }
  return 'rgb(' + [r, g, b].map(c => Math.round(clamp(c, 0, 255) * f)).join(',') + ')';
}
function setT(T, quiet) {
  eng.T = Math.max(0, T); store.set('T', eng.T); tField.render(); paintT();
  if (!quiet) renderTPop();
}
function paintT() {
  const c = blackbody(eng.T), d = $('bbDot');
  d.style.background = c; d.style.color = c;
}
const fmtT = v => v >= 10000 ? (v / 1000).toFixed(1) + 'k' : v >= 100 ? Math.round(v).toString() : v >= 10 ? v.toFixed(1) : v.toFixed(2);
const tField = scrub($('tField'), {
  get: () => eng.T, set: v => setT(v, true), min: 0, max: 6000, off: 40, map: 'log', unit: 'K', hardMin: 0, hardMax: 50000,
  fmt: fmtT, edit: v => String(+v.toFixed(2)), parse: txt => parseTemp(txt)
});
function parseTemp(txt) {
  const m = txt.trim().toLowerCase().replace(',', '.').match(/^(-?[\d.]+)\s*(k|°?c|°?f)?$/);
  if (!m) return NaN;
  const v = parseFloat(m[1]), u = (m[2] || 'k').replace('°', '');
  return u === 'c' ? v + 273.15 : u === 'f' ? (v - 32) * 5 / 9 + 273.15 : v;
}
paintT();
const T_PRESETS = [[0, 'Absolute zero', 'all motion freezes'], [77, 'Liquid nitrogen', ''], [STP_T, 'Room', '25 °C'], [373.15, 'Boiling water', '100 °C'], [1000, 'Red heat', ''], [3000, 'Flame', 'bonds start breaking'], [6000, 'Sun surface', 'most bonds break']];
const tHist = [], pHist = [];

/* ======================= popovers ======================= */
let openPop = null;
function showPop(pop, anchor, place) {
  closePop();
  pop.classList.add('open'); openPop = { pop, anchor };
  const a = anchor.getBoundingClientRect(), w = pop.offsetWidth, h = pop.offsetHeight;
  let x, y;
  if (place === 'below') { x = clamp(a.right - w, 8, innerWidth - w - 8); y = a.bottom + 8; pop.style.setProperty('--oy', '0%'); }
  else if (place === 'right') { x = a.right + 10; y = clamp(a.top + a.height / 2 - h / 2, 8, innerHeight - h - 8); pop.style.setProperty('--ox', '0%'); pop.style.setProperty('--oy', '50%'); }
  else { x = clamp(a.left + a.width / 2 - w / 2, 8, innerWidth - w - 8); y = a.top - h - 10; pop.style.setProperty('--oy', '100%'); }
  pop.style.left = x + 'px'; pop.style.top = y + 'px';
}
function closePop() { if (openPop) { openPop.pop.classList.remove('open'); openPop = null; } }
document.addEventListener('pointerdown', e => {
  if (openPop && !openPop.pop.contains(e.target) && !openPop.anchor.contains(e.target)) closePop();
}, true);

function renderTPop() {
  const pop = $('tPop'); if (!pop.classList.contains('open')) return;
  pop.querySelector('.head b').textContent = eng.temperature().toFixed(1) + ' K';
  pop.querySelectorAll('.row[data-t]').forEach(r => r.classList.toggle('on', Math.abs(+r.dataset.t - eng.T) < 0.01));
  pop.querySelector('.toggle').classList.toggle('on', eng.thermostat);
}
function openTPop() {
  const pop = $('tPop');
  pop.innerHTML = '<div class="head">Measured <b></b></div><canvas class="spark" id="tSpark"></canvas><div class="sep"></div>' +
    T_PRESETS.map((p, i) => '<button class="row" style="--i:' + i + '" data-t="' + p[0] + '"><span class="sw" style="background:' + blackbody(p[0]) + ';box-shadow:0 0 8px ' + blackbody(p[0]) + '"></span><span class="k">' + fmtT(p[0]) + ' K</span><span class="lbl">' + p[1] + '</span><span class="note">' + p[2] + '</span></button>').join('') +
    '<div class="sep"></div><button class="toggle" title="On: a heat bath holds the target temperature. Off: the box is isolated and energy is conserved.">Heat bath (thermostat)<i></i></button>';
  pop.querySelectorAll('.row').forEach(r => r.onclick = () => { setT(+r.dataset.t); toast('Target ' + fmtT(+r.dataset.t) + ' K — ' + r.querySelector('.lbl').textContent); });
  pop.querySelector('.toggle').onclick = () => { eng.thermostat = !eng.thermostat; store.set('thermostat', eng.thermostat); renderTPop(); toast(eng.thermostat ? 'Heat bath on — temperature is held at the target' : 'Heat bath off — isolated box, energy conserved'); };
  showPop(pop, $('gT'), 'below'); renderTPop();
}
$('gT').addEventListener('click', e => { if (e.target.closest('.scrub')) return; openPop && openPop.pop.id === 'tPop' ? closePop() : openTPop(); });
$('bbDot').title = 'Temperature presets';

function openPPop() {
  const pop = $('pPop');
  pop.innerHTML = '<div class="head">Wall pressure <b id="pNow"></b></div><canvas class="spark" id="pSpark"></canvas>' +
    '<div class="kv"><span>Volume</span><span id="pVol"></span></div><div class="kv"><span>Ideal gas N·k·T/V</span><span id="pIdeal"></span></div><div class="kv"><span>Atoms</span><span id="pAtoms"></span></div>' +
    '<div class="sep"></div><div class="kv"><span>Box width</span><span class="scrub" id="bW"></span></div><div class="kv"><span>Box height</span><span class="scrub" id="bH"></span></div><div class="kv"><span>Slab depth</span><span class="scrub" id="bD"></span></div>';
  const mk = (id, get, set, min, max) => scrub($(id), { get, set: (v, c) => { set(v); if (c) { saveBox(); } }, min, max, unit: 'nm', fmt: v => v.toFixed(2), hardMin: min, hardMax: max });
  mk('bW', () => (eng.box.x1 - eng.box.x0) / 10, v => { eng.box.x1 = eng.box.x0 + v * 10; }, 1, 5);
  mk('bH', () => (eng.box.y1 - eng.box.y0) / 10, v => { eng.box.y1 = eng.box.y0 + v * 10; }, 1, 5);
  mk('bD', () => (eng.box.z1 - eng.box.z0) / 10, v => { eng.box.z0 = -v * 5; eng.box.z1 = v * 5; }, 0.4, 3);
  showPop(pop, $('gP'), 'below'); renderPPop();
}
function renderPPop() {
  const pop = $('pPop'); if (!pop.classList.contains('open')) return;
  const b = eng.box, V = (b.x1 - b.x0) * (b.y1 - b.y0) * (b.z1 - b.z0);
  $('pNow').textContent = fmtP(eng.pressureEMA);
  $('pVol').textContent = (V / 1000).toFixed(1) + ' nm³';
  const nMol = species.count || 0;
  $('pIdeal').textContent = fmtP(nMol * KB * eng.temperature() / V * 16605.39);
  $('pAtoms').textContent = eng.N;
}
function saveBox() { const b = eng.box; store.set('box', { w: b.x1 - b.x0, h: b.y1 - b.y0, d: b.z1 - b.z0 }); }
function fmtP(bar) {
  const a = Math.abs(bar);
  if (a >= 1e4) return (bar / 1000).toFixed(1) + ' kbar';
  if (a >= 100) return bar.toFixed(0) + ' bar';
  if (a >= 1) return bar.toFixed(1) + ' bar';
  return bar.toFixed(2) + ' bar';
}
$('gP').addEventListener('click', () => { openPop && openPop.pop.id === 'pPop' ? closePop() : openPPop(); });
function spark(id, data, color, lo) {
  const c = $(id); if (!c) return;
  const w = c.clientWidth, h = c.clientHeight, dpr = Math.min(2, devicePixelRatio || 1);
  if (c.width !== w * dpr) { c.width = w * dpr; c.height = h * dpr; }
  const x = c.getContext('2d'); x.setTransform(dpr, 0, 0, dpr, 0, 0); x.clearRect(0, 0, w, h);
  if (data.length < 2) return;
  let mn = Math.min(...data), mx = Math.max(...data); if (lo != null) mn = Math.min(mn, lo);
  if (mx - mn < 1e-6) { mx += 1; mn -= 1; }
  x.beginPath();
  data.forEach((v, i) => { const px = i / (data.length - 1) * w, py = h - 3 - (v - mn) / (mx - mn) * (h - 6); i ? x.lineTo(px, py) : x.moveTo(px, py); });
  x.strokeStyle = color; x.lineWidth = 1.2; x.stroke();
  x.lineTo(w, h); x.lineTo(0, h); x.closePath();
  const g = x.createLinearGradient(0, 0, 0, h); g.addColorStop(0, color.replace('rgb', 'rgba').replace(')', ',.18)')); g.addColorStop(1, 'rgba(0,0,0,0)');
  x.fillStyle = g; x.fill();
}

/* ======================= speed ======================= */
const SPEEDS = [[0.1, 'Crawl'], [0.25, 'Slow'], [0.5, 'Unhurried'], [1, 'Normal'], [1.5, 'Fast'], [2, 'Very fast'], [4, 'Flat out']];
const speedField = scrub($('speedField'), {
  get: () => time.speed, set: (v) => { time.speed = v; store.set('speed', v); renderSpeedPop(); }, min: 0.1, max: 4, off: 0, map: 'log', unit: '×', hardMin: 0.01, hardMax: 20,
  fmt: v => v < 1 ? v.toFixed(v < 0.1 ? 2 : 2).replace(/0$/, '') : v.toFixed(1), parse: t => parseFloat(t.replace(/[x×]/g, '').replace(',', '.'))
});
function renderSpeedPop() {
  const pop = $('speedPop'); if (!pop.classList.contains('open')) return;
  const p = clamp(Math.log(time.speed / 0.1) / Math.log(40), 0, 1);
  pop.querySelector('.head b').textContent = speedField && (time.speed < 1 ? +time.speed.toFixed(2) : time.speed.toFixed(1)) + '×';
  pop.querySelector('.bar i').style.width = (p * 100) + '%';
  pop.querySelectorAll('.row').forEach(r => r.classList.toggle('on', Math.abs(+r.dataset.s - time.speed) < 1e-6));
}
function openSpeedPop() {
  const pop = $('speedPop');
  pop.innerHTML = '<div class="head"><b></b><span class="bar"><i></i></span></div><div class="sep"></div>' +
    SPEEDS.map((s, i) => '<button class="row" style="--i:' + i + '" data-s="' + s[0] + '"><span class="k">' + s[0] + '×</span><span class="lbl">' + s[1] + '</span><span class="note">' + fmtRate(stepsPerSecond(s[0])) + '</span><svg class="chk" viewBox="0 0 16 16"><path d="M3 8.5l3.2 3L13 4.5"/></svg></button>').join('');
  pop.querySelectorAll('.row').forEach(r => r.onclick = () => { setSpeed(+r.dataset.s); closePop(); });
  showPop(pop, $('speed'), 'above'); renderSpeedPop();
}
function setSpeed(v) { time.speed = v; store.set('speed', v); speedField.render(); renderSpeedPop(); }
$('speedMenuBtn').onclick = () => openPop && openPop.pop.id === 'speedPop' ? closePop() : openSpeedPop();

/* ======================= play / step ======================= */
let playMorph = { t: 0, target: 0, raf: 0 };
const PLAY_A = [[6, 4], [10.5, 6.6], [10.5, 13.4], [6, 16]], PLAY_B = [[10.5, 6.6], [15, 9.2], [15, 10.8], [10.5, 13.4]];
const PAUSE_A = [[5.5, 4.5], [8.6, 4.5], [8.6, 15.5], [5.5, 15.5]], PAUSE_B = [[11.4, 4.5], [14.5, 4.5], [14.5, 15.5], [11.4, 15.5]];
function morphIcon() {
  const k = playMorph.t, e = k < 0.5 ? 2 * k * k : 1 - Math.pow(-2 * k + 2, 2) / 2;
  const lerp = (A, B) => 'M' + A.map((p, i) => (p[0] + (B[i][0] - p[0]) * e).toFixed(2) + ' ' + (p[1] + (B[i][1] - p[1]) * e).toFixed(2)).join(' L') + ' Z';
  $('playA').setAttribute('d', lerp(PLAY_A, PAUSE_A)); $('playB').setAttribute('d', lerp(PLAY_B, PAUSE_B));
}
function animatePlayIcon() {
  cancelAnimationFrame(playMorph.raf);
  const t0 = performance.now(), from = playMorph.t, to = playMorph.target;
  const tick = now => { const k = Math.min(1, (now - t0) / 220); playMorph.t = from + (to - from) * k; morphIcon(); if (k < 1) playMorph.raf = requestAnimationFrame(tick); };
  playMorph.raf = requestAnimationFrame(tick);
}
function setPlaying(p) {
  if (time.playing === p) return;
  time.playing = p; time.acc = 0; time.limited = false;
  const b = $('playBtn');
  b.classList.toggle('running', p); b.setAttribute('aria-label', p ? 'Pause' : 'Play');
  playMorph.target = p ? 1 : 0; animatePlayIcon();
}
function togglePlay() { setPlaying(!time.playing); }
function stepOnce(dir) {
  setPlaying(false);
  if (dir < 0) {
    if (!eng.stepBack()) { flashBtn('backBtn'); return false; }
  } else eng.step();
  time.manual = { t0: performance.now(), dur: 140 };
  return true;
}
function flashBtn(id) { const b = $(id); b.animate([{ transform: 'translateX(0)' }, { transform: 'translateX(-4px)' }, { transform: 'translateX(3px)' }, { transform: 'translateX(0)' }], { duration: 260 }); }
$('playBtn').onclick = togglePlay;
function holdButton(id, dir) {
  const b = $(id); let timer = null, rate = 0, stop = false;
  const run = () => {
    if (stop) return;
    const n = Math.max(1, Math.round(rate));
    for (let k = 0; k < n; k++) if (!stepOnce(dir)) { stop = true; break; }
    rate = Math.min(rate * 1.18 + 0.4, 400);
    timer = setTimeout(run, 45);
  };
  b.addEventListener('pointerdown', e => { if (e.button) return; stop = false; rate = 1; stepOnce(dir); timer = setTimeout(run, 380); });
  const end = () => { stop = true; clearTimeout(timer); };
  b.addEventListener('pointerup', end); b.addEventListener('pointerleave', end); b.addEventListener('pointercancel', end);
}
holdButton('fwdBtn', 1); holdButton('backBtn', -1);

/* ======================= tools & elements ======================= */
const TOOL_ICONS = {
  grab: '<svg viewBox="0 0 20 20"><path d="M7 9V4.8a1.3 1.3 0 012.6 0V9M9.6 8.5V3.6a1.3 1.3 0 012.6 0v5M12.2 8.8V5.2a1.3 1.3 0 012.6 0v6.3c0 3.4-2.2 5.8-5.3 5.8-2.2 0-3.6-1-4.7-2.8L3 11.2a1.3 1.3 0 012.1-1.4L7 12"/></svg>',
  erase: '<svg viewBox="0 0 20 20"><path d="M8.2 16.5h8.3M3.9 12.3l7.4-7.6a1.6 1.6 0 012.3 0l2.2 2.2a1.6 1.6 0 010 2.3l-6.9 7.3H7.6z"/><path d="M7.4 8.8l4.3 4.3"/></svg>',
  heat: '<svg viewBox="0 0 20 20"><path d="M10 17.5c-3 0-5-2.1-5-4.8 0-3.3 3.2-4.5 3.2-8.2 2.9 1.4 4 3.6 3.6 6 1-.4 1.6-1.4 1.8-2.4 1.3 1.4 1.4 3 1.4 4.6 0 2.7-2 4.8-5 4.8z"/></svg>'
};
const TOOLS = [['grab', 'Grab & pan', 'Drag atoms (pulls while running, moves the molecule while paused). Drag empty space to pan. Shift-drag to select.'],
  ['erase', 'Erase', 'Click or drag over atoms to remove them. Alt removes whole molecules.'],
  ['heat', 'Heat brush', 'Drag to heat atoms under the brush; Shift or right-drag cools. Alt+scroll resizes.']];
let tool = 'grab', armed = null; // armed element symbol for placing
const dockEls = store.get('dockEls', ['H', 'C', 'N', 'O', 'F', 'S', 'P', 'Cl']);
function buildDock() {
  $('tools').innerHTML = TOOLS.map(t => '<button class="tool" data-tool="' + t[0] + '" aria-label="' + t[1] + '">' + TOOL_ICONS[t[0]] + '<sup></sup></button>').join('');
  $('elements').innerHTML = dockEls.map(s => '<button class="el" data-el="' + s + '" style="--c:' + BY_SYM[s].color + '" aria-label="' + BY_SYM[s].name + '">' + s + '<i></i><sup></sup></button>').join('');
  $('tools').querySelectorAll('.tool').forEach(b => {
    b.onclick = () => setTool(b.dataset.tool);
    tipOn(b, () => { const t = TOOLS.find(x => x[0] === b.dataset.tool); return '<b>' + t[1] + '</b> <span class="m">' + keyHint('tool.' + t[0]) + '</span><div style="margin-top:4px">' + t[2] + '</div>'; }, 'right');
  });
  $('elements').querySelectorAll('.el').forEach(b => {
    b.onclick = () => arm(armed === b.dataset.el ? null : b.dataset.el);
    tipOn(b, () => elTip(b.dataset.el), 'right');
  });
  refreshDock(); refreshKeyHints();
}
function elTip(sym) {
  const e = BY_SYM[sym];
  return '<b>' + e.name + '</b> <span class="m">' + keyHint('el.' + sym) + '</span>' +
    '<div class="row"><span>Mass</span><span class="m">' + e.mass.toFixed(3) + ' u</span></div>' +
    '<div class="row"><span>van der Waals radius</span><span class="m">' + (e.rvdw / 10).toFixed(3) + ' nm</span></div>' +
    '<div class="row"><span>Valence</span><span class="m">' + e.valences.join(', ') + '</span></div>' +
    '<div class="row"><span>Electronegativity</span><span class="m">' + (e.chi || '—') + '</span></div>' +
    '<div style="margin-top:5px">Click the field to place · drag to throw</div>';
}
function refreshDock() {
  document.querySelectorAll('.tool').forEach(b => b.classList.toggle('on', !armed && !placing && tool === b.dataset.tool));
  document.querySelectorAll('.el').forEach(b => b.classList.toggle('on', armed === b.dataset.el));
  canvas.className = tool === 'erase' && !armed ? 'erase' : (tool === 'grab' && !armed && !placing) ? 'grab' : '';
}
function setTool(t) { tool = t; armed = null; placing = null; brush.active = false; refreshDock(); }
function arm(sym) {
  if (sym && !dockEls.includes(sym)) {
    if (dockEls.length >= 9) dockEls.pop();
    dockEls.push(sym); store.set('dockEls', dockEls); buildDock();
  }
  armed = sym; placing = null; refreshDock();
}
// periodic flyout
const PT54 = 'H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe'.split(' ');
function ptPos(z) {
  if (z === 1) return [1, 1]; if (z === 2) return [1, 18];
  if (z <= 10) return [2, z <= 4 ? z - 2 : z + 8];
  if (z <= 18) return [3, z <= 12 ? z - 10 : z];
  if (z <= 36) return [4, z - 18];
  return [5, z - 36];
}
function openElPop() {
  const pop = $('elPop');
  pop.innerHTML = '<div class="cap">Elements with bonding parameters. Pick one to place it; it joins the dock.</div><div class="grid">' + PT54.map((s, i) => {
    const [r, c] = ptPos(i + 1), e = BY_SYM[s];
    return '<button class="cell ' + (e ? '' : 'off') + (dockEls.includes(s) ? ' pinned' : '') + '" style="grid-row:' + r + ';grid-column:' + c + ';--c:' + (e ? e.color : '#333') + '" data-el="' + s + '">' + s + (e ? '<i></i>' : '') + '</button>';
  }).join('') + '</div>';
  pop.querySelectorAll('.cell:not(.off)').forEach(b => { b.onclick = () => { arm(b.dataset.el); closePop(); }; tipOn(b, () => elTip(b.dataset.el), 'above'); });
  showPop(pop, $('moreEl'), 'right');
}
$('moreEl').onclick = () => openPop && openPop.pop.id === 'elPop' ? closePop() : openElPop();

/* tooltips */
const tip = $('tip'); let tipTimer = 0;
function tipOn(el, html, place) {
  el.addEventListener('pointerenter', () => { clearTimeout(tipTimer); tipTimer = setTimeout(() => showTip(el, html(), place), 380); });
  el.addEventListener('pointerleave', hideTip); el.addEventListener('pointerdown', hideTip);
}
function showTip(el, html, place) {
  tip.innerHTML = html; tip.classList.add('show');
  const a = el.getBoundingClientRect(), w = tip.offsetWidth, h = tip.offsetHeight;
  let x = a.right + 12, y = a.top + a.height / 2 - h / 2;
  if (place === 'above') { x = a.left + a.width / 2 - w / 2; y = a.top - h - 8; }
  tip.style.left = clamp(x, 6, innerWidth - w - 6) + 'px'; tip.style.top = clamp(y, 6, innerHeight - h - 6) + 'px';
}
function showTipAt(x, y, html) {
  tip.innerHTML = html; tip.classList.add('show');
  const w = tip.offsetWidth, h = tip.offsetHeight;
  tip.style.left = clamp(x + 16, 6, innerWidth - w - 6) + 'px'; tip.style.top = clamp(y + 16, 6, innerHeight - h - 6) + 'px';
}
function hideTip() { clearTimeout(tipTimer); tip.classList.remove('show'); }

/* ======================= placing atoms ======================= */
function thermalVel(sym, T) { const m = BY_SYM[sym].mass, s = Math.sqrt(KB * T / (m * 1e4)); return [s * eng.gauss(), s * eng.gauss(), s * eng.gauss()]; }
function tooClose(x, y, z, minD) {
  for (let i = 0; i < eng.N; i++) { const dx = eng.pos[3 * i] - x, dy = eng.pos[3 * i + 1] - y, dz = eng.pos[3 * i + 2] - z; if (dx * dx + dy * dy + dz * dz < minD * minD) return true; }
  return false;
}
function placeAtom(sym, x, y, v) {
  const z = (eng.rand() - 0.5) * 0.6;
  if (tooClose(x, y, z, 1.1)) { toast('Too close to another atom'); return; }
  pushUndo();
  const tv = thermalVel(sym, eng.T);
  eng.addAtom(sym, x, y, z, { v: v ? [v[0] + tv[0] * 0.3, v[1] + tv[1] * 0.3, tv[2] * 0.3] : tv });
  edited();
  bondTrack.pending.push({ x, y, t0: performance.now(), dur: 420, kind: 'form', big: true });
}
// fling: screen drag → velocity (Å/fs). 10 Å of drag ≈ 0.02 Å/fs (2 km/s)
function flingVel(dx, dy) { const k = 0.002, v = [dx * k, dy * k]; const s = Math.hypot(v[0], v[1]), cap = 0.12; if (s > cap) { v[0] *= cap / s; v[1] *= cap / s; } return v; }
function flingLabel(v, mass) { const ke = 0.5 * mass * (v[0] * v[0] + v[1] * v[1]) * 1e4; return (Math.hypot(v[0], v[1]) * 1e5 / 1000).toFixed(1) + ' km/s · ' + ke.toFixed(ke < 10 ? 1 : 0) + ' kJ/mol'; }

/* ======================= pointer on the field ======================= */
const pointers = new Map();
let gesture = null, hoverAtom = -1, lastMouse = { x: 0, y: 0, wx: 0, wy: 0 };
const selection = new Set();
const brush = { r: store.get('brushR', 4), active: false, cool: false, x: 0, y: 0 };
function hitAtom(wx, wy, extra = 0) {
  let best = -1, bd = Infinity;
  const k = R.mode === 'space' ? 1 : R.mode === 'balls' ? 0.45 : 0.5;
  for (let i = 0; i < eng.N; i++) {
    const dx = rp[3 * i] - wx, dy = rp[3 * i + 1] - wy, d = dx * dx + dy * dy;
    const r = Math.max(ELEMENTS[eng.type[i]].rvdw * k, 7 / R.scale) + extra;
    if (d < r * r && d < bd) { bd = d; best = i; }
  }
  return best;
}
function boxEdgeAt(sx, sy) {
  const b = eng.box, [x0, y0] = R.toScreen(b.x0, b.y0), [x1, y1] = R.toScreen(b.x1, b.y1), t = 7;
  const nearX = Math.abs(sx - x1) < t && sy > y0 - t && sy < y1 + t, nearY = Math.abs(sy - y1) < t && sx > x0 - t && sx < x1 + t;
  return nearX && nearY ? 'xy' : nearX ? 'x' : nearY ? 'y' : null;
}
function fragmentOf(i) { const fr = eng.fragments(0.5); return fr.list[fr.comp[i]]; }

canvas.addEventListener('contextmenu', e => e.preventDefault());
canvas.addEventListener('pointerdown', e => {
  canvas.setPointerCapture(e.pointerId);
  pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
  closePop(); hideTip();
  if (pointers.size === 2) { // pinch
    const [a, b] = [...pointers.values()];
    gesture = { type: 'pinch', d: Math.hypot(a.x - b.x, a.y - b.y), mx: (a.x + b.x) / 2, my: (a.y + b.y) / 2 };
    return;
  }
  const [wx, wy] = R.toWorld(e.clientX, e.clientY);
  const hit = hitAtom(wx, wy);
  if (e.button === 1 || (e.button === 0 && spaceHeld)) { gesture = { type: 'pan', sx: e.clientX, sy: e.clientY, cx: R.cam.cx, cy: R.cam.cy }; canvas.className = 'grabbing'; return; }
  if (placing) { gesture = { type: 'placeMol', sx: e.clientX, sy: e.clientY, wx, wy }; return; }
  if (armed && e.button === 0) { gesture = { type: 'placeAtom', sx: e.clientX, sy: e.clientY, wx, wy }; return; }
  if (e.button === 2) {
    if (tool === 'heat') { brush.active = true; brush.cool = true; gesture = { type: 'brush' }; return; }
    if (hit >= 0) { pushUndo(); eng.removeAtoms(e.altKey ? fragmentOf(hit) : [hit]); edited(); return; }
    gesture = { type: 'pan', sx: e.clientX, sy: e.clientY, cx: R.cam.cx, cy: R.cam.cy }; canvas.className = 'grabbing'; return;
  }
  if (tool === 'erase') { pushUndo(); gesture = { type: 'erase', alt: e.altKey }; eraseAt(wx, wy, e.altKey); return; }
  if (tool === 'heat') { brush.active = true; brush.cool = e.shiftKey; gesture = { type: 'brush' }; return; }
  const edge = boxEdgeAt(e.clientX, e.clientY);
  if (edge && hit < 0) { pushUndo(); gesture = { type: 'box', edge }; return; }
  if (e.shiftKey) { gesture = { type: 'marquee', x0: wx, y0: wy, x1: wx, y1: wy }; return; }
  if (hit >= 0) {
    if (time.playing) { gesture = { type: 'tweezer', i: hit }; eng.tweezer = { i: hit, x: wx, y: wy, k: 30 }; eng.checkpoints.length = 0; }
    else {
      pushUndo();
      const set = e.altKey ? [hit] : (selection.has(hit) ? [...selection] : fragmentOf(hit));
      gesture = { type: 'move', set, wx, wy, moved: false };
    }
    canvas.className = 'grabbing';
    return;
  }
  selection.clear();
  gesture = { type: 'pan', sx: e.clientX, sy: e.clientY, cx: R.cam.cx, cy: R.cam.cy };
  canvas.className = 'grabbing';
});
canvas.addEventListener('pointermove', e => {
  if (pointers.has(e.pointerId)) pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
  const [wx, wy] = R.toWorld(e.clientX, e.clientY);
  lastMouse = { x: e.clientX, y: e.clientY, wx, wy };
  brush.x = wx; brush.y = wy;
  if (placing) { placing.x = wx; placing.y = wy; }
  if (!gesture) {
    hoverAtom = (armed || placing) ? -1 : hitAtom(wx, wy);
    const edge = !armed && !placing && tool === 'grab' ? boxEdgeAt(e.clientX, e.clientY) : null;
    boxHot = edge && hoverAtom < 0 ? edge : null;
    if (boxHot) canvas.className = 'resize-' + boxHot; else refreshDock();
    if (hoverAtom >= 0 && tool === 'grab') showTipAt(e.clientX, e.clientY, atomTip(hoverAtom)); else hideTip();
    return;
  }
  hideTip();
  const g = gesture;
  if (g.type === 'pinch' && pointers.size === 2) {
    const [a, b] = [...pointers.values()], d = Math.hypot(a.x - b.x, a.y - b.y), mx = (a.x + b.x) / 2, my = (a.y + b.y) / 2;
    zoomAt(mx, my, g.d / d);
    R.cam.cx -= (mx - g.mx) / R.scale; R.cam.cy -= (my - g.my) / R.scale; clampCam();
    g.d = d; g.mx = mx; g.my = my; return;
  }
  if (g.type === 'pan') { R.cam.cx = g.cx - (e.clientX - g.sx) / R.scale; R.cam.cy = g.cy - (e.clientY - g.sy) / R.scale; camAnim = null; clampCam(); }
  else if (g.type === 'tweezer') { if (eng.tweezer) { eng.tweezer.x = wx; eng.tweezer.y = wy; } }
  else if (g.type === 'move') {
    const dx = wx - g.wx, dy = wy - g.wy; g.wx = wx; g.wy = wy; g.moved = true;
    for (const i of g.set) { eng.pos[3 * i] += dx; eng.pos[3 * i + 1] += dy; eng.prev[3 * i] = eng.pos[3 * i]; eng.prev[3 * i + 1] = eng.pos[3 * i + 1]; }
    eng.needRebuild = true; eng.needForces = true;
  } else if (g.type === 'erase') eraseAt(wx, wy, g.alt);
  else if (g.type === 'marquee') { g.x1 = wx; g.y1 = wy; }
  else if (g.type === 'box') {
    const b = eng.box;
    if (g.edge.includes('x')) b.x1 = clamp(wx, b.x0 + 10, b.x0 + 50);
    if (g.edge.includes('y')) b.y1 = clamp(wy, b.y0 + 10, b.y0 + 50);
  } else if (g.type === 'placeAtom' || g.type === 'placeMol') { g.cx = wx; g.cy = wy; }
});
function endPointer(e) {
  pointers.delete(e.pointerId);
  const g = gesture; if (!g) return;
  if (g.type === 'pinch') { if (pointers.size < 2) gesture = null; return; }
  gesture = null;
  const [wx, wy] = R.toWorld(e.clientX, e.clientY);
  if (g.type === 'tweezer') eng.tweezer = null;
  else if (g.type === 'move') { if (g.moved) edited(); else { undoStack.pop(); if (!e.shiftKey) { selection.clear(); } } }
  else if (g.type === 'erase') edited();
  else if (g.type === 'brush') brush.active = false;
  else if (g.type === 'box') { saveBox(); edited(); }
  else if (g.type === 'marquee') {
    const x0 = Math.min(g.x0, g.x1), x1 = Math.max(g.x0, g.x1), y0 = Math.min(g.y0, g.y1), y1 = Math.max(g.y0, g.y1);
    if (!e.ctrlKey && !e.metaKey) selection.clear();
    for (let i = 0; i < eng.N; i++) if (rp[3 * i] >= x0 && rp[3 * i] <= x1 && rp[3 * i + 1] >= y0 && rp[3 * i + 1] <= y1) selection.add(i);
    if (selection.size) toast(selection.size + ' atoms selected — Delete removes, drag moves (paused)');
  } else if (g.type === 'placeAtom') {
    const dragPx = Math.hypot(e.clientX - g.sx, e.clientY - g.sy);
    placeAtom(armed, g.wx, g.wy, dragPx > 8 ? flingVel(wx - g.wx, wy - g.wy) : null);
  } else if (g.type === 'placeMol') {
    const dragPx = Math.hypot(e.clientX - g.sx, e.clientY - g.sy);
    dropMolecule(g.wx, g.wy, dragPx > 8 ? flingVel(wx - g.wx, wy - g.wy) : null, e.shiftKey);
  }
  refreshDock();
}
canvas.addEventListener('pointerup', endPointer);
canvas.addEventListener('pointercancel', endPointer);
canvas.addEventListener('pointerleave', () => { hoverAtom = -1; hideTip(); });
canvas.addEventListener('wheel', e => {
  e.preventDefault();
  if (placing && (e.shiftKey || e.altKey)) { placing.rot += Math.sign(e.deltaY) * Math.PI / 12; return; }
  if (tool === 'heat' && e.altKey) { brush.r = clamp(brush.r * Math.exp(-e.deltaY * 0.002), 1, 20); store.set('brushR', brush.r); return; }
  const d = e.deltaMode === 1 ? e.deltaY * 16 : e.deltaY;
  zoomAt(e.clientX, e.clientY, Math.exp(d * (e.ctrlKey ? 0.01 : 0.0015)));
}, { passive: false });
function eraseAt(wx, wy, alt) {
  const i = hitAtom(wx, wy, 0.3); if (i < 0) return;
  const list = alt ? fragmentOf(i) : [i];
  bondTrack.pending.push({ x: rp[3 * i], y: rp[3 * i + 1], t0: performance.now(), dur: 380, kind: 'break', big: true });
  eng.removeAtoms(list); rp = new Float64Array(eng.pos.subarray(0, 3 * eng.N)); eng.touch(); eng.refresh(); bondTrack.reset = true; species.reset = true;
}
function atomTip(i) {
  const e = ELEMENTS[eng.type[i]], frag = fragmentOf(i), f = eng.formulaOf(frag) + (eng.isRadical(frag) ? '·' : '');
  const v2 = eng.vel[3 * i] ** 2 + eng.vel[3 * i + 1] ** 2 + eng.vel[3 * i + 2] ** 2, Tloc = e.mass * v2 * 1e4 / (3 * KB);
  const q = eng.q[i];
  return '<b>' + e.name + '</b> in <span class="m">' + esc(pretty(f)) + '</span>' +
    '<div class="row"><span>Partial charge</span><span class="m">' + (q >= 0 ? '+' : '−') + Math.abs(q).toFixed(2) + ' e</span></div>' +
    '<div class="row"><span>Bonded neighbours</span><span class="m">' + eng.Z[i].toFixed(1) + ' / ' + eng.val[i] + '</span></div>' +
    '<div class="row"><span>Speed</span><span class="m">' + (Math.sqrt(v2) * 100).toFixed(2) + ' km/s</span></div>' +
    '<div class="row"><span>Kinetic temperature</span><span class="m">' + Tloc.toFixed(0) + ' K</span></div>' +
    '<div class="row"><span>Depth z</span><span class="m">' + (eng.pos[3 * i + 2] / 10).toFixed(3) + ' nm</span></div>';
}
let spaceHeld = false, boxHot = null;

/* heat brush: rescale velocities of atoms under the brush a little every frame */
function applyBrush(dt) {
  if (!brush.active) return;
  const f = brush.cool ? Math.exp(-dt * 4) : Math.exp(dt * 2.2), r2 = brush.r * brush.r, list = [];
  for (let i = 0; i < eng.N; i++) { const dx = eng.pos[3 * i] - brush.x, dy = eng.pos[3 * i + 1] - brush.y; if (dx * dx + dy * dy < r2) list.push(i); }
  if (!list.length) return;
  // give frozen atoms a seed velocity so heating works from absolute zero
  for (const i of list) { const v2 = eng.vel[3 * i] ** 2 + eng.vel[3 * i + 1] ** 2 + eng.vel[3 * i + 2] ** 2; if (!brush.cool && v2 < 1e-8) { const tv = thermalVel(ELEMENTS[eng.type[i]].sym, 30); eng.vel.set(tv, 3 * i); } }
  eng.scaleVelocities(list, f); eng.checkpoints.length = 0;
}

/* ======================= species, reactions, bond events ======================= */
const species = { map: new Map(), reset: true, last: 0, count: 0 };
const bondTrack = { set: new Map(), reset: true, pending: [] };
const flashes = [];
const feedItems = [];
function speciesNow() {
  const fr = eng.fragments(0.5), map = new Map(), frags = [], byId = new Map();
  fr.list.forEach((g, k) => {
    const f = eng.formulaOf(g) + (eng.isRadical(g) ? '·' : '');
    map.set(f, (map.get(f) || 0) + 1);
    let x = 0, y = 0; for (const i of g) { x += eng.pos[3 * i]; y += eng.pos[3 * i + 1]; byId.set(eng.ids[i], k); }
    frags.push({ f, x: x / g.length, y: y / g.length, n: g.length });
  });
  return { map, frags, byId, n: fr.list.length };
}
/* A reaction = a connected cluster of old and new fragments that share atoms. */
function reactionsBetween(prev, cur) {
  const parent = new Map(), find = a => { while (parent.get(a) !== a) { parent.set(a, parent.get(parent.get(a))); a = parent.get(a); } return a; };
  const node = k => { if (!parent.has(k)) parent.set(k, k); return k; };
  for (const [id, o] of prev.byId) { const n = cur.byId.get(id); if (n === undefined) continue; const a = find(node('o' + o)), b = find(node('n' + n)); if (a !== b) parent.set(a, b); }
  const groups = new Map();
  for (const k of parent.keys()) { const r = find(k); if (!groups.has(r)) groups.set(r, { o: [], n: [] }); groups.get(r)[k[0]].push(+k.slice(1)); }
  const out = [];
  for (const g of groups.values()) {
    if (g.o.length === 1 && g.n.length === 1 && prev.frags[g.o[0]].f === cur.frags[g.n[0]].f) continue;
    const side = (list, frs) => { const c = new Map(); list.forEach(k => c.set(frs[k].f, (c.get(frs[k].f) || 0) + 1)); return [...c.entries()].sort((a, b) => a[0].localeCompare(b[0])).map(([f, n]) => (n > 1 ? n + ' ' : '') + pretty(f)).join(' + '); };
    const lhs = side(g.o, prev.frags), rhs = side(g.n, cur.frags);
    if (lhs === rhs) continue;
    const big = g.n.reduce((a, k) => cur.frags[k].n > cur.frags[a].n ? k : a, g.n[0]);
    out.push({ lhs, rhs, loc: [cur.frags[big].x, cur.frags[big].y] });
  }
  return out;
}
function updateSpecies(now) {
  if (now - species.last < 200) return;
  species.last = now;
  const cur = speciesNow();
  species.count = cur.n;
  if (!species.reset && species.prev) {
    const rx = reactionsBetween(species.prev, cur);
    rx.slice(0, 3).forEach(r => logReaction(r.lhs, r.rhs, r.loc));
    if (rx.length > 3) logReaction('', '+' + (rx.length - 3) + ' more at once', null);
  }
  species.reset = false; species.map = cur.map; species.prev = cur;
  renderInventory(cur.map);
}
function renderInventory(map) {
  const items = [...map.entries()].sort((a, b) => b[1] - a[1] || a[0].localeCompare(b[0])).slice(0, 10);
  const html = items.map(([k, n]) => '<span' + (k.endsWith('·') ? ' class="rad" title="radical: has an unpaired valence"' : '') + '><b>' + n + '</b>' + esc(pretty(k)) + '</span>').join('');
  const inv = $('inventory');
  if (inv.dataset.h !== html) { inv.dataset.h = html; inv.innerHTML = html; }
}
function logReaction(lhs, rhs, loc) {
  const li = document.createElement('li');
  li.innerHTML = '<time>' + fmtTime(eng.time) + '</time><span class="eq">' + (lhs ? esc(lhs) + ' <span class="arrow">→</span> ' : '') + esc(rhs) + '</span>';
  li.title = 'Show where it happened';
  li.onclick = () => { if (loc) camTo(loc[0], loc[1], Math.min(R.cam.span, 16), true); };
  $('feed').appendChild(li);
  feedItems.push({ li, t: performance.now() });
  while (feedItems.length > 5) { const o = feedItems.shift(); o.li.remove(); }
  if (loc) flashes.push({ x: loc[0], y: loc[1], t0: performance.now(), dur: 900, kind: 'form', big: true });
}
function ageFeed(now) {
  for (const f of feedItems) { const a = now - f.t; f.li.classList.toggle('old', a > 9000); f.li.classList.toggle('gone', a > 25000); }
}
function updateBondEvents(bonds, now) {
  const cur = new Map();
  for (const b of bonds) {
    const key = eng._key(b.i, b.j), prev = bondTrack.set.get(key);
    if (b.strength > 0.6 || (prev && b.strength > 0.35)) cur.set(key, [b.i, b.j]);
  }
  if (!bondTrack.reset) {
    let n = 0;
    for (const [k, ij] of cur) if (!bondTrack.set.has(k) && n++ < 12) flashes.push({ x: (rp[3 * ij[0]] + rp[3 * ij[1]]) / 2, y: (rp[3 * ij[0] + 1] + rp[3 * ij[1] + 1]) / 2, t0: now, dur: 650, kind: 'form' });
    for (const [k, ij] of bondTrack.set) if (!cur.has(k) && n++ < 12 && ij[0] < eng.N && ij[1] < eng.N) flashes.push({ x: (rp[3 * ij[0]] + rp[3 * ij[1]]) / 2, y: (rp[3 * ij[0] + 1] + rp[3 * ij[1] + 1]) / 2, t0: now, dur: 650, kind: 'break' });
  }
  bondTrack.reset = false; bondTrack.set = cur;
  while (bondTrack.pending.length) flashes.push(bondTrack.pending.shift());
  for (let k = flashes.length - 1; k >= 0; k--) if (now - flashes[k].t0 > flashes[k].dur) flashes.splice(k, 1);
}

/* ======================= molecules: inbox ======================= */
const INBOX_KEY = 'chemPlayground.inbox';
function readInbox() { try { return JSON.parse(localStorage.getItem(INBOX_KEY) || '[]').filter(validMol); } catch (e) { return []; } }
function writeInbox(list) { try { localStorage.setItem(INBOX_KEY, JSON.stringify(list.slice(0, 40))); } catch (e) { } }
function validMol(m) { return m && m.format === 'chem-playground/molecule@1' && Array.isArray(m.atoms) && m.atoms.length && m.atoms.every(a => BY_SYM[a.el]); }
let inbox = readInbox();
function receive(mol, open) {
  if (!validMol(mol)) {
    const bad = mol && Array.isArray(mol.atoms) ? mol.atoms.find(a => !BY_SYM[a.el]) : null;
    toast(bad ? bad.el + ' has no bonding parameters in the Playground yet' : 'That is not a Nomenclature molecule');
    return;
  }
  if (!inbox.some(m => m.id === mol.id)) { inbox.unshift(mol); writeInbox(inbox); }
  renderTray(mol.id);
  const pill = $('trayPill'); pill.classList.remove('ping'); void pill.offsetWidth; pill.classList.add('ping');
  if (open) startConditioning(mol);
  else { $('tray').classList.add('open'); toast((mol.name || pretty(mol.formula)) + ' arrived from Nomenclature', 'ok'); }
}
function thumbSVG(m) {
  const A = m.atoms; let x0 = Infinity, y0 = Infinity, x1 = -Infinity, y1 = -Infinity;
  A.forEach(a => { x0 = Math.min(x0, a.x); y0 = Math.min(y0, a.y); x1 = Math.max(x1, a.x); y1 = Math.max(y1, a.y); });
  const pad = 1.2, w = Math.max(3, x1 - x0 + 2 * pad), h = Math.max(2, y1 - y0 + 2 * pad);
  let s = '<svg viewBox="' + (x0 - pad) + ' ' + (y0 - pad) + ' ' + w + ' ' + h + '" preserveAspectRatio="xMidYMid meet">';
  for (const b of m.bonds || []) {
    const p = A[b.a], q = A[b.b]; if (!p || !q) continue;
    const dx = q.x - p.x, dy = q.y - p.y, L = Math.hypot(dx, dy) || 1, nx = -dy / L * 0.13, ny = dx / L * 0.13;
    const offs = b.order >= 3 ? [-1, 0, 1] : b.order === 2 ? [-0.5, 0.5] : [0];
    offs.forEach(o => { s += '<line x1="' + (p.x + nx * o * 2) + '" y1="' + (p.y + ny * o * 2) + '" x2="' + (q.x + nx * o * 2) + '" y2="' + (q.y + ny * o * 2) + '" stroke="rgba(233,228,216,.55)" stroke-width="0.09"/>'; });
  }
  for (const a of A) { const e = BY_SYM[a.el]; s += '<circle cx="' + a.x + '" cy="' + a.y + '" r="' + (a.el === 'H' ? 0.18 : 0.32) + '" fill="' + e.color + '"/>'; }
  return s + '</svg>';
}
function renderTray(freshId) {
  $('trayCount').textContent = inbox.length;
  const list = $('trayList');
  if (!inbox.length) {
    list.innerHTML = '<div class="tray-empty">Draw a molecule in <b>Nomenclature</b> and press <b>Playground</b>. It lands here with every hydrogen filled in, ready to be conditioned and placed.</div>';
    return;
  }
  list.innerHTML = inbox.map(m => '<button class="mol' + (m.id === freshId ? ' fresh' : '') + '" data-id="' + m.id + '">' + thumbSVG(m) +
    '<span class="nm">' + esc(m.name || pretty(m.formula)) + '</span><span class="fm">' + esc(pretty(m.formula)) + ' · ' + m.atoms.length + ' atoms</span><span class="x" data-x="1" title="Remove">×</span></button>').join('');
  list.querySelectorAll('.mol').forEach(b => b.onclick = e => {
    const m = inbox.find(x => x.id === b.dataset.id);
    if (e.target.dataset.x) { inbox = inbox.filter(x => x !== m); writeInbox(inbox); renderTray(); return; }
    $('tray').classList.remove('open'); startConditioning(m);
  });
}
$('trayPill').onclick = () => $('tray').classList.toggle('open');
function toggleTray() { $('tray').classList.toggle('open'); }
// same-origin delivery from Nomenclature
if ('BroadcastChannel' in window) {
  const ch = new BroadcastChannel('chem-playground');
  ch.onmessage = e => { const d = e.data; if (d && d.type === 'molecule' && d.mol) { ch.postMessage({ type: 'ack', id: d.mol.id }); receive(d.mol, false); } };
}
window.addEventListener('storage', e => { if (e.key === INBOX_KEY) { inbox = readInbox(); renderTray(); } });
function readHash() {
  const m = location.hash.match(/#m=([A-Za-z0-9_-]+)/); if (!m) return;
  try {
    const json = decodeURIComponent(escape(atob(m[1].replace(/-/g, '+').replace(/_/g, '/'))));
    receive(JSON.parse(json), true);
  } catch (e) { toast('Could not read the molecule link'); }
  history.replaceState(null, '', location.pathname + location.search);
}
window.addEventListener('hashchange', readHash);
document.addEventListener('paste', e => {
  const t = e.clipboardData && e.clipboardData.getData('text'); if (!t || t[0] !== '{') return;
  try { const m = JSON.parse(t); if (m.format && m.format.startsWith('chem-playground/molecule')) { receive(m, true); e.preventDefault(); } } catch (err) { }
});

/* ======================= conditioning ======================= */
let cond = null, condR = null, placing = null;
const condPrefs = store.get('cond', { preset: 'normal', T: 500, P: 1 });
function condCanvasSize() { const c = $('sheetCanvas'); return [c.clientWidth || 300, c.clientHeight || 220]; }
function condParams() {
  if (condPrefs.preset === 'zero') return { T: 0, P: 0 };
  if (condPrefs.preset === 'custom') return { T: condPrefs.T, P: condPrefs.P };
  return { T: STP_T, P: ATM };
}
function startConditioning(mol) {
  placing = null; armed = null; refreshDock();
  $('sheet').hidden = false; $('sheet').classList.remove('leaving');
  if (!condR) { condR = new FieldRenderer($('sheetCanvas'), { showBox: false, cell: 4 }); condR.mode = 'density'; }
  const [w, h] = condCanvasSize(); condR.resize(w, h);
  $('sheetName').textContent = mol.name || pretty(mol.formula);
  $('sheetFormula').textContent = pretty(mol.formula) + ' · ' + mol.atoms.length + ' atoms' + (mol.smiles ? ' · ' + mol.smiles : '');
  cond = { mol };
  paintSeg(); buildCondEngine();
}
function buildCondEngine() {
  const mol = cond.mol, { T, P } = condParams();
  const e = new Engine({ T, thermostat: true, tau: 60, seed: 7 });
  e.recording = false;
  e.box = { x0: -1e3, x1: 1e3, y0: -1e3, y1: 1e3, z0: -1e3, z1: 1e3 };
  const sum = new Array(mol.atoms.length).fill(0);
  for (const b of mol.bonds || []) { sum[b.a] += b.order; sum[b.b] += b.order; }
  mol.atoms.forEach((a, i) => {
    const z = (a.z || 0) + (Math.sin(i * 12.9898) * 43758.5453 % 1) * 0.35; // deterministic tiny out-of-plane nudge
    e.addAtom(a.el, a.x, a.y, z, { thermal: false, charge: a.charge || 0, V: Engine.valenceForBonds(a.el, a.charge || 0, sum[i]) });
  });
  e.refresh();
  for (const b of mol.bonds || []) if (b.order > 1) e.setBondOrder(b.a, b.b, b.order);
  let rmax = 0; for (let i = 0; i < e.N; i++) rmax = Math.max(rmax, Math.hypot(e.pos[3 * i], e.pos[3 * i + 1], e.pos[3 * i + 2]) + ELEMENTS[e.type[i]].rvdw);
  let Rs = rmax + 6;
  if (T > 0 && P > 0) {
    const V = KB * T / (P / 16605.39); // volume per molecule (Å³) at this T and P
    Rs = clamp(Math.cbrt(3 * V / (4 * Math.PI)), rmax * 0.55, 80);
  }
  e.sphere = { x: 0, y: 0, z: 0, R: Rs };
  Object.assign(cond, { e, T, P, phase: 'relax', relaxIt: 0, steps: 0, target: 2000, Rs, rmax });
  condR.cam = { cx: 0, cy: 0, span: Math.max(8, rmax * 2.6) };
  $('condPlace').disabled = true;
}
function tickConditioning(budgetMs) {
  if (!cond) return;
  const c = cond, e = c.e, t0 = performance.now();
  if (c.phase === 'relax') {
    const it = e.minimize(40, 0.8);
    c.relaxIt += 40;
    if (it < 40 || c.relaxIt > 2400) {
      if (c.T > 0) { e.thermalize(c.T); e.zeroMomentum(); e.needForces = true; c.phase = 'equil'; }
      else c.phase = 'done';
      e.prev.set(e.pos.subarray(0, 3 * e.N));
    }
  } else {
    if (c.T === 0) { e.minimize(10, 0.2); }
    else {
      while (performance.now() - t0 < budgetMs) { for (let k = 0; k < 20; k++) e.step(); c.steps += 20; }
      if (c.phase === 'equil' && c.steps >= c.target) c.phase = 'done';
    }
  }
  // status
  const fr = e.fragments(0.5), parts = fr.list.length;
  const st = $('condStatus'), bar = $('condBar');
  st.classList.toggle('warn', parts > 1 && c.phase !== 'relax');
  let text;
  if (c.phase === 'relax') { text = 'Relaxing geometry in 3D · E ' + e.Epot.toFixed(0) + ' kJ/mol'; bar.style.width = Math.min(100, c.relaxIt / 24) + '%'; }
  else if (c.phase === 'equil') { text = 'Equilibrating at ' + fmtT(c.T) + ' K, ' + fmtP(c.P) + ' · ' + (c.steps / 1000).toFixed(1) + ' / ' + (c.target / 1000).toFixed(1) + ' ps'; bar.style.width = (c.steps / c.target * 100) + '%'; }
  else {
    bar.style.width = '100%';
    if (parts > 1) text = 'Falls apart at ' + fmtT(c.T) + ' K → ' + fr.list.map(g => pretty(e.formulaOf(g) + (e.isRadical(g) ? '·' : ''))).slice(0, 5).join(' + ') + ' · placeable anyway';
    else text = c.T === 0 ? 'Ready · relaxed to its energy minimum, motionless · E ' + e.Epot.toFixed(0) + ' kJ/mol' : 'Ready · intact at ' + e.temperature().toFixed(0) + ' K · E ' + e.Epot.toFixed(0) + ' kJ/mol';
  }
  $('condText').textContent = text;
  $('condPlace').disabled = c.phase === 'relax';
  // camera follows the molecule
  let cx = 0, cy = 0; for (let i = 0; i < e.N; i++) { cx += e.pos[3 * i]; cy += e.pos[3 * i + 1]; } cx /= e.N; cy /= e.N;
  condR.cam.cx += (cx - condR.cam.cx) * 0.2; condR.cam.cy += (cy - condR.cam.cy) * 0.2;
  const n3 = 3 * e.N, bonds = e.bonds(0.3);
  condR.draw({ N: e.N, type: e.type, pos: e.pos.subarray(0, n3), bonds, now: performance.now() });
  // container outline
  const ctx = condR.ctx, [sx, sy] = condR.toScreen(e.sphere.x, e.sphere.y), rr = e.sphere.R * condR.scale;
  if (rr < 400) { ctx.strokeStyle = 'rgba(127,167,201,.25)'; ctx.setLineDash([2, 5]); ctx.beginPath(); ctx.arc(sx, sy, rr, 0, Math.PI * 2); ctx.stroke(); ctx.setLineDash([]); }
}
function paintSeg() {
  $('condSeg').querySelectorAll('button').forEach(b => b.setAttribute('aria-checked', String(b.dataset.v === condPrefs.preset)));
  $('customRow').classList.toggle('show', condPrefs.preset === 'custom');
  $('customSmall').textContent = fmtT(condPrefs.T) + ' K · ' + fmtP(condPrefs.P);
  cTField.render(); cPField.render();
}
$('condSeg').querySelectorAll('button').forEach(b => b.onclick = () => { condPrefs.preset = b.dataset.v; store.set('cond', condPrefs); paintSeg(); if (cond) buildCondEngine(); });
const cTField = scrub($('cT'), { get: () => condPrefs.T, set: (v, c) => { condPrefs.T = v; store.set('cond', condPrefs); $('customSmall').textContent = fmtT(v) + ' K · ' + fmtP(condPrefs.P); if (c && cond) buildCondEngine(); }, min: 0, max: 6000, off: 40, map: 'log', unit: 'K', hardMin: 0, hardMax: 50000, fmt: fmtT, parse: parseTemp });
const cPField = scrub($('cP'), { get: () => condPrefs.P, set: (v, c) => { condPrefs.P = v; store.set('cond', condPrefs); $('customSmall').textContent = fmtT(condPrefs.T) + ' K · ' + fmtP(v); if (c && cond) buildCondEngine(); }, min: 0.01, max: 1e5, off: 0, map: 'log', unit: 'bar', hardMin: 0, hardMax: 1e7, fmt: v => v >= 1000 ? (v / 1000).toFixed(v >= 1e4 ? 0 : 1) + 'k' : v >= 10 ? v.toFixed(0) : v.toFixed(2) });
function closeSheet() {
  if (!cond && $('sheet').hidden) return;
  cond = null; const s = $('sheet'); s.classList.add('leaving');
  setTimeout(() => { if (!cond) s.hidden = true; s.classList.remove('leaving'); }, 240);
}
$('condCancel').onclick = closeSheet;
$('condPlace').onclick = beginPlacing;
function beginPlacing() {
  if (!cond || cond.phase === 'relax') return;
  const e = cond.e;
  let cx = 0, cy = 0, cz = 0, M = 0, vx = 0, vy = 0, vz = 0;
  for (let i = 0; i < e.N; i++) { const m = e.mass[i]; cx += m * e.pos[3 * i]; cy += m * e.pos[3 * i + 1]; cz += m * e.pos[3 * i + 2]; vx += m * e.vel[3 * i]; vy += m * e.vel[3 * i + 1]; vz += m * e.vel[3 * i + 2]; M += m; }
  cx /= M; cy /= M; cz /= M; vx /= M; vy /= M; vz /= M;
  const atoms = [];
  for (let i = 0; i < e.N; i++) atoms.push({ t: e.type[i], sym: ELEMENTS[e.type[i]].sym, x: e.pos[3 * i] - cx, y: e.pos[3 * i + 1] - cy, z: e.pos[3 * i + 2] - cz, v: [e.vel[3 * i] - vx, e.vel[3 * i + 1] - vy, e.vel[3 * i + 2] - vz], q: e.formal[i], V: e.val[i] });
  const bonds = e.bonds(0.5).map(b => ({ a: b.i, b: b.j, order: b.order }));
  placing = { name: cond.mol.name || pretty(cond.mol.formula), atoms, bonds, mass: M, rot: 0, x: lastMouse.wx, y: lastMouse.wy, T: cond.T };
  closeSheet(); refreshDock();
  toast('Click to place · drag to throw · Q/E rotate · Shift-click places several · Esc stops');
}
function ghostAtoms(x, y) {
  const c = Math.cos(placing.rot), s = Math.sin(placing.rot);
  return placing.atoms.map(a => ({ t: a.t, x: x + a.x * c - a.y * s, y: y + a.x * s + a.y * c, z: a.z }));
}
function ghostOK(g) {
  const b = eng.box;
  for (const a of g) {
    if (a.x < b.x0 || a.x > b.x1 || a.y < b.y0 || a.y > b.y1) return false;
    if (tooClose(a.x, a.y, a.z, 1.5)) return false;
  }
  return true;
}
function dropMolecule(x, y, fling, keep) {
  const g = ghostAtoms(x, y);
  if (!ghostOK(g)) { toast('No room here — move the ghost to a free spot inside the box'); return; }
  pushUndo();
  const c = Math.cos(placing.rot), s = Math.sin(placing.rot), base = eng.N;
  placing.atoms.forEach((a, k) => {
    const v = [a.v[0] * c - a.v[1] * s + (fling ? fling[0] : 0), a.v[0] * s + a.v[1] * c + (fling ? fling[1] : 0), a.v[2]];
    eng.addAtom(a.sym, g[k].x, g[k].y, g[k].z, { v, charge: a.q, V: a.V });
  });
  eng.refresh();
  for (const b of placing.bonds) if (b.order > 1.05) eng.setBondOrder(base + b.a, base + b.b, b.order);
  edited();
  bondTrack.pending.push({ x, y, t0: performance.now(), dur: 600, kind: 'form', big: true });
  if (keep) { placing.rot += Math.PI * 0.37; } else { placing = null; refreshDock(); }
}

/* ======================= actions & keybinds ======================= */
const ACTIONS = [];
const act = (id, group, name, keys, run, o = {}) => ACTIONS.push({ id, group, name, keys, run, ...o });
act('time.play', 'Time', 'Play / pause', ['Space'], togglePlay);
act('time.fwd', 'Time', 'Step forward 1 fs', ['.'], () => stepOnce(1), { repeat: true });
act('time.back', 'Time', 'Step back 1 fs', [','], () => stepOnce(-1), { repeat: true });
act('time.faster', 'Time', 'Faster', [']'], () => setSpeed(nextSpeed(1)));
act('time.slower', 'Time', 'Slower', ['['], () => setSpeed(nextSpeed(-1)));
act('time.normal', 'Time', 'Normal speed (1×)', ['\\'], () => setSpeed(1));
act('tool.grab', 'Tools', 'Grab & pan', ['V'], () => setTool('grab'));
act('tool.erase', 'Tools', 'Erase', ['X'], () => setTool('erase'));
act('tool.heat', 'Tools', 'Heat brush', ['B'], () => setTool('heat'));
act('el.more', 'Tools', 'All elements…', ['E'], () => openElPop());
for (const [sym, key] of [['H', 'H'], ['C', 'C'], ['N', 'N'], ['O', 'O'], ['F', 'F'], ['S', 'S'], ['P', 'P'], ['Cl', 'L'], ['Na', 'A'], ['Br', ''], ['I', ''], ['He', ''], ['Ar', '']])
  act('el.' + sym, 'Elements', BY_SYM[sym].name, key ? [key] : [], () => arm(armed === sym ? null : sym));
act('cond.zero', 'Conditions', 'Absolute zero (0 K)', ['Z'], () => { setT(0); toast('Target 0 K — the heat bath drains all motion'); });
act('cond.room', 'Conditions', 'Room temperature', ['R'], () => { setT(STP_T); toast('Target 298 K'); });
act('cond.flame', 'Conditions', 'Flame (3000 K)', ['Shift+R'], () => { setT(3000); toast('Target 3000 K'); });
act('cond.hotter', 'Conditions', 'Hotter ×1.25', ['Shift+ArrowUp'], () => setT(eng.T < 5 ? 25 : eng.T * 1.25));
act('cond.colder', 'Conditions', 'Colder ÷1.25', ['Shift+ArrowDown'], () => setT(eng.T < 5 ? 0 : eng.T / 1.25));
act('cond.typeT', 'Conditions', 'Type a temperature', ['T'], () => tField.startEdit());
act('cond.bath', 'Conditions', 'Heat bath on / off', [], () => { eng.thermostat = !eng.thermostat; store.set('thermostat', eng.thermostat); toast(eng.thermostat ? 'Heat bath on' : 'Heat bath off — energy conserved'); });
act('view.in', 'View', 'Zoom in', ['='], () => zoomAt(innerWidth / 2, innerHeight / 2, 1 / 1.35));
act('view.out', 'View', 'Zoom out', ['-'], () => zoomAt(innerWidth / 2, innerHeight / 2, 1.35));
act('view.fit', 'View', 'Fit the box', ['0'], () => fitBox(true));
act('view.mode', 'View', 'Cycle look: density / ball-and-stick / space-filling', ['D'], cycleMode);
act('edit.undo', 'Edit', 'Undo', ['Ctrl+Z'], undo);
act('edit.redo', 'Edit', 'Redo', ['Ctrl+Shift+Z', 'Ctrl+Y'], redo);
act('edit.delete', 'Edit', 'Delete selection', ['Delete', 'Backspace'], deleteSelection);
act('edit.all', 'Edit', 'Select all', ['Ctrl+A'], () => { for (let i = 0; i < eng.N; i++) selection.add(i); });
act('edit.pin', 'Edit', 'Pin / unpin selection', ['K'], pinSelection);
act('edit.clear', 'Edit', 'Clear the field', ['Shift+Delete'], () => { if (!eng.N) return; pushUndo(); eng.clear(); selection.clear(); edited(); toast('Field cleared — Ctrl+Z brings it back'); });
act('edit.stopMotion', 'Edit', 'Stop all motion now', [], () => { eng.vel.fill(0); eng.checkpoints.length = 0; toast('All velocities set to zero'); });
act('ui.molecules', 'Panels', 'Molecules tray', ['M'], toggleTray);
act('ui.keys', 'Panels', 'Shortcuts', ['?'], openKeys);
act('ui.palette', 'Panels', 'Command palette', ['Ctrl+K', '/'], openPalette);
function nextSpeed(d) {
  const list = SPEEDS.map(s => s[0]);
  if (d > 0) return list.find(v => v > time.speed + 1e-9) ?? list[list.length - 1];
  return [...list].reverse().find(v => v < time.speed - 1e-9) ?? list[0];
}
function cycleMode() {
  const modes = ['density', 'balls', 'space'], names = { density: 'Electron density', balls: 'Ball and stick', space: 'Space-filling (van der Waals)' };
  R.mode = modes[(modes.indexOf(R.mode) + 1) % 3]; store.set('mode', R.mode); toast(names[R.mode]);
}
function deleteSelection() {
  if (!selection.size) { if (hoverAtom >= 0) { pushUndo(); eng.removeAtoms([hoverAtom]); edited(); } return; }
  pushUndo(); eng.removeAtoms([...selection]); selection.clear(); edited();
}
function pinSelection() {
  const list = selection.size ? [...selection] : hoverAtom >= 0 ? [hoverAtom] : [];
  if (!list.length) return toast('Select atoms (Shift-drag) or hover one to pin it');
  const pin = list.some(i => !eng.pinned[i]);
  for (const i of list) { eng.pinned[i] = pin ? 1 : 0; if (pin) { eng.vel[3 * i] = eng.vel[3 * i + 1] = eng.vel[3 * i + 2] = 0; } }
  eng.checkpoints.length = 0;
  toast(pin ? list.length + ' atoms pinned in place' : 'Unpinned');
}

// key handling
const defaultKeys = Object.fromEntries(ACTIONS.map(a => [a.id, a.keys.slice()]));
let keymap = Object.assign({}, defaultKeys, store.get('keys', {}));
/* Letters and digits bind to the physical key (works on Cyrillic, German, … layouts);
   symbols use the typed character, falling back to the US position. */
const CODE_SYM = { Slash: ['/', '?'], Period: ['.', '>'], Comma: [',', '<'], BracketLeft: ['[', '{'], BracketRight: [']', '}'], Backslash: ['\\', '|'], Equal: ['=', '+'], Minus: ['-', '_'], Semicolon: [';', ':'], Quote: ["'", '"'], Backquote: ['`', '~'] };
function comboOf(e) {
  if (['Control', 'Shift', 'Alt', 'Meta', 'AltGraph'].includes(e.key)) return null;
  const c = e.code || '';
  let k, sym = false;
  if (/^Key[A-Z]$/.test(c)) k = c.slice(3);
  else if (/^Digit\d$/.test(c)) k = c.slice(5);
  else if (e.key === ' ' || c === 'Space') k = 'Space';
  else {
    k = e.key;
    const ascii = k && k.length === 1 && k.charCodeAt(0) < 128;
    if ((!k || k === 'Unidentified' || (k.length === 1 && !ascii)) && CODE_SYM[c]) k = CODE_SYM[c][e.shiftKey ? 1 : 0];
    if (k && k.length === 1) { sym = true; k = k.toUpperCase(); }
  }
  if (!k) return null;
  const mods = [];
  if (e.ctrlKey || e.metaKey) mods.push('Ctrl');
  if (e.altKey) mods.push('Alt');
  if (e.shiftKey && !sym) mods.push('Shift');
  return [...mods, k].join('+');
}
function keyHint(id) { const k = (keymap[id] || [])[0]; return k ? prettyKey(k) : ''; }
function prettyKey(k) { return k.replace('Ctrl+', navigator.platform.includes('Mac') ? '⌘' : 'Ctrl+').replace('Shift+', '⇧').replace('ArrowUp', '↑').replace('ArrowDown', '↓').replace('ArrowLeft', '←').replace('ArrowRight', '→').replace('Space', '␣ Space').replace('Delete', 'Del').replace('Backspace', '⌫'); }
function refreshKeyHints() {
  document.querySelectorAll('.tool').forEach(b => b.querySelector('sup').textContent = keyHint('tool.' + b.dataset.tool).replace('⇧', ''));
  document.querySelectorAll('.el').forEach(b => b.querySelector('sup').textContent = keyHint('el.' + b.dataset.el));
  $('playBtn').title = 'Play / pause (' + keyHint('time.play') + ')';
  $('fwdBtn').title = 'Step forward 1 fs (' + keyHint('time.fwd') + ') — hold to keep stepping';
  $('backBtn').title = 'Step back 1 fs (' + keyHint('time.back') + ') — hold to rewind';
}
let listening = null;
window.addEventListener('keydown', e => {
  if (listening) { e.preventDefault(); e.stopPropagation(); captureKey(e); return; }
  const t = e.target;
  if (t && (t.isContentEditable || t.tagName === 'INPUT' || t.tagName === 'TEXTAREA')) return;
  if (e.key === ' ') spaceHeld = true;
  if (e.key === 'Escape') {
    if (!$('keysVeil').hidden) { closeKeys(); return; }
    if (!$('cmdVeil').hidden) { closePalette(); return; }
    if (openPop) { closePop(); return; }
    if (placing) { placing = null; refreshDock(); return; }
    if (cond) { closeSheet(); return; }
    if (armed) { arm(null); return; }
    if (selection.size) { selection.clear(); return; }
    if ($('tray').classList.contains('open')) { $('tray').classList.remove('open'); return; }
    setTool('grab'); return;
  }
  if (e.key === 'Enter' && cond && !placing) { e.preventDefault(); beginPlacing(); return; }
  if (placing && (e.code === 'KeyQ' || e.code === 'KeyE') && !e.ctrlKey) { placing.rot += (e.code === 'KeyQ' ? -1 : 1) * Math.PI / 12; e.preventDefault(); return; }
  const combo = comboOf(e); if (!combo) return;
  const a = ACTIONS.find(x => (keymap[x.id] || []).includes(combo));
  if (!a) return;
  if (e.repeat && !a.repeat) { e.preventDefault(); return; }
  e.preventDefault();
  a.run();
});
window.addEventListener('keyup', e => { if (e.key === ' ') spaceHeld = false; });

function openKeys() {
  closePop();
  const groups = [...new Set(ACTIONS.map(a => a.group))];
  $('keysGrid').innerHTML = groups.map(g => '<div class="kgroup"><h3>' + g + '</h3>' + ACTIONS.filter(a => a.group === g).map(a =>
    '<div class="krow"><span>' + esc(a.name) + '</span><button class="kcap' + ((keymap[a.id] || []).length ? '' : ' empty') + '" data-id="' + a.id + '">' + ((keymap[a.id] || []).length ? esc(prettyKey(keymap[a.id][0])) : 'add') + '</button></div>').join('') + '</div>').join('') +
    '<div class="kgroup"><h3>Mouse</h3>' + [['Scroll', 'Zoom (0.05 – 5 nm)'], ['Drag empty space', 'Pan'], ['Drag an atom', 'Pull it (running) / move its molecule (paused)'], ['Alt + drag atom', 'Move one atom (paused)'], ['Shift + drag', 'Select a region'], ['Right-click atom', 'Delete (Alt: whole molecule)'], ['Drag while placing', 'Throw with that velocity'], ['Drag box edge', 'Resize the container'], ['Q / E', 'Rotate a molecule before placing']].map(r => '<div class="krow"><span>' + r[1] + '</span><span class="kcap" style="border:0;background:none;color:var(--faint)">' + r[0] + '</span></div>').join('') + '</div>';
  $('keysGrid').querySelectorAll('button.kcap').forEach(b => b.onclick = () => {
    document.querySelectorAll('.kcap.listening').forEach(x => x.classList.remove('listening'));
    listening = b; b.classList.add('listening'); b.textContent = 'press a key';
  });
  $('keysVeil').hidden = false;
}
function captureKey(e) {
  const b = listening; listening = null; b.classList.remove('listening');
  const id = b.dataset.id;
  if (e.key === 'Escape') { openKeys(); return; }
  if (e.key === 'Backspace' && !keymap[id]?.length) { openKeys(); return; }
  const combo = comboOf(e); if (!combo) { openKeys(); return; }
  let moved = null;
  for (const a of ACTIONS) if (a.id !== id && (keymap[a.id] || []).includes(combo)) { keymap[a.id] = keymap[a.id].filter(k => k !== combo); moved = a.name; }
  keymap[id] = [combo];
  store.set('keys', keymap);
  openKeys(); refreshKeyHints();
  const nb = $('keysGrid').querySelector('[data-id="' + id + '"]'); nb && nb.classList.add('flash');
  if (moved) toast(prettyKey(combo) + ' moved here from “' + moved + '”');
}
function closeKeys() { listening = null; $('keysVeil').hidden = true; }
$('keysBtn').onclick = openKeys; $('keysClose').onclick = closeKeys;
$('keysReset').onclick = () => { keymap = Object.assign({}, defaultKeys); store.set('keys', keymap); openKeys(); refreshKeyHints(); toast('Shortcuts reset'); };
$('keysVeil').addEventListener('pointerdown', e => { if (e.target === $('keysVeil')) closeKeys(); });

/* ======================= command palette ======================= */
let cmdItems = [], cmdSel = 0;
function openPalette() { closePop(); $('cmdVeil').hidden = false; $('cmdInput').value = ''; renderPalette(); $('cmdInput').focus(); }
function closePalette() { $('cmdVeil').hidden = true; }
function renderPalette() {
  const q = $('cmdInput').value.trim().toLowerCase();
  const items = [];
  const tm = parseTemp(q); if (/k$|°c$|°f$|c$/.test(q) && Number.isFinite(tm)) items.push({ grp: 'Set', name: 'Temperature ' + fmtT(tm) + ' K', run: () => setT(tm) });
  const sm = q.match(/^([\d.]+)\s*[x×]$/); if (sm) items.push({ grp: 'Set', name: 'Speed ' + sm[1] + '×', run: () => setSpeed(+sm[1]) });
  for (const m of inbox) if (!q || (m.name || '').toLowerCase().includes(q) || m.formula.toLowerCase().includes(q)) items.push({ grp: 'Molecule', name: (m.name || pretty(m.formula)), run: () => startConditioning(m) });
  for (const a of ACTIONS) if (!q || a.name.toLowerCase().includes(q) || a.group.toLowerCase().includes(q)) items.push({ grp: a.group, name: a.name, key: keyHint(a.id), run: a.run });
  for (const e of ELEMENTS) if (q && (e.name.toLowerCase().startsWith(q) || e.sym.toLowerCase() === q) && !ACTIONS.some(a => a.id === 'el.' + e.sym)) items.push({ grp: 'Element', name: e.name, run: () => arm(e.sym) });
  cmdItems = items.slice(0, 40); cmdSel = 0;
  $('cmdList').innerHTML = cmdItems.map((it, i) => '<li data-i="' + i + '" class="' + (i ? '' : 'on') + '"><span><span class="grp">' + esc(it.grp) + '</span>' + esc(it.name) + '</span>' + (it.key ? '<kbd>' + esc(it.key) + '</kbd>' : '') + '</li>').join('') || '<li>No matching action</li>';
  $('cmdList').querySelectorAll('li[data-i]').forEach(li => { li.onclick = () => runPalette(+li.dataset.i); li.onpointermove = () => selPalette(+li.dataset.i); });
}
function selPalette(i) { cmdSel = clamp(i, 0, cmdItems.length - 1); $('cmdList').querySelectorAll('li').forEach((li, k) => li.classList.toggle('on', k === cmdSel)); const on = $('cmdList').children[cmdSel]; on && on.scrollIntoView({ block: 'nearest' }); }
function runPalette(i) { const it = cmdItems[i]; closePalette(); if (it) it.run(); }
$('cmdInput').addEventListener('input', renderPalette);
$('cmdInput').addEventListener('keydown', e => {
  if (e.key === 'ArrowDown') { e.preventDefault(); selPalette(cmdSel + 1); }
  else if (e.key === 'ArrowUp') { e.preventDefault(); selPalette(cmdSel - 1); }
  else if (e.key === 'Enter') { e.preventDefault(); runPalette(cmdSel); }
  else if (e.key === 'Escape') { e.preventDefault(); closePalette(); }
});
$('cmdVeil').addEventListener('pointerdown', e => { if (e.target === $('cmdVeil')) closePalette(); });

/* ======================= persistence ======================= */
let saveTimer = 0;
function scheduleSave() { clearTimeout(saveTimer); saveTimer = setTimeout(saveScene, 800); }
function saveScene() {
  try { localStorage.setItem('cp.scene', JSON.stringify(eng.toJSON())); } catch (e) { }
}
function loadScene() {
  try {
    const s = JSON.parse(localStorage.getItem('cp.scene') || 'null');
    if (!s || !Array.isArray(s.atoms)) return;
    for (const a of s.atoms) if (BY_SYM[a[0]]) eng.addAtom(a[0], a[1], a[2], a[3], { v: [a[4], a[5], a[6]], charge: a[7], V: a[8] });
    eng.time = s.time || 0;
  } catch (e) { }
}
window.addEventListener('beforeunload', saveScene);
setInterval(() => { if (time.playing) saveScene(); }, 5000);
function updateEmpty() { $('emptyHint').classList.toggle('gone', eng.N > 0 || !!cond || !!placing); }

/* ======================= main loop ======================= */
let lastUI = 0, lastSpark = 0;
function frame(now) {
  const frameMs = now - time.last, dt = Math.min(0.1, frameMs / 1000); time.last = now;
  time.frameEMA = (time.frameEMA || 16) * 0.9 + Math.min(100, frameMs) * 0.1;
  // camera animation
  if (camAnim) {
    const k = Math.min(1, (now - camAnim.t0) / camAnim.dur), e = 1 - Math.pow(1 - k, 3);
    for (const key of ['cx', 'cy', 'span']) R.cam[key] = camAnim.from[key] + (camAnim.to[key] - camAnim.from[key]) * e;
    if (k >= 1) camAnim = null;
  }
  // physics: fixed 1 fs steps inside a render-decoupled accumulator
  // physics gets ~55 % of the real frame interval, so slow displays still simulate at full rate
  const frameBudget = clamp(time.frameEMA * 0.55, 6, 30) - (cond ? 3 : 0);
  if (time.playing) {
    time.acc += dt * stepsPerSecond(time.speed);
    const t0 = performance.now(); let n = 0; time.limited = false;
    while (time.acc >= 1) {
      eng.step(); time.acc -= 1; n++;
      if ((n & 7) === 0 && performance.now() - t0 > frameBudget) { if (time.acc >= 1) { time.limited = true; time.acc %= 1; } break; }
    }
    time.stepsWindow += n;
    applyBrush(dt);
  } else if (brush.active) { applyBrush(dt); }
  if (now - time.windowStart > 500) { const r = time.stepsWindow / ((now - time.windowStart) / 1000); time.rateEMA = time.rateEMA * 0.4 + r * 0.6; time.stepsWindow = 0; time.windowStart = now; }
  if (cond) tickConditioning(4);
  // interpolation alpha between the previous and the current step
  let alpha = 1;
  if (time.playing) alpha = time.acc;
  else if (time.manual) { const k = Math.min(1, (now - time.manual.t0) / time.manual.dur); alpha = 1 - Math.pow(1 - k, 3); if (k >= 1) time.manual = null; }
  const n3 = 3 * eng.N;
  if (rp.length < n3 || rp.length > n3 + 3000) rp = new Float64Array(n3 + 300);
  const P = eng.pos, Q = eng.prev;
  for (let k = 0; k < n3; k++) rp[k] = Q[k] + (P[k] - Q[k]) * alpha;
  if (eng.needForces) eng.refresh(); // edits while paused
  const bonds = eng.bonds(0.3);
  updateBondEvents(bonds, now);
  if (time.playing || species.reset) updateSpecies(now);
  if ((eng.N > 0) === !$('emptyHint').classList.contains('gone')) updateEmpty();
  // ghost
  let ghost = null;
  if (placing) {
    const g = gesture && gesture.type === 'placeMol' ? gesture : null;
    const gx = g ? g.wx : placing.x, gy = g ? g.wy : placing.y, atoms = ghostAtoms(gx, gy);
    ghost = { atoms, bonds: placing.bonds, ok: ghostOK(atoms) };
    if (g && g.cx != null && Math.hypot(g.cx - g.wx, g.cy - g.wy) * R.scale > 8) {
      const v = flingVel(g.cx - g.wx, g.cy - g.wy);
      ghost.fling = { x0: gx, y0: gy, x1: g.cx, y1: g.cy, label: flingLabel(v, placing.mass) };
    }
  } else if (armed && gesture && gesture.type === 'placeAtom' && gesture.cx != null && Math.hypot(gesture.cx - gesture.wx, gesture.cy - gesture.wy) * R.scale > 8) {
    const v = flingVel(gesture.cx - gesture.wx, gesture.cy - gesture.wy);
    ghost = { atoms: [{ t: BY_SYM[armed].t, x: gesture.wx, y: gesture.wy, z: 0 }], bonds: [], ok: true, fling: { x0: gesture.wx, y0: gesture.wy, x1: gesture.cx, y1: gesture.cy, label: flingLabel(v, BY_SYM[armed].mass) } };
  } else if (armed && !gesture) {
    ghost = { atoms: [{ t: BY_SYM[armed].t, x: lastMouse.wx, y: lastMouse.wy, z: 0 }], bonds: [], ok: true };
  }
  R.draw({
    N: eng.N, type: eng.type, pos: rp, bonds, box: eng.box, boxHot: gesture && gesture.type === 'box' ? gesture.edge : boxHot,
    hover: gesture ? -1 : hoverAtom, eraseHover: tool === 'erase' && !armed, selected: selection, pinned: eng.pinned.subarray(0, eng.N), ghost, flashes, now,
    tweezer: eng.tweezer, brush: tool === 'heat' && !armed && !placing ? brush : null,
    marquee: gesture && gesture.type === 'marquee' ? gesture : null
  });
  // UI text, ~10 Hz
  if (now - lastUI > 100) {
    lastUI = now;
    $('clock').textContent = fmtTime(eng.time);
    const rate = $('rate');
    if (time.playing) {
      rate.textContent = (time.limited ? 'CPU-bound · ' : '') + fmtRate(time.rateEMA) + (time.limited ? ' · ' + (time.rateEMA / stepsPerSecond(time.speed)).toFixed(2) + '×' : '');
      rate.classList.toggle('limited', time.limited);
    } else { rate.textContent = eng.canStepBack() ? 'paused · ' + fmtTime(Math.min(eng.historySpan(), eng.time)) + ' rewindable' : 'paused'; rate.classList.remove('limited'); }
    $('tMeasured').textContent = eng.N ? eng.temperature().toFixed(eng.temperature() < 10 ? 2 : 0) + ' K now' : '';
    $('pVal').textContent = eng.N ? fmtP(eng.pressureEMA) : '—';
    $('backBtn').disabled = !eng.canStepBack();
    // scale bar
    const s = R.scale, cands = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2], nm = cands.find(v => v * 10 * s >= 60) || 2;
    $('scale').querySelector('i').style.width = (nm * 10 * s) + 'px';
    $('scale').querySelector('span').textContent = (nm >= 1 ? nm + ' nm' : (nm * 10).toFixed(nm < 0.1 ? 1 : 0) + ' Å') + '  ·  view ' + (R.cam.span / 10).toFixed(2) + ' nm';
    ageFeed(now);
    renderTPop(); renderPPop();
  }
  if (now - lastSpark > 250) {
    lastSpark = now;
    if (time.playing) { tHist.push(eng.temperature()); pHist.push(eng.pressureEMA); if (tHist.length > 120) tHist.shift(); if (pHist.length > 120) pHist.shift(); }
    if ($('tPop').classList.contains('open')) spark('tSpark', tHist, 'rgb(255,178,63)', 0);
    if ($('pPop').classList.contains('open')) spark('pSpark', pHist, 'rgb(127,167,201)', 0);
  }
  requestAnimationFrame(frame);
}

/* ======================= boot ======================= */
resize();
loadScene();
eng.refresh();
fitBox(false);
buildDock();
renderTray();
paintSeg();
updateEmpty();
readHash();
requestAnimationFrame(t => { time.last = t; frame(t); });
window.chemPlayground = { eng, R, receive, setT, setSpeed, togglePlay, startConditioning, beginPlacing, dropMolecule, tickConditioning, get cond() { return cond; }, get placing() { return placing; } }; // console / test hooks
})();
