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
const SPAN_MIN = 0.5, SPAN_MAX = 150; // manual visible width 0.05–15 nm; Fit may go wider

/* Three ways to hold a temperature, and the difference matters:
   wall   — a heater at the boundary; the interior warms through collisions, so the reading
            lags the setpoint and never quite reaches it against an open wall.
   kelvin — every atom, every step, at exactly the setpoint. Paired with a void wall it is a
            steady state by construction: what the wall takes, the stat returns the same tick.
   off    — no stat at all; setting a temperature is a one-off edit of the velocities. */
const bathMode = () => !eng.thermostat ? 'off' : eng.thermostatMode === 'kelvin' ? 'kelvin' : 'wall';
const wallHeater = () => bathMode() === 'wall';
function setBathMode(mode) {
  if (bathMode() === mode) return;
  pushUndo();
  eng.thermostat = mode !== 'off';
  if (mode !== 'off') eng.thermostatMode = mode === 'kelvin' ? 'kelvin' : 'wall';
  if (eng.thermostat) eng.setTemperature(eng.T);
  eng.checkpoints.length = 0;
  store.set('thermostat', eng.thermostat); store.set('bathMode', eng.thermostatMode);
  tField.render(); paintT(); scheduleSave();
  if ($('tPop').classList.contains('open')) openTPop();
}
function toggleBath() { setBathMode(eng.thermostat ? 'off' : (store.get('bathMode', 'wall') === 'kelvin' ? 'kelvin' : 'wall')); }
function resampleMotion() {
  if (!eng.N) return toast('Place atoms or molecules first');
  pushUndo(); eng.thermalize(eng.T); edited();
  toast('Fresh thermal velocities at ' + fmtT(eng.T) + ' K · press play to run');
}
function safeStep(engine) {
  try { engine.step(); return true; }
  catch (error) {
    if (engine === eng) setPlaying(false);
    toast(error.message);
    return false;
  }
}
function paintViews() {} // contours are the only field representation
let hintTimer = 0;
function updateLab() {
  $('statusDot').classList.toggle('running', time.playing);
  $('statusDot').title = time.playing ? 'Running' : 'Paused';
  $('atomCount').textContent = eng.N + (eng.N === 1 ? ' atom' : ' atoms');
  $('liveTemperature').textContent = eng.N ? eng.temperature().toFixed(1) + ' K' : '—';
  $('liveEnergy').textContent = eng.N ? (eng.Epot + eng.kinetic()).toFixed(1) + ' kJ/mol' : '—';
  const mode = bathMode();
  $('bathBtn').textContent = mode === 'wall' ? (eng.wallT - 273.15).toFixed(1) + ' °C wall'
    : mode === 'kelvin' ? fmtT(eng.T) + ' K held' : 'Off · direct setting';
  $('bathBtn').setAttribute('aria-pressed', eng.thermostat);
  $('physicsNotice').hidden = !eng.clamped;
  $('physicsNotice').textContent = eng.clamped ? eng.clamped + ' safety speed clamps · energy affected' : '';
  const hint = placing ? 'Click to place · Q/E rotate · right-click for the hand' : armed ? armed + ' selected · click to place · right-click for the hand' : tool === 'heat' ? 'Drag to heat · Shift-drag to cool' : tool === 'erase' ? 'Click an atom to erase · Alt: molecule' : viewTilted() ? 'Turning · release to face the chamber again' : '';
  const ctx = $('toolContext');
  if (ctx.textContent !== hint) {
    ctx.textContent = hint; ctx.classList.toggle('show',!!hint);
    clearTimeout(hintTimer); hintTimer = setTimeout(() => ctx.classList.remove('show'), 2600);
  }
}

/* ======================= engine ======================= */
const savedBox = store.get('box', { w: 36, h: 20, d: 12 });
const eng = new Engine({ width: savedBox.w, height: savedBox.h, depth: savedBox.d, T: store.get('T', STP_T) });
eng.thermostat = store.get('thermostat', true);
eng.thermostatMode = store.get('bathMode', 'wall') === 'kelvin' ? 'kelvin' : 'wall';
const savedBounds = store.get('bounds', null);
if (savedBounds) {
  if (savedBounds.mode === 'forcefield' || savedBounds.mode === 'solid') eng.boundsMode = savedBounds.mode;
  // version 3 reworked the channels; older saves only carried the two that still exist
  eng.voidTemperature = !!savedBounds.vt;
  eng.voidPressure = savedBounds.version >= 3 && !!savedBounds.vp;
  eng.voidVelocity = savedBounds.version >= 3 && !!savedBounds.vv;
  eng.pressureControl = !!savedBounds.pc;
  if (Number.isFinite(savedBounds.pt)) eng.pressureTarget = savedBounds.pt;
}
eng.recording = true;

const canvas = $('field');
const inspector = new AtomInspector({ engine: eng, pause: () => setPlaying(false), focus: () => canvas.focus(),
  screenOf: i => R.toScreen(rp[3 * i], rp[3 * i + 1]),
  remove: id => {
  const i = eng.ids.subarray(0, eng.N).indexOf(id);
  if (i < 0) return;
  pushUndo(); eng.removeAtoms([i]); edited();
} });
const R = new FieldRenderer(canvas);
R.mode = 'density';
const appearance = Object.assign({quality:'balanced',motion:'system',grid:true,labels:true,budget:'balanced'},store.get('appearance',{}));
const motionQuery=matchMedia('(prefers-reduced-motion: reduce)');
function applyAppearance(){
  document.documentElement.dataset.motion=appearance.motion==='reduced'||(appearance.motion==='system'&&motionQuery.matches)?'reduced':'full';
  R.showGrid=!!appearance.grid; R.showLabels=!!appearance.labels;
  const cell={fine:3,balanced:5,fast:8}[appearance.quality]||5;
  if(R.cell!==cell){R.cell=cell;R.resize();}
  document.querySelectorAll('svg').forEach(svg=>{if(svg.pauseAnimations){if(document.documentElement.dataset.motion==='reduced')svg.pauseAnimations();else svg.unpauseAnimations();}});
  store.set('appearance',appearance);
}
motionQuery.addEventListener('change',applyAppearance);
applyAppearance();
$('fitBtn').onclick = () => fitBox(true);
$('measurementsBtn').onclick = () => {
  const panel = $('measurementsPanel'), open = !panel.classList.contains('open');
  panel.classList.toggle('open', open); panel.inert = !open;
  $('measurementsBtn').setAttribute('aria-expanded', open);
};
$('bathBtn').onclick = toggleBath;
$('resampleBtn').onclick = resampleMotion;
document.addEventListener('visibilitychange', () => {
  time.last = performance.now(); time.acc = 0;
  if (document.hidden && time.playing) { setPlaying(false); toast('Paused while this tab is hidden'); }
});
let rp = new Float64Array(0);        // interpolated render positions

/* ======================= time ======================= */
const time = { playing: false, speed: store.get('speed', 1), acc: 0, last: performance.now(), rateEMA: 0, stepsWindow: 0, windowStart: performance.now(), limited: false, manual: null };
const stepsPerSecond = ChemProtocol.stepsPerSecond;
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
  tField.render(); paintT();
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
  scheduleSave();
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
  const mobile=W<=760, left=mobile?30:104, right=40, top=mobile?160:150, bottom=mobile?230:136;
  const aw=Math.max(100,W-left-right), ah=Math.max(100,H-top-bottom);
  const span=clamp(Math.max((b.x1-b.x0)*W/aw,(b.y1-b.y0)*W/ah),SPAN_MIN,5000);
  camTo((b.x0+b.x1)/2-(left-right)*span/(2*W), (b.y0+b.y1)/2-(top-bottom)*span/(2*W),span,animate);
}
let camAnim = null;
function camTo(cx, cy, span, animate) {
  if (!animate || document.documentElement.dataset.motion==='reduced') { R.cam.cx = cx; R.cam.cy = cy; R.cam.span = span; camAnim = null; return; }
  camAnim = { from: { ...R.cam }, to: { cx, cy, span }, t0: performance.now(), dur: 420 };
}
function zoomAt(sx, sy, factor) {
  const [wx, wy] = R.toWorld(sx, sy);
  const span = clamp(R.cam.span * factor, SPAN_MIN, Math.max(SPAN_MAX, R.cam.span));
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
  el.setAttribute('role', 'spinbutton');
  el.setAttribute('aria-label', el.getAttribute('aria-label') || el.title || ({tField:'Target temperature', speedField:'Playback speed', cT:'Conditioning temperature', cP:'Nominal confinement pressure', bW:'Box width', bH:'Box height', bD:'Slab depth'}[el.id] || 'Value'));
  const api = {
    render() {
      if (el.classList.contains('editing')) return;
      const v = o.get(), p = toP(v);
      el.innerHTML = esc(o.fmt(v)) + (o.unit ? '<span class="u">' + o.unit + '</span>' : '');
      el.style.setProperty('--p', clamp(p, 0, 1));
      el.setAttribute('aria-valuenow', o.ariaValue ? o.ariaValue(v) : v);
      el.setAttribute('aria-valuetext', o.fmt(v) + (o.unit || ''));
      el.classList.toggle('over', p > 1.0001 || p < -0.0001);
    },
    set(v, commit) { if (!Number.isFinite(v)) return; v = o.snap ? o.snap(v) : v; if (o.hardMin != null) v = Math.max(o.hardMin, v); if (o.hardMax != null) v = Math.min(o.hardMax, v); o.set(v, commit); api.render(); }
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
    // inside a scrolling panel the wheel belongs to the panel, not to the value under the pointer
    if (o.wheel === false) return;
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
  if (!Number.isFinite(T)) return;
  eng.setTemperature(clamp(T, 0, 50000));
  store.set('T', eng.T); tField.render(); paintT(); scheduleSave();
  if (!quiet) renderTPop();
}
function paintT() { $('bbDot').style.color = blackbody(eng.T); }
/* The kinetic temperature of a handful of atoms genuinely wanders: its relative spread is
   about sqrt(2/3N), so twenty atoms swing by tens of kelvin between one instant and the next.
   That is the sample being small, not friction. The gauge reads a running mean and shows the
   measured spread beside it, so the number is stable while the physics stays honest. */
const tStat = { mean: 0, var: 0, live: false };
function readTemperature() {
  const val = $('tMeasuredVal'), band = $('tMeasuredBand'), gauge = $('tMeasured');
  if (!eng.N) { tStat.live = false; val.textContent = ''; band.style.opacity = '0'; return; }
  const T = eng.temperature();
  if (!tStat.live) { tStat.mean = T; tStat.var = 0; tStat.live = true; }
  const a = 0.09, d = T - tStat.mean;
  tStat.mean += a * d;
  tStat.var = (1 - a) * (tStat.var + a * d * d);
  const sd = Math.sqrt(tStat.var), mean = tStat.mean;
  val.textContent = mean.toFixed(mean < 10 ? 2 : 0) + ' K';
  band.style.opacity = '1';
  const half = clamp(sd / Math.max(1, mean) * 130, 1.5, 46);
  band.style.setProperty('--lo', (50 - half) + '%');
  band.style.setProperty('--hi', (50 - half) + '%');
  gauge.title = 'Measured kinetic temperature ' + mean.toFixed(1) + ' K ± ' + sd.toFixed(1) +
    ' K (now ' + T.toFixed(1) + ' K). A sample of ' + eng.N + ' atoms fluctuates by about ' +
    (mean * Math.sqrt(2 / (3 * Math.max(1, eng.N)))).toFixed(1) + ' K — that is its size, not friction.';
}
const fmtT = v => v >= 10000 ? (v / 1000).toFixed(1) + 'k' : v >= 100 ? Math.round(v).toString() : v >= 10 ? v.toFixed(1) : v.toFixed(2);
const tField = scrub($('tField'), {
  get: () => eng.T, set: v => setT(v, true), get min() { return wallHeater() ? 288.15 : 0; }, get max() { return wallHeater() ? 623.15 : 6000; }, off: 40, map: 'log', get unit() { return wallHeater() ? '°C' : 'K'; }, hardMin: 0, hardMax: 50000,
  fmt: v => wallHeater() ? (v - 273.15).toFixed(1) : fmtT(v),
  ariaValue: v => wallHeater() ? v - 273.15 : v,
  edit: v => String(+(wallHeater() ? v - 273.15 : v).toFixed(2)),
  parse: txt => parseTemp(wallHeater() && /^-?[\d.,]+$/.test(txt.trim()) ? txt + ' C' : txt)
});
function parseTemp(txt) {
  const m = txt.trim().toLowerCase().replace(',', '.').match(/^(-?[\d.]+)\s*(k|°?c|°?f)?$/);
  if (!m) return NaN;
  const v = parseFloat(m[1]), u = (m[2] || 'k').replace('°', '');
  return u === 'c' ? v + 273.15 : u === 'f' ? (v - 32) * 5 / 9 + 273.15 : v;
}
paintT();
const T_PRESETS = [[0, 'Absolute zero', 'classical cooling'], [77, 'Liquid nitrogen', ''], [STP_T, 'Room', '25 °C'], [373.15, 'Boiling water', '100 °C'], [1000, 'Red heat', ''], [3000, 'Flame', 'high-temperature model'], [6000, 'Sun surface', 'extreme conditions']];
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
  pop.style.left = clamp(x, 8, Math.max(8, innerWidth - w - 8)) + 'px'; pop.style.top = clamp(y, 8, Math.max(8, innerHeight - h - 8)) + 'px';
}
function closePop() { if (openPop) { openPop.pop.classList.remove('open'); openPop = null; } }
document.addEventListener('pointerdown', e => {
  if (openPop && !openPop.pop.contains(e.target) && !openPop.anchor.contains(e.target)) closePop();
}, true);

function renderTPop() {
  const pop = $('tPop'); if (!pop.classList.contains('open')) return;
  pop.querySelector('.head b').textContent = eng.temperature().toFixed(1) + ' K';
  pop.querySelectorAll('.row[data-t]').forEach(r => r.classList.toggle('on', Math.abs(+r.dataset.t - eng.T) < 0.01));
  pop.querySelectorAll('#tPopBath button').forEach(b => b.setAttribute('aria-pressed', b.dataset.v === bathMode()));
  const wall = $('wallNow'); if (wall) wall.textContent = (eng.wallT - 273.15).toFixed(2) + ' °C';
  const heat = $('wallHeat'); if (heat) heat.textContent = eng.heatToSample.toFixed(2) + ' kJ/mol';
  const sw = $('statWork'); if (sw) sw.textContent = eng.kelvinWork.toFixed(1) + ' kJ/mol';
}
function openTPop() {
  const pop = $('tPop');
  const presets = wallHeater() ? [[288.15, 'Cool', '15 °C'], [298.15, 'Room', '25 °C'], [373.15, 'Warm', '100 °C'], [473.15, 'Hot', '200 °C'], [623.15, 'Maximum', '350 °C']] : T_PRESETS;
  pop.innerHTML = '<div class="head">Measured <b></b></div><canvas class="spark" id="tSpark"></canvas><div class="sep"></div>' +
    presets.map((p, i) => '<button class="row" style="--i:' + i + '" data-t="' + p[0] + '"><span class="sw" style="background:' + blackbody(p[0]) + '"></span><span class="k">' + (wallHeater() ? (p[0] - 273.15).toFixed(0) + ' °C' : fmtT(p[0]) + ' K') + '</span><span class="lbl">' + p[1] + '</span></button>').join('') +
    '<div class="sep"></div><div class="seg-mini bath-seg" id="tPopBath"><button data-v="wall">Wall</button><button data-v="kelvin">Kelvin</button><button data-v="off">Off</button></div>' +
    (wallHeater() ? '<div class="kv"><span>Wall now</span><span id="wallNow"></span></div><div class="kv"><span>Heat into sample</span><span id="wallHeat"></span></div><div class="kv"><span>Heater response τ</span><span class="scrub" id="wallTau"></span></div><p class="thermal-note">The wall reaches 63% of a temperature change in τ. Time is simulated, not playback time. Only the 0.2 nm boundary layer contacts the reservoir; the interior warms through interactions. This is a specified nanoscale heater, not a calibrated household thermostat.</p>' : bathMode() === 'kelvin' ? '<div class="kv"><span>Stat work</span><span id="statWork"></span></div><p class="thermal-note">Every atom is rescaled to this temperature on every step, so the measured value is the set value. Against a void wall that is a steady state: what the wall removes, the stat returns the same tick. It fixes the total kinetic energy, not the distribution.</p>' : '<p class="thermal-note">Editing temperature immediately sets the sample’s kinetic temperature. Then it evolves without a heat bath. This changes the velocities, not the geometry or phase.</p>');
  pop.querySelectorAll('.row[data-t]').forEach(r => r.onclick = () => { setT(+r.dataset.t); });
  pop.querySelectorAll('#tPopBath button').forEach(b => {
    b.setAttribute('aria-pressed', b.dataset.v === bathMode());
    b.onclick = () => { setBathMode(b.dataset.v); openTPop(); };
  });
  if (wallHeater()) scrub($('wallTau'), { get: () => eng.wallTau / 1000, set: v => { eng.wallTau = v * 1000; eng.checkpoints.length = 0; scheduleSave(); }, min: 1, max: 100, hardMin: .01, hardMax: 1e15, unit: 'ps', fmt: v => String(+v.toPrecision(4)) });
  showPop(pop, $('gT'), 'below'); renderTPop();
}
$('gT').addEventListener('click', e => { if (e.target.closest('.scrub')) return; openPop && openPop.pop.id === 'tPop' ? closePop() : openTPop(); });
$('bbDot').title = 'Temperature presets';

function openPPop() {
  const pop = $('pPop');
  pop.innerHTML = '<div class="head">Wall pressure <b id="pNow"></b></div><canvas class="spark" id="pSpark"></canvas>' +
    '<div class="kv"><span>Volume</span><span id="pVol"></span></div><div class="kv"><span>Ideal-gas estimate (fragments)</span><span id="pIdeal"></span></div><div class="kv"><span>Atoms</span><span id="pAtoms"></span></div>' +
    '<div class="sep"></div><div class="kv"><span>Box width</span><span class="scrub" id="bW"></span></div><div class="kv"><span>Box height</span><span class="scrub" id="bH"></span></div><div class="kv"><span>Slab depth</span><span class="scrub" id="bD"></span></div>';
  const mk = (id, get, set, min, max) => scrub($(id), { get, set: (v, c) => { set(v); eng.touch(); if (c) { saveBox(); scheduleSave(); } }, min, max, unit: 'nm', fmt: v => v.toFixed(2), hardMin: min, hardMax: max });
  mk('bW', () => (eng.box.x1 - eng.box.x0) / 10, v => { eng.box.x1 = eng.box.x0 + v * 10; }, 1, 50);
  mk('bH', () => (eng.box.y1 - eng.box.y0) / 10, v => { eng.box.y1 = eng.box.y0 + v * 10; }, 1, 50);
  mk('bD', () => (eng.box.z1 - eng.box.z0) / 10, v => { eng.box.z0 = -v * 5; eng.box.z1 = v * 5; }, 0.1, 50);
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
$('gP').addEventListener('keydown', e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); e.stopPropagation(); $('gP').click(); } });
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
function setSpeed(v) { if (!Number.isFinite(v)) return; v = clamp(v, 0.01, 20); time.speed = v; store.set('speed', v); speedField.render(); renderSpeedPop(); }
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
  } else if (!safeStep(eng)) return false;
  bondTrack.reset = true; species.reset = true;
  $('feed').replaceChildren();
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
  b.addEventListener('click', e => { if (e.detail === 0) stepOnce(dir); });
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
function staggerDock() { // the dock deals itself out on first paint
  document.querySelectorAll('.dock .tool, .dock .el').forEach((b, i) => b.style.animationDelay = (0.16 + i * 0.022).toFixed(3) + 's');
}
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
  refreshDock(); refreshKeyHints(); staggerDock();
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
  document.querySelectorAll('.tool').forEach(b => (b.classList.toggle('on', !armed && !placing && tool === b.dataset.tool), b.setAttribute('aria-pressed', !armed && !placing && tool === b.dataset.tool)));
  document.querySelectorAll('.el').forEach(b => (b.classList.toggle('on', armed === b.dataset.el), b.setAttribute('aria-pressed', armed === b.dataset.el)));
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
/* (x, y) arrive in the frame the pointer works in; the atom itself is born in world coordinates,
   a little off the mid-plane so nothing lands exactly on top of anything else. */
function placeAtom(sym, x, y, v) {
  const jitter = (eng.rand() - 0.5) * 0.6;
  const w = viewTilted() ? fromView(x, y, viewPivot()[2] + jitter) : [x, y, jitter];
  if (tooClose(w[0], w[1], w[2], 1.1)) { toast('Too close to another atom'); return; }
  pushUndo();
  const tv = thermalVel(sym, eng.T);
  const fl = v ? viewDelta(v[0], v[1]) : null;
  eng.addAtom(sym, w[0], w[1], w[2], { v: fl ? [fl[0] + tv[0] * 0.3, fl[1] + tv[1] * 0.3, fl[2] + tv[2] * 0.3] : tv });
  edited();
  bondTrack.pending.push({ x, y, t0: performance.now(), dur: 420, kind: 'form', big: true });
}
// fling: screen drag → velocity (Å/fs). 10 Å of drag ≈ 0.02 Å/fs (2 km/s)
function flingVel(dx, dy) { const k = 0.002, v = [dx * k, dy * k]; const s = Math.hypot(v[0], v[1]), cap = 0.12; if (s > cap) { v[0] *= cap / s; v[1] *= cap / s; } return v; }
function flingLabel(v, mass) { const ke = 0.5 * mass * (v[0] * v[0] + v[1] * v[1]) * 1e4; return (Math.hypot(v[0], v[1]) * 1e5 / 1000).toFixed(1) + ' km/s · ' + ke.toFixed(ke < 10 ? 1 : 0) + ' kJ/mol'; }

/* ======================= view rotation ======================= */
/* The chamber is a real 3D slab; the default view looks straight down the z axis at it.
   Right-dragging empty space turns the slab about its own centre. Everything the pointer does
   keeps working because the render positions, the hit test and the pointer all share one frame:
   world coordinates go through `toView` once per frame, and anything that writes back into the
   engine — dragging an atom, placing one, the tweezer — comes back through `fromView`. */
const view = { yaw: 0, pitch: 0 };
R.view = view;
const viewTilted = () => Math.abs(view.yaw) > 1e-6 || Math.abs(view.pitch) > 1e-6;
function viewPivot() { const b = eng.box; return [(b.x0 + b.x1) / 2, (b.y0 + b.y1) / 2, (b.z0 + b.z1) / 2]; }
function toView(x, y, z) {
  const p = viewPivot(), cy = Math.cos(view.yaw), sy = Math.sin(view.yaw), cp = Math.cos(view.pitch), sp = Math.sin(view.pitch);
  const dx = x - p[0], dy = y - p[1], dz = z - p[2];
  const x1 = cy * dx + sy * dz, z1 = cy * dz - sy * dx;
  return [p[0] + x1, p[1] + cp * dy - sp * z1, p[2] + sp * dy + cp * z1];
}
function fromView(x, y, z) {
  const p = viewPivot(), cy = Math.cos(view.yaw), sy = Math.sin(view.yaw), cp = Math.cos(view.pitch), sp = Math.sin(view.pitch);
  const dx = x - p[0], dy = y - p[1], dz = z - p[2];
  const y1 = cp * dy + sp * dz, z1 = cp * dz - sp * dy;
  return [p[0] + cy * dx - sy * z1, p[1] + y1, p[2] + sy * dx + cy * z1];
}
/* A pointer names a point on the plane through the chamber centre that faces the viewer. */
function pointerWorld(sx, sy) {
  const [vx, vy] = R.toWorld(sx, sy);
  if (!viewTilted()) return [vx, vy, viewPivot()[2]];
  return fromView(vx, vy, viewPivot()[2]);
}
function viewDelta(dvx, dvy, dvz) { return viewDelta3([dvx, dvy, dvz || 0]); }
function viewDelta3(v) {         // a vector seen on screen, expressed in world axes
  if (!viewTilted()) return [v[0], v[1], v[2] || 0];
  const p = viewPivot();
  const a = fromView(p[0], p[1], p[2]), b = fromView(p[0] + v[0], p[1] + v[1], p[2] + (v[2] || 0));
  return [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
}
function boxCorners() {
  const b = eng.box, out = [];
  for (const z of [b.z0, b.z1]) for (const y of [b.y0, b.y1]) for (const x of [b.x0, b.x1]) out.push(toView(x, y, z));
  return out;
}
function setTilt(yaw, pitch, animate) {
  view.yaw = yaw;
  view.pitch = clamp(pitch, -Math.PI / 2 + 0.02, Math.PI / 2 - 0.02);
  if (animate) { /* the caller animates by stepping this each frame */ }
}
function faceChamber() {
  const from = { yaw: view.yaw, pitch: view.pitch }, t0 = performance.now();
  if (!viewTilted()) return;
  const spin = () => {
    const k = Math.min(1, (performance.now() - t0) / 380), e = 1 - Math.pow(1 - k, 3);
    view.yaw = from.yaw * (1 - e); view.pitch = from.pitch * (1 - e);
    if (k < 1) requestAnimationFrame(spin);
  };
  spin();
}

/* ======================= pointer on the field ======================= */
const pointers = new Map();
let gesture = null, hoverAtom = -1, lastMouse = { x: 0, y: 0, wx: 0, wy: 0 };
const selection = new Set();
const brush = { r: store.get('brushR', 4), active: false, cool: false, x: 0, y: 0 };
function hitAtom(wx, wy, extra = 0) {
  let best = -1, bd = Infinity;
  const k = 0.5;
  for (let i = 0; i < eng.N; i++) {
    const dx = rp[3 * i] - wx, dy = rp[3 * i + 1] - wy, d = dx * dx + dy * dy;
    const r = Math.max(ELEMENTS[eng.type[i]].rvdw * k, 7 / R.scale) + extra;
    if (d < r * r && d < bd) { bd = d; best = i; }
  }
  return best;
}
function boxEdgeAt(sx, sy) {
  if (viewTilted()) return null;                 // the faces are no longer screen-aligned
  const b = eng.box, [x0, y0] = R.toScreen(b.x0, b.y0), [x1, y1] = R.toScreen(b.x1, b.y1), t = 7;
  const nearX = Math.abs(sx - x1) < t && sy > y0 - t && sy < y1 + t, nearY = Math.abs(sy - y1) < t && sx > x0 - t && sx < x1 + t;
  return nearX && nearY ? 'xy' : nearX ? 'x' : nearY ? 'y' : null;
}
function fragmentOf(i) { const fr = eng.fragments(); return fr.list[fr.comp[i]]; }

canvas.addEventListener('contextmenu', e => e.preventDefault());
canvas.addEventListener('dblclick', e => {
  if(armed||placing||tool!=='grab')return;
  const hit=hitAtom(...R.toWorld(e.clientX,e.clientY));
  if(hit>=0){eng.tweezer=null;gesture=null;inspector.open(hit);}
});
canvas.addEventListener('pointerdown', e => {
  try { canvas.setPointerCapture(e.pointerId); } catch (err) { }
  pointers.set(e.pointerId, { x: e.clientX, y: e.clientY });
  closePop(); hideTip();
  if (pointers.size === 2) { // pinch
    const [a, b] = [...pointers.values()];
    gesture = { type: 'pinch', d: Math.hypot(a.x - b.x, a.y - b.y), mx: (a.x + b.x) / 2, my: (a.y + b.y) / 2 };
    return;
  }
  const [wx, wy] = R.toWorld(e.clientX, e.clientY);
  const hit = hitAtom(wx, wy);
  if (e.button === 2 && (placing || armed)) {   // put it down and take the hand back
    const what = placing ? 'Molecule' : BY_SYM[armed].name;
    placing = null; arm(null); setTool('grab');
    toast(what + ' put down · hand tool');
    return;
  }
  if (e.button === 2 && hit >= 0) { inspector.open(hit); return; }
  if (e.button === 1 || (e.button === 0 && spaceHeld)) { gesture = { type: 'pan', sx: e.clientX, sy: e.clientY, cx: R.cam.cx, cy: R.cam.cy }; canvas.className = 'grabbing'; return; }
  if (placing) { gesture = { type: 'placeMol', sx: e.clientX, sy: e.clientY, wx, wy }; return; }
  if (armed && e.button === 0) { gesture = { type: 'placeAtom', sx: e.clientX, sy: e.clientY, wx, wy }; return; }
  if (e.button === 2) {
    if (tool === 'heat') { brush.active = true; brush.cool = true; gesture = { type: 'brush' }; return; }
    // empty scene: turn the chamber in its third dimension for as long as the button is held
    gesture = { type: 'orbit', sx: e.clientX, sy: e.clientY, yaw: view.yaw, pitch: view.pitch };
    canvas.className = 'orbiting'; return;
  }
  if (tool === 'erase') { pushUndo(); gesture = { type: 'erase', alt: e.altKey }; eraseAt(wx, wy, e.altKey); return; }
  if (tool === 'heat') { brush.active = true; brush.cool = e.shiftKey; gesture = { type: 'brush' }; return; }
  const edge = boxEdgeAt(e.clientX, e.clientY);
  if (edge && hit < 0) { pushUndo(); gesture = { type: 'box', edge }; return; }
  if (e.shiftKey) { gesture = { type: 'marquee', x0: wx, y0: wy, x1: wx, y1: wy }; return; }
  if (hit >= 0) {
    if (time.playing) { const w = pointerWorld(e.clientX, e.clientY); gesture = { type: 'tweezer', i: hit }; eng.tweezer = { i: hit, x: w[0], y: w[1], k: 30 }; eng.checkpoints.length = 0; }
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
    if (hoverAtom >= 0 && tool === 'grab') showTipAt(e.clientX, e.clientY, '<b>' + ELEMENTS[eng.type[hoverAtom]].name + '</b><div class="row">Right-click to inspect</div>'); else hideTip();
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
  else if (g.type === 'orbit') {
    view.yaw = g.yaw + (e.clientX - g.sx) * 0.006;
    view.pitch = clamp(g.pitch - (e.clientY - g.sy) * 0.006, -Math.PI / 2 + 0.02, Math.PI / 2 - 0.02);
  }
  else if (g.type === 'tweezer') { if (eng.tweezer) { const w = pointerWorld(e.clientX, e.clientY); eng.tweezer.x = w[0]; eng.tweezer.y = w[1]; eng.needForces = true; } }
  else if (g.type === 'move') {
    const d = viewDelta(wx - g.wx, wy - g.wy); g.wx = wx; g.wy = wy; g.moved = true;
    for (const i of g.set) for (let c = 0; c < 3; c++) { eng.pos[3 * i + c] += d[c]; eng.prev[3 * i + c] = eng.pos[3 * i + c]; }
    eng.needRebuild = true; eng.needForces = true;
  } else if (g.type === 'erase') eraseAt(wx, wy, g.alt);
  else if (g.type === 'marquee') { g.x1 = wx; g.y1 = wy; }
  else if (g.type === 'box') {
    const b = eng.box;
    if (g.edge.includes('x')) b.x1 = clamp(wx, b.x0 + 10, b.x0 + 500);
    if (g.edge.includes('y')) b.y1 = clamp(wy, b.y0 + 10, b.y0 + 500);
    eng.touch();
  } else if (g.type === 'placeAtom' || g.type === 'placeMol') { g.cx = wx; g.cy = wy; }
});
function endPointer(e) {
  pointers.delete(e.pointerId);
  const g = gesture; if (!g) return;
  if (g.type === 'pinch') { if (pointers.size < 2) gesture = null; return; }
  gesture = null;
  const [wx, wy] = R.toWorld(e.clientX, e.clientY);
  if (g.type === 'tweezer') { eng.tweezer = null; eng.touch(); }
  else if (g.type === 'move') { if (g.moved) edited(); else { undoStack.pop(); if (!e.shiftKey) { selection.clear(); } } }
  else if (g.type === 'erase') edited();
  else if (g.type === 'brush') brush.active = false;
  else if (g.type === 'orbit') faceChamber();   // the turn lasts as long as the button is held
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
  for (let i = 0; i < eng.N; i++) {
    const v = viewTilted() ? toView(eng.pos[3 * i], eng.pos[3 * i + 1], eng.pos[3 * i + 2]) : eng.pos.subarray(3 * i, 3 * i + 3);
    const dx = v[0] - brush.x, dy = v[1] - brush.y;
    if (dx * dx + dy * dy < r2) list.push(i);
  }
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
  const fr = eng.fragments(), map = new Map(), frags = [], byId = new Map();
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
function validMol(m) { return ChemProtocol.validMolecule(m, BY_SYM); }
let inbox = readInbox();
function receive(mol, open) {
  if (!validMol(mol)) {
    const bad = mol && Array.isArray(mol.atoms) ? mol.atoms.find(a => a && !Object.hasOwn(BY_SYM, a.el)) : null;
    toast(bad ? bad.el + ' has no bonding parameters in the Playground yet' : 'Invalid molecule: check coordinates, elements, units, and bond indices');
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
  if (!validMol(mol)) { toast('Invalid molecule'); return; }
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
  const e = new Engine({ T, thermostat: true, thermostatMode: 'csvr', tau: 60, seed: 7 });
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
  } else if (c.phase === 'equil') {
    if (c.T === 0) { e.minimize(10, 0.2); }
    else {
      while (performance.now() - t0 < budgetMs) { for (let k = 0; k < 20 && c.steps < c.target; k++) { if (!safeStep(e)) { c.phase = 'error'; break; } c.steps++; } if (c.phase === 'error' || c.steps >= c.target) break; }
      if (c.phase === 'equil' && c.steps >= c.target) c.phase = 'done';
    }
  }
  // status
  const fr = e.fragments(), parts = fr.list.length;
  const st = $('condStatus'), bar = $('condBar');
  st.classList.toggle('warn', parts > 1 && c.phase !== 'relax');
  let text;
  if (c.phase === 'error') { text = 'Preview stopped: extreme collision. Try a lower temperature or larger cell.'; bar.style.width = '0%'; }
  else if (c.phase === 'relax') { text = 'Relaxing geometry in 3D · E ' + e.Epot.toFixed(0) + ' kJ/mol'; bar.style.width = Math.min(100, c.relaxIt / 24) + '%'; }
  else if (c.phase === 'equil') { text = 'Sampling at ' + fmtT(c.T) + ' K' + ' · ' + (c.steps / 1000).toFixed(1) + ' / ' + (c.target / 1000).toFixed(1) + ' ps'; bar.style.width = (c.steps / c.target * 100) + '%'; }
  else {
    bar.style.width = '100%';
    if (parts > 1) text = 'Model fragmented at ' + fmtT(c.T) + ' K → ' + fr.list.map(g => pretty(e.formulaOf(g) + (e.isRadical(g) ? '·' : ''))).slice(0, 5).join(' + ') + ' · placeable anyway';
    else text = c.T === 0 ? 'Ready · geometry relaxation complete, motionless · E ' + e.Epot.toFixed(0) + ' kJ/mol' : 'Ready · intact at ' + e.temperature().toFixed(0) + ' K · E ' + e.Epot.toFixed(0) + ' kJ/mol';
  }
  $('condText').textContent = text;
  $('condPlace').disabled = c.phase !== 'done';
  // camera follows the molecule
  let cx = 0, cy = 0; for (let i = 0; i < e.N; i++) { cx += e.pos[3 * i]; cy += e.pos[3 * i + 1]; } cx /= e.N; cy /= e.N;
  condR.cam.cx += (cx - condR.cam.cx) * 0.2; condR.cam.cy += (cy - condR.cam.cy) * 0.2;
  const n3 = 3 * e.N, bonds = e.bonds();
  condR.draw({ N: e.N, type: e.type, pos: e.pos.subarray(0, n3), bonds, now: performance.now() });
  // container outline
  const ctx = condR.ctx, [sx, sy] = condR.toScreen(e.sphere.x, e.sphere.y), rr = e.sphere.R * condR.scale;
  if (rr < 400) { ctx.strokeStyle = 'rgba(127,167,201,.25)'; ctx.setLineDash([2, 5]); ctx.beginPath(); ctx.arc(sx, sy, rr, 0, Math.PI * 2); ctx.stroke(); ctx.setLineDash([]); }
}
function paintSeg() {
  $('condSeg').querySelectorAll('button').forEach(b => b.setAttribute('aria-checked', String(b.dataset.v === condPrefs.preset)));
  $('customRow').classList.toggle('show', condPrefs.preset === 'custom');
  $('customRow').inert = condPrefs.preset !== 'custom';
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
  if (!cond || cond.phase !== 'done') return;
  const e = cond.e;
  let cx = 0, cy = 0, cz = 0, M = 0, vx = 0, vy = 0, vz = 0;
  for (let i = 0; i < e.N; i++) { const m = e.mass[i]; cx += m * e.pos[3 * i]; cy += m * e.pos[3 * i + 1]; cz += m * e.pos[3 * i + 2]; vx += m * e.vel[3 * i]; vy += m * e.vel[3 * i + 1]; vz += m * e.vel[3 * i + 2]; M += m; }
  cx /= M; cy /= M; cz /= M; vx /= M; vy /= M; vz /= M;
  const atoms = [];
  for (let i = 0; i < e.N; i++) atoms.push({ t: e.type[i], sym: ELEMENTS[e.type[i]].sym, x: e.pos[3 * i] - cx, y: e.pos[3 * i + 1] - cy, z: e.pos[3 * i + 2] - cz, v: $('placeAtRest').checked ? [0, 0, 0] : [e.vel[3 * i] - vx, e.vel[3 * i + 1] - vy, e.vel[3 * i + 2] - vz], q: e.formal[i], V: e.val[i] });
  const bonds = e.bonds(0.35).map(b => ({ a: b.i, b: b.j, order: b.order }));
  placing = { name: cond.mol.name || pretty(cond.mol.formula), atoms, bonds, mass: M, rot: 0, x: lastMouse.wx, y: lastMouse.wy, T: cond.T, atRest: $('placeAtRest').checked };
  closeSheet(); refreshDock();
  toast('Click to place · drag to throw · Q/E rotate · Shift-click places several · Esc stops');
}
function ghostAtoms(x, y) {
  const c = Math.cos(placing.rot), s = Math.sin(placing.rot);
  return placing.atoms.map(a => ({ t: a.t, x: x + a.x * c - a.y * s, y: y + a.x * s + a.y * c, z: a.z }));
}
// a ghost atom as the engine would receive it: the preview is drawn in the frame you see
function ghostWorld(a) { return viewTilted() ? fromView(a.x, a.y, viewPivot()[2] + a.z) : [a.x, a.y, a.z]; }
function ghostOK(g) {
  const b = eng.box;
  for (const a of g) {
    const w = ghostWorld(a);
    if (w[0] < b.x0 || w[0] > b.x1 || w[1] < b.y0 || w[1] > b.y1 || w[2] < b.z0 || w[2] > b.z1) return false;
    if (tooClose(w[0], w[1], w[2], 1.5)) return false;
  }
  return true;
}
function dropMolecule(x, y, fling, keep) {
  const g = ghostAtoms(x, y);
  if (!ghostOK(g)) { toast('No room here — move the ghost to a free spot inside the box'); return; }
  pushUndo();
  const c = Math.cos(placing.rot), s = Math.sin(placing.rot), base = eng.N;
  // Conditioning removes center-of-mass motion. Restore thermal translation
  // once per molecule, otherwise isolated imports vibrate forever in place.
  const thrown = fling ? viewDelta(fling[0], fling[1]) : null;
  const drift = placing.atRest ? [0, 0, 0] : thrown || ChemProtocol.thermalTranslation(placing.mass, placing.T, () => eng.gauss());
  placing.atoms.forEach((a, k) => {
    const w = ghostWorld(g[k]);
    const vv = [a.v[0] * c - a.v[1] * s, a.v[0] * s + a.v[1] * c, a.v[2]];
    const vw = viewTilted() ? viewDelta3(vv) : vv;
    eng.addAtom(a.sym, w[0], w[1], w[2], { v: [vw[0] + drift[0], vw[1] + drift[1], vw[2] + drift[2]], charge: a.q, V: a.V });
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
function directTemperature(T) {
  pushUndo(); eng.thermostat = false; store.set('thermostat', false); setT(T);
  if ($('tPop').classList.contains('open')) openTPop();
  toast('Sample set to ' + fmtT(T) + ' K · wall thermostat off');
}
act('cond.zero', 'Conditions', 'Set sample to 0 K; thermostat off', ['Z'], () => directTemperature(0));
act('cond.room', 'Conditions', 'Room temperature', ['R'], () => { setT(STP_T); toast('Target 298 K'); });
act('cond.flame', 'Conditions', 'Set sample to 3000 K; thermostat off', ['Shift+R'], () => directTemperature(3000));
act('cond.hotter', 'Conditions', 'Hotter ×1.25', ['Shift+ArrowUp'], () => setT(eng.T < 5 ? 25 : eng.T * 1.25));
act('cond.colder', 'Conditions', 'Colder ÷1.25', ['Shift+ArrowDown'], () => setT(eng.T < 5 ? 0 : eng.T / 1.25));
act('cond.typeT', 'Conditions', 'Type a temperature', ['T'], () => tField.startEdit());
act('cond.bath', 'Conditions', 'Thermostat on / off', [], () => { toggleBath(); toast(eng.thermostat ? (bathMode() === 'kelvin' ? 'Kelvin stat on' : 'Wall heater on') : 'Thermostat off'); });
act('cond.kelvin', 'Conditions', 'Hold every atom at this temperature', [], () => { setBathMode('kelvin'); toast('Kelvin stat holding ' + fmtT(eng.T) + ' K'); });
act('cond.resample', 'Conditions', 'Resample thermal motion', ['Shift+T'], resampleMotion);
act('view.in', 'View', 'Zoom in', ['='], () => zoomAt(innerWidth / 2, innerHeight / 2, 1 / 1.35));
act('view.out', 'View', 'Zoom out', ['-'], () => zoomAt(innerWidth / 2, innerHeight / 2, 1.35));
act('view.fit', 'View', 'Fit the box', ['0'], () => fitBox(true));
act('edit.undo', 'Edit', 'Undo', ['Ctrl+Z'], undo);
act('edit.redo', 'Edit', 'Redo', ['Ctrl+Shift+Z', 'Ctrl+Y'], redo);
act('edit.delete', 'Edit', 'Delete selection', ['Delete', 'Backspace'], deleteSelection);
act('edit.all', 'Edit', 'Select all', ['Ctrl+A'], () => { for (let i = 0; i < eng.N; i++) selection.add(i); });
act('edit.pin', 'Edit', 'Pin / unpin selection', ['K'], pinSelection);
act('edit.clear', 'Edit', 'Clear the field', ['Shift+Delete'], () => { if (!eng.N) return; pushUndo(); eng.clear(); selection.clear(); edited(); toast('Field cleared — Ctrl+Z brings it back'); });
act('edit.stopMotion', 'Edit', 'Stop all motion now', [], () => { eng.vel.fill(0); eng.checkpoints.length = 0; toast('All velocities set to zero'); });
act('ui.molecules', 'Panels', 'Molecules tray', ['M'], toggleTray);
act('ui.keys', 'Panels', 'Shortcuts', ['?'], () => openConsole('keys'));
act('ui.console', 'Panels', 'Console', ['Shift+M'], () => toggleConsole());
act('ui.palette', 'Panels', 'Command palette', ['Ctrl+K', '/'], openPalette);
function nextSpeed(d) {
  const list = SPEEDS.map(s => s[0]);
  if (d > 0) return list.find(v => v > time.speed + 1e-9) ?? list[list.length - 1];
  return [...list].reverse().find(v => v < time.speed - 1e-9) ?? list[0];
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
  if (t && t.closest('select')) return; // native picker owns its keyboard navigation
  if (t && t.closest('button, a') && [' ', 'Enter'].includes(e.key)) return;
  if (t && (t.isContentEditable || t.tagName === 'INPUT' || t.tagName === 'TEXTAREA')) return;
  if (e.key === ' ') spaceHeld = true;
  if (e.key === 'Escape') {
    if (consoleIsOpen()) { closeConsole(); return; }
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

/* The keybind grid is an app inside the console; it paints into whatever stage is mounted. */
function paintKeys() {
  const grid = $('keysGrid'); if (!grid) return;
  const groups = [...new Set(ACTIONS.map(a => a.group))];
  grid.innerHTML = groups.map(g => '<div class="kgroup"><h3>' + g + '</h3>' + ACTIONS.filter(a => a.group === g).map(a =>
    '<div class="krow"><span>' + esc(a.name) + '</span><button class="kcap' + ((keymap[a.id] || []).length ? '' : ' empty') + '" data-id="' + a.id + '">' + ((keymap[a.id] || []).length ? esc(prettyKey(keymap[a.id][0])) : 'add') + '</button></div>').join('') + '</div>').join('') +
    '<div class="kgroup"><h3>Mouse</h3>' + [['Scroll', 'Zoom (0.05–15 nm; Fit can show larger boxes)'], ['Drag empty space', 'Pan'], ['Right-drag empty space', 'Turn the chamber in 3D while held; it faces you again on release'], ['Drag an atom', 'Pull it (running) / move its molecule (paused)'], ['Alt + drag atom', 'Move one atom (paused)'], ['Shift + drag', 'Select a region'], ['Right-click atom', 'Atom inspector with its calculated electron cloud'], ['Right-click while placing', 'Put the atom or molecule down and take the hand'], ['Drag while placing', 'Throw with that velocity'], ['Drag box edge', 'Resize the chamber'], ['Q / E', 'Rotate a molecule before placing']].map(r => '<div class="krow"><span>' + r[1] + '</span><span class="kcap" style="border:0;background:none;color:var(--faint)">' + r[0] + '</span></div>').join('') + '</div>';
  grid.querySelectorAll('button.kcap').forEach(b => b.onclick = () => {
    document.querySelectorAll('.kcap.listening').forEach(x => x.classList.remove('listening'));
    listening = b; b.classList.add('listening'); b.textContent = 'press a key';
  });
}
function captureKey(e) {
  const b = listening; listening = null; b.classList.remove('listening');
  const id = b.dataset.id;
  if (e.key === 'Escape') { paintKeys(); return; }
  if (e.key === 'Backspace' && !keymap[id]?.length) { paintKeys(); return; }
  const combo = comboOf(e); if (!combo) { paintKeys(); return; }
  let moved = null;
  for (const a of ACTIONS) if (a.id !== id && (keymap[a.id] || []).includes(combo)) { keymap[a.id] = keymap[a.id].filter(k => k !== combo); moved = a.name; }
  keymap[id] = [combo];
  store.set('keys', keymap);
  paintKeys(); refreshKeyHints();
  const nb = $('keysGrid') && $('keysGrid').querySelector('[data-id="' + id + '"]'); nb && nb.classList.add('flash');
  if (moved) toast(prettyKey(combo) + ' moved here from “' + moved + '”');
}


/* ======================= console ======================= */
/* Apps live on the scene rather than over it: opening one holds the clock where it stands,
   the field keeps drawing behind, and closing hands time back exactly where it was left.
   Every control here writes straight through to the engine, so a change shows in the chamber
   while it is still being made. */
const ICONS = {
  environment: '<svg viewBox="0 0 16 16"><rect x="2" y="2.6" width="12" height="10.8" rx="2"/><path d="M4.6 5.4v5.2M11.4 5.4v5.2"/><circle cx="8" cy="8" r="1.5"/></svg>',
  display: '<svg viewBox="0 0 16 16"><rect x="2" y="2" width="12" height="9" rx="2"/><path d="M5 14h6M8 11v3"/></svg>',
  about: '<svg viewBox="0 0 16 16"><circle cx="8" cy="8" r="6.2"/><path d="M8 7.3v3.9"/><path d="M8 4.9h.01"/></svg>',
  keys: '<svg viewBox="0 0 16 16"><rect x="1.6" y="4.2" width="12.8" height="8" rx="1.7"/><path d="M4.2 6.9h.01M6.6 6.9h.01M9 6.9h.01M11.4 6.9h.01M5.2 9.7h5.6"/></svg>'
};
const APPS = [
  { id: 'environment', name: 'Conditions', render: renderEnvironmentApp },
  { id: 'display', name: 'Display', render: renderDisplayApp },
  { id: 'about', name: 'Model', render: renderAboutApp },
  { id: 'keys', name: 'Keys', render: renderKeysApp }
];
let consoleApp = store.get('consoleApp', 'environment');
let consoleResume = false, consoleTick = null, consoleCloseTimer = null;
const consoleIsOpen = () => !$('console').hidden;

function openConsole(id) {
  if (id) consoleApp = id;
  if (consoleIsOpen() && !$('console').classList.contains('closing')) { selectApp(consoleApp); return; }
  clearTimeout(consoleCloseTimer);
  closePop(); hideTip();
  consoleResume = time.playing;
  setPlaying(false);
  const win = $('consoleWin');
  win.style.removeProperty('--panel-x'); win.style.removeProperty('--panel-y');
  $('console').hidden = false;
  $('console').classList.remove('closing');
  $('menuBtn').setAttribute('aria-expanded', 'true');
  buildRail();
  selectApp(consoleApp);
  consoleTick = setInterval(() => { $('consoleClock').textContent = fmtTime(eng.time); refreshApp(); }, 250);
  $('consoleClock').textContent = fmtTime(eng.time);
  $('consoleClose').focus();
}
function closeConsole() {
  if (!consoleIsOpen()) return;
  closeDropdowns();
  clearInterval(consoleTick); consoleTick = null;
  clearTimeout(consoleCloseTimer);
  const veil = $('console');
  veil.classList.add('closing');
  $('menuBtn').setAttribute('aria-expanded', 'false');
  consoleCloseTimer = setTimeout(() => { veil.hidden = true; veil.classList.remove('closing'); }, document.documentElement.dataset.motion === 'reduced' ? 0 : 180);
  if (consoleResume) setPlaying(true);
  consoleResume = false;
  $('menuBtn').focus();
}
function toggleConsole(id) {
  if (consoleIsOpen() && !$('console').classList.contains('closing') && (!id || id === consoleApp)) closeConsole();
  else openConsole(id);
}
function buildRail() {
  $('consoleRail').innerHTML = APPS.map(a =>
    '<button class="app-tab" role="tab" data-app="' + a.id + '" aria-selected="false">' + ICONS[a.id] + '<span>' + a.name + '</span></button>').join('');
  $('consoleRail').querySelectorAll('.app-tab').forEach(b => b.onclick = () => selectApp(b.dataset.app));
}
function selectApp(id) {
  closeDropdowns();
  const app = APPS.find(a => a.id === id) || APPS[0];
  consoleApp = app.id; store.set('consoleApp', app.id);
  $('consoleTitle').textContent = app.name;
  $('consoleRail').querySelectorAll('.app-tab').forEach(b => b.setAttribute('aria-selected', b.dataset.app === app.id));
  const stage = $('consoleStage');
  stage.scrollTop = 0; stage.innerHTML = '';
  app.render(stage);
}
function refreshApp() { if (appRefresh) appRefresh(); }
let appRefresh = null;

$('menuBtn').onclick = () => toggleConsole();
$('consoleClose').onclick = closeConsole;
document.addEventListener('keydown',e=>{
  if(e.key!=='Tab')return;
  const panel=consoleIsOpen()?$('consoleWin'):inspector.panel.classList.contains('open')?inspector.panel:null;
  if(!panel)return;
  const nodes=[...panel.querySelectorAll('button,select,input,summary,a[href],[tabindex="0"]')].filter(x=>!x.disabled&&x.getClientRects().length);
  if(!nodes.length)return;
  const first=nodes[0],last=nodes[nodes.length-1];
  if(e.shiftKey&&(document.activeElement===first||!panel.contains(document.activeElement))){e.preventDefault();last.focus();}
  else if(!e.shiftKey&&(document.activeElement===last||!panel.contains(document.activeElement))){e.preventDefault();first.focus();}
});
// Keep the panel attached to its default anchor when the viewport changes.
window.addEventListener('resize', () => {
  $('consoleWin').style.removeProperty('--panel-x');
  $('consoleWin').style.removeProperty('--panel-y');
});
// Drag the panel within the visible viewport.
$('consoleBar').addEventListener('pointerdown', e => {
  if (e.target.closest('button')) return;
  const win = $('consoleWin'), r = win.getBoundingClientRect();
  const ox = e.clientX - r.left, oy = e.clientY - r.top;
  const bar = $('consoleBar');
  bar.setPointerCapture(e.pointerId);
  const move = ev => {
    const x = clamp(ev.clientX - ox, 12, Math.max(12,innerWidth-r.width-12));
    const y = clamp(ev.clientY - oy, 12, Math.max(12,innerHeight-r.height-12));
    win.style.setProperty('--panel-x', x + 'px');
    win.style.setProperty('--panel-y', y + 'px');
  };
  const up = () => { bar.removeEventListener('pointermove', move); bar.removeEventListener('pointerup', up); };
  bar.addEventListener('pointermove', move); bar.addEventListener('pointerup', up);
});

/* ---- small builders shared by the apps ---- */
function el(tag, cls, html) {
  const n = document.createElement(tag);
  if (cls) n.className = cls;
  if (html != null) n.innerHTML = html;
  return n;
}
const TICK = '<svg viewBox="0 0 12 12"><path d="M2.4 6.2l2.4 2.4 4.8-5"/></svg>';
/* One dropdown component, two behaviours: `multi` turns the marks into checkboxes and keeps
   the sheet open, because those options are not alternatives to each other. */
function closeDropdowns(){document.querySelectorAll('.drop.open .drop-btn').forEach(b=>b.click());}
function dropdown({ options, multi, get, set, summary, label = 'Boundary type' }) {
  const wrap = el('div', 'drop');
  const btn = el('button', 'drop-btn',
    '<i class="drop-lead"></i><span class="drop-val"></span>' +
    '<svg class="caret" viewBox="0 0 10 10"><path d="M2 3.8l3 3 3-3"/></svg>');
  btn.setAttribute('aria-haspopup', multi ? 'true' : 'listbox');
  btn.setAttribute('aria-expanded', 'false');
  wrap.appendChild(btn);
  let sheet = null;
  const paint = () => {
    const v = get(), text = summary(v), val = btn.querySelector('.drop-val');
    val.textContent = text.label;
    val.classList.toggle('set', !!text.accent);
    btn.querySelector('.drop-lead').innerHTML = text.icon || '';
    if (sheet) sheet.querySelectorAll('.drop-opt').forEach(o => {
      o.setAttribute(multi ? 'aria-checked' : 'aria-selected', multi ? v.includes(o.dataset.v) : v === o.dataset.v);
    });
  };
  // the sheet lives outside `wrap`, so both have to count as inside
  const outside = e => { if (!wrap.contains(e.target) && !(sheet && sheet.contains(e.target))) close(); };
  const close = (focus = false) => {
    if (!sheet) return;
    sheet.remove(); sheet = null; wrap.classList.remove('open');
    btn.setAttribute('aria-expanded', 'false');
    document.removeEventListener('pointerdown', outside);
    const stage = $('consoleStage');
    if (stage) stage.removeEventListener('scroll', close);
    removeEventListener('resize', close);
    if (focus === true) btn.focus();
  };
  /* The sheet floats above everything instead of living inside the scrolling panel: an
     absolutely positioned child is clipped by that panel, so options past its bottom edge
     could not be clicked at all. It flips upward when there is no room below. */
  const placeSheet = () => {
    const r = btn.getBoundingClientRect();
    sheet.style.width = r.width + 'px';
    sheet.style.maxHeight = Math.max(100,Math.max(r.top-15,innerHeight-r.bottom-15))+'px';
    sheet.style.left = r.left + 'px';
    const h = sheet.offsetHeight, below = innerHeight - r.bottom - 10;
    sheet.style.top = (below >= h || r.top < h + 10 ? r.bottom + 5 : r.top - h - 5) + 'px';
  };
  const open = () => {
    if (sheet) return;
    document.querySelectorAll('.drop.open .drop-btn').forEach(b => b !== btn && b.click());
    sheet = el('div', 'drop-sheet');
    sheet.setAttribute('role', multi ? 'group' : 'listbox');
    sheet.setAttribute('aria-label', label);
    sheet.innerHTML = options.map(o =>
      '<button class="drop-opt" role="' + (multi ? 'checkbox' : 'option') + '" data-v="' + o.v + '" ' +
      (multi ? 'aria-checked="false"' : 'aria-selected="false"') + '>' +
      '<span class="drop-mark' + (multi ? ' box' : '') + '">' + TICK + '</span>' +
      '<i class="opt-icon">' + (o.icon || '') + '</i><strong>' + esc(o.name) + '</strong></button>').join('');
    sheet.querySelectorAll('.drop-opt').forEach(o => o.onclick = ev => {
      ev.stopPropagation(); set(o.dataset.v); paint();
      if (!multi) close(true);
    });
    sheet.addEventListener('keydown', navigate);
    sheet.addEventListener('focusout', leave);
    document.body.appendChild(sheet); wrap.classList.add('open');
    placeSheet();
    btn.setAttribute('aria-expanded', 'true'); paint();
    document.addEventListener('pointerdown', outside);
    const stage = $('consoleStage');
    if (stage) stage.addEventListener('scroll', close);
    addEventListener('resize', close);
    (sheet.querySelector('[aria-selected="true"],[aria-checked="true"]') || sheet.firstElementChild).focus();
  };
  btn.onclick = e => { e.stopPropagation(); sheet ? close() : open(); };
  function navigate(e) {
    if (e.key === 'Escape' && sheet) { e.preventDefault(); e.stopPropagation(); close(true); }
    else if (['ArrowDown','ArrowUp','Home','End'].includes(e.key)) {
      e.preventDefault(); e.stopPropagation();
      if (!sheet) { open(); return; }
      const items = [...sheet.querySelectorAll('.drop-opt')], i = items.indexOf(document.activeElement);
      const next = e.key === 'Home' ? 0 : e.key === 'End' ? items.length - 1 : (i + (e.key === 'ArrowDown' ? 1 : -1) + items.length) % items.length;
      items[next].focus();
    } else if (e.key === 'Tab' && sheet) close(true);
  }
  const leave=()=>queueMicrotask(()=>{if(!wrap.contains(document.activeElement) && !(sheet && sheet.contains(document.activeElement)))close();});
  wrap.addEventListener('keydown',navigate);
  wrap.addEventListener('focusout',leave);
  paint(); wrap.repaint = paint;
  return wrap;
}

/* ---- environment ---- */
function boundarySettingsChanged(note) {
  eng.checkpoints.length = 0; eng.needForces = true; eng.touch();
  store.set('bounds', { version: 3, mode: eng.boundsMode, vt: eng.voidTemperature, vp: eng.voidPressure, vv: eng.voidVelocity,
    pc: eng.pressureControl, pt: eng.pressureTarget });
  scheduleSave();
  if (note) toast(note);
}
/* Glyphs carry the meaning here; words only name a thing once. */
const GLYPH = {
  // a hard bar turns the path in a sharp V; a field is arcs, and the path curves inside them
  solid: '<svg class="gl" viewBox="0 0 24 24"><path class="face" d="M18 2.5v19"/><path class="path" d="M4 19.5L16.5 12L4 4.5"/><circle class="ball" cx="16.5" cy="12" r="2"/></svg>',
  forcefield: '<svg class="gl" viewBox="0 0 24 24"><path class="halo" d="M19 3.5a13 13 0 0 1 0 17"/><path class="face soft" d="M15 5.5a10 10 0 0 1 0 13"/><path class="path" d="M3 19.5C9 17 12.5 14.5 13.5 12C12.5 9.5 9 7 3 4.5"/><circle class="ball" cx="13.5" cy="12" r="2"/></svg>',
  heat: '<svg class="gl" viewBox="0 0 24 24"><path class="stem" d="M9.4 13V5a2.6 2.6 0 0 1 5.2 0v8a4.4 4.4 0 1 1-5.2 0z"/><circle class="ball" cx="12" cy="16.3" r="2.6"/></svg>',
  press: '<svg class="gl" viewBox="0 0 24 24"><path class="face" d="M4 17a8 8 0 1 1 16 0"/><path class="path" d="M12 17l4-5"/></svg>',
  vel: '<svg class="gl" viewBox="0 0 24 24"><path class="face" d="M19 3v18"/><path class="path" d="M3 12h11M11 8.5l3.5 3.5L11 15.5"/><circle class="ball" cx="16.5" cy="12" r="1.8"/></svg>',
  setP: '<svg class="gl" viewBox="0 0 24 24"><path class="face" d="M4 17a8 8 0 1 1 16 0"/><path class="path" d="M12 17l4-5M7 6.5l1.2 1.6M17 6.5l-1.2 1.6"/></svg>',
  none: '<svg class="gl" viewBox="0 0 24 24"><circle class="face" cx="12" cy="12" r="8" stroke-dasharray="3 3"/></svg>',
  w: '<svg class="gl" viewBox="0 0 24 24"><path class="path" d="M3 12h18M6 9l-3 3 3 3M18 9l3 3-3 3"/></svg>',
  h: '<svg class="gl" viewBox="0 0 24 24"><path class="path" d="M12 3v18M9 6l3-3 3 3M9 18l3 3 3-3"/></svg>',
  d: '<svg class="gl" viewBox="0 0 24 24"><path class="path" d="M6 18L18 6M18 11V6h-5M6 13v5h5"/></svg>',
  wall: '<svg class="gl" viewBox="0 0 24 24"><path class="face" d="M6 3v18"/><path class="path" d="M10 7c2 1.4 2-1.4 4 0M10 12c2 1.4 2-1.4 4 0M10 17c2 1.4 2-1.4 4 0"/></svg>',
  off: '<svg class="gl" viewBox="0 0 24 24"><path class="path" d="M12 4v7"/><path class="face" d="M7.4 7.4a7 7 0 1 0 9.2 0"/></svg>',
  kelvin: '<svg class="gl" viewBox="0 0 24 24"><path class="face" d="M4 3v18M20 3v18"/><circle class="ball" cx="9" cy="8" r="1.6"/><circle class="ball" cx="15" cy="12.5" r="1.6"/><circle class="ball" cx="10" cy="17" r="1.6"/></svg>',
  drain: '<svg class="gl" viewBox="0 0 24 24"><path class="path" d="M4 12h12M12 8l4 4-4 4"/><path class="face" d="M20 5v14" stroke-dasharray="3 3"/></svg>',
  gauge: '<svg class="gl" viewBox="0 0 24 24"><path class="face" d="M4 17a8 8 0 1 1 16 0"/><path class="path" d="M12 17l5-4"/></svg>'
};
function renderEnvironmentApp(stage) {
  const group = (title, node) => { const g = el('div', 'app-group', '<h3>' + title + '</h3>'); g.appendChild(node); stage.appendChild(g); return g; };

  group('Boundary', dropdown({
    options: [
      { v: 'solid', name: 'Solid', icon: GLYPH.solid },
      { v: 'forcefield', name: 'Forcefield', icon: GLYPH.forcefield }
    ],
    get: () => eng.boundsMode,
    set: v => { eng.boundsMode = v; boundarySettingsChanged(); },
    summary: v => ({ label: v === 'solid' ? 'Solid walls' : 'Soft forcefield', icon: GLYPH[v] })
  }));

  group('Void walls', dropdown({
    multi: true, label: 'Void wall channels',
    options: [
      { v: 'temperature', name: 'Temperature', icon: GLYPH.heat },
      { v: 'pressure', name: 'Pressure', icon: GLYPH.press },
      { v: 'velocity', name: 'Velocity', icon: GLYPH.vel }
    ],
    get: () => [eng.voidTemperature && 'temperature', eng.voidPressure && 'pressure', eng.voidVelocity && 'velocity'].filter(Boolean),
    set: v => {
      if (v === 'temperature') eng.voidTemperature = !eng.voidTemperature;
      else if (v === 'pressure') eng.voidPressure = !eng.voidPressure;
      else eng.voidVelocity = !eng.voidVelocity;
      boundarySettingsChanged();
    },
    summary: v => ({
      label: v.length ? v.map(x=>x[0].toUpperCase()+x.slice(1)).join(' · ') : 'None',
      accent: v.length > 0,
      icon: v.length ? GLYPH[v[0] === 'temperature' ? 'heat' : v[0] === 'pressure' ? 'press' : 'vel'] : GLYPH.none
    })
  }));

  const dims = el('div', 'dim-row',
    ['envW', 'envH', 'envD'].map((id, k) =>
      '<label class="dim" aria-label="' + ['Width','Height','Depth'][k] + '"><i>' + GLYPH[['w', 'h', 'd'][k]] + '</i><span class="scrub" id="' + id + '" data-unit="nm"></span></label>').join(''));
  group('Chamber', dims);

  const baro = el('div', 'baro-row',
    '<span class="seg-mini" id="envBaro"><button data-v="off" aria-label="Passive pressure" title="Passive — the gauge only reports">Free</button>' +
    '<button data-v="on" aria-label="Hold a pressure" title="Hold this pressure — the chamber breathes toward it">Hold</button></span>' +
    '<label class="dim" id="baroTargetWrap" title="Target pressure"><span class="scrub" id="baroTarget" title="Target pressure in bar" data-unit="bar"></span></label>');
  group('Pressure', baro);

  const bath = el('span', 'seg-mini',
    '<button data-v="wall" aria-label="Wall heater" title="Wall heater — warms the boundary; the interior follows">' + GLYPH.wall + '<span>Wall</span></button>' +
    '<button data-v="kelvin" aria-label="Kelvin stat" title="Kelvin stat — rescales total kinetic energy to the target each step">' + GLYPH.kelvin + '<span>Kelvin</span></button>' +
    '<button data-v="off" aria-label="Thermostat off" title="Off — setting a temperature is a one-off edit">' + GLYPH.off + '<span>Off</span></button>');
  bath.id = 'envBath';
  group('Thermostat', bath);

  const stat = el('div', 'stat-strip compact-stat', '<span>Energy removed</span>' +
    [['void', GLYPH.drain]]
      .map(([k, g]) => '<div class="stat" data-k="' + k + '">' + g + '<b>—</b></div>').join(''));
  stage.appendChild(stat);

  const mk = (id, get, set, min, max) => { $(id).title=({envW:'Chamber width',envH:'Chamber height',envD:'Chamber depth'})[id]+' in nanometres'; return scrub($(id), {
    get, set: (v, commit) => { set(v); eng.touch(); if (commit) { saveBox(); scheduleSave(); } },
    min, max, unit: 'nm', fmt: v => v.toFixed(2), hardMin: min, hardMax: max, wheel: false
  }); };
  mk('envW', () => (eng.box.x1 - eng.box.x0) / 10, v => { eng.box.x1 = eng.box.x0 + v * 10; }, 1, 50);
  mk('envH', () => (eng.box.y1 - eng.box.y0) / 10, v => { eng.box.y1 = eng.box.y0 + v * 10; }, 1, 50);
  mk('envD', () => (eng.box.z1 - eng.box.z0) / 10, v => { eng.box.z0 = -v * 5; eng.box.z1 = v * 5; }, 0.1, 50);
  bath.querySelectorAll('button').forEach(b => b.onclick = () => { setBathMode(b.dataset.v); paintBath(); });
  scrub($('baroTarget'), {
    get: () => eng.pressureTarget, set: (v, commit) => { eng.pressureTarget = v; if (commit) boundarySettingsChanged(); },
    min: 0, max: 200, off: 1, map: 'log', unit: 'bar', fmt: v => v < 10 ? v.toFixed(2) : v.toFixed(0),
    hardMin: 0, hardMax: 100000, wheel: false
  });
  $('envBaro').querySelectorAll('button').forEach(b => b.onclick = () => {
    eng.pressureControl = b.dataset.v === 'on';
    boundarySettingsChanged(); paintBaro();
  });
  function paintBaro() {
    $('envBaro').querySelectorAll('button').forEach(b => b.setAttribute('aria-pressed', (b.dataset.v === 'on') === !!eng.pressureControl));
    $('baroTargetWrap').classList.toggle('inactive', !eng.pressureControl);
    $('baroTargetWrap').inert = !eng.pressureControl;
  }
  function paintBath() {
    bath.querySelectorAll('button').forEach(b => b.setAttribute('aria-pressed', b.dataset.v === bathMode()));
  }
  paintBath(); paintBaro();
  const cells = stat.querySelectorAll('.stat');
  appRefresh = () => {
    if (!stage.isConnected) return;
    cells[0].lastElementChild.textContent = (eng.voidTemperature || eng.voidPressure || eng.voidVelocity) ? eng.voidHeat.toFixed(1) + ' kJ/mol' : '—';
    cells[0].title = 'Energy removed by void walls';
    paintBaro();
    paintBath();
  };
  appRefresh();
}

function renderDisplayApp(stage){
  appRefresh=null;

  const select=(title,note,key,options)=>{
    const row=el('label','preference-row','<span>'+title+'</span>');
    row.title=note;
    const input=document.createElement('select');input.setAttribute('aria-label',title);
    for(const [value,label] of options){const option=document.createElement('option');option.value=value;option.textContent=label;input.append(option);}
    input.value=appearance[key];input.onchange=()=>{appearance[key]=input.value;applyAppearance();};row.append(input);stage.append(row);
  };
  select('Contour detail','Higher detail uses more rendering time.','quality',[['fine','Fine'],['balanced','Balanced'],['fast','Fast']]);
  select('Interface motion','Follows your system by default.','motion',[['system','System'],['full','Fluid'],['reduced','Reduced']]);
  select('Compute priority','Reserve more time for interaction or simulation.','budget',[['responsive','Interaction'],['balanced','Balanced'],['throughput','Simulation']]);
  for(const [key,title,note] of [['grid','Reference grid','Keep a subtle spatial reference.'],['labels','Atom labels','Show symbols alongside atomic centres.']]){
    const row=el('label','preference-row','<span>'+title+'</span>');
    row.title=note;
    const input=document.createElement('input');input.type='checkbox';input.checked=appearance[key];input.setAttribute('aria-label',title);input.onchange=()=>{appearance[key]=input.checked;applyAppearance();};row.append(input);stage.append(row);
  }
  const reset=el('button','ghost preference-reset','Reset display');reset.onclick=()=>{Object.assign(appearance,{quality:'balanced',motion:'system',grid:true,labels:true,budget:'balanced'});applyAppearance();selectApp('display');};stage.append(reset);
}

/* ---- about ---- */
function renderAboutApp(stage) {
  appRefresh = null;
  stage.appendChild(el('p', 'app-lede',
    'An experimental reactive force field with fitted molecular examples. The structures and energies below are checked; rates and mechanisms are qualitative.'));
  const facts = [
    ['Physical clock', 'Each step advances 1 fs. At 1×, the target is 20,000 steps/s: <b>20 ps per real second</b>. At 0.1×, one step takes one second. × is a playback setting, never a change to the physics.'],
    ['Why a mixture may look inert', 'A second at the target rate covers 20 picoseconds. A reaction may need activation, a solvent, or chemistry this model cannot represent. For stationary imports, use <b>Resample thermal motion</b>.'],
    ['What is checked', 'H₂ 0.741 Å / 436 kJ/mol, water 104.2°, the water dimer at −26 kJ/mol, 2 H₂ + O₂ → 2 H₂O at −486 (lit. −484), H + H₂ barrier 43 (lit. 40), ethene + Cl· barrierless. Na + Cl₂ → NaCl runs at 300 K; CH₄ + O₂ at 3500 K passes through CH₂O and OH·.'],
    ['Known to be wrong', 'Spin is absent (O₂ is patched to act as the triplet diradical it is). Cl + H₂ and OH + H₂ have barriers that are too high, H + O₂ → OH + O too low. CO gets a double bond instead of a triple. Water needs below ~150 K to freeze, because lone pairs have no direction. No tunnelling, excited states or solvent.'],
    ['Three dimensions', 'A 3D chamber seen face on. Solid walls reflect atomic centres at all six faces; the optional soft field allows penetration. Electron contours can extend past a face. Right-drag empty space to tilt the view; zoom changes the view, never atom sizes or the chamber volume.'],
    ['Temperature', 'Wall mode exchanges heat at the boundary. Kelvin mode rescales total kinetic energy each step: it suppresses fluctuations and can alter reaction dynamics. Off leaves motion unthermostatted.'],
    ['Energy &amp; pressure', 'Force checks hold internal bond orders fixed; changing those heuristic orders can cause drift, and extreme-collision speed clamps remove energy. Wall pressure is averaged normal stress. Optional pressure control adjusts chamber width and height; it is a heuristic controller, not validated NPT sampling.'],
    ['What you see', 'Field contours show tabulated atomic sizes through Gaussian surfaces. The atom inspector runs a real Hartree–Fock calculation; the chamber view does not.'],
    ['Zero kelvin', 'A classical geometry minimum with zero initial velocities. Quantum zero-point motion is not represented.']
  ];
  stage.appendChild(el('dl', 'app-facts', facts.map(f => '<div><dt>' + f[0] + '</dt><dd>' + f[1] + '</dd></div>').join('')));
  stage.appendChild(el('p', 'app-source',
    'Methods: <a href="https://manual.gromacs.org/current/reference-manual/algorithms/molecular-dynamics.html" target="_blank" rel="noopener">MD integration ↗</a> · ' +
    '<a href="https://docs.lammps.org/pair_reaxff.html" target="_blank" rel="noopener">Validated reactive force fields ↗</a>. This engine is not ReaxFF.'));
}

/* ---- keybinds ---- */
function renderKeysApp(stage) {
  appRefresh = null;
  const head = el('div', 'keys-hint', '<span>Select a key to rebind.</span>');
  const reset = el('button', 'ghost', 'Reset all');
  reset.onclick = () => { keymap = Object.assign({}, defaultKeys); store.set('keys', keymap); refreshKeyHints(); selectApp('keys'); };
  head.appendChild(reset);
  stage.appendChild(head);
  const grid = el('div', 'keys-grid');
  grid.id = 'keysGrid';
  stage.appendChild(grid);
  paintKeys();
}

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
    // the chamber travels with the scene, so restored atoms are never left outside their walls
    const box = s.box;
    if (box && ['x0', 'x1', 'y0', 'y1', 'z0', 'z1'].every(k => Number.isFinite(box[k])) && box.x1 > box.x0 && box.y1 > box.y0 && box.z1 > box.z0) {
      eng.box = { ...box }; saveBox();
    }
    for (const a of s.atoms) if (BY_SYM[a[0]]) eng.addAtom(a[0], a[1], a[2], a[3], { v: [a[4], a[5], a[6]], charge: a[7], V: a[8] });
    eng.time = s.time || 0;
    for (const key of ['wallT', 'wallTarget', 'wallTau', 'heatToSample', 'heaterWork', 'kelvinWork']) if (Number.isFinite(s[key])) eng[key] = s[key];
    if (s.boundsMode === 'forcefield' || s.boundsMode === 'solid') eng.boundsMode = s.boundsMode;
    if (typeof s.voidTemperature === 'boolean') eng.voidTemperature = s.voidTemperature;
    eng.voidPressure = s.dampingVersion >= 3 && !!s.voidPressure;
    eng.voidVelocity = s.dampingVersion >= 3 && !!s.voidVelocity;
    if (typeof s.pressureControl === 'boolean') eng.pressureControl = s.pressureControl;
    if (Number.isFinite(s.pressureTarget)) eng.pressureTarget = s.pressureTarget;
    for(const key of ['voidTau','voidSkin']) if(Number.isFinite(s[key]) && s[key]>0) eng[key]=s[key];
    if(Number.isFinite(s.voidHeat) && s.voidHeat>=0) eng.voidHeat=s.voidHeat;
  } catch (e) { }
}
window.addEventListener('beforeunload', saveScene);
setInterval(() => { if (time.playing) saveScene(); }, 5000);

/* The inspector card rides with its atom: it is repositioned every frame, and the field draws
   the leader that says which atom it belongs to. */
function inspectLeader() {
  const i = inspector.i;
  if (i < 0 || i >= eng.N) return null;
  const [sx, sy] = R.toScreen(rp[3 * i], rp[3 * i + 1]);
  inspector.place(sx, sy);
  return inspector.anchor ? { x: rp[3 * i], y: rp[3 * i + 1], ax: inspector.anchor[0], ay: inspector.anchor[1], r: ELEMENTS[eng.type[i]].rvdw } : null;
}

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
  const frameBudget = clamp(time.frameEMA * ({responsive:.3,balanced:.55,throughput:.75}[appearance.budget]||.55), 3, 24) - (cond ? 3 : 0);
  if (time.playing) {
    time.acc += dt * stepsPerSecond(time.speed);
    const t0 = performance.now(); let n = 0; time.limited = false;
    while (time.acc >= 1) {
      if (!safeStep(eng)) { time.acc = 0; break; } time.acc -= 1; n++;
      if ((n & 7) === 0 && performance.now() - t0 > frameBudget) { if (time.acc >= 1) { time.limited = true; time.acc %= 1; } break; }
    }
    time.stepsWindow += n;
    applyBrush(dt);
  } else if (brush.active) { applyBrush(dt); }
  if (now - time.windowStart > 500) { const r = time.stepsWindow / ((now - time.windowStart) / 1000); time.rateEMA = time.rateEMA * 0.4 + r * 0.6; time.stepsWindow = 0; time.windowStart = now; }
  if (cond) tickConditioning(4);
  if (eng.tweezer) eng.checkpoints.length = 0;
  // interpolation alpha between the previous and the current step
  let alpha = 1;
  if (time.playing) alpha = time.acc;
  else if (time.manual) { const k = Math.min(1, (now - time.manual.t0) / time.manual.dur); alpha = 1 - Math.pow(1 - k, 3); if (k >= 1) time.manual = null; }
  const n3 = 3 * eng.N;
  if (rp.length < n3 || rp.length > n3 + 3000) rp = new Float64Array(n3 + 300);
  const P = eng.pos, Q = eng.prev;
  for (let k = 0; k < n3; k++) rp[k] = Q[k] + (P[k] - Q[k]) * alpha;
  if (viewTilted()) for (let i = 0; i < eng.N; i++) {
    const v = toView(rp[3 * i], rp[3 * i + 1], rp[3 * i + 2]);
    rp[3 * i] = v[0]; rp[3 * i + 1] = v[1]; rp[3 * i + 2] = v[2];
  }
  if (eng.needForces) eng.refresh(); // edits while paused
  const bonds = eng.bonds();
  updateBondEvents(bonds, now);
  if (time.playing || species.reset) updateSpecies(now);
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
    N: eng.N, type: eng.type, pos: rp, bonds, box: eng.box, box3: viewTilted() ? boxCorners() : null,
    bounds: { mode: eng.boundsMode, range: eng.fieldRange, voidT: eng.voidTemperature, voidP: eng.voidPressure },
    boxHot: gesture && gesture.type === 'box' ? gesture.edge : boxHot,
    hover: gesture ? -1 : hoverAtom, eraseHover: tool === 'erase' && !armed, selected: selection, pinned: eng.pinned.subarray(0, eng.N), ghost, flashes, now,
    tweezer: eng.tweezer, brush: tool === 'heat' && !armed && !placing ? brush : null,
    marquee: gesture && gesture.type === 'marquee' ? gesture : null,
    inspect: inspectLeader()
  });
  // UI text, ~10 Hz
  if (now - lastUI > 100) {
    lastUI = now;
    updateLab();
    $('clock').textContent = fmtTime(eng.time);
    const rate = $('rate');
    if (time.playing) {
      rate.textContent = (time.limited ? 'CPU-bound · ' : '') + fmtRate(time.rateEMA) + (time.limited ? ' · ' + Math.round(time.rateEMA / stepsPerSecond(time.speed) * 100) + '% target' : '');
      rate.classList.toggle('limited', time.limited);
    } else { rate.textContent = eng.canStepBack() ? 'paused · ' + fmtTime(Math.min(eng.historySpan(), eng.time)) + ' rewindable' : 'paused'; rate.classList.remove('limited'); }
    readTemperature();
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
paintViews();
loadScene();
eng.refresh();
fitBox(false);
buildDock();
renderTray();
paintSeg();
readHash();
requestAnimationFrame(t => { time.last = t; frame(t); });
window.chemPlayground = { eng, R, receive, setT, setSpeed, togglePlay, pause: () => setPlaying(false), startConditioning, beginPlacing, dropMolecule, tickConditioning, get cond() { return cond; }, get placing() { return placing; } }; // console / test hooks
})();
