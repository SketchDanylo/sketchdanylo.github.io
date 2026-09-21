/* Chem Playground — field renderer.
 * Coordinates: world Å, x right, y down, +z toward the viewer. Camera {cx, cy, span}: span is the
 * visible width in Å. Density mode sums one Gaussian per atom whose 0.135 isoline sits on the atom's
 * van der Waals radius, then draws marching-squares isolines, so the outer contour is the molecular
 * surface at true scale and inner rings map the nuclei and bond ridges. */
(function (root) {
'use strict';
const { ELEMENTS } = root.ChemEngine;
const ISO = [0.135, 0.36, 0.72, 1.08];           // outermost ≈ van der Waals surface
const SIGMA = ELEMENTS.map(e => e.rvdw / Math.sqrt(2 * Math.log(1 / ISO[0])));
const RGB = ELEMENTS.map(e => { const h = e.color.slice(1); return [parseInt(h.slice(0, 2), 16), parseInt(h.slice(2, 4), 16), parseInt(h.slice(4, 6), 16)]; });
const LUT_N = 2048, LUT_MAX = 4.6, LUT = new Float32Array(LUT_N + 1);
for (let k = 0; k <= LUT_N; k++) LUT[k] = Math.exp(-k / LUT_N * LUT_MAX);

function shade(rgb, f) { return 'rgb(' + rgb.map(c => Math.max(0, Math.min(255, Math.round(c * f)))).join(',') + ')'; }

class FieldRenderer {
  constructor(canvas, opts = {}) {
    this.canvas = canvas; this.ctx = canvas.getContext('2d');
    this.cam = { cx: 25, cy: 15, span: 40 };
    this.mode = 'density';
    this.cell = opts.cell || 5;
    this.showBox = opts.showBox !== false;
    this.showGrid = opts.showGrid !== false;
    this.sprites = new Map();
    this.glow = document.createElement('canvas'); this.gctx = this.glow.getContext('2d');
    this.resize();
  }
  resize(w, h) {
    const dpr = Math.min(2, window.devicePixelRatio || 1);
    const W = w || this.canvas.clientWidth || this.canvas.width, H = h || this.canvas.clientHeight || this.canvas.height;
    this.W = W; this.H = H; this.dpr = dpr;
    this.canvas.width = Math.round(W * dpr); this.canvas.height = Math.round(H * dpr);
    const cs = this.cell;
    this.gw = Math.ceil(W / cs) + 2; this.gh = Math.ceil(H / cs) + 2;
    const n = this.gw * this.gh;
    this.rho = new Float32Array(n); this.cr = new Float32Array(n); this.cg = new Float32Array(n); this.cb = new Float32Array(n);
    this.glow.width = this.gw; this.glow.height = this.gh;
    this.img = this.gctx.createImageData(this.gw, this.gh);
  }
  get scale() { return this.W / this.cam.span; }
  toScreen(x, y) { const s = this.scale; return [(x - this.cam.cx) * s + this.W / 2, (y - this.cam.cy) * s + this.H / 2]; }
  toWorld(sx, sy) { const s = this.scale; return [(sx - this.W / 2) / s + this.cam.cx, (sy - this.H / 2) / s + this.cam.cy]; }

  /* scene: {N, type, pos (render positions, 3N), bonds:[{i,j,order,strength}], box, hover, selected:Set, ghost, flashes, tweezer, brush, time} */
  draw(sc) {
    const ctx = this.ctx, dpr = this.dpr, W = this.W, H = this.H;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    // field
    const g = ctx.createRadialGradient(W * 0.5, H * 0.45, 0, W * 0.5, H * 0.45, Math.hypot(W, H) * 0.6);
    g.addColorStop(0, '#242424'); g.addColorStop(1, '#111111');
    ctx.fillStyle = g; ctx.fillRect(0, 0, W, H);
    if (this.showGrid) this._grid(sc);
    if (this.showBox && sc.box && !sc.box3) { this._boxUnder(sc.box); this._boundary(sc.box, sc.bounds); }
    if (this.showBox && sc.box3) this._box3(sc.box3, 'under');
    const atoms = this._atomsView(sc);
    this._density(sc, atoms);
    this._bonds(sc, 'hair');
    this._nuclei(sc, atoms, 0.17);
    if (sc.ghost) this._ghost(sc.ghost);
    this._flashes(sc);
    this._overlay(sc);
    if (sc.inspect) this._leader(sc.inspect, sc.now);
    if (this.showBox && sc.box3) this._box3(sc.box3, 'over');
    else if (this.showBox && sc.box) this._boxOver(sc.box, sc.boxHot, sc.bounds);
  }

  _atomsView(sc) {
    // visible atoms sorted back-to-front
    const out = [], s = this.scale, P = sc.pos;
    const x0 = this.cam.cx - this.W / 2 / s - 4, x1 = this.cam.cx + this.W / 2 / s + 4;
    const y0 = this.cam.cy - this.H / 2 / s - 4, y1 = this.cam.cy + this.H / 2 / s + 4;
    for (let i = 0; i < sc.N; i++) {
      const x = P[3 * i], y = P[3 * i + 1];
      if (x < x0 || x > x1 || y < y0 || y > y1) continue;
      out.push(i);
    }
    out.sort((a, b) => P[3 * a + 2] - P[3 * b + 2]);
    return out;
  }

  _grid(sc) {
    const ctx = this.ctx, s = this.scale;
    const steps = [0.1, 0.2, 0.5, 1, 2, 5, 10];
    let st = steps.find(v => v * s >= 26) || 10;
    const [wx0, wy0] = this.toWorld(0, 0), [wx1, wy1] = this.toWorld(this.W, this.H);
    const r = Math.min(1.3, 0.6 + st * s / 120);
    ctx.fillStyle = 'rgba(179,179,185,0.16)';
    for (let x = Math.floor(wx0 / st) * st; x <= wx1; x += st) {
      const sx = (x - this.cam.cx) * s + this.W / 2;
      const major = Math.abs(x / 10 - Math.round(x / 10)) < 1e-6;
      for (let y = Math.floor(wy0 / st) * st; y <= wy1; y += st) {
        const sy = (y - this.cam.cy) * s + this.H / 2;
        const majorY = Math.abs(y / 10 - Math.round(y / 10)) < 1e-6;
        if (major && majorY) { ctx.fillRect(sx - 3, sy - 0.5, 6, 1); ctx.fillRect(sx - 0.5, sy - 3, 1, 6); }
        else ctx.fillRect(sx - r / 2, sy - r / 2, r, r);
      }
    }
  }
  _boxUnder(b) {
    const ctx = this.ctx, [x0, y0] = this.toScreen(b.x0, b.y0), [x1, y1] = this.toScreen(b.x1, b.y1);
    ctx.save();
    ctx.fillStyle = 'rgba(7,7,7,0.5)';
    ctx.beginPath(); ctx.rect(0, 0, this.W, this.H); ctx.rect(x0, y0, x1 - x0, y1 - y0); ctx.fill('evenodd');
    ctx.restore();
  }
  /* What the chamber face is, drawn on the chamber itself: a solid bound is a crisp hard line,
     a forcefield is a cushion that glows inward, and a void wall fades the outside away. */
  _boundary(b, cfg) {
    if (!cfg) return;
    const ctx = this.ctx, s = this.scale;
    const [x0, y0] = this.toScreen(b.x0, b.y0), [x1, y1] = this.toScreen(b.x1, b.y1);
    ctx.save();
    if (cfg.mode === 'forcefield') {
      const w = Math.min(26, Math.max(3, cfg.range * s));
      const band = (gx0, gy0, gx1, gy1, rx, ry, rw, rh) => {
        const g = ctx.createLinearGradient(gx0, gy0, gx1, gy1);
        g.addColorStop(0, 'rgba(179,179,185,0.16)'); g.addColorStop(1, 'rgba(179,179,185,0)');
        ctx.fillStyle = g; ctx.fillRect(rx, ry, rw, rh);
      };
      band(x0, 0, x0 + w, 0, x0, y0, w, y1 - y0);
      band(x1, 0, x1 - w, 0, x1 - w, y0, w, y1 - y0);
      band(0, y0, 0, y0 + w, x0, y0, x1 - x0, w);
      band(0, y1, 0, y1 - w, x0, y1 - w, x1 - x0, w);
    }
    if (cfg.voidT || cfg.voidP) {
      // the outside: warm where energy leaves, cool where impulse is swallowed
      const w = Math.min(16, Math.max(5, 1.6 * s));
      const tint = cfg.voidT ? '255,150,90' : '179,179,185';
      const band = (gx0, gy0, gx1, gy1, rx, ry, rw, rh) => {
        const g = ctx.createLinearGradient(gx0, gy0, gx1, gy1);
        g.addColorStop(0, 'rgba(' + tint + ',0.10)'); g.addColorStop(1, 'rgba(' + tint + ',0)');
        ctx.fillStyle = g; ctx.fillRect(rx, ry, rw, rh);
      };
      const inset = Math.min(18, (y1 - y0) * 0.18), insetX = Math.min(18, (x1 - x0) * 0.18);
      band(x0, 0, x0 - w, 0, x0 - w, y0 + inset, w, y1 - y0 - 2 * inset);
      band(x1, 0, x1 + w, 0, x1, y0 + inset, w, y1 - y0 - 2 * inset);
      band(0, y0, 0, y0 - w, x0 + insetX, y0 - w, x1 - x0 - 2 * insetX, w);
      band(0, y1, 0, y1 + w, x0 + insetX, y1, x1 - x0 - 2 * insetX, w);
    }
    ctx.restore();
  }
  _boxOver(b, hot, cfg) {
    const ctx = this.ctx, [x0, y0] = this.toScreen(b.x0, b.y0), [x1, y1] = this.toScreen(b.x1, b.y1);
    ctx.save();
    const field = cfg && cfg.mode === 'forcefield';
    ctx.strokeStyle = field ? 'rgba(179,179,185,0.5)' : 'rgba(179,179,185,0.38)';
    ctx.lineWidth = field ? 1.4 : 1;
    if (field) ctx.setLineDash([5, 4]);
    ctx.strokeRect(Math.round(x0) + 0.5, Math.round(y0) + 0.5, Math.round(x1 - x0), Math.round(y1 - y0));
    ctx.setLineDash([]);
    // corner brackets
    ctx.strokeStyle = 'rgba(233,228,216,0.75)'; ctx.lineWidth = 1.5;
    const L = 12;
    for (const [cx, cy, dx, dy] of [[x0, y0, 1, 1], [x1, y0, -1, 1], [x0, y1, 1, -1], [x1, y1, -1, -1]]) {
      ctx.beginPath(); ctx.moveTo(cx + dx * L, cy); ctx.lineTo(cx, cy); ctx.lineTo(cx, cy + dy * L); ctx.stroke();
    }
    if (hot) {
      ctx.strokeStyle = '#ffb23f'; ctx.lineWidth = 2;
      ctx.beginPath();
      if (hot.includes('x')) { ctx.moveTo(x1, y0); ctx.lineTo(x1, y1); }
      if (hot.includes('y')) { ctx.moveTo(x0, y1); ctx.lineTo(x1, y1); }
      ctx.stroke();
    }
    // dimension lines, drawn like a technical drawing: width over the top edge, height beside the right edge
    ctx.font = '400 10px "Martian Mono", monospace';
    const dimCol = 'rgba(162,162,167,0.8)', hotCol = '#ffb23f';
    const wl = ((b.x1 - b.x0) / 10).toFixed(2) + ' nm', hl = ((b.y1 - b.y0) / 10).toFixed(2) + ' nm';
    const ty = Math.max(y0 - 9, 10), mx = (Math.max(x0, 0) + Math.min(x1, this.W)) / 2;
    ctx.strokeStyle = 'rgba(179,179,185,0.35)'; ctx.lineWidth = 1;
    if (y0 - 9 > 4) {
      const tw = ctx.measureText(wl).width / 2 + 8;
      ctx.beginPath(); ctx.moveTo(x0, ty); ctx.lineTo(mx - tw, ty); ctx.moveTo(mx + tw, ty); ctx.lineTo(x1, ty);
      ctx.moveTo(x0, ty - 3); ctx.lineTo(x0, ty + 3); ctx.moveTo(x1, ty - 3); ctx.lineTo(x1, ty + 3); ctx.stroke();
    }
    ctx.fillStyle = hot && hot.includes('x') ? hotCol : dimCol; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    ctx.fillText(wl, mx, ty);
    const rx = Math.min(x1 + 10, this.W - 10), my = (Math.max(y0, 0) + Math.min(y1, this.H)) / 2;
    ctx.save(); ctx.translate(rx, my); ctx.rotate(-Math.PI / 2);
    ctx.fillStyle = hot && hot.includes('y') ? hotCol : dimCol; ctx.fillText(hl, 0, 0); ctx.restore();
    ctx.restore();
  }

  /* Tilted chamber: the caller hands over the eight corners already in view coordinates, so the
     wireframe, the atoms and the contours all live in the same rotated frame. Edges that run
     behind the sample are drawn faint, which is the whole cue for which way the box is turned. */
  _box3(corners, pass) {
    const ctx = this.ctx;
    const P = corners.map(c => this.toScreen(c[0], c[1]));
    const zs = corners.map(c => c[2]);
    const zmin = Math.min(...zs), zmax = Math.max(...zs), span = Math.max(1e-6, zmax - zmin);
    const EDGES = [[0, 1], [1, 3], [3, 2], [2, 0], [4, 5], [5, 7], [7, 6], [6, 4], [0, 4], [1, 5], [2, 6], [3, 7]];
    ctx.save();
    if (pass === 'under') {
      for (const [a, b] of EDGES) {
        const depth = ((zs[a] + zs[b]) / 2 - zmin) / span;
        if (depth > 0.5) continue;                       // far edges only
        ctx.strokeStyle = 'rgba(179,179,185,' + (0.12 + 0.12 * depth) + ')';
        ctx.lineWidth = 1;
        ctx.beginPath(); ctx.moveTo(P[a][0], P[a][1]); ctx.lineTo(P[b][0], P[b][1]); ctx.stroke();
      }
    } else {
      for (const [a, b] of EDGES) {
        const depth = ((zs[a] + zs[b]) / 2 - zmin) / span;
        if (depth <= 0.5) continue;                      // near edges over the sample
        ctx.strokeStyle = 'rgba(200,200,205,' + (0.3 + 0.35 * depth) + ')';
        ctx.lineWidth = 1.1;
        ctx.beginPath(); ctx.moveTo(P[a][0], P[a][1]); ctx.lineTo(P[b][0], P[b][1]); ctx.stroke();
      }
      ctx.fillStyle = 'rgba(233,228,216,.6)';
      for (let i = 0; i < 8; i++) { const near = (zs[i] - zmin) / span > 0.5; if (near) ctx.fillRect(P[i][0] - 1.5, P[i][1] - 1.5, 3, 3); }
    }
    ctx.restore();
  }

  /* Ties the inspector card to its atom: a ring on the atom and a hairline to the card edge. */
  _leader(ins, now) {
    const ctx = this.ctx, [x, y] = this.toScreen(ins.x, ins.y);
    const r = Math.max(9, ins.r * 0.5 * this.scale);
    const pulse = 0.5 + 0.5 * Math.sin((now || 0) / 520);
    ctx.save();
    ctx.strokeStyle = 'rgba(255,178,63,' + (0.5 + 0.3 * pulse) + ')'; ctx.lineWidth = 1.2;
    ctx.beginPath(); ctx.arc(x, y, r, 0, 7); ctx.stroke();
    ctx.strokeStyle = 'rgba(255,178,63,.3)'; ctx.lineWidth = 1; ctx.setLineDash([3, 4]);
    const dx = ins.ax - x, dy = ins.ay - y, d = Math.hypot(dx, dy) || 1;
    ctx.beginPath(); ctx.moveTo(x + dx / d * (r + 3), y + dy / d * (r + 3)); ctx.lineTo(ins.ax, ins.ay); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = 'rgba(255,178,63,.75)';
    ctx.beginPath(); ctx.arc(ins.ax, ins.ay, 2.2, 0, 7); ctx.fill();
    ctx.restore();
  }

  _density(sc, atoms) {
    const cs = this.cell, gw = this.gw, gh = this.gh, rho = this.rho, cr = this.cr, cg = this.cg, cb = this.cb;
    rho.fill(0); cr.fill(0); cg.fill(0); cb.fill(0);
    const s = this.scale, P = sc.pos, T = sc.type;
    const ox = this.W / 2 - this.cam.cx * s, oy = this.H / 2 - this.cam.cy * s;
    for (const i of atoms) {
      const t = T[i], sig = SIGMA[t] * s / cs, rad = Math.ceil(3.03 * sig);
      const px = (P[3 * i] * s + ox) / cs + 1, py = (P[3 * i + 1] * s + oy) / cs + 1;
      const gx0 = Math.max(0, Math.floor(px - rad)), gx1 = Math.min(gw - 1, Math.ceil(px + rad));
      const gy0 = Math.max(0, Math.floor(py - rad)), gy1 = Math.min(gh - 1, Math.ceil(py + rad));
      if (gx0 > gx1 || gy0 > gy1) continue;
      const inv = LUT_N / (2 * sig * sig * LUT_MAX), c = RGB[t];
      // depth: atoms behind the mid-plane are slightly fainter
      const depth = 1 + Math.max(-0.35, Math.min(0.2, P[3 * i + 2] * 0.04));
      for (let gy = gy0; gy <= gy1; gy++) {
        const dy = gy - py, dy2 = dy * dy, row = gy * gw;
        for (let gx = gx0; gx <= gx1; gx++) {
          const dx = gx - px, k = ((dx * dx + dy2) * inv) | 0;
          if (k >= LUT_N) continue;
          const w = LUT[k] * depth, idx = row + gx;
          rho[idx] += w; cr[idx] += w * c[0]; cg[idx] += w * c[1]; cb[idx] += w * c[2];
        }
      }
    }
    // tinted glow image
    const d = this.img.data;
    for (let k = 0, n = gw * gh; k < n; k++) {
      const r = rho[k], o = 4 * k;
      if (r < 0.01) { d[o + 3] = 0; continue; }
      const a = Math.min(1, r / 1.6);
      d[o] = cr[k] / r; d[o + 1] = cg[k] / r; d[o + 2] = cb[k] / r;
      d[o + 3] = 255 * (0.05 + 0.3 * a * a) * Math.min(1, r / ISO[0]);
    }
    this.gctx.putImageData(this.img, 0, 0);
    const ctx = this.ctx;
    ctx.save();
    ctx.imageSmoothingEnabled = true; ctx.imageSmoothingQuality = 'high';
    ctx.globalCompositeOperation = 'lighter';
    ctx.drawImage(this.glow, -cs, -cs, gw * cs, gh * cs);
    ctx.restore();
    // isolines (marching squares)
    ctx.save();
    ctx.lineJoin = 'round'; ctx.lineCap = 'round';
    for (let L = ISO.length - 1; L >= 0; L--) {
      const lev = ISO[L];
      ctx.beginPath();
      for (let gy = 0; gy < gh - 1; gy++) {
        const r0 = gy * gw, r1 = r0 + gw;
        for (let gx = 0; gx < gw - 1; gx++) {
          const a = rho[r0 + gx], b = rho[r0 + gx + 1], c = rho[r1 + gx + 1], dd = rho[r1 + gx];
          const code = (a > lev ? 8 : 0) | (b > lev ? 4 : 0) | (c > lev ? 2 : 0) | (dd > lev ? 1 : 0);
          if (code === 0 || code === 15) continue;
          const x = (gx - 1) * cs, y = (gy - 1) * cs;
          const top = x + cs * (lev - a) / (b - a), right = y + cs * (lev - b) / (c - b);
          const bottom = x + cs * (lev - dd) / (c - dd), left = y + cs * (lev - a) / (dd - a);
          switch (code) {
            case 1: case 14: ctx.moveTo(x, left); ctx.lineTo(bottom, y + cs); break;
            case 2: case 13: ctx.moveTo(bottom, y + cs); ctx.lineTo(x + cs, right); break;
            case 3: case 12: ctx.moveTo(x, left); ctx.lineTo(x + cs, right); break;
            case 4: case 11: ctx.moveTo(top, y); ctx.lineTo(x + cs, right); break;
            case 5: ctx.moveTo(x, left); ctx.lineTo(top, y); ctx.moveTo(bottom, y + cs); ctx.lineTo(x + cs, right); break;
            case 6: case 9: ctx.moveTo(top, y); ctx.lineTo(bottom, y + cs); break;
            case 7: case 8: ctx.moveTo(x, left); ctx.lineTo(top, y); break;
            case 10: ctx.moveTo(top, y); ctx.lineTo(x + cs, right); ctx.moveTo(x, left); ctx.lineTo(bottom, y + cs); break;
          }
        }
      }
      ctx.strokeStyle = L === 0 ? 'rgba(165,199,228,0.62)' : 'rgba(179,179,185,' + (0.34 - L * 0.06) + ')';
      ctx.lineWidth = L === 0 ? 1.15 : 0.8;
      ctx.stroke();
    }
    ctx.restore();
  }

  _bonds(sc, style) {
    const ctx = this.ctx, s = this.scale, P = sc.pos;
    ctx.save(); ctx.lineCap = 'round';
    const base = style === 'stick' ? Math.max(1.2, 0.16 * s) : Math.max(0.9, Math.min(2.2, 0.045 * s));
    for (const b of sc.bonds) {
      const i = b.i, j = b.j;
      const [x1, y1] = this.toScreen(P[3 * i], P[3 * i + 1]), [x2, y2] = this.toScreen(P[3 * j], P[3 * j + 1]);
      const dx = x2 - x1, dy = y2 - y1, L = Math.hypot(dx, dy) || 1, nx = -dy / L, ny = dx / L;
      const a = Math.min(1, b.strength) * (style === 'stick' ? 0.95 : 0.8);
      const n = b.order, gap = Math.max(2.2, Math.min(0.2 * s, 7));
      const lines = n > 2.6 ? [-1, 0, 1] : n > 1.7 ? [-0.5, 0.5] : n > 1.25 ? [-0.5, 0.5] : [0];
      const dashed = n > 1.25 && n <= 1.7;
      lines.forEach((o, k) => {
        ctx.beginPath();
        ctx.setLineDash(dashed && k === 1 ? [Math.max(2, 0.12 * s), Math.max(2, 0.1 * s)] : []);
        ctx.moveTo(x1 + nx * o * gap, y1 + ny * o * gap); ctx.lineTo(x2 + nx * o * gap, y2 + ny * o * gap);
        if (style === 'stick') {
          ctx.strokeStyle = 'rgba(8,12,18,' + a + ')'; ctx.lineWidth = (lines.length > 1 ? base * 0.62 : base) + 2; ctx.stroke();
          ctx.strokeStyle = 'rgba(214,208,196,' + a + ')'; ctx.lineWidth = lines.length > 1 ? base * 0.62 : base; ctx.stroke();
        } else { ctx.strokeStyle = 'rgba(233,228,216,' + a + ')'; ctx.lineWidth = lines.length > 1 ? base * 0.8 : base; ctx.stroke(); }
      });
    }
    ctx.setLineDash([]);
    ctx.restore();
  }
  _sprite(t, r, dim) {
    const R = Math.max(1, Math.round(r * 2) / 2), key = t * 100000 + R * 10 + dim;
    let sp = this.sprites.get(key);
    if (sp) return sp;
    if (this.sprites.size > 900) this.sprites.clear();
    const dpr = this.dpr, size = Math.ceil(2 * R * dpr + 4), c = document.createElement('canvas');
    c.width = c.height = size;
    const x = c.getContext('2d'), m = size / 2, rr = R * dpr;
    const rgb = RGB[t], f = [0.84, 0.93, 1][dim]; // gentle depth cue: never enough to change an element's apparent colour
    const gr = x.createRadialGradient(m - rr * 0.35, m - rr * 0.4, rr * 0.05, m, m, rr);
    gr.addColorStop(0, shade(rgb, 1.35 * f)); gr.addColorStop(0.45, shade(rgb, 0.95 * f)); gr.addColorStop(1, shade(rgb, 0.38 * f));
    x.fillStyle = gr; x.beginPath(); x.arc(m, m, rr, 0, Math.PI * 2); x.fill();
    x.strokeStyle = 'rgba(0,0,0,0.35)'; x.lineWidth = Math.max(0.5, dpr * 0.6); x.stroke();
    sp = { c, size: size / dpr }; this.sprites.set(key, sp);
    return sp;
  }
  _spheres(sc, atoms, k) {
    const ctx = this.ctx, s = this.scale, P = sc.pos, T = sc.type;
    for (const i of atoms) {
      const t = T[i], z = P[3 * i + 2];
      const r = ELEMENTS[t].rvdw * k * s; // Orthographic view: depth changes shading, never physical radius.
      const [x, y] = this.toScreen(P[3 * i], P[3 * i + 1]);
      const dim = z < -2 ? 0 : z < 1 ? 1 : 2;
      if (r > 220) {
        const rgb = RGB[t], gr = ctx.createRadialGradient(x - r * 0.35, y - r * 0.4, r * 0.05, x, y, r);
        gr.addColorStop(0, shade(rgb, 1.3)); gr.addColorStop(1, shade(rgb, 0.4));
        ctx.fillStyle = gr; ctx.beginPath(); ctx.arc(x, y, r, 0, Math.PI * 2); ctx.fill();
      } else {
        const sp = this._sprite(t, r, dim);
        ctx.drawImage(sp.c, x - sp.size / 2, y - sp.size / 2, sp.size, sp.size);
      }
      if (k < 1 && s > 55) this._label(t, x, y, r);
    }
  }
  _nuclei(sc, atoms, k) {
    const ctx = this.ctx, s = this.scale, P = sc.pos, T = sc.type;
    for (const i of atoms) {
      const t = T[i], r = Math.max(1.8, ELEMENTS[t].rvdw * k * s);
      const [x, y] = this.toScreen(P[3 * i], P[3 * i + 1]);
      const z = P[3 * i + 2], dim = z < -2 ? 0 : z < 1 ? 1 : 2;
      if (r < 3) { ctx.fillStyle = shade(RGB[t], [0.84, 0.94, 1.05][dim]); ctx.beginPath(); ctx.arc(x, y, r, 0, Math.PI * 2); ctx.fill(); }
      else { const sp = this._sprite(t, r, dim); ctx.drawImage(sp.c, x - sp.size / 2, y - sp.size / 2, sp.size, sp.size); }
      if (s > 80) this._label(t, x, y, r, true);
    }
  }
  _label(t, x, y, r, below) {
    const ctx = this.ctx;
    ctx.font = '500 ' + Math.min(14, Math.max(9, r * 0.8)) + 'px "Martian Mono", monospace';
    ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
    ctx.fillStyle = below ? 'rgba(233,228,216,0.75)' : 'rgba(10,14,20,0.8)';
    if(this.showLabels!==false) ctx.fillText(ELEMENTS[t].sym, x, below ? y + r + 9 : y);
  }
  _ghost(gh) {
    const ctx = this.ctx, s = this.scale;
    ctx.save(); ctx.globalAlpha = 0.55;
    for (const b of gh.bonds || []) {
      const A = gh.atoms[b.a], B = gh.atoms[b.b];
      const [x1, y1] = this.toScreen(A.x, A.y), [x2, y2] = this.toScreen(B.x, B.y);
      ctx.strokeStyle = gh.ok ? 'rgba(233,228,216,.7)' : 'rgba(255,107,91,.8)'; ctx.lineWidth = Math.max(1, 0.06 * s);
      ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x2, y2); ctx.stroke();
    }
    for (const a of gh.atoms) {
      const [x, y] = this.toScreen(a.x, a.y), r = Math.max(2.5, ELEMENTS[a.t].rvdw * 0.3 * s);
      ctx.fillStyle = gh.ok ? ELEMENTS[a.t].color : '#ff6b5b';
      ctx.beginPath(); ctx.arc(x, y, r, 0, Math.PI * 2); ctx.fill();
      ctx.strokeStyle = 'rgba(255,255,255,.35)'; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.arc(x, y, ELEMENTS[a.t].rvdw * s, 0, Math.PI * 2); ctx.setLineDash([3, 4]); ctx.stroke(); ctx.setLineDash([]);
    }
    ctx.restore();
    if (gh.fling) {
      const [x1, y1] = this.toScreen(gh.fling.x0, gh.fling.y0), [x2, y2] = this.toScreen(gh.fling.x1, gh.fling.y1);
      const ang = Math.atan2(y2 - y1, x2 - x1);
      ctx.save(); ctx.strokeStyle = '#ffb23f'; ctx.fillStyle = '#ffb23f'; ctx.lineWidth = 1.6;
      ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x2, y2); ctx.stroke();
      ctx.beginPath(); ctx.moveTo(x2, y2); ctx.lineTo(x2 - 10 * Math.cos(ang - 0.4), y2 - 10 * Math.sin(ang - 0.4)); ctx.lineTo(x2 - 10 * Math.cos(ang + 0.4), y2 - 10 * Math.sin(ang + 0.4)); ctx.fill();
      ctx.font = '400 10.5px "Martian Mono", monospace'; ctx.textBaseline = 'bottom';
      ctx.fillText(gh.fling.label, x2 + 10, y2 - 6);
      ctx.restore();
    }
  }
  _flashes(sc) {
    if (!sc.flashes || !sc.flashes.length) return;
    const ctx = this.ctx, s = this.scale, now = sc.now;
    ctx.save();
    for (const f of sc.flashes) {
      const t = (now - f.t0) / f.dur; if (t < 0 || t > 1) continue;
      const [x, y] = this.toScreen(f.x, f.y), e = 1 - Math.pow(1 - t, 3);
      const r = Math.max(6, (0.4 + e * 1.4) * s * (f.big ? 1.8 : 1));
      ctx.strokeStyle = f.kind === 'form' ? 'rgba(255,178,63,' + (1 - t) * 0.9 + ')' : 'rgba(143,211,255,' + (1 - t) * 0.7 + ')';
      ctx.lineWidth = f.kind === 'form' ? 1.6 : 1;
      if (f.kind !== 'form') ctx.setLineDash([2, 4]);
      ctx.beginPath(); ctx.arc(x, y, r, 0, Math.PI * 2); ctx.stroke(); ctx.setLineDash([]);
    }
    ctx.restore();
  }
  _overlay(sc) {
    const ctx = this.ctx, s = this.scale, P = sc.pos, T = sc.type;
    ctx.save();
    if (sc.selected && sc.selected.size) {
      ctx.strokeStyle = 'rgba(255,178,63,.85)'; ctx.lineWidth = 1.2; ctx.setLineDash([3, 3]);
      for (const i of sc.selected) {
        if (i >= sc.N) continue;
        const [x, y] = this.toScreen(P[3 * i], P[3 * i + 1]);
        ctx.beginPath(); ctx.arc(x, y, Math.max(6, ELEMENTS[T[i]].rvdw * 0.55 * s), 0, Math.PI * 2); ctx.stroke();
      }
      ctx.setLineDash([]);
    }
    if (sc.pinned) {
      ctx.fillStyle = 'rgba(143,211,255,.9)';
      for (let i = 0; i < sc.N; i++) if (sc.pinned[i]) { const [x, y] = this.toScreen(P[3 * i], P[3 * i + 1]); ctx.fillRect(x - 2, y - 2, 4, 4); }
    }
    if (sc.hover >= 0 && sc.hover < sc.N) {
      const i = sc.hover, [x, y] = this.toScreen(P[3 * i], P[3 * i + 1]);
      ctx.strokeStyle = sc.eraseHover ? '#ff6b5b' : 'rgba(255,178,63,.95)'; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.arc(x, y, Math.max(7, ELEMENTS[T[i]].rvdw * (this.mode === 'space' ? 1.05 : 0.5) * s), 0, Math.PI * 2); ctx.stroke();
    }
    if (sc.tweezer) {
      const tw = sc.tweezer, i = tw.i;
      if (i < sc.N) {
        const [x1, y1] = this.toScreen(P[3 * i], P[3 * i + 1]), [x2, y2] = this.toScreen(tw.x, tw.y);
        ctx.strokeStyle = 'rgba(255,178,63,.8)'; ctx.lineWidth = 1.2; ctx.setLineDash([4, 3]);
        ctx.beginPath(); ctx.moveTo(x1, y1); ctx.lineTo(x2, y2); ctx.stroke(); ctx.setLineDash([]);
        ctx.fillStyle = '#ffb23f'; ctx.beginPath(); ctx.arc(x2, y2, 3.5, 0, Math.PI * 2); ctx.fill();
      }
    }
    if (sc.brush) {
      const b = sc.brush, [x, y] = this.toScreen(b.x, b.y), r = b.r * s;
      ctx.strokeStyle = b.cool ? 'rgba(143,211,255,.8)' : 'rgba(255,178,63,.8)'; ctx.lineWidth = 1.2;
      ctx.setLineDash([2, 5]); ctx.beginPath(); ctx.arc(x, y, r, 0, Math.PI * 2); ctx.stroke(); ctx.setLineDash([]);
      if (b.active) { ctx.fillStyle = b.cool ? 'rgba(143,211,255,.07)' : 'rgba(255,178,63,.08)'; ctx.fill(); }
    }
    if (sc.marquee) {
      const m = sc.marquee, [x1, y1] = this.toScreen(m.x0, m.y0), [x2, y2] = this.toScreen(m.x1, m.y1);
      ctx.fillStyle = 'rgba(255,178,63,.06)'; ctx.strokeStyle = 'rgba(255,178,63,.6)'; ctx.lineWidth = 1;
      ctx.fillRect(x1, y1, x2 - x1, y2 - y1); ctx.strokeRect(x1 + 0.5, y1 + 0.5, x2 - x1, y2 - y1);
    }
    ctx.restore();
  }
}
root.FieldRenderer = FieldRenderer;
root.FieldRenderer.ISO = ISO;
})(typeof self !== 'undefined' ? self : this);
