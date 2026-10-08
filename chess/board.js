// The board: pieces that glide between positions, drag or click to move,
// and arrows that spring from one move to the next instead of blinking.
import { PIECES } from './pieces.js?v=2026-10-08-a';

const FILES = 'abcdefgh';
const ARROW_UNIT = 100; // one square in arrow-space

function squareXY(square, orientation) {
  const file = FILES.indexOf(square[0]);
  const rank = Number(square[1]) - 1;
  return orientation === 'w' ? { x: file, y: 7 - rank } : { x: 7 - file, y: rank };
}

export function parseFen(fen) {
  const map = new Map();
  const rows = fen.split(' ')[0].split('/');
  rows.forEach((row, i) => {
    let file = 0;
    for (const ch of row) {
      if (/\d/.test(ch)) { file += Number(ch); continue; }
      const color = ch === ch.toUpperCase() ? 'w' : 'b';
      map.set(FILES[file] + (8 - i), color + ch.toUpperCase());
      file++;
    }
  });
  return map;
}

export function pieceImage(code) {
  return `url("data:image/svg+xml;utf8,${encodeURIComponent(PIECES[code])}")`;
}

export class Board {
  constructor(root, handlers = {}) {
    this.root = root;
    this.handlers = handlers;
    this.orientation = 'w';
    this.pieces = [];             // { el, code, square }
    this.position = new Map();
    this.movable = null;          // Map from -> [to]
    this.selected = null;
    this.edit = false;
    this.brush = null;
    this.arrows = new Map();
    this.hovered = null;
    this.drag = null;
    this.raf = 0;
    this.build();
  }

  build() {
    this.root.classList.add('board');
    this.root.innerHTML = `
      <div class="b-squares"></div>
      <div class="b-coords" aria-hidden="true"></div>
      <div class="b-marks"></div>
      <div class="b-pieces"></div>
      <svg class="b-arrows" viewBox="0 0 800 800" aria-hidden="true"></svg>
      <div class="b-promo" hidden></div>`;
    this.squaresEl = this.root.querySelector('.b-squares');
    this.coordsEl = this.root.querySelector('.b-coords');
    this.marksEl = this.root.querySelector('.b-marks');
    this.piecesEl = this.root.querySelector('.b-pieces');
    this.svg = this.root.querySelector('.b-arrows');
    this.promoEl = this.root.querySelector('.b-promo');
    this.squareEls = new Map();
    for (let rank = 1; rank <= 8; rank++) {
      for (const file of FILES) {
        const square = file + rank;
        const el = document.createElement('div');
        el.className = `sq ${(FILES.indexOf(file) + rank) % 2 === 0 ? 'dark' : 'light'}`;
        el.dataset.square = square;
        this.squaresEl.appendChild(el);
        this.squareEls.set(square, el);
      }
    }
    this.layoutSquares();
    this.renderCoords();

    this.root.addEventListener('pointerdown', event => this.pointerDown(event));
    this.root.addEventListener('pointermove', event => this.pointerHover(event));
    this.root.addEventListener('pointerleave', () => this.setHovered(null));
    this.root.addEventListener('contextmenu', event => event.preventDefault());
    window.addEventListener('pointermove', event => this.pointerMove(event));
    window.addEventListener('pointerup', event => this.pointerUp(event));
    window.addEventListener('pointercancel', () => this.cancelDrag());
  }

  layoutSquares() {
    this.squareEls.forEach((el, square) => {
      const { x, y } = squareXY(square, this.orientation);
      el.style.left = `${x * 12.5}%`;
      el.style.top = `${y * 12.5}%`;
    });
  }

  renderCoords() {
    const files = this.orientation === 'w' ? FILES : [...FILES].reverse().join('');
    const ranks = this.orientation === 'w' ? '87654321' : '12345678';
    this.coordsEl.innerHTML =
      [...files].map((f, i) => `<span class="cf ${(i + (this.orientation === 'w' ? 0 : 1)) % 2 ? 'on-light' : 'on-dark'}" style="left:${i * 12.5}%">${f}</span>`).join('') +
      [...ranks].map((r, i) => `<span class="cr ${(i + (this.orientation === 'w' ? 0 : 1)) % 2 ? 'on-dark' : 'on-light'}" style="top:${i * 12.5}%">${r}</span>`).join('');
  }

  setOrientation(color) {
    if (color === this.orientation) return;
    this.orientation = color;
    this.root.classList.add('flipping');
    this.layoutSquares();
    this.renderCoords();
    this.pieces.forEach(piece => this.place(piece.el, piece.square));
    this.renderMarks(this.marks || []);
    this.arrows.forEach(arrow => this.retarget(arrow));
    this.kick();
    clearTimeout(this.flipTimer);
    this.flipTimer = setTimeout(() => this.root.classList.remove('flipping'), 520);
  }

  place(el, square) {
    const { x, y } = squareXY(square, this.orientation);
    el.style.transform = `translate(${x * 100}%, ${y * 100}%)`;
  }

  // Morph from the current pieces to a new position: keep, slide, fade.
  setPosition(fen, { instant = false } = {}) {
    const next = parseFen(fen);
    this.position = next;
    const unmatched = new Map(next);
    const old = [];
    for (const piece of this.pieces) {
      if (unmatched.get(piece.square) === piece.code) {
        unmatched.delete(piece.square);
        old.push(piece);
      } else {
        piece.free = true;
      }
    }
    const free = this.pieces.filter(piece => piece.free);
    const kept = [...old];
    // Slide the nearest piece of the same kind into each new square.
    const wanted = [...unmatched];
    wanted.sort((a, b) => this.nearest(a, free).distance - this.nearest(b, free).distance);
    for (const [square, code] of wanted) {
      const { piece } = this.nearest([square, code], free);
      if (piece) {
        free.splice(free.indexOf(piece), 1);
        piece.free = false;
        piece.square = square;
        if (instant) piece.el.classList.add('no-anim');
        this.place(piece.el, square);
        piece.el.classList.add('moving');
        clearTimeout(piece.moveTimer);
        piece.moveTimer = setTimeout(() => piece.el.classList.remove('moving'), 360);
        kept.push(piece);
      } else {
        const el = document.createElement('div');
        el.className = 'pc entering';
        el.style.backgroundImage = pieceImage(code);
        this.place(el, square);
        this.piecesEl.appendChild(el);
        afterPaint(() => el.classList.remove('entering'));
        kept.push({ el, code, square });
      }
    }
    for (const piece of free) {
      piece.el.classList.add('leaving');
      setTimeout(() => piece.el.remove(), 320);
    }
    this.pieces = kept;
    this.pieces.forEach(piece => {
      piece.el.dataset.square = piece.square;
      if (instant) requestAnimationFrame(() => piece.el.classList.remove('no-anim'));
    });
    this.select(null);
  }

  nearest([square, code], pool) {
    let best = null;
    let distance = Infinity;
    const a = squareXY(square, 'w');
    for (const piece of pool) {
      if (piece.code !== code) continue;
      const b = squareXY(piece.square, 'w');
      const d = Math.hypot(a.x - b.x, a.y - b.y);
      if (d < distance) { distance = d; best = piece; }
    }
    return { piece: best, distance };
  }

  pieceAt(square) { return this.pieces.find(piece => piece.square === square) || null; }

  setMovable(map) { this.movable = map; if (!map) this.select(null); }

  setHighlights({ lastMove = null, check = null } = {}) {
    this.squareEls.forEach((el, square) => {
      el.classList.toggle('last', !!lastMove && (lastMove.from === square || lastMove.to === square));
      el.classList.toggle('check', check === square);
    });
  }

  select(square) {
    this.selected = square;
    const targets = square && this.movable ? this.movable.get(square) || [] : [];
    this.squareEls.forEach((el, sq) => {
      el.classList.toggle('sel', sq === square);
      el.classList.toggle('dot', targets.includes(sq) && !this.position.has(sq));
      el.classList.toggle('ring', targets.includes(sq) && this.position.has(sq));
    });
  }

  // ---------- pointer ----------

  squareFromPoint(clientX, clientY) {
    const rect = this.root.getBoundingClientRect();
    const x = Math.floor(((clientX - rect.left) / rect.width) * 8);
    const y = Math.floor(((clientY - rect.top) / rect.height) * 8);
    if (x < 0 || x > 7 || y < 0 || y > 7) return null;
    const file = this.orientation === 'w' ? x : 7 - x;
    const rank = this.orientation === 'w' ? 7 - y : y;
    return FILES[file] + (rank + 1);
  }

  canDrag(square) {
    if (this.edit) return this.position.has(square);
    return !!(this.movable && this.movable.has(square));
  }

  pointerDown(event) {
    if (event.button !== 0 && event.pointerType === 'mouse') {
      if (event.button === 2 && this.edit) {
        const square = this.squareFromPoint(event.clientX, event.clientY);
        if (square) this.handlers.onEdit && this.handlers.onEdit({ type: 'remove', square });
      }
      return;
    }
    const square = this.squareFromPoint(event.clientX, event.clientY);
    if (!square) return;
    event.preventDefault();

    if (this.edit && this.brush) {
      this.handlers.onEdit && this.handlers.onEdit({ type: 'paint', square, code: this.brush });
      this.painting = true;
      this.lastPainted = square;
      return;
    }

    // Second click of a click-click move.
    if (this.selected && this.selected !== square && this.movable && (this.movable.get(this.selected) || []).includes(square)) {
      const from = this.selected;
      this.select(null);
      this.handlers.onMove && this.handlers.onMove(from, square);
      return;
    }

    if (!this.canDrag(square)) {
      this.select(null);
      return;
    }
    const piece = this.pieceAt(square);
    if (!piece) return;
    const rect = this.root.getBoundingClientRect();
    this.drag = { from: square, piece, startX: event.clientX, startY: event.clientY, rect, active: false, wasSelected: this.selected === square };
    if (!this.edit) this.select(square);
  }

  pointerMove(event) {
    if (this.painting && this.edit && this.brush) {
      const square = this.squareFromPoint(event.clientX, event.clientY);
      if (square && square !== this.lastPainted) {
        this.lastPainted = square;
        this.handlers.onEdit && this.handlers.onEdit({ type: 'put', square, code: this.brush });
      }
      return;
    }
    const drag = this.drag;
    if (!drag) return;
    if (!drag.active) {
      if (Math.hypot(event.clientX - drag.startX, event.clientY - drag.startY) < 4) return;
      drag.active = true;
      drag.piece.el.classList.add('dragging');
      this.root.classList.add('is-dragging');
      this.setHovered(null);
    }
    const size = drag.rect.width / 8;
    const x = event.clientX - drag.rect.left - size / 2;
    const y = event.clientY - drag.rect.top - size / 2;
    drag.piece.el.style.transform = `translate(${x}px, ${y}px) scale(1.12)`;
    const over = this.squareFromPoint(event.clientX, event.clientY);
    this.squareEls.forEach((el, sq) => el.classList.toggle('over', sq === over && sq !== drag.from));
  }

  pointerUp(event) {
    if (this.painting) { this.painting = false; return; }
    const drag = this.drag;
    if (!drag) return;
    this.drag = null;
    this.root.classList.remove('is-dragging');
    this.squareEls.forEach(el => el.classList.remove('over'));
    const target = this.squareFromPoint(event.clientX, event.clientY);
    const el = drag.piece.el;
    if (!drag.active) {
      // Plain click: toggle selection.
      if (drag.wasSelected) this.select(null);
      return;
    }
    el.classList.remove('dragging');
    if (this.edit) {
      el.classList.add('no-anim');
      this.place(el, target || drag.from);
      requestAnimationFrame(() => el.classList.remove('no-anim'));
      this.handlers.onEdit && this.handlers.onEdit(target ? { type: 'move', from: drag.from, to: target } : { type: 'remove', square: drag.from });
      return;
    }
    const legal = target && this.movable && (this.movable.get(drag.from) || []).includes(target);
    if (legal) {
      // Drop exactly where the pointer let go, no slide.
      el.classList.add('no-anim');
      this.place(el, target);
      requestAnimationFrame(() => requestAnimationFrame(() => el.classList.remove('no-anim')));
      this.select(null);
      this.handlers.onMove && this.handlers.onMove(drag.from, target, { dropped: true });
    } else {
      this.place(el, drag.from);
    }
  }

  cancelDrag() {
    if (!this.drag) return;
    const { piece, from } = this.drag;
    piece.el.classList.remove('dragging');
    this.place(piece.el, from);
    this.drag = null;
    this.root.classList.remove('is-dragging');
  }

  // Drag a piece in from the editor palette.
  startExternalDrag(code, event) {
    const ghost = document.createElement('div');
    ghost.className = 'pc ghost';
    ghost.style.backgroundImage = pieceImage(code);
    const size = this.root.getBoundingClientRect().width / 8;
    ghost.style.width = ghost.style.height = `${size}px`;
    document.body.appendChild(ghost);
    const move = e => { ghost.style.transform = `translate(${e.clientX - size / 2}px, ${e.clientY - size / 2}px) scale(1.12)`; };
    move(event);
    const up = e => {
      window.removeEventListener('pointermove', move);
      window.removeEventListener('pointerup', up);
      ghost.remove();
      const square = this.squareFromPoint(e.clientX, e.clientY);
      if (square) this.handlers.onEdit && this.handlers.onEdit({ type: 'put', square, code });
      else if (Math.hypot(e.clientX - event.clientX, e.clientY - event.clientY) < 6) this.handlers.onPaletteClick && this.handlers.onPaletteClick(code);
    };
    window.addEventListener('pointermove', move);
    window.addEventListener('pointerup', up);
  }

  askPromotion(square, color) {
    return new Promise(resolve => {
      const { x, y } = squareXY(square, this.orientation);
      const down = y === 0;
      this.promoEl.hidden = false;
      this.promoEl.innerHTML = '';
      ['q', 'n', 'r', 'b'].forEach((type, i) => {
        const button = document.createElement('button');
        button.className = 'promo-pick';
        button.style.left = `${x * 12.5}%`;
        button.style.top = `${(down ? y + i : y - i) * 12.5}%`;
        button.style.backgroundImage = pieceImage(color + type.toUpperCase());
        button.style.transitionDelay = `${i * 30}ms`;
        button.setAttribute('aria-label', type);
        button.addEventListener('pointerdown', event => { event.stopPropagation(); finish(type); });
        this.promoEl.appendChild(button);
      });
      const finish = type => {
        this.promoEl.classList.remove('open');
        setTimeout(() => { this.promoEl.hidden = true; }, 200);
        this.promoEl.removeEventListener('pointerdown', cancel);
        resolve(type);
      };
      const cancel = () => finish(null);
      this.promoEl.addEventListener('pointerdown', cancel);
      requestAnimationFrame(() => this.promoEl.classList.add('open'));
    });
  }

  // ---------- marks (attacked, shielded, targets) ----------

  setMarks(marks) {
    this.marks = marks;
    this.renderMarks(marks);
  }

  renderMarks(marks) {
    const existing = new Map([...this.marksEl.children].map(el => [el.dataset.key, el]));
    const seen = new Set();
    for (const mark of marks) {
      const key = `${mark.type}:${mark.square}`;
      seen.add(key);
      let el = existing.get(key);
      if (!el) {
        el = document.createElement('div');
        el.className = `mark ${mark.type} entering`;
        el.dataset.key = key;
        this.marksEl.appendChild(el);
        afterPaint(() => el.classList.remove('entering'));
      }
      el.classList.remove('leaving');
      el.innerHTML = mark.count ? `<i>${mark.count}</i>` : '';
      const { x, y } = squareXY(mark.square, this.orientation);
      el.style.transform = `translate(${x * 100}%, ${y * 100}%)`;
      el.title = mark.title || '';
    }
    existing.forEach((el, key) => {
      if (seen.has(key)) return;
      el.classList.add('leaving');
      setTimeout(() => { if (el.classList.contains('leaving')) el.remove(); }, 300);
    });
  }

  // ---------- arrows ----------

  // arrows: [{ id, from, to, color, width, opacity, dashed, hover: data }]
  setArrows(list) {
    const seen = new Set();
    for (const spec of list) {
      seen.add(spec.id);
      let arrow = this.arrows.get(spec.id);
      const target = this.arrowTarget(spec);
      if (!arrow) {
        const el = document.createElementNS('http://www.w3.org/2000/svg', 'path');
        el.setAttribute('class', `arrow ${spec.className || ''}`);
        this.svg.appendChild(el);
        // New arrows grow out of their origin square.
        arrow = { id: spec.id, el, cur: { ...target, tx: target.fx, ty: target.fy, opacity: 0 }, vel: { fx: 0, fy: 0, tx: 0, ty: 0 } };
        this.arrows.set(spec.id, arrow);
      }
      arrow.spec = spec;
      arrow.target = target;
      arrow.dying = false;
      arrow.el.setAttribute('class', `arrow ${spec.className || ''}${this.hovered === spec.id ? ' hovered' : ''}`);
      arrow.el.style.fill = spec.color;
    }
    this.arrows.forEach(arrow => {
      if (!seen.has(arrow.id)) {
        arrow.dying = true;
        arrow.target = { ...arrow.target, opacity: 0 };
      }
    });
    if (this.hovered && !seen.has(this.hovered)) this.setHovered(null);
    this.svg.append(...[...this.arrows.values()].sort((a, b) => (a.spec.z || 0) - (b.spec.z || 0)).map(arrow => arrow.el));
    this.kick();
  }

  arrowTarget(spec) {
    const a = squareXY(spec.from, this.orientation);
    const b = squareXY(spec.to, this.orientation);
    return {
      fx: (a.x + 0.5) * ARROW_UNIT,
      fy: (a.y + 0.5) * ARROW_UNIT,
      tx: (b.x + 0.5) * ARROW_UNIT,
      ty: (b.y + 0.5) * ARROW_UNIT,
      width: spec.width || 14,
      opacity: spec.opacity ?? 0.85
    };
  }

  retarget(arrow) {
    arrow.target = { ...this.arrowTarget(arrow.spec), opacity: arrow.dying ? 0 : (arrow.spec.opacity ?? 0.85) };
  }

  kick() {
    if (!this.raf) this.raf = requestAnimationFrame(time => this.step(time));
    this.armFallback();
  }

  // A tab that stops painting never fires the next frame; settle the arrows anyway.
  armFallback() {
    clearTimeout(this.fallback);
    this.fallback = setTimeout(() => {
      if (this.raf) cancelAnimationFrame(this.raf);
      this.raf = 0;
      this.lastTime = 0;
      this.settle();
    }, 300);
  }

  settle() {
    for (const arrow of [...this.arrows.values()]) {
      Object.assign(arrow.cur, arrow.target);
      if (arrow.dying) { arrow.el.remove(); this.arrows.delete(arrow.id); continue; }
      arrow.el.setAttribute('d', arrowPath(arrow.cur, arrow.spec));
      arrow.el.style.opacity = arrow.cur.opacity.toFixed(3);
    }
  }

  step(time) {
    const dt = Math.min(0.034, this.lastTime ? (time - this.lastTime) / 1000 : 0.016);
    this.lastTime = time;
    let moving = false;
    const stiffness = 210;
    const damping = 26;
    for (const arrow of [...this.arrows.values()]) {
      const cur = arrow.cur;
      const target = arrow.target;
      for (const key of ['fx', 'fy', 'tx', 'ty']) {
        const force = (target[key] - cur[key]) * stiffness - arrow.vel[key] * damping;
        arrow.vel[key] += force * dt;
        cur[key] += arrow.vel[key] * dt;
        if (Math.abs(target[key] - cur[key]) > 0.1 || Math.abs(arrow.vel[key]) > 0.1) moving = true;
        else { cur[key] = target[key]; arrow.vel[key] = 0; }
      }
      const hoverBoost = this.hovered === arrow.id ? 1.22 : 1;
      const widthTarget = target.width * hoverBoost;
      cur.width += (widthTarget - cur.width) * Math.min(1, dt * 14);
      cur.opacity += (target.opacity - cur.opacity) * Math.min(1, dt * 10);
      if (Math.abs(widthTarget - cur.width) > 0.05 || Math.abs(target.opacity - cur.opacity) > 0.005) moving = true;
      if (arrow.dying && cur.opacity < 0.01) {
        arrow.el.remove();
        this.arrows.delete(arrow.id);
        continue;
      }
      arrow.el.setAttribute('d', arrowPath(cur, arrow.spec));
      arrow.el.style.opacity = cur.opacity.toFixed(3);
    }
    if (moving) {
      this.raf = requestAnimationFrame(next => this.step(next));
      this.armFallback();
    } else {
      this.raf = 0;
      this.lastTime = 0;
      clearTimeout(this.fallback);
    }
  }

  // Book arrows can be hovered: a hit test on the shafts, because the arrow layer lets clicks through.
  pointerHover(event) {
    if (this.drag && this.drag.active) return;
    const rect = this.root.getBoundingClientRect();
    const px = ((event.clientX - rect.left) / rect.width) * 800;
    const py = ((event.clientY - rect.top) / rect.height) * 800;
    let best = null;
    let bestDistance = Infinity;
    this.arrows.forEach(arrow => {
      if (arrow.dying || !arrow.spec.hover) return;
      const d = segmentDistance(px, py, arrow.target);
      const reach = arrow.target.width / 2 + 12;
      if (d < reach && d < bestDistance) { bestDistance = d; best = arrow; }
    });
    this.setHovered(best ? best.id : null, event);
  }

  setHovered(id, event) {
    if (id === this.hovered) {
      if (id && event) this.handlers.onArrowHover && this.handlers.onArrowHover(this.arrows.get(id).spec, event);
      return;
    }
    if (this.hovered && this.arrows.get(this.hovered)) this.arrows.get(this.hovered).el.classList.remove('hovered');
    this.hovered = id;
    this.root.classList.toggle('arrow-hover', !!id);
    if (id) this.arrows.get(id).el.classList.add('hovered');
    this.handlers.onArrowHover && this.handlers.onArrowHover(id ? this.arrows.get(id).spec : null, event);
    this.kick();
  }
}

// Run once the browser has painted the current state (so a transition can start from it),
// with a timer for tabs that are not painting at all.
function afterPaint(fn) {
  let done = false;
  const run = () => { if (!done) { done = true; fn(); } };
  requestAnimationFrame(() => requestAnimationFrame(run));
  setTimeout(run, 60);
}

function segmentDistance(px, py, t) {
  const dx = t.tx - t.fx;
  const dy = t.ty - t.fy;
  const length2 = dx * dx + dy * dy || 1;
  const k = Math.max(0, Math.min(1, ((px - t.fx) * dx + (py - t.fy) * dy) / length2));
  return Math.hypot(px - (t.fx + k * dx), py - (t.fy + k * dy));
}

function arrowPath(c, spec) {
  const dx = c.tx - c.fx;
  const dy = c.ty - c.fy;
  const length = Math.hypot(dx, dy);
  if (length < 1) return '';
  const ux = dx / length;
  const uy = dy / length;
  const nx = -uy;
  const ny = ux;
  const w = c.width / 2;
  const head = Math.min(c.width * 2.3, length * 0.55);
  const headW = c.width * 1.45;
  const start = Math.min(spec.inset ?? 22, length * 0.25);
  const tip = spec.tipInset ?? 10;
  const sx = c.fx + ux * start;
  const sy = c.fy + uy * start;
  const ex = c.tx - ux * tip;
  const ey = c.ty - uy * tip;
  const bx = ex - ux * head;
  const by = ey - uy * head;
  const p = (x, y) => `${x.toFixed(1)} ${y.toFixed(1)}`;
  return [
    `M${p(sx + nx * w, sy + ny * w)}`,
    `L${p(bx + nx * w, by + ny * w)}`,
    `L${p(bx + nx * headW, by + ny * headW)}`,
    `L${p(ex, ey)}`,
    `L${p(bx - nx * headW, by - ny * headW)}`,
    `L${p(bx - nx * w, by - ny * w)}`,
    `L${p(sx - nx * w, sy - ny * w)}`,
    // rounded tail
    `A${w.toFixed(1)} ${w.toFixed(1)} 0 0 1 ${p(sx + nx * w, sy + ny * w)}`,
    'Z'
  ].join(' ');
}
