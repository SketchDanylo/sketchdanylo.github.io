import { Chess, validateFen } from './lib/chess.js?v=2026-10-08-a';
import { Board, parseFen, pieceImage } from './board.js?v=2026-10-08-a';
import { Engine, winPercent, formatScore } from './engine.js?v=2026-10-08-a';
import { lookup, bookMoves, family, ecoGroup } from './book.js?v=2026-10-08-a';
import {
  classify, classLabel, explainMove, playLine, plyLabels, lineText, stepNote, pressure, other, moveFromUci
} from './explain.js?v=2026-10-08-a';

const START = 'rnbqkbnr/pppppppp/8/8/8/8/PPPPPPPP/RNBQKBNR w KQkq - 0 1';
const VIEW_DEPTH = 22;
const COLORS = {
  book: '#c4955f',
  best: '#7fd6a0',
  mine: '#9db8ff',
  theirs: '#f0827a',
  threat: '#f0827a',
  guard: '#8fdcb0'
};
const CLASS_ICON = { book: '♜', best: '★', excellent: '!', good: '✓', inaccuracy: '?!', mistake: '?', blunder: '??' };

const $ = id => document.getElementById(id);
const state = {
  nodes: [{ fen: START }],
  index: 0,
  mode: 'both',
  orientation: 'w',
  show: { book: true, best: true, threats: false },
  explore: null,
  editing: null
};

// ---------- engine ----------

const engine = new Engine();
const engineState = $('engine-state');
function setEngineState(kind, text) {
  engineState.className = `engine-state ${kind}`;
  engineState.querySelector('span').textContent = text;
}
if (engine.failed) setEngineState('error', 'Engine unavailable');

engine.on(message => {
  if (message.type === 'ready') {
    setEngineState('idle', 'Engine ready');
    analyzeView();
  } else if (message.type === 'error') {
    setEngineState('error', 'Engine unavailable');
  } else if (message.type === 'info') {
    if (message.fen === viewFen()) setEngineState('thinking', `Thinking · depth ${message.entry.depth}`);
    scheduleSoft();
    if (state.waitingEngine === message.fen) tryEngineReply();
  } else if (message.type === 'done') {
    if (message.fen === viewFen()) setEngineState('idle', `Depth ${engine.entry(message.fen)?.depth ?? '–'}`);
    scheduleSoft();
    if (state.waitingEngine === message.fen) tryEngineReply(true);
  }
});

function analyzeView() {
  const fen = viewFen();
  if (terminal(fen) || state.editing) return;
  engine.analyze(fen, { depth: VIEW_DEPTH, multipv: 3 });
}

// ---------- position helpers ----------

function chessAt(fen) {
  const chess = new Chess();
  chess.load(fen, { skipValidation: true });
  return chess;
}

const terminalCache = new Map();
function terminal(fen) {
  if (terminalCache.has(fen)) return terminalCache.get(fen);
  const chess = chessAt(fen);
  let result = null;
  if (chess.isCheckmate()) result = { cp: null, mate: 0, mated: chess.turn(), text: chess.turn() === 'w' ? 'Black wins by checkmate' : 'White wins by checkmate' };
  else if (chess.isStalemate()) result = { cp: 0, mate: null, text: 'Stalemate' };
  else if (chess.isInsufficientMaterial()) result = { cp: 0, mate: null, text: 'Draw · insufficient material' };
  terminalCache.set(fen, result);
  return result;
}

function evalOf(fen) {
  const end = terminal(fen);
  if (end) return { score: end, pv: [], depth: 99, final: true };
  const entry = engine.entry(fen);
  if (!entry || !entry.lines[0]) return null;
  const line = entry.lines[0];
  return { score: { cp: line.cp, mate: line.mate }, pv: line.pv, depth: entry.depth };
}

function viewFen() {
  if (state.editing) return state.editing.fen;
  if (state.explore) return currentLine()[state.explore.step].fen;
  return state.nodes[state.index].fen;
}

function turnOf(fen) { return fen.split(' ')[1]; }

function legalMap(fen) {
  const map = new Map();
  for (const move of chessAt(fen).moves({ verbose: true })) {
    if (!map.has(move.from)) map.set(move.from, []);
    map.get(move.from).push(move.to);
  }
  return map;
}

function sanOf(fen, uci) {
  const move = moveFromUci(chessAt(fen), uci);
  return move ? move.san : uci;
}

function moveLabel(fen, uci) {
  const [, color, , , , full] = fen.split(' ');
  return `${full || 1}${color === 'w' ? '.' : '…'}${sanOf(fen, uci)}`;
}

function userColor() { return state.mode === 'both' ? null : state.mode; }
function isUserMove(node) { return !!node.move && (state.mode === 'both' || node.move.color === state.mode); }

// ---------- board ----------

const board = new Board($('board'), {
  onMove: (from, to) => handleBoardMove(from, to),
  onEdit: action => editAction(action),
  onPaletteClick: code => setBrush(state.editing && state.editing.brush === code ? null : code),
  onArrowHover: (spec, event) => {
    if (spec && spec.hover) showBookCard(spec.hover.fen, spec.hover.uci, event);
    else hideCard();
  }
});

async function handleBoardMove(from, to) {
  const fen = state.nodes[state.index].fen;
  const chess = chessAt(fen);
  const piece = chess.get(from);
  let promotion;
  if (piece && piece.type === 'p' && (to[1] === '8' || to[1] === '1')) {
    promotion = await board.askPromotion(to, piece.color);
    if (!promotion) { board.setPosition(fen); return; }
  }
  play(from + to + (promotion || ''));
}

function play(uci, { byEngine = false } = {}) {
  const before = state.nodes[state.index].fen;
  const chess = chessAt(before);
  const move = moveFromUci(chess, uci);
  if (!move) return false;
  hideCard(true);
  state.nodes = state.nodes.slice(0, state.index + 1);
  state.nodes.push({ fen: chess.fen(), move, uci: move.from + move.to + (move.promotion || ''), byEngine });
  state.index++;
  // Feedback needs a baseline for the position the move was played from.
  engine.prioritize(before, 12);
  refresh();
  maybeEngineMove();
  return true;
}

function goTo(index) {
  index = Math.max(0, Math.min(state.nodes.length - 1, index));
  if (index === state.index) return;
  state.index = index;
  clearTimeout(state.engineTimer);
  state.waitingEngine = null;
  hideCard(true);
  refresh();
}

// Full update: board position plus everything around it.
function refresh() {
  const fen = viewFen();
  board.setPosition(fen);
  const lastMove = state.explore
    ? (currentLine()[state.explore.step].move || null)
    : state.nodes[state.index].move || null;
  const chess = chessAt(fen);
  let check = null;
  if (chess.inCheck()) {
    for (const row of chess.board()) for (const p of row) if (p && p.type === 'k' && p.color === chess.turn()) check = p.square;
  }
  board.setHighlights({ lastMove, check });
  board.setMovable(canMove(fen) ? legalMap(fen) : null);
  analyzeView();
  softUpdate();
  renderMoveList();
  updateNavButtons();
}

function canMove(fen) {
  if (state.explore || state.editing || terminal(fen)) return false;
  if (state.mode === 'both') return true;
  return turnOf(fen) === state.mode;
}

// Cheap update after engine output: arrows, eval bar, lines, feedback.
// Batched to one update per frame; the timer covers tabs that are not painting.
let softPending = false;
function scheduleSoft() {
  if (softPending) return;
  softPending = true;
  const run = () => {
    if (!softPending) return;
    softPending = false;
    softUpdate();
  };
  requestAnimationFrame(run);
  setTimeout(run, 120);
}
function softUpdate() {
  if (state.editing) return;
  renderArrows();
  renderEvalBar();
  renderLines();
  if (state.explore) renderExplore();
  else {
    renderOpening();
    renderFeedback();
    updateMoveDots();
  }
}

// ---------- arrows & marks ----------

function renderArrows() {
  const fen = viewFen();
  const arrows = [];
  let marks = [];
  if (state.explore) {
    const line = currentLine();
    const next = line[state.explore.step + 1];
    if (next) {
      const mine = next.move.color === state.explore.color;
      arrows.push({ id: 'line', from: next.move.from, to: next.move.to, color: mine ? COLORS.mine : COLORS.theirs, width: 15, opacity: 0.88, z: 3 });
    }
    const overlay = threatOverlay(fen, state.explore.color, true);
    arrows.push(...overlay.arrows);
    marks = overlay.marks;
  } else {
    const entry = engine.entry(fen);
    const best = state.show.best && !terminal(fen) && entry && entry.lines[0] ? entry.lines[0].pv[0] : null;
    if (state.show.book) {
      const all = bookMoves(fen);
      const total = all.reduce((sum, move) => sum + move.count, 0) || 1;
      const moves = all.slice(0, 6);
      moves.forEach(move => {
        const share = move.count / total;
        const same = best === move.uci;
        arrows.push({
          id: `book:${move.uci}`,
          from: move.uci.slice(0, 2), to: move.uci.slice(2, 4),
          color: COLORS.book,
          width: same ? 26 : 9 + 10 * Math.sqrt(share),
          opacity: same ? 0.7 : 0.5 + 0.38 * Math.sqrt(share),
          className: 'book',
          tipInset: same ? 2 : 10,
          hover: { fen, uci: move.uci },
          z: 1
        });
      });
    }
    if (best) {
      // Opacity grows with depth: a shallow guess is faint, a deep answer is solid.
      const confidence = Math.min(1, entry.depth / 18);
      arrows.push({ id: 'best', from: best.slice(0, 2), to: best.slice(2, 4), color: COLORS.best, width: 14, opacity: 0.45 + 0.45 * confidence, className: 'best', z: 2 });
    }
    if (state.show.threats && !terminal(fen)) {
      const overlay = threatOverlay(fen, turnOf(fen), false);
      arrows.push(...overlay.arrows);
      marks = overlay.marks;
    }
  }
  board.setArrows(arrows);
  board.setMarks(marks);
}

// Attacks on `color`'s pieces (red), the pieces shielding them (green), and loose enemy pieces.
function threatOverlay(fen, color, withTargets) {
  const arrows = [];
  const marks = [];
  const mine = pressure(fen, color);
  let threatCount = 0;
  let guardCount = 0;
  const VAL = { p: 1, n: 3, b: 3, r: 5, q: 9, k: 100 };
  mine.sort((a, b) => VAL[b.piece.type] - VAL[a.piece.type]);
  for (const item of mine) {
    if (item.hanging) {
      marks.push({ type: 'hanging', square: item.square, title: item.check ? 'In check' : `Your ${item.piece.type} is under-protected` });
      for (const from of item.attackers.slice(0, 2)) {
        if (threatCount++ >= 6) break;
        arrows.push({ id: `threat:${from}${item.square}`, from, to: item.square, color: COLORS.threat, width: 7, opacity: 0.8, className: 'threat', inset: 18, tipInset: 22, z: 4 });
      }
    } else {
      marks.push({ type: 'shield', square: item.square, count: item.defenders.length, title: `Attacked ${item.attackers.length}×, defended ${item.defenders.length}×` });
      for (const from of item.attackers.slice(0, 1)) {
        if (threatCount++ >= 6) break;
        arrows.push({ id: `threat:${from}${item.square}`, from, to: item.square, color: COLORS.threat, width: 5, opacity: 0.45, className: 'threat', inset: 18, tipInset: 22, z: 4 });
      }
      for (const from of item.defenders.slice(0, 2)) {
        if (guardCount++ >= 6) break;
        arrows.push({ id: `guard:${from}${item.square}`, from, to: item.square, color: COLORS.guard, width: 5, opacity: 0.55, className: 'guard', inset: 18, tipInset: 24, z: 4 });
      }
    }
  }
  if (withTargets) {
    for (const item of pressure(fen, other(color))) {
      if (item.hanging && !item.check) marks.push({ type: 'target', square: item.square, title: 'Loose enemy piece' });
    }
  }
  return { arrows, marks };
}

// ---------- eval bar ----------

const evalBar = $('evalbar');
function renderEvalBar() {
  const fen = viewFen();
  const result = evalOf(fen);
  const fill = $('eval-fill');
  const top = $('eval-top');
  const bottom = $('eval-bottom');
  evalBar.classList.toggle('flipped', state.orientation === 'b');
  if (!result || state.editing) {
    fill.style.height = '50%';
    top.textContent = bottom.textContent = '';
    return;
  }
  const white = winPercent(result.score);
  const shown = result.score.mate !== null && result.score.mate !== undefined
    ? (scoreWhiteWins(result.score) ? 100 : 0)
    : Math.max(4, Math.min(96, white));
  fill.style.height = `${shown}%`;
  const text = formatScore(result.score, { signed: false }).replace('−', '').replace('-', '');
  const whiteBetter = scoreWhiteWins(result.score) || (result.score.mate == null && result.score.cp >= 0);
  // The label sits at the end of the side that is better.
  const whiteAtBottom = state.orientation === 'w';
  const atBottom = whiteBetter === whiteAtBottom;
  const label = result.score.mate === 0 ? (result.score.mated === 'w' ? '0–1' : '1–0') : text;
  top.textContent = atBottom ? '' : label;
  bottom.textContent = atBottom ? label : '';
  top.className = `eval-label top ${whiteBetter ? 'dark' : 'light'}`;
  bottom.className = `eval-label bottom ${whiteBetter ? 'dark' : 'light'}`;
}
function scoreWhiteWins(score) {
  if (score.mate === 0) return score.mated === 'b';
  return score.mate != null && score.mate > 0;
}

// ---------- engine lines ----------

function renderLines() {
  const fen = viewFen();
  const list = $('lines');
  const entry = engine.entry(fen);
  const end = terminal(fen);
  const depthText = $('depth-text');
  const depthFill = $('depth-fill');
  if (end) {
    list.innerHTML = `<li class="placeholder">${end.text}.</li>`;
    depthText.textContent = 'final';
    depthFill.style.width = '100%';
    return;
  }
  if (!entry || !entry.lines.length) {
    list.innerHTML = `<li class="placeholder">${engine.failed ? 'The engine could not start in this browser.' : 'Analysing…'}</li>`;
    depthText.textContent = 'depth –';
    depthFill.style.width = '0%';
    return;
  }
  depthText.textContent = `depth ${entry.depth}`;
  depthFill.style.width = `${Math.min(100, (entry.depth / VIEW_DEPTH) * 100)}%`;
  const seen = new Set();
  const rows = [];
  for (const line of entry.lines) {
    if (!line || seen.has(line.pv[0])) continue;
    seen.add(line.pv[0]);
    const steps = playLine(fen, line.pv, 9);
    if (steps.length < 2) continue;
    const text = lineText(fen, steps.slice(1));
    const [first, ...rest] = text.split(' ');
    const score = { cp: line.cp, mate: line.mate };
    const negative = (score.mate != null ? score.mate < 0 : score.cp < 0);
    rows.push(`<li data-uci="${line.pv[0]}" class="${line.depth < entry.depth - 2 ? 'stale' : ''}"><span class="score ${negative ? 'neg' : ''}">${formatScore(score)}</span><span class="pv"><b>${first}</b> ${rest.join(' ')}</span></li>`);
  }
  const html = rows.join('');
  if (list._html !== html) { list.innerHTML = html; list._html = html; }
}

$('lines').addEventListener('click', event => {
  const row = event.target.closest('li[data-uci]');
  if (!row || state.explore) return;
  if (canMove(viewFen())) play(row.dataset.uci);
});

// ---------- opening card ----------

function openingContext(index) {
  let named = null;
  let leftAt = null;
  for (let i = 0; i <= index; i++) {
    const record = lookup(state.nodes[i].fen);
    if (record && record.exact) named = record.exact;
    if (!record && leftAt === null && i > 0) leftAt = i;
    if (record) leftAt = null;
  }
  const here = lookup(state.nodes[index].fen);
  return { named, here, leftAt };
}

function renderOpening() {
  const body = $('opening-body');
  const { named, here, leftAt } = openingContext(state.index);
  const fen = state.nodes[state.index].fen;
  const moves = here ? here.moves.slice(0, 6) : [];
  const total = here ? here.moves.reduce((sum, move) => sum + move.count, 0) : 1;
  const chips = moves.length
    ? `<div class="next-book"><span class="next-label">Theory continues</span>${moves.map(move => {
      const share = Math.round((move.count / total) * 100);
      return `<button class="bk" data-uci="${move.uci}">${sanOf(fen, move.uci)}<small>${share < 1 ? '<1' : share}%</small></button>`;
    }).join('')}</div>`
    : '';
  let html;
  let key;
  if (!named) {
    const fresh = state.index === 0 && fen === START;
    key = fresh ? 'start' : 'custom';
    html = fresh
      ? `<div class="op-top"><span class="eco">—</span><span class="op-status inbook">Starting position</span></div>
         <p class="op-name">Pick an opening</p>
         <p class="op-summary">Brown arrows are book moves: thicker means more theory behind it. Hover one to see where it leads. The green arrow is the engine’s choice and firms up as it searches deeper.</p>${chips}`
      : `<div class="op-top"><span class="eco">—</span><span class="op-status">${here ? 'In book' : 'Custom position'}</span></div>
         <p class="op-name">${here ? 'Known position' : 'No opening name'}</p>
         <p class="op-summary">${here ? 'This position appears in opening theory.' : 'This position is not in the opening book. The engine still analyses it.'}</p>${chips}`;
  } else {
    const info = family(named.family);
    const group = ecoGroup(named.eco);
    const status = here
      ? '<span class="op-status inbook">In book</span>'
      : `<span class="op-status">Out of book since ${leftAt ? moveLabel(state.nodes[leftAt - 1].fen, state.nodes[leftAt].uci) : 'here'}</span>`;
    key = `op:${named.index}`;
    html = `<div class="op-top"><span class="eco">${named.eco}</span>${status}</div>
      <p class="op-name">${named.family}</p>
      ${named.variation ? `<p class="op-var">${named.variation}</p>` : ''}
      ${info ? `<div class="tags">${info.tags.map(tag => `<span class="tag">${tag}</span>`).join('')}</div>` : ''}
      <p class="op-summary">${info ? info.summary : group || ''}</p>
      ${info ? `<ul class="ideas"><li><span class="side-dot w"></span><span>${info.white}</span></li><li><span class="side-dot b"></span><span>${info.black}</span></li></ul>` : ''}
      ${chips}`;
  }
  morph(body, html, key);
}

$('opening-body').addEventListener('click', event => {
  const chip = event.target.closest('.bk');
  if (!chip) return;
  hideCard(true);
  if (canMove(viewFen())) play(chip.dataset.uci);
});
$('opening-body').addEventListener('pointerover', event => {
  const chip = event.target.closest('.bk');
  if (!chip || event.pointerType === 'touch') return;
  const rect = chip.getBoundingClientRect();
  showBookCard(state.nodes[state.index].fen, chip.dataset.uci, { clientX: rect.left, clientY: rect.bottom + 4 }, { anchor: 'below' });
});
$('opening-body').addEventListener('pointerout', event => {
  const chip = event.target.closest('.bk');
  if (chip && !chip.contains(event.relatedTarget)) hideCard();
});

// ---------- hover card ----------

const card = $('hovercard');
const cardBody = $('hovercard-body');
let cardKey = null;
let cardHideTimer = 0;

function bookCardHtml(fen, uci) {
  const chess = chessAt(fen);
  const move = moveFromUci(chess, uci);
  if (!move) return null;
  const next = lookup(chess.fen());
  const target = next ? (next.exact || next.representative) : null;
  const moves = bookMoves(fen);
  const total = moves.reduce((sum, item) => sum + item.count, 0) || 1;
  const mine = moves.find(item => item.uci === uci);
  const share = mine ? mine.count / total : 0;
  const info = target ? family(target.family) : null;
  const lineSteps = target ? playLine(START, target.moves, 30) : [];
  const lineStr = lineSteps.length > 1 ? lineText(START, lineSteps.slice(1)) : '';
  const entry = engine.entry(fen);
  const engineLine = entry && entry.lines.find(line => line && line.pv[0] === uci);
  const isBest = entry && entry.lines[0] && entry.lines[0].pv[0] === uci;
  const descendants = mine ? mine.count : 0;
  return `<div class="hc-kicker"><span class="move">${moveLabel(fen, uci)}</span>Book move${target ? ` · ${target.eco}` : ''}</div>
    <p class="hc-name">${target ? target.family : 'Opening theory'}</p>
    ${target && target.variation ? `<p class="hc-var">${target.variation}</p>` : ''}
    ${lineStr ? `<p class="hc-line">${lineStr}</p>` : ''}
    ${info ? `<p class="hc-summary">${info.summary}</p>` : (target ? `<p class="hc-summary">${ecoGroup(target.eco) || ''}</p>` : '')}
    <div class="hc-bar"><div class="hc-bar-label"><span>Share of book lines here</span><span>${Math.round(share * 100)}%</span></div>
      <div class="hc-bar-track"><div class="hc-bar-fill" style="width:${Math.max(3, share * 100)}%"></div></div></div>
    <div class="hc-engine">${descendants} named line${descendants === 1 ? '' : 's'} go through it${engineLine ? ` · engine ${formatScore({ cp: engineLine.cp, mate: engineLine.mate })}` : ''}${isBest ? ' · <span style="color:var(--best)">engine’s pick</span>' : ''}</div>`;
}

function showBookCard(fen, uci, event, { anchor = 'pointer' } = {}) {
  clearTimeout(cardHideTimer);
  const key = `${fen}|${uci}`;
  const wasShown = card.classList.contains('show');
  if (key !== cardKey) {
    const html = bookCardHtml(fen, uci);
    if (!html) return;
    cardKey = key;
    swapCardContent(html, wasShown);
  }
  positionCard(event, anchor, !wasShown);
  card.classList.add('show');
}

function swapCardContent(html, animate) {
  const old = cardBody.firstElementChild;
  const startHeight = cardBody.offsetHeight;
  const content = document.createElement('div');
  content.className = 'hc-content';
  content.innerHTML = html;
  if (!animate || !old) {
    cardBody.replaceChildren(content);
    cardBody.style.height = '';
    return;
  }
  // Cross-fade and let the height glide to the new content.
  old.classList.add('swap');
  clearTimeout(cardBody._swap);
  cardBody._swap = setTimeout(() => {
    content.classList.add('swap');
    cardBody.replaceChildren(content);
    cardBody.style.height = `${startHeight}px`;
    const target = content.offsetHeight + 32;
    requestAnimationFrame(() => {
      cardBody.style.height = `${target}px`;
      content.classList.remove('swap');
    });
    setTimeout(() => { cardBody.style.height = ''; }, 360);
  }, 110);
}

function positionCard(event, anchor, snap) {
  const width = card.offsetWidth || 320;
  const height = cardBody.offsetHeight || 220;
  let x = event.clientX + 22;
  let y = event.clientY + 18;
  if (anchor === 'below') { x = event.clientX; y = event.clientY; }
  if (x + width > window.innerWidth - 12) x = event.clientX - width - 22;
  if (x < 12) x = 12;
  if (y + height > window.innerHeight - 12) y = Math.max(12, window.innerHeight - height - 12);
  card.classList.toggle('snap', snap);
  card.style.setProperty('--x', `${Math.round(x)}px`);
  card.style.setProperty('--y', `${Math.round(y)}px`);
  if (snap) requestAnimationFrame(() => card.classList.remove('snap'));
}

function hideCard(now = false) {
  clearTimeout(cardHideTimer);
  const hide = () => { card.classList.remove('show'); cardKey = null; };
  if (now) hide();
  else cardHideTimer = setTimeout(hide, 140);
}

// ---------- feedback ----------

function judge(i) {
  const node = state.nodes[i];
  const prev = state.nodes[i - 1];
  if (!node || !prev || !node.move) return null;
  const inBook = bookMoves(prev.fen).some(move => move.uci === node.uci);
  const beforeEval = evalOf(prev.fen);
  const afterEval = evalOf(node.fen);
  const before = beforeEval && beforeEval.depth >= 10 ? beforeEval : null;
  const after = afterEval && afterEval.depth >= 10 ? afterEval : null;
  // A named line is only "book" while it is also sound: the database names plenty of dubious tries.
  const engineVerdict = classify({ inBook: false, playedUci: node.uci, before, after, color: node.move.color });
  const extra = { before: beforeEval, after: afterEval, inBook };
  if (inBook && (!engineVerdict || engineVerdict.drop < 8)) return { kind: 'book', drop: 0, ...extra };
  return engineVerdict ? { ...engineVerdict, ...extra } : { kind: 'pending', ...extra };
}

function feedbackIndex() {
  for (let i = state.index; i > 0; i--) {
    if (isUserMove(state.nodes[i])) return i;
  }
  return 0;
}

function renderFeedback() {
  const body = $('feedback-body');
  const cardEl = $('feedback-card');
  const i = feedbackIndex();
  const end = terminal(state.nodes[state.index].fen);
  if (i === 0) {
    cardEl.dataset.kind = '';
    morph(body, `<div class="fb"><div class="fb-icon k-book">♞</div><div>
      <p class="fb-title">Your move</p>
      <p class="fb-text">Play on the board. After every move you get a verdict; when it isn’t the engine’s top choice, <b>Why?</b> walks you through what goes wrong.</p></div></div>`, 'intro');
    return;
  }
  const node = state.nodes[i];
  const verdict = judge(i);
  const kind = verdict.kind;
  cardEl.dataset.kind = kind;
  const label = moveLabel(state.nodes[i - 1].fen, node.uci);
  const bestUci = verdict.before && verdict.before.pv[0];
  const bestSan = bestUci ? sanOf(state.nodes[i - 1].fen, bestUci) : null;
  const color = node.move.color;
  const forMover = score => score ? formatScore(color === 'w' ? score : flip(score)) : '…';
  let title;
  let text;
  let actions = '';
  if (kind === 'pending') {
    title = `${label}`;
    text = 'Checking the move…';
  } else if (kind === 'book') {
    const { named } = openingContext(i);
    title = `${label} is book`;
    text = named ? `Theory: ${named.name}.` : 'A known opening move.';
  } else if (kind === 'best') {
    title = `${label} is the best move`;
    text = `The engine’s top choice${verdict.after ? ` (${forMover(verdict.after.score)} for you)` : ''}.`;
  } else {
    const verb = { excellent: 'is excellent', good: 'is good', inaccuracy: 'is an inaccuracy', mistake: 'is a mistake', blunder: 'is a blunder' }[kind];
    title = `${label} ${verb}`;
    const { named } = verdict.inBook ? openingContext(i) : { named: null };
    text = (named ? `It has a name — ${named.name} — but it isn’t sound. ` : '') + (bestSan
      ? `The engine preferred <b>${bestSan}</b>. ${forMover(verdict.before.score)} → ${forMover(verdict.after.score)} for you${verdict.drop >= 1 ? ` — ${Math.round(verdict.drop)}% less winning chance` : ''}.`
      : '');
    const serious = ['inaccuracy', 'mistake', 'blunder'].includes(kind);
    actions = `<div class="fb-actions">
      <button class="pill why" data-act="why" data-index="${i}">${serious ? 'Why?' : 'Compare with best'}</button>
      ${serious ? `<button class="pill" data-act="retry" data-index="${i}">Retry</button>` : ''}
    </div>`;
  }
  if (end) text += ` <b>${end.text}.</b>`;
  const icon = kind === 'pending' ? '' : CLASS_ICON[kind];
  const html = `<div class="fb"><div class="fb-icon k-${kind}">${icon}</div><div>
    <p class="fb-title">${title}</p><p class="fb-text">${text}</p>${actions}</div></div>`;
  morph(body, html, `fb:${i}:${kind}`);
}

function flip(score) {
  if (!score) return score;
  if (score.mate != null) return { ...score, mate: score.mate === 0 ? 0 : -score.mate };
  return { cp: -score.cp, mate: null };
}

$('feedback-body').addEventListener('click', event => {
  const button = event.target.closest('button[data-act]');
  if (!button) return;
  const i = Number(button.dataset.index);
  if (button.dataset.act === 'why') openExplore(i);
  if (button.dataset.act === 'retry') retry(i);
});

function retry(i) {
  state.nodes = state.nodes.slice(0, i);
  state.index = i - 1;
  closeExplore(false);
  toast('Try another move');
  refresh();
}

// ---------- move list ----------

function renderMoveList() {
  const list = $('movelist');
  const nodes = state.nodes;
  const signature = nodes.map(node => node.uci || '').join(' ') + (nodes[0].fen === START ? '' : nodes[0].fen);
  if (list._sig !== signature) {
    list._sig = signature;
    if (nodes.length === 1) {
      list.innerHTML = '<li class="empty">No moves yet.</li>';
    } else {
      const first = nodes[0].fen.split(' ');
      let number = Number(first[5] || 1);
      let html = '';
      let color = first[1];
      if (color === 'b') html += `<li class="num">${number}.</li><li></li>`;
      for (let i = 1; i < nodes.length; i++) {
        if (color === 'w') html += `<li class="num">${number}.</li>`;
        html += `<li class="mv" data-i="${i}"><span class="dot"></span>${nodes[i].move.san}</li>`;
        if (color === 'b') number++;
        color = other(color);
      }
      // Keep existing rows still; only the new ones animate in.
      const old = list.querySelectorAll('.mv').length;
      list.innerHTML = html;
      list.querySelectorAll('.mv').forEach((el, k) => { if (k < old) el.style.animation = 'none'; });
    }
  }
  list.querySelectorAll('.mv').forEach(el => el.classList.toggle('cur', Number(el.dataset.i) === state.index));
  scrollWithin(list, list.querySelector('.mv.cur'));
  updateMoveDots();
}

function updateMoveDots() {
  $('movelist').querySelectorAll('.mv').forEach(el => {
    const i = Number(el.dataset.i);
    const dot = el.querySelector('.dot');
    const node = state.nodes[i];
    if (!node || !isUserMove(node)) { dot.classList.remove('shown'); return; }
    const verdict = judge(i);
    const kind = verdict && verdict.kind !== 'pending' ? verdict.kind : null;
    dot.classList.toggle('shown', !!kind);
    dot.style.background = kind ? `var(--${kind})` : '';
    el.title = kind ? classLabel(kind) : '';
  });
}

$('movelist').addEventListener('click', event => {
  const row = event.target.closest('.mv');
  if (!row) return;
  if (state.explore) closeExplore(false);
  goTo(Number(row.dataset.i));
});

// ---------- explore ("why?") ----------

function currentLine() {
  return state.explore.lines[state.explore.tab];
}

function openExplore(i) {
  const node = state.nodes[i];
  const prev = state.nodes[i - 1];
  const verdict = judge(i);
  const before = verdict.before;
  const after = verdict.after;
  const playedPv = [node.uci, ...((after && after.pv) || []).slice(0, 11)];
  const bestPv = ((before && before.pv) || []).slice(0, 12);
  const history = state.nodes.slice(1, i).map(n => n.move);
  const info = explainMove({ fenBefore: prev.fen, move: node.move, best: before, after, history });
  state.explore = {
    index: i,
    color: node.move.color,
    kind: verdict.kind,
    info,
    tab: 'played',
    step: 1,
    lines: { played: playLine(prev.fen, playedPv, 12), best: playLine(prev.fen, bestPv, 12) }
  };
  // Quick evaluations for every step so the curve fills in.
  for (const line of Object.values(state.explore.lines)) for (const step of line) engine.prioritize(step.fen, 11);
  stopAutoplay();
  showView('explore');
  setSegmented($('explore-tabs'), 'played');
  $('explore-tabs').querySelector('[data-line="best"]').textContent = state.explore.lines.best[1] ? `Best: ${state.explore.lines.best[1].move.san}` : 'Best move';
  renderExploreHeader();
  refresh();
}

function closeExplore(doRefresh = true) {
  if (!state.explore) return;
  stopAutoplay();
  state.explore = null;
  $('explore-badge').classList.remove('show');
  showView('game');
  if (doRefresh) refresh();
}

function renderExploreHeader() {
  const ex = state.explore;
  const chip = $('explore-chip');
  chip.textContent = classLabel(ex.kind);
  chip.className = `chip k-${ex.kind}`;
  $('explore-title').textContent = ex.info.headline;
  $('explore-reasons').innerHTML = ex.info.points.map((point, k) => `<li style="animation-delay:${k * 60}ms">${point}</li>`).join('');
}

function renderExplore() {
  const ex = state.explore;
  const line = currentLine();
  const labels = plyLabels(line[0].fen, line.slice(1));
  const rows = line.slice(1).map((step, k) => {
    const i = k + 1;
    const result = evalOf(step.fen);
    const ev = result ? formatScore(ex.color === 'w' ? result.score : flip(result.score)) : '';
    const note = stepNote(line[i - 1].fen, step.move);
    const mine = step.move.color === ex.color;
    const cls = [i === ex.step ? 'cur' : '', i > ex.step ? 'future' : ''].join(' ');
    return `<li data-step="${i}" class="${cls}"><span class="san ${mine ? 'yours' : 'theirs'}">${labels[k]}</span><span class="note">${note}</span><span class="ev">${ev}</span></li>`;
  }).join('');
  const steps = $('steps');
  if (steps._html !== rows) {
    steps.innerHTML = rows;
    steps._html = rows;
    scrollWithin(steps, steps.querySelector('li.cur'));
  }
  renderSpark();
  const badge = $('explore-badge');
  const stepMove = line[ex.step].move;
  badge.textContent = ex.step === 0
    ? 'Before your move'
    : `${ex.tab === 'played' ? 'Your line' : 'Best line'} · ${labels[ex.step - 1]} · step ${ex.step}/${line.length - 1}`;
  badge.classList.add('show');
  $('step-back').disabled = ex.step === 0;
  $('step-forward').disabled = ex.step >= line.length - 1;
  void stepMove;
}

function renderSpark() {
  const ex = state.explore;
  const svg = $('spark');
  const width = 300;
  const height = 64;
  const maxSteps = Math.max(ex.lines.played.length, ex.lines.best.length) - 1 || 1;
  const pointsFor = line => line.map((step, i) => {
    const result = evalOf(step.fen);
    if (!result) return null;
    const win = winPercent(result.score);
    const mine = ex.color === 'w' ? win : 100 - win;
    return [(i / maxSteps) * width, height - (mine / 100) * height];
  });
  const path = points => {
    const known = points.filter(Boolean);
    if (known.length < 2) return '';
    return known.map((p, k) => `${k ? 'L' : 'M'}${p[0].toFixed(1)} ${p[1].toFixed(1)}`).join(' ');
  };
  const played = pointsFor(ex.lines.played);
  const best = pointsFor(ex.lines.best);
  const playedColor = `var(--${['blunder', 'mistake', 'inaccuracy'].includes(ex.kind) ? ex.kind : 'line'})`;
  const active = ex.tab === 'played' ? played : best;
  const cursor = active[ex.step];
  svg.innerHTML = `
    <line class="mid" x1="0" x2="${width}" y1="${height / 2}" y2="${height / 2}"/>
    <path class="curve ${ex.tab === 'played' ? 'faded' : ''}" d="${path(best)}" stroke="var(--best)"/>
    <path class="curve ${ex.tab === 'best' ? 'faded' : ''}" d="${path(played)}" stroke="${playedColor}"/>
    ${cursor ? `<circle class="cursor" r="4" cx="${cursor[0].toFixed(1)}" cy="${cursor[1].toFixed(1)}"/>` : ''}
    <text x="0" y="10" fill="var(--faint)" font-size="9">you win</text>
    <text x="0" y="${height - 2}" fill="var(--faint)" font-size="9">you lose</text>`;
}

function setStep(step) {
  const line = currentLine();
  state.explore.step = Math.max(0, Math.min(line.length - 1, step));
  refresh();
}

$('steps').addEventListener('click', event => {
  const row = event.target.closest('li[data-step]');
  if (row) { stopAutoplay(); setStep(Number(row.dataset.step)); }
});
$('step-back').addEventListener('click', () => { stopAutoplay(); setStep(state.explore.step - 1); });
$('step-forward').addEventListener('click', () => { stopAutoplay(); setStep(state.explore.step + 1); });
$('explore-exit').addEventListener('click', () => closeExplore());

let autoplay = 0;
$('step-play').addEventListener('click', () => {
  if (autoplay) { stopAutoplay(); return; }
  const line = currentLine();
  if (state.explore.step >= line.length - 1) setStep(0);
  $('step-play').textContent = 'Pause';
  autoplay = setInterval(() => {
    if (!state.explore || state.explore.step >= currentLine().length - 1) { stopAutoplay(); return; }
    setStep(state.explore.step + 1);
  }, 1100);
});
function stopAutoplay() {
  clearInterval(autoplay);
  autoplay = 0;
  $('step-play').textContent = 'Play line';
}

setupSegmented($('explore-tabs'), value => {
  stopAutoplay();
  state.explore.tab = value;
  state.explore.step = Math.min(state.explore.step, currentLine().length - 1);
  if (value === 'best' && state.explore.step === 0) state.explore.step = 1;
  refresh();
}, 'line');

// ---------- engine opponent ----------

function maybeEngineMove() {
  clearTimeout(state.engineTimer);
  state.waitingEngine = null;
  if (state.mode === 'both' || state.explore || state.editing) return;
  if (state.index !== state.nodes.length - 1) return;
  const fen = state.nodes[state.index].fen;
  if (turnOf(fen) === state.mode || terminal(fen)) return;
  const all = bookMoves(fen);
  // Stay in mainstream theory so the user meets real lines: skip the rare sidelines,
  // and favour moves with more theory behind them.
  const top = all.length ? all[0].count : 0;
  const book = all.filter(move => move.count >= top * 0.12);
  if (book.length) {
    const weight = move => Math.pow(move.count, 1.3);
    const total = book.reduce((sum, move) => sum + weight(move), 0);
    let pick = Math.random() * total;
    const choice = book.find(move => (pick -= weight(move)) < 0) || book[0];
    state.engineTimer = setTimeout(() => play(choice.uci, { byEngine: true }), 650);
    return;
  }
  state.waitingEngine = fen;
  state.engineAsked = performance.now();
  tryEngineReply();
}

function tryEngineReply(force = false) {
  const fen = state.waitingEngine;
  if (!fen || fen !== state.nodes[state.index].fen) return;
  const entry = engine.entry(fen);
  if (!entry || !entry.lines[0]) return;
  if (entry.depth < 14 && !force) return;
  state.waitingEngine = null;
  const wait = Math.max(0, 600 - (performance.now() - state.engineAsked));
  state.engineTimer = setTimeout(() => {
    if (state.nodes[state.index].fen === fen) play(entry.lines[0].pv[0], { byEngine: true });
  }, wait);
}

// ---------- editor ----------

function openEditor() {
  closeExplore(false);
  hideCard(true);
  clearTimeout(state.engineTimer);
  state.waitingEngine = null;
  const fen = state.nodes[state.index].fen;
  const parts = fen.split(' ');
  state.editing = { map: parseFen(fen), turn: parts[1], castling: parts[2], brush: null, fen, original: fen };
  board.edit = true;
  board.setMovable(null);
  board.setArrows([]);
  board.setMarks([]);
  board.setHighlights({});
  showView('edit');
  renderEvalBar();
  setSegmented($('edit-turn'), parts[1]);
  document.querySelectorAll('[data-castle]').forEach(box => { box.checked = parts[2].includes(box.dataset.castle); });
  $('edit-toggle').classList.add('on');
  setBrush(null);
  syncEditor();
}

function closeEditor() {
  state.editing = null;
  board.edit = false;
  board.brush = null;
  $('edit-toggle').classList.remove('on');
  showView('game');
}

function buildEditFen() {
  const ed = state.editing;
  let placement = '';
  for (let rank = 8; rank >= 1; rank--) {
    let empty = 0;
    for (const file of 'abcdefgh') {
      const code = ed.map.get(file + rank);
      if (!code) { empty++; continue; }
      if (empty) { placement += empty; empty = 0; }
      placement += code[0] === 'w' ? code[1] : code[1].toLowerCase();
    }
    if (empty) placement += empty;
    if (rank > 1) placement += '/';
  }
  // Only keep castling rights the pieces still allow.
  const has = (square, code) => ed.map.get(square) === code;
  const rights = [];
  document.querySelectorAll('[data-castle]').forEach(box => {
    const right = box.dataset.castle;
    const possible = {
      K: has('e1', 'wK') && has('h1', 'wR'),
      Q: has('e1', 'wK') && has('a1', 'wR'),
      k: has('e8', 'bK') && has('h8', 'bR'),
      q: has('e8', 'bK') && has('a8', 'bR')
    }[right];
    box.disabled = !possible;
    box.closest('.tg').style.opacity = possible ? '' : '0.4';
    if (possible && box.checked) rights.push(right);
  });
  return `${placement} ${ed.turn} ${rights.join('') || '-'} - 0 1`;
}

function editorProblem(fen) {
  const map = parseFen(fen);
  const count = code => [...map.values()].filter(value => value === code).length;
  if (count('wK') !== 1 || count('bK') !== 1) return 'Each side needs exactly one king.';
  for (const [square, code] of map) {
    if (code[1] === 'P' && (square[1] === '1' || square[1] === '8')) return 'Pawns can’t stand on the first or last rank.';
  }
  const check = validateFen(fen);
  if (!check.ok) return check.error;
  const parts = fen.split(' ');
  parts[1] = other(parts[1]);
  const swapped = chessAt(parts.join(' '));
  if (swapped.inCheck()) return 'The side that isn’t to move is in check.';
  return '';
}

function syncEditor({ fromInput = false } = {}) {
  const fen = buildEditFen();
  state.editing.fen = fen;
  board.setPosition(fen);
  if (!fromInput) $('fen-input').value = fen;
  const problem = editorProblem(fen);
  $('edit-error').textContent = problem;
  $('edit-done').disabled = !!problem;
  $('edit-done').style.opacity = problem ? '0.45' : '';
}

function editAction(action) {
  const ed = state.editing;
  if (!ed) return;
  if (action.type === 'put') ed.map.set(action.square, action.code);
  if (action.type === 'paint') {
    if (ed.map.get(action.square) === action.code) ed.map.delete(action.square);
    else ed.map.set(action.square, action.code);
  }
  if (action.type === 'remove') ed.map.delete(action.square);
  if (action.type === 'move') {
    const code = ed.map.get(action.from);
    if (code && action.from !== action.to) { ed.map.delete(action.from); ed.map.set(action.to, code); }
  }
  syncEditor();
}

function setBrush(code) {
  if (!state.editing) return;
  state.editing.brush = code;
  board.brush = code;
  document.querySelectorAll('#palette button').forEach(button => button.classList.toggle('on', button.dataset.code === code));
}

function buildPalette() {
  const palette = $('palette');
  for (const color of 'wb') {
    for (const type of 'KQRBNP') {
      const code = color + type;
      const button = document.createElement('button');
      button.dataset.code = code;
      button.style.backgroundImage = pieceImage(code);
      button.setAttribute('aria-label', code);
      button.addEventListener('pointerdown', event => {
        event.preventDefault();
        board.startExternalDrag(code, event);
      });
      palette.appendChild(button);
    }
  }
}

$('fen-input').addEventListener('input', event => {
  const value = event.target.value.trim();
  const parts = value.split(/\s+/);
  if (parts[0].split('/').length !== 8) { $('edit-error').textContent = 'That doesn’t look like a FEN yet.'; return; }
  try {
    state.editing.map = parseFen(value);
    state.editing.turn = parts[1] === 'b' ? 'b' : 'w';
    setSegmented($('edit-turn'), state.editing.turn);
    const rights = parts[2] || '-';
    document.querySelectorAll('[data-castle]').forEach(box => { box.checked = rights.includes(box.dataset.castle); });
    syncEditor({ fromInput: true });
  } catch (error) {
    $('edit-error').textContent = 'Couldn’t read that FEN.';
  }
});
setupSegmented($('edit-turn'), value => { state.editing.turn = value; syncEditor(); }, 'turn');
document.querySelectorAll('[data-castle]').forEach(box => box.addEventListener('change', () => syncEditor()));
$('edit-start').addEventListener('click', () => {
  state.editing.map = parseFen(START);
  state.editing.turn = 'w';
  setSegmented($('edit-turn'), 'w');
  document.querySelectorAll('[data-castle]').forEach(box => { box.checked = true; });
  syncEditor();
});
$('edit-clear').addEventListener('click', () => {
  state.editing.map = new Map([['e1', 'wK'], ['e8', 'bK']]);
  syncEditor();
});
$('edit-cancel').addEventListener('click', () => {
  closeEditor();
  refresh();
});
$('edit-done').addEventListener('click', () => {
  const fen = state.editing.fen;
  if (editorProblem(fen)) return;
  closeEditor();
  state.nodes = [{ fen }];
  state.index = 0;
  refresh();
  maybeEngineMove();
});
$('edit-toggle').addEventListener('click', () => {
  if (state.editing) { closeEditor(); refresh(); } else openEditor();
});

// ---------- views, controls, segmented ----------

function showView(name) {
  for (const view of ['game', 'explore', 'edit']) {
    $(`view-${view}`).hidden = view !== name;
  }
  document.querySelectorAll('.segmented').forEach(placeGlider);
}

function setupSegmented(el, onChange, attr) {
  el.addEventListener('click', event => {
    const button = event.target.closest('button');
    if (!button) return;
    const value = button.dataset[attr];
    setSegmented(el, value);
    onChange(value);
  });
}

function setSegmented(el, value) {
  el.querySelectorAll('button').forEach(button => {
    const on = Object.values(button.dataset).includes(value);
    button.classList.toggle('on', on);
    button.setAttribute('aria-checked', on ? 'true' : 'false');
  });
  placeGlider(el);
}

function placeGlider(el) {
  const glider = el.querySelector('.seg-glider');
  const on = el.querySelector('button.on');
  if (!glider || !on || !on.offsetWidth) return;
  glider.style.width = `${on.offsetWidth}px`;
  glider.style.transform = `translateX(${on.offsetLeft}px)`;
}
// Re-measure whenever a control changes size (fonts loading, view shown, phone rotated).
const gliderObserver = new ResizeObserver(entries => entries.forEach(entry => placeGlider(entry.target)));
document.querySelectorAll('.segmented').forEach(el => gliderObserver.observe(el));

setupSegmented($('mode'), value => {
  state.mode = value;
  if (value !== 'both') setOrientation(value);
  refresh();
  maybeEngineMove();
}, 'mode');

function setOrientation(color) {
  state.orientation = color;
  board.setOrientation(color);
  renderEvalBar();
}

$('flip').addEventListener('click', () => setOrientation(state.orientation === 'w' ? 'b' : 'w'));
$('nav-start').addEventListener('click', () => navigate('start'));
$('nav-back').addEventListener('click', () => navigate('back'));
$('nav-forward').addEventListener('click', () => navigate('forward'));
$('nav-end').addEventListener('click', () => navigate('end'));
$('new-game').addEventListener('click', () => {
  closeExplore(false);
  if (state.editing) closeEditor();
  state.nodes = [{ fen: START }];
  state.index = 0;
  refresh();
  maybeEngineMove();
});
$('copy-pgn').addEventListener('click', async () => {
  const chess = new Chess(state.nodes[0].fen);
  state.nodes.slice(1).forEach(node => chess.move(node.move.san));
  try {
    await navigator.clipboard.writeText(chess.pgn());
    toast('PGN copied');
  } catch (error) {
    toast('Couldn’t reach the clipboard');
  }
});

for (const key of ['book', 'best', 'threats']) {
  $(`tg-${key}`).addEventListener('change', event => {
    state.show[key] = event.target.checked;
    if (!state.show.book) hideCard(true);
    renderArrows();
  });
}

function navigate(where) {
  if (state.explore) {
    const line = currentLine();
    stopAutoplay();
    const step = { start: 0, back: state.explore.step - 1, forward: state.explore.step + 1, end: line.length - 1 }[where];
    setStep(step);
    return;
  }
  if (state.editing) return;
  const index = { start: 0, back: state.index - 1, forward: state.index + 1, end: state.nodes.length - 1 }[where];
  goTo(index);
}

function updateNavButtons() {
  const atStart = state.explore ? state.explore.step === 0 : state.index === 0;
  const atEnd = state.explore ? state.explore.step >= currentLine().length - 1 : state.index >= state.nodes.length - 1;
  $('nav-start').disabled = $('nav-back').disabled = atStart || !!state.editing;
  $('nav-end').disabled = $('nav-forward').disabled = atEnd || !!state.editing;
}

document.addEventListener('keydown', event => {
  if (event.target.closest('input, textarea')) return;
  if (event.key === 'ArrowLeft') { event.preventDefault(); navigate('back'); }
  else if (event.key === 'ArrowRight') { event.preventDefault(); navigate('forward'); }
  else if (event.key === 'Home') navigate('start');
  else if (event.key === 'End') navigate('end');
  else if (event.key === 'f' || event.key === 'F') setOrientation(state.orientation === 'w' ? 'b' : 'w');
  else if (event.key === 'Escape') {
    if (state.explore) closeExplore();
    else if (state.editing) { closeEditor(); refresh(); }
  }
});

// ---------- small helpers ----------

// Swap a panel's content with a cross-fade while its height glides to fit.
function morph(el, html, key) {
  let inner = el.firstElementChild;
  if (!inner) {
    inner = document.createElement('div');
    inner.className = 'morph-inner';
    el.appendChild(inner);
    inner.innerHTML = html;
    el._html = html;
    el.dataset.key = key;
    return;
  }
  if (el.dataset.key === key) {
    if (el._html === html || el._pending) return;
    const from = el.offsetHeight;
    inner.innerHTML = html;
    el._html = html;
    glideHeight(el, from);
    return;
  }
  el.dataset.key = key;
  el._html = html;
  const token = (el._token || 0) + 1;
  el._token = token;
  el._pending = true;
  el.style.height = `${el.offsetHeight}px`;
  el.classList.add('swapping');
  setTimeout(() => {
    if (el._token !== token) return;
    inner.innerHTML = el._html;
    el._pending = false;
    el.style.height = `${inner.offsetHeight}px`;
    el.classList.remove('swapping');
    setTimeout(() => { if (el._token === token) el.style.height = ''; }, 420);
  }, 170);
}

function glideHeight(el, from) {
  const to = el.firstElementChild.offsetHeight;
  if (Math.abs(to - from) < 1) return;
  el.style.height = `${from}px`;
  void el.offsetHeight;
  el.style.height = `${to}px`;
  clearTimeout(el._glide);
  el._glide = setTimeout(() => { el.style.height = ''; }, 420);
}

// Keep a row visible inside its own scroll box without moving the page.
function scrollWithin(box, row) {
  if (!row) return;
  const top = row.offsetTop - box.offsetTop;
  if (top < box.scrollTop) box.scrollTo({ top: top - 4, behavior: 'smooth' });
  else if (top + row.offsetHeight > box.scrollTop + box.clientHeight) box.scrollTo({ top: top + row.offsetHeight - box.clientHeight + 4, behavior: 'smooth' });
}

let toastTimer = 0;
function toast(text) {
  const el = $('toast');
  el.textContent = text;
  el.classList.add('show');
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => el.classList.remove('show'), 1800);
}

// Quiet starfield, same family as the home page.
function starfield() {
  const canvas = $('space-bg');
  const context = canvas.getContext('2d');
  const reduce = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  const stars = [];
  const pointer = { x: 0, y: 0, tx: 0, ty: 0 };
  let width = 0;
  let height = 0;
  function resize() {
    const ratio = Math.min(window.devicePixelRatio || 1, 2);
    width = window.innerWidth;
    height = window.innerHeight;
    canvas.width = width * ratio;
    canvas.height = height * ratio;
    context.setTransform(ratio, 0, 0, ratio, 0, 0);
    stars.length = 0;
    const count = Math.min(150, Math.max(70, Math.floor((width * height) / 11000)));
    for (let i = 0; i < count; i++) {
      stars.push({ x: Math.random() * width, y: Math.random() * height, r: 0.4 + Math.random() * 1.1, depth: 0.15 + Math.random() * 0.85, alpha: 0.12 + Math.random() * 0.45, twinkle: Math.random() * 6.28 });
    }
  }
  function draw(time) {
    pointer.x += (pointer.tx - pointer.x) * 0.05;
    pointer.y += (pointer.ty - pointer.y) * 0.05;
    context.clearRect(0, 0, width, height);
    for (const star of stars) {
      context.beginPath();
      context.fillStyle = `rgba(255,255,255,${Math.max(0.03, star.alpha + Math.sin(time * 0.0012 + star.twinkle) * 0.1)})`;
      context.arc(star.x + pointer.x * star.depth * 14, star.y + pointer.y * star.depth * 10, star.r, 0, 6.283);
      context.fill();
    }
    if (!reduce && !document.hidden) requestAnimationFrame(draw);
  }
  window.addEventListener('resize', resize);
  window.addEventListener('pointermove', event => {
    pointer.tx = (event.clientX / width - 0.5) * 2;
    pointer.ty = (event.clientY / height - 0.5) * 2;
  }, { passive: true });
  document.addEventListener('visibilitychange', () => { if (!document.hidden && !reduce) requestAnimationFrame(draw); });
  resize();
  requestAnimationFrame(draw);
}

// ---------- start ----------

buildPalette();
starfield();
board.setPosition(START, { instant: true });
refresh();
const placeAllGliders = () => document.querySelectorAll('.segmented').forEach(placeGlider);
setTimeout(placeAllGliders, 0);
if (document.fonts) document.fonts.ready.then(placeAllGliders);
