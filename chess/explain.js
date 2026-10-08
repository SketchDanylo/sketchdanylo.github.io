// Turns engine lines into reasons a person can follow:
// what the move left hanging, what the reply wins, which defender walked away.
import { Chess } from './lib/chess.js?v=2026-10-08-a';
import { winPercent, formatScore } from './engine.js?v=2026-10-08-a';

export const VALUE = { p: 1, n: 3, b: 3, r: 5, q: 9, k: 100 };
export const NAME = { p: 'pawn', n: 'knight', b: 'bishop', r: 'rook', q: 'queen', k: 'king' };
const SQUARES = [];
for (let r = 8; r >= 1; r--) for (const f of 'abcdefgh') SQUARES.push(f + r);

export function other(color) { return color === 'w' ? 'b' : 'w'; }

function load(fen) {
  const chess = new Chess();
  chess.load(fen, { skipValidation: true });
  return chess;
}

export function moveFromUci(chess, uci) {
  try {
    return chess.move({ from: uci.slice(0, 2), to: uci.slice(2, 4), promotion: uci[4] || 'q' });
  } catch (error) {
    return null;
  }
}

// Replays an engine line, keeping each position and the verbose move that led there.
export function playLine(fen, uciMoves, maxPlies = 12) {
  const chess = load(fen);
  const steps = [{ fen, move: null }];
  for (const uci of uciMoves.slice(0, maxPlies)) {
    const move = moveFromUci(chess, uci);
    if (!move) break;
    steps.push({ fen: chess.fen(), move });
  }
  return steps;
}

export function sanLine(fen, uciMoves, maxPlies = 10) {
  return playLine(fen, uciMoves, maxPlies).slice(1).map(step => step.move.san);
}

function material(chess) {
  const sum = { w: 0, b: 0 };
  for (const row of chess.board()) for (const piece of row) {
    if (piece && piece.type !== 'k') sum[piece.color] += VALUE[piece.type];
  }
  return sum;
}

// Every piece of `color` that the other side attacks, with attackers and defenders.
export function pressure(fen, color) {
  const chess = load(fen);
  const result = [];
  for (const square of SQUARES) {
    const piece = chess.get(square);
    if (!piece || piece.color !== color) continue;
    const attackers = chess.attackers(square, other(color));
    if (!attackers.length) continue;
    if (piece.type === 'k') {
      result.push({ square, piece, attackers, defenders: [], hanging: true, check: true });
      continue;
    }
    const defenders = chess.attackers(square, color);
    const cheapest = Math.min(...attackers.map(sq => VALUE[chess.get(sq).type]));
    const hanging = defenders.length === 0 || cheapest < VALUE[piece.type] || attackers.length > defenders.length;
    result.push({ square, piece, attackers, defenders, hanging, check: false });
  }
  return result;
}

function piecesAttackedBy(chess, square) {
  // Squares holding enemy pieces that the piece on `square` hits.
  const piece = chess.get(square);
  if (!piece) return [];
  const hits = [];
  for (const target of SQUARES) {
    const victim = chess.get(target);
    if (!victim || victim.color === piece.color) continue;
    if (chess.attackers(target, piece.color).includes(square)) hits.push({ square: target, piece: victim });
  }
  return hits;
}

function forkOf(chess, move) {
  const hits = piecesAttackedBy(chess, move.to);
  const mover = VALUE[move.piece];
  const juicy = hits.filter(hit => {
    if (hit.piece.type === 'k') return true;
    if (hit.piece.type === 'p') return false;
    if (VALUE[hit.piece.type] > mover) return true;
    return chess.attackers(hit.square, hit.piece.color).length === 0;
  });
  return juicy.length >= 2 ? juicy : null;
}

function numbering(fen) {
  const parts = fen.split(' ');
  return { color: parts[1], number: Number(parts[5] || 1) };
}

export function plyLabels(fen, steps) {
  let { color, number } = numbering(fen);
  return steps.map(step => {
    const label = `${number}${color === 'w' ? '.' : '…'}`;
    if (color === 'b') number++;
    color = other(color);
    return label + step.move.san;
  });
}

// "2…Nf6 3.Bc4 Nc6": move numbers only where a reader needs them.
export function lineText(fen, steps) {
  let { color, number } = numbering(fen);
  return steps.map((step, i) => {
    let text = step.move.san;
    if (color === 'w') text = `${number}.${text}`;
    else if (i === 0) text = `${number}…${text}`;
    if (color === 'b') number++;
    color = other(color);
    return text;
  }).join(' ');
}

// What a single move in a line does, in a few words.
export function stepNote(fenBefore, move) {
  const chess = load(fenBefore);
  const notes = [];
  const played = moveFromUci(chess, move.from + move.to + (move.promotion || ''));
  if (!played) return '';
  if (played.isKingsideCastle() || played.isQueensideCastle()) notes.push('castles');
  if (played.captured) notes.push(`takes the ${NAME[played.captured]}`);
  if (played.promotion) notes.push(`promotes to a ${NAME[played.promotion]}`);
  if (chess.isCheckmate()) notes.push('checkmate');
  else if (chess.inCheck()) notes.push('check');
  const fork = forkOf(chess, played);
  if (fork) notes.push(`forks ${fork.map(hit => NAME[hit.piece.type] === 'king' ? 'the king' : `the ${NAME[hit.piece.type]}`).join(' and ')}`);
  const loose = pressure(chess.fen(), played.color).filter(item => item.hanging && !item.check);
  if (loose.some(item => item.square === played.to) && !played.captured) notes.push('can be taken');
  return notes.join(', ');
}

// Material swing for `color` along a line, played until it goes quiet.
function swing(fen, uciMoves, color) {
  const chess = load(fen);
  const start = material(chess);
  const startDiff = start[color] - start[other(color)];
  let plies = 0;
  for (const uci of uciMoves) {
    const move = moveFromUci(chess, uci);
    if (!move) break;
    plies++;
    if (plies >= 6 && !move.captured && !chess.inCheck()) break;
    if (plies >= 12) break;
  }
  const end = material(chess);
  return end[color] - end[other(color)] - startDiff;
}

// The opponent capture that actually wins something (not one that is simply recaptured).
function winningCapture(fen, uciMoves, color) {
  const steps = playLine(fen, uciMoves, 10);
  for (let i = 1; i < steps.length; i++) {
    const move = steps[i].move;
    if (move.color === color || !move.captured) continue;
    const reply = steps[i + 1] && steps[i + 1].move;
    const back = reply && reply.captured ? VALUE[reply.captured] : 0;
    if (back < VALUE[move.captured]) return { move, index: i, steps };
  }
  return null;
}

function scoreFor(score, color) {
  return color === 'w' ? winPercent(score) : 100 - winPercent(score);
}

function forMover(score, color) {
  if (!score) return null;
  if (score.mate !== null && score.mate !== undefined) {
    return { ...score, mate: color === 'w' ? score.mate : -score.mate, mated: score.mated };
  }
  return { cp: color === 'w' ? score.cp : -score.cp, mate: null };
}

const CLASS_TEXT = {
  best: 'Best move',
  excellent: 'Excellent',
  good: 'Good',
  inaccuracy: 'Inaccuracy',
  mistake: 'Mistake',
  blunder: 'Blunder',
  book: 'Book move'
};
export function classLabel(kind) { return CLASS_TEXT[kind] || ''; }

export function classify({ inBook, playedUci, before, after, color }) {
  if (inBook) return { kind: 'book', drop: 0 };
  if (!before || !after) return null;
  const best = before.pv && before.pv[0];
  const drop = Math.max(0, scoreFor(before.score, color) - scoreFor(after.score, color));
  if (best === playedUci) return { kind: 'best', drop: 0 };
  if (before.alternatives && before.alternatives.includes(playedUci) && drop < 2) return { kind: 'excellent', drop };
  if (drop < 2) return { kind: 'excellent', drop };
  if (drop < 5) return { kind: 'good', drop };
  if (drop < 10) return { kind: 'inaccuracy', drop };
  if (drop < 20) return { kind: 'mistake', drop };
  return { kind: 'blunder', drop };
}

// The headline and reasons shown when the user asks "why?".
export function explainMove({ fenBefore, move, best, after, history }) {
  const color = move.color;
  const enemy = other(color);
  const points = [];
  const bestUci = best && best.pv && best.pv[0];
  const chessBefore = load(fenBefore);
  const bestMove = bestUci ? moveFromUci(load(fenBefore), bestUci) : null;
  const bestSan = bestMove ? bestMove.san : null;
  const replyPv = (after && after.pv) || [];
  const afterFen = (() => { const c = load(fenBefore); moveFromUci(c, move.from + move.to + (move.promotion || '')); return c.fen(); })();
  const replyLine = playLine(afterFen, replyPv, 10);
  const replyLabels = plyLabels(afterFen, replyLine.slice(1));
  const moverBefore = forMover(best && best.score, color);
  const moverAfter = forMover(after && after.score, color);
  let headline = '';

  // 1. Mates, either allowed or missed.
  if (moverAfter && moverAfter.mate !== null && moverAfter.mate < 0) {
    const n = Math.abs(moverAfter.mate);
    headline = n === 1 ? `${move.san} allows mate in one.` : `${move.san} allows a forced mate in ${n}.`;
    points.push(`The opponent mates with <b>${lineText(afterFen, replyLine.slice(1, n * 2))}</b>.`);
  } else if (moverBefore && moverBefore.mate !== null && moverBefore.mate > 0 && !(moverAfter && moverAfter.mate > 0)) {
    headline = `You had a forced mate in ${moverBefore.mate}.`;
    points.push(`It starts with <b>${bestSan}</b>${moverBefore.mate === 1 ? ' — mate right away' : ''}.`);
  }

  // 2. Material the reply wins.
  const lost = -swing(afterFen, replyPv, color);
  const capture = winningCapture(afterFen, replyPv, color);
  if (lost >= 2 && capture) {
    const victimSquare = capture.move.to;
    const victim = NAME[capture.move.captured];
    const label = replyLabels[capture.index - 1];
    if (!headline) headline = `${move.san} loses ${lost >= 3 && capture.move.captured !== 'p' ? `the ${victim}` : `${lost === 1 ? 'a pawn' : `material (${lost} pawns’ worth)`}`}.`;
    const afterPressure = pressure(afterFen, color);
    const exposed = afterPressure.find(item => item.square === victimSquare && item.hanging);
    if (capture.index === 1 && victimSquare === move.to) {
      points.push(`Your ${victim} on ${victimSquare} can simply be taken: <b>${label}</b>.`);
    } else if (capture.index === 1 && exposed) {
      const before = pressure(fenBefore, color).find(item => item.square === victimSquare);
      const guardLeft = before && before.defenders.includes(move.from) && !exposed.defenders.includes(move.to);
      if (guardLeft) points.push(`Your ${NAME[move.piece]} on ${move.from} was guarding the ${victim} on ${victimSquare}. Once it moved, <b>${label}</b> wins it.`);
      else points.push(`It leaves the ${victim} on ${victimSquare} ${exposed.defenders.length ? 'under-protected' : 'undefended'}, and <b>${label}</b> wins it.`);
    } else {
      points.push(`After <b>${lineText(afterFen, replyLine.slice(1, capture.index + 1))}</b> the ${victim} on ${victimSquare} falls.`);
    }
  }

  // 3. Forks in the reply.
  if (replyLine[1]) {
    const replyMove = replyLine[1].move;
    const chess = load(replyLine[1].fen);
    const fork = forkOf(chess, replyMove);
    if (fork) {
      const targets = fork.map(hit => hit.piece.type === 'k' ? 'king' : NAME[hit.piece.type]);
      points.push(`<b>${replyLabels[0]}</b> is a fork: it hits your ${targets.join(' and ')} at once.`);
      if (!headline) headline = `${move.san} walks into a fork.`;
    }
  }

  // 4. What the best move would have won.
  if (bestUci && best.pv) {
    const gain = swing(fenBefore, best.pv, color);
    if (gain >= 2 && !headline.startsWith('You had')) {
      const bestLine = playLine(fenBefore, best.pv, 10);
      const grab = bestLine.slice(1).find(step => step.move.color === color && step.move.captured);
      points.push(`You missed <b>${bestSan}</b>, which wins ${grab ? `the ${NAME[grab.move.captured]}` : 'material'}.`);
      if (!headline) headline = `${move.san} misses a chance to win material.`;
    }
  }

  // 5. Opening principles when no tactic explains it.
  const fullmove = Number(fenBefore.split(' ')[5] || 1);
  if (points.length === 0 && fullmove <= 14) {
    const mine = (history || []).filter(item => item.color === color);
    const minorsHome = countUndeveloped(chessBefore, color);
    if (move.piece === 'q' && fullmove <= 7 && !move.captured) {
      points.push('Bringing the queen out this early lets your opponent develop while attacking it, gaining time.');
    }
    if (move.piece !== 'p' && move.piece !== 'k' && mine.some(item => item.to === move.from) && minorsHome >= 2) {
      points.push(`That piece already moved, and ${minorsHome} of your minor pieces are still at home. In the opening, every move should bring a new piece into play.`);
    }
    if (move.piece === 'p' && 'ah'.includes(move.from[0]) && fullmove <= 10) {
      points.push('Edge-pawn moves do little for the center or your development.');
    }
    if (move.piece === 'p' && move.from[0] === 'f' && fullmove <= 10) {
      points.push(`Moving the f-pawn early opens the diagonal toward your king.`);
    }
    if (move.piece === 'k' && !move.isKingsideCastle() && !move.isQueensideCastle()) {
      points.push('The king gives up its right to castle and stays in the middle.');
    }
    if (bestMove) {
      if (bestMove.isKingsideCastle() || bestMove.isQueensideCastle()) {
        points.push(`<b>${bestSan}</b> was the moment to castle: king safe, rooks connected.`);
      } else if ((bestMove.piece === 'n' || bestMove.piece === 'b') && (bestMove.from[1] === '1' || bestMove.from[1] === '8')) {
        points.push(`<b>${bestSan}</b> develops another piece toward the center.`);
      } else if (bestMove.piece === 'p' && 'cde'.includes(bestMove.from[0])) {
        points.push(`<b>${bestSan}</b> claims the center with a pawn.`);
      }
    }
  }

  // 6. Nothing concrete: say what the opponent gets.
  if (!headline) {
    const dropText = moverBefore && moverAfter ? `${formatScore(moverBefore)} → ${formatScore(moverAfter)}` : '';
    headline = bestSan ? `${bestSan} was stronger.` : `${move.san} is not the most accurate.`;
    if (replyLabels.length) {
      points.push(`The engine expects <b>${lineText(afterFen, replyLine.slice(1, 4))}</b>, and your position gets worse${dropText ? ` (${dropText} for you)` : ''}.`);
    }
  }
  if (bestSan && !points.some(point => point.includes(`<b>${bestSan}</b>`))) {
    points.push(`Better was <b>${bestSan}</b>${moverBefore ? ` (${formatScore(moverBefore)} for you)` : ''}.`);
  }
  return { headline, points };
}

function countUndeveloped(chess, color) {
  const home = color === 'w' ? ['b1', 'g1', 'c1', 'f1'] : ['b8', 'g8', 'c8', 'f8'];
  const kinds = ['n', 'n', 'b', 'b'];
  return home.filter((square, i) => {
    const piece = chess.get(square);
    return piece && piece.color === color && piece.type === kinds[i];
  }).length;
}
