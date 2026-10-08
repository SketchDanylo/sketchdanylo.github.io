// Stockfish in a worker, with a small scheduler.
// Every result is cached per position (scores stored from White's point of view),
// so the move classifier and the explainer can reuse whatever depth was reached.

const ENGINE_URL = new URL('./engine/stockfish-19-lite-single.js', import.meta.url);

export class Engine {
  constructor() {
    this.cache = new Map();      // fen -> { depth, lines: [{ cp, mate, pv, depth }] }
    this.listeners = new Set();
    this.ready = false;
    this.current = null;         // job being searched
    this.pending = null;         // job waiting for the current search to stop
    this.urgent = [];            // quick searches that jump the queue (feedback needs them)
    this.background = [];        // queue of { fen, depth, multipv }
    this.failed = false;
    try {
      this.worker = new Worker(ENGINE_URL);
    } catch (error) {
      this.failed = true;
      return;
    }
    this.worker.onmessage = event => this.receive(String(event.data));
    this.worker.onerror = () => { this.failed = true; this.emit({ type: 'error' }); };
    this.send('uci');
  }

  send(command) { this.worker && this.worker.postMessage(command); }
  on(listener) { this.listeners.add(listener); return () => this.listeners.delete(listener); }
  emit(message) { this.listeners.forEach(listener => listener(message)); }

  entry(fen) { return this.cache.get(fen) || null; }

  receive(line) {
    if (line === 'uciok') {
      this.send('setoption name Hash value 64');
      this.send('setoption name UCI_ShowWDL value false');
      this.send('isready');
      return;
    }
    if (line === 'readyok') {
      if (!this.ready) {
        this.ready = true;
        this.emit({ type: 'ready' });
        this.pump();
      }
      return;
    }
    if (line.startsWith('info ') && this.current) {
      this.readInfo(line);
      return;
    }
    if (line.startsWith('bestmove')) {
      const finished = this.current;
      this.current = null;
      if (finished && !finished.stopped) {
        const entry = this.cache.get(finished.fen);
        if (entry) entry.complete = Math.max(entry.complete || 0, finished.depth);
        this.emit({ type: 'done', fen: finished.fen, job: finished });
      }
      this.pump();
    }
  }

  readInfo(line) {
    const job = this.current;
    if (line.includes(' lowerbound') || line.includes(' upperbound')) return;
    const depthMatch = line.match(/ depth (\d+)/);
    const pvMatch = line.match(/ pv (.+)$/);
    const scoreMatch = line.match(/ score (cp|mate) (-?\d+)/);
    if (!depthMatch || !pvMatch || !scoreMatch) return;
    const depth = Number(depthMatch[1]);
    const multipv = Number((line.match(/ multipv (\d+)/) || [0, 1])[1]);
    const whiteToMove = job.fen.split(' ')[1] === 'w';
    const sign = whiteToMove ? 1 : -1;
    const value = Number(scoreMatch[2]) * sign;
    // "mate 0" means the side to move is already mated; the page scores finished games itself.
    if (scoreMatch[1] === 'mate' && value === 0) return;
    const scored = scoreMatch[1] === 'cp' ? { cp: value, mate: null } : { cp: null, mate: value };

    let entry = this.cache.get(job.fen);
    if (!entry) {
      entry = { depth: 0, lines: [], complete: 0 };
      this.cache.set(job.fen, entry);
    }
    // Keep a deeper result from an earlier search rather than overwrite it with a shallow one.
    if (depth < entry.depth && multipv === 1 && !job.fresh) return;
    if (multipv === 1) {
      job.fresh = true;
      entry.depth = depth;
    }
    if (multipv > job.multipv) return;
    entry.lines[multipv - 1] = { ...scored, pv: pvMatch[1].trim().split(' '), depth };
    this.emit({ type: 'info', fen: job.fen, entry, job });
  }

  // The foreground search: whatever position is on the board right now.
  analyze(fen, { depth = 22, multipv = 3 } = {}) {
    const job = { fen, depth, multipv, kind: 'view' };
    const entry = this.cache.get(fen);
    if (entry && (entry.complete || 0) >= depth && entry.lines.length >= multipv) {
      this.pending = null;
      if (this.current && this.current.kind === 'view') this.stopCurrent();
      this.emit({ type: 'info', fen, entry, job });
      this.emit({ type: 'done', fen, job });
      return;
    }
    this.pending = job;
    if (this.current) this.stopCurrent();
    else this.pump();
  }

  // Low-priority depth for positions the user is not looking at.
  queue(fen, depth = 14) {
    const entry = this.cache.get(fen);
    if (entry && (entry.complete || 0) >= depth) return;
    if (this.background.some(job => job.fen === fen && job.depth >= depth)) return;
    this.background.push({ fen, depth, multipv: 1, kind: 'background' });
    if (!this.current && !this.pending) this.pump();
  }

  // A quick search that runs before anything else (but never interrupts the view search mid-way
  // unless the view is about to change anyway).
  prioritize(fen, depth = 12) {
    const entry = this.cache.get(fen);
    if (entry && Math.max(entry.complete || 0, entry.depth) >= depth) return;
    if (this.urgent.some(job => job.fen === fen)) return;
    this.urgent.push({ fen, depth, multipv: 1, kind: 'urgent' });
    if (!this.current) this.pump();
  }

  stopCurrent() {
    if (!this.current || this.current.stopped) return;
    this.current.stopped = true;
    this.send('stop');
  }

  pump() {
    if (!this.ready || this.current) return;
    let job = null;
    while (!job && this.urgent.length) {
      const candidate = this.urgent.shift();
      const entry = this.cache.get(candidate.fen);
      if (!entry || Math.max(entry.complete || 0, entry.depth) < candidate.depth) job = candidate;
    }
    if (!job) {
      job = this.pending;
      this.pending = null;
    }
    if (!job) {
      while (this.background.length) {
        const candidate = this.background.shift();
        const entry = this.cache.get(candidate.fen);
        if (!entry || (entry.complete || 0) < candidate.depth) { job = candidate; break; }
      }
    }
    if (!job) { this.emit({ type: 'idle' }); return; }
    this.current = job;
    this.send(`setoption name MultiPV value ${job.multipv}`);
    this.send(`position fen ${job.fen}`);
    this.send(`go depth ${job.depth}`);
  }
}

// ---- score helpers (all from White's point of view) ----

export function scoreOf(line) {
  if (!line) return null;
  return { cp: line.cp, mate: line.mate };
}

export function scoreToCp(score) {
  if (!score) return 0;
  if (score.mate !== null && score.mate !== undefined) {
    if (score.mate === 0) return score.mated === 'w' ? -10000 : 10000;
    return Math.sign(score.mate) * (10000 - Math.abs(score.mate) * 10);
  }
  return score.cp;
}

// Lichess' win-probability curve.
export function winPercent(score) {
  const cp = Math.max(-1000, Math.min(1000, scoreToCp(score)));
  return 50 + 50 * (2 / (1 + Math.exp(-0.00368208 * cp)) - 1);
}

export function formatScore(score, { signed = true } = {}) {
  if (!score) return '…';
  if (score.mate !== null && score.mate !== undefined) {
    if (score.mate === 0) return '#';
    return `${score.mate > 0 ? '' : '−'}M${Math.abs(score.mate)}`;
  }
  const pawns = score.cp / 100;
  const text = Math.abs(pawns) >= 10 ? Math.abs(pawns).toFixed(0) : Math.abs(pawns).toFixed(1);
  if (Math.abs(pawns) < 0.05) return '0.0';
  return `${pawns > 0 ? (signed ? '+' : '') : '−'}${text}`;
}
