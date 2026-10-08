// Opening book: which moves are theory from a position, and what that theory is called.
import { LINES, POSITIONS } from './openings.js?v=2026-10-08-a';

// Same key the build script used: a hash of the first four FEN fields.
export function bookKey(fen) {
  const epd = fen.split(' ').slice(0, 4).join(' ');
  let h1 = 0x811c9dc5;
  let h2 = 0x01000193 ^ 0x5bd1e995;
  for (let i = 0; i < epd.length; i++) {
    const c = epd.charCodeAt(i);
    h1 = Math.imul(h1 ^ c, 16777619) >>> 0;
    h2 = Math.imul(h2 ^ c, 2246822519) >>> 0;
  }
  return h1.toString(36) + (h2 & 0xfffff).toString(36);
}

function line(index) {
  if (index < 0) return null;
  const [eco, name, moves] = LINES[index];
  const split = name.indexOf(':');
  return {
    index,
    eco,
    name,
    family: split === -1 ? name : name.slice(0, split),
    variation: split === -1 ? '' : name.slice(split + 1).trim(),
    moves: moves.split(' ')
  };
}

export function lookup(fen) {
  const record = POSITIONS[bookKey(fen)];
  if (!record) return null;
  const [exact, representative, list] = record;
  const moves = list
    ? list.split(',').map(item => {
      const [uci, count] = item.split(':');
      return { uci, count: Number(count || 1) };
    })
    : [];
  return { exact: line(exact), representative: line(representative), moves };
}

export function bookMoves(fen) {
  const record = lookup(fen);
  return record ? record.moves : [];
}

export function family(name) {
  const base = name.replace(/ (Accepted|Declined|Refused)$/, '');
  return FAMILIES[name] || FAMILIES[base] || null;
}

export function ecoGroup(eco) {
  return ECO_GROUPS[eco ? eco[0] : ''] || null;
}

const ECO_GROUPS = {
  A: 'Flank openings and unusual first moves: fights for the center from a distance.',
  B: 'Semi-open games: Black answers 1.e4 with something other than 1…e5.',
  C: 'Open games after 1.e4 e5, plus the French Defense.',
  D: 'Closed games: 1.d4 d5 and the Queen’s Gambit complex.',
  E: 'Indian defenses: Black meets 1.d4 with …Nf6 and controls the center with pieces.'
};

// Short, practical notes for the families a learner meets most.
const FAMILIES = {
  "Sicilian Defense": {
    tags: ['Sharp', 'Asymmetric'],
    summary: 'The most popular reply to 1.e4. Black trades a wing pawn for a center pawn and plays for an unbalanced, double-edged game.',
    white: 'Open the d-file with d4, develop fast and attack on the kingside.',
    black: 'Use the half-open c-file, counterattack on the queenside and keep the extra central pawn.'
  },
  "Ruy Lopez": {
    tags: ['Classical', 'Strategic'],
    summary: 'One of the oldest openings. The bishop on b5 pressures the knight that defends e5, building slow, lasting pressure.',
    white: 'Prepare d4 with c3, keep the tension and reroute the queen’s knight to g3 or e3.',
    black: 'Chase the bishop with …a6 and …b5, hold e5 firmly and look for …d5 at the right moment.'
  },
  "Italian Game": {
    tags: ['Classical', 'Open'],
    summary: 'The bishop goes to c4 and eyes f7, the weakest point in Black’s camp. Ranges from quiet maneuvering to wild gambits.',
    white: 'Build the center with c3 and d4, or play slowly with d3 and a later d4.',
    black: 'Develop symmetrically, castle quickly and strike back with …d5 when it is safe.'
  },
  "Queen's Gambit Declined": {
    tags: ['Solid', 'Classical'],
    summary: 'Black keeps a firm pawn on d5 instead of taking on c4. A rock-solid choice trusted at every level.',
    white: 'Press the d5-pawn, use the c-file and aim for a minority attack with b4–b5.',
    black: 'Finish development, then free the position with …c5 or …e5 and solve the c8-bishop.'
  },
  "Queen's Gambit Accepted": {
    tags: ['Classical'],
    summary: 'Black takes the c4-pawn, not to keep it but to gain time and free the game.',
    white: 'Win the pawn back with e3 and Bxc4, then use the strong center.',
    black: 'Hit the center with …c5 and …e6, developing smoothly.'
  },
  "Queen's Gambit": {
    tags: ['Classical'],
    summary: 'White offers the c-pawn to pull Black’s d-pawn away from the center.',
    white: 'Gain central control with e4 if Black ever lets go of d5.',
    black: 'Decide: hold d5 (Declined, Slav) or take and return the pawn (Accepted).'
  },
  "Slav Defense": {
    tags: ['Solid'],
    summary: 'Black supports d5 with the c-pawn, keeping the light-squared bishop free to develop.',
    white: 'Play for e4 and a space advantage.',
    black: 'Develop the bishop to f5 or g4 before playing …e6.'
  },
  "Semi-Slav Defense": {
    tags: ['Rich', 'Sharp'],
    summary: 'A combination of …c6 and …e6: very solid yet full of dynamic, deeply analysed lines.',
    white: 'Choose between the sharp Botvinnik/Moscow lines and quieter setups.',
    black: 'Prepare …dxc4 followed by …b5, or the freeing …e5.'
  },
  "French Defense": {
    tags: ['Solid', 'Counterattacking'],
    summary: 'After …e6 and …d5 Black challenges e4 at once. The price is a sometimes passive light-squared bishop.',
    white: 'Gain space with e5 and attack on the kingside.',
    black: 'Undermine White’s pawn chain with …c5 and …f6.'
  },
  "Caro-Kann Defense": {
    tags: ['Solid'],
    summary: '…c6 and …d5: like the French, but the light-squared bishop gets out before …e6.',
    white: 'Use the extra space; the Advance and Panov lines are the sharpest tries.',
    black: 'Develop the bishop to f5 or g4, then build a healthy, compact position.'
  },
  "English Opening": {
    tags: ['Flexible', 'Positional'],
    summary: '1.c4 controls d5 from the flank. Often a reversed Sicilian with an extra tempo.',
    white: 'Fianchetto the king’s bishop and play on the long diagonal and queenside.',
    black: 'Take central space with …e5 or …d5, or transpose into Indian setups.'
  },
  "King's Gambit Accepted": {
    tags: ['Romantic', 'Gambit'],
    summary: 'White gives the f-pawn to rip open the f-file and seize the center. Attack-or-die chess from the 1800s.',
    white: 'Play d4, Bc4 and O-O; aim everything at f7.',
    black: 'Hold the extra pawn with …g5 or give it back with …d5 for easy development.'
  },
  "King's Gambit Declined": {
    tags: ['Gambit'],
    summary: 'Black refuses the f-pawn, usually with …Bc5 to keep White’s king from castling short.',
    white: 'Develop and castle, keeping the f-file in mind.',
    black: 'Use the a7–g1 diagonal and develop normally.'
  },
  "King's Gambit": {
    tags: ['Romantic', 'Gambit'],
    summary: 'White offers the f-pawn on move two for a fast attack.',
    white: 'Open the f-file and build a big center.',
    black: 'Take and defend, or decline with a solid setup.'
  },
  "King's Indian Defense": {
    tags: ['Dynamic', 'Hypermodern'],
    summary: 'Black lets White build a big center, then attacks it with …e5 or …c5 and often storms the kingside.',
    white: 'Use the space advantage and open the queenside with c5.',
    black: 'Lock the center and launch …f5, …f4 and …g5 at White’s king.'
  },
  "Nimzo-Indian Defense": {
    tags: ['Positional', 'Hypermodern'],
    summary: '…Bb4 pins the knight to stop e4. Black is often willing to trade the bishop to double White’s pawns.',
    white: 'Keep the bishop pair and play for e4 or a big center.',
    black: 'Control e4, provoke weaknesses and blockade on the light squares.'
  },
  "Queen's Indian Defense": {
    tags: ['Solid', 'Hypermodern'],
    summary: 'Black fianchettoes the queen’s bishop to watch e4 from afar.',
    white: 'Fight for e4 with Nc3 and Qc2, or fianchetto too.',
    black: 'Keep a grip on e4 and the long diagonal.'
  },
  "Grünfeld Defense": {
    tags: ['Dynamic', 'Hypermodern'],
    summary: 'Black invites a big pawn center, then attacks it with pieces and the c-pawn.',
    white: 'Hold the d4–e4 center and use it to attack.',
    black: 'Pressure d4 with the g7-bishop, …c5 and …Nc6.'
  },
  "Benoni Defense": {
    tags: ['Sharp', 'Asymmetric'],
    summary: 'Black accepts less space for active play: a queenside majority and a powerful g7-bishop.',
    white: 'Use the space and push e5 when possible.',
    black: 'Push the queenside pawns with …b5 and use the e-file.'
  },
  "Dutch Defense": {
    tags: ['Aggressive', 'Asymmetric'],
    summary: '1…f5 grabs e4 and promises kingside play, at some cost to king safety.',
    white: 'Fianchetto and aim for e4; the king’s diagonal can become weak.',
    black: 'Set up a Stonewall or Leningrad and attack on the kingside.'
  },
  "Queen's Pawn Game": {
    tags: ['Flexible'],
    summary: 'White opens 1.d4 without an early c4. Includes system openings like the London and Colle.',
    white: 'Build a reliable setup and play for a slow kingside initiative.',
    black: 'Take the center with …d5 and …c5, and develop freely.'
  },
  "Indian Defense": {
    tags: ['Flexible', 'Hypermodern'],
    summary: '1…Nf6 keeps options open: King’s Indian, Nimzo, Grünfeld and more can still arise.',
    white: 'Choose a setup: c4 and Nc3 for the main lines, or a system.',
    black: 'Control e4 with pieces before committing the pawns.'
  },
  "Alekhine Defense": {
    tags: ['Provocative', 'Hypermodern'],
    summary: 'The knight invites White’s pawns forward, hoping they become targets.',
    white: 'Gain big space but avoid overextending.',
    black: 'Undermine the pawn center with …d6, …c5 and pieces.'
  },
  "Scotch Game": {
    tags: ['Open'],
    summary: 'White opens the center at once with d4, trading off the e5-pawn.',
    white: 'Use the open position and quick development.',
    black: 'Counter with …Bc5 or …Nf6 and hit the e4-pawn.'
  },
  "Petrov's Defense": {
    tags: ['Solid', 'Symmetrical'],
    summary: 'Black counterattacks e4 instead of defending e5: famously drawish and reliable.',
    white: 'Look for a small but long-lasting initiative.',
    black: 'Equalize with symmetric, harmonious development; beware the Nxf7 tricks early on.'
  },
  "Four Knights Game": {
    tags: ['Classical', 'Symmetrical'],
    summary: 'Both sides develop the knights first — calm and principled.',
    white: 'Bb5 or d4 creates small, lasting pressure.',
    black: 'Mirror carefully and strike in the center.'
  },
  "Scandinavian Defense": {
    tags: ['Direct'],
    summary: '1…d5 hits e4 immediately. The queen often recaptures early and has to find a safe square.',
    white: 'Gain time chasing the queen with Nc3.',
    black: 'Place the queen on a5 or d6 and develop quickly.'
  },
  "Philidor Defense": {
    tags: ['Solid', 'Passive'],
    summary: '…d6 supports e5 with a pawn: solid, but a little cramped.',
    white: 'Use the extra space and develop actively.',
    black: 'Keep e5 and aim for a compact Hanham setup.'
  },
  "Nimzowitsch Defense": {
    tags: ['Offbeat'],
    summary: '1…Nc6 against 1.e4: a rare, flexible try.',
    white: 'Take the center with d4.',
    black: 'Hit d4 with pieces and …e5 or …d5.'
  },
  "Réti Opening": {
    tags: ['Hypermodern', 'Flexible'],
    summary: '1.Nf3 with c4 and a fianchetto: control the center from the flanks.',
    white: 'Pressure d5 and the long diagonal.',
    black: 'Occupy the center and see if White over-finesses.'
  },
  "Zukertort Opening": {
    tags: ['Flexible'],
    summary: '1.Nf3: develop first, decide on the pawn structure later.',
    white: 'Transpose into a favorite setup.',
    black: 'Take the center with …d5 or mirror with …Nf6.'
  },
  "Vienna Game": {
    tags: ['Open', 'Aggressive'],
    summary: 'Nc3 before Nf3 keeps the f-pawn free for a later f4.',
    white: 'Play f4 and attack on the kingside.',
    black: 'Strike in the center with …Nf6 and …d5.'
  },
  "Bishop's Opening": {
    tags: ['Open'],
    summary: 'Bc4 straight away aims at f7 while keeping the f-pawn flexible.',
    white: 'Can transpose into the Italian or play a Vienna setup.',
    black: 'Develop the knight to f6 and challenge with …d5.'
  },
  "Modern Defense": {
    tags: ['Hypermodern'],
    summary: '…g6 and …Bg7 without an early …Nf6: flexible but space-conceding.',
    white: 'Take the full center and develop.',
    black: 'Undermine with …c5, …e5 or …d6 setups.'
  },
  "Pirc Defense": {
    tags: ['Hypermodern'],
    summary: '…d6, …Nf6 and …g6: Black lets White build a center and attacks it later.',
    white: 'Sharp tries like the Austrian Attack with f4.',
    black: 'Castle and counter with …c5 or …e5.'
  },
  "Catalan Opening": {
    tags: ['Positional'],
    summary: 'Queen’s Gambit plus a g2-fianchetto: long-term pressure on the light squares.',
    white: 'Pressure the queenside with the g2-bishop.',
    black: 'Take on c4 and hold the pawn or free the game with …c5.'
  },
  "Bird Opening": {
    tags: ['Offbeat', 'Aggressive'],
    summary: '1.f4 controls e5 and starts a reversed Dutch.',
    white: 'Kingside attack after Nf3, e3 and Be2/d3.',
    black: 'Strike back with …e5 (From Gambit) or develop solidly.'
  },
  "Tarrasch Defense": {
    tags: ['Active'],
    summary: 'Black accepts an isolated d-pawn in return for free, active pieces.',
    white: 'Blockade and attack the isolated pawn.',
    black: 'Use the activity before the endgame arrives.'
  },
  "London System": {
    tags: ['System', 'Solid'],
    summary: 'd4, Bf4, e3, c3: a reliable setup that can be played against almost anything.',
    white: 'Build the pyramid and play on the kingside.',
    black: 'Challenge with …c5 and …Qb6, hitting b2.'
  },
  "Trompowsky Attack": {
    tags: ['Offbeat'],
    summary: '2.Bg5 attacks the knight at once and avoids a lot of theory.',
    white: 'Trade on f6 to give Black doubled pawns, or keep it tense.',
    black: '…Ne4 and …c5 are the most direct replies.'
  },
  "Benko Gambit": {
    tags: ['Gambit', 'Positional'],
    summary: 'Black gives a queenside pawn for lasting pressure on the a- and b-files.',
    white: 'Keep the extra pawn and consolidate.',
    black: 'Rooks to a8 and b8, bishop on g7: long-term queenside pressure.'
  },
  "Owen Defense": {
    tags: ['Offbeat'],
    summary: '1…b6: a queenside fianchetto against 1.e4.',
    white: 'Take the full center.',
    black: 'Pressure e4 from b7.'
  },
  "Center Game": {
    tags: ['Open'],
    summary: 'White takes back on d4 with the queen — direct but it loses time.',
    white: 'Castle long and attack.',
    black: 'Develop with tempo against the queen.'
  },
  "Danish Gambit": {
    tags: ['Gambit', 'Romantic'],
    summary: 'White sacrifices one or two pawns for two raking bishops.',
    white: 'Aim the bishops at f7 and g7.',
    black: 'Return material with …d5 to finish development.'
  },
  "Englund Gambit": {
    tags: ['Gambit', 'Trappy'],
    summary: '1.d4 e5?!: a dubious but trap-filled gambit.',
    white: 'Hold the pawn, watch out for …Qe7 tricks.',
    black: 'Set traps with …Nc6 and …Qe7.'
  },
  "Latvian Gambit": {
    tags: ['Gambit', 'Risky'],
    summary: '2…f5!?: an ultra-sharp, objectively shaky counter-gambit.',
    white: 'Nxe5 and Qh5+ ideas punish the weakened king.',
    black: 'Play fast and hope for complications.'
  },
  "Elephant Gambit": {
    tags: ['Gambit', 'Risky'],
    summary: '2…d5!? hits back in the center at the cost of a pawn.',
    white: 'Take and develop calmly.',
    black: 'Rapid development for the pawn.'
  },
  "Ponziani Opening": {
    tags: ['Classical'],
    summary: '3.c3 prepares d4: an old, slightly forgotten line.',
    white: 'Build a big center with d4.',
    black: '…d5 or …Nf6 hit back in the center.'
  },
  "Blackmar-Diemer Gambit": {
    tags: ['Gambit', 'Aggressive'],
    summary: 'White gambits the e-pawn for the open f-file and fast development.',
    white: 'Open lines toward Black’s king.',
    black: 'Keep the pawn and stay calm.'
  },
  "Bogo-Indian Defense": {
    tags: ['Solid'],
    summary: '…Bb4+ against Nf3 setups: an easy, solid way to develop.',
    white: 'Block with Bd2 or Nbd2 and keep the bishop pair.',
    black: 'Trade a bishop and play with simple development.'
  },
  "King's Indian Attack": {
    tags: ['System'],
    summary: 'A reversed King’s Indian for White: Nf3, g3, Bg2, O-O, d3 and e4.',
    white: 'Push e5 and attack on the kingside.',
    black: 'Expand on the queenside.'
  },
  "Three Knights Opening": {
    tags: ['Classical'],
    summary: 'A close cousin of the Four Knights.',
    white: 'Develop and decide on d4 or Bb5.',
    black: 'Develop the last knight or try …g6.'
  },
  "Old Indian Defense": {
    tags: ['Solid'],
    summary: '…d6 and …e5 without a fianchetto: compact and solid.',
    white: 'Take space.',
    black: 'Hold e5 and develop calmly.'
  },
  "Neo-Grünfeld Defense": {
    tags: ['Hypermodern'],
    summary: 'A Grünfeld against a fianchetto setup.',
    white: 'Use the g2-bishop on the long diagonal.',
    black: 'Pressure d4 just like in the Grünfeld.'
  },
  "Polish Opening": {
    tags: ['Offbeat'],
    summary: '1.b4: a queenside expansion aimed at the long diagonal.',
    white: 'Bb2 pressures e5.',
    black: 'Attack the b4-pawn and build a center.'
  },
  "Van Geet Opening": {
    tags: ['Offbeat'],
    summary: '1.Nc3: flexible, rare and often transposing.',
    white: 'Aim for e4 and a quick initiative.',
    black: 'Take the center with …d5 or …e5.'
  },
  "Hungarian Opening": {
    tags: ['Offbeat'],
    summary: '1.g3: a quiet fianchetto start.',
    white: 'Fianchetto and stay flexible.',
    black: 'Occupy the center.'
  },
  "Grob Opening": {
    tags: ['Offbeat', 'Risky'],
    summary: '1.g4: an eccentric, loosening first move.',
    white: 'Bg2 pressure on b7 and d5.',
    black: 'Take the center and punish the weakened kingside.'
  },
  "Pterodactyl Defense": {
    tags: ['Offbeat'],
    summary: '…g6, …Bg7 and an early …c5 and …Qa5: a hypermodern oddity.',
    white: 'Develop and take space.',
    black: 'Pressure along the long diagonal and the a5–e1 diagonal.'
  },
  "Nimzo-Larsen Attack": {
    tags: ['Flexible'],
    summary: '1.b3: the queen’s bishop goes to b2 and watches e5.',
    white: 'Control the long diagonal.',
    black: 'Build a broad center.'
  },
  "Richter-Veresov Attack": {
    tags: ['Offbeat'],
    summary: 'd4, Nc3 and Bg5 — an aggressive setup aiming for f3 and e4.',
    white: 'Play f3 and e4.',
    black: 'Challenge with …c5 or …Bf5.'
  },
  "Torre Attack": {
    tags: ['System'],
    summary: 'd4, Nf3 and Bg5: a sound, low-theory system.',
    white: 'e3, Bd3 and a later e4.',
    black: '…c5 and …Qb6 create early pressure.'
  },
  "Colle System": {
    tags: ['System'],
    summary: 'd4, Nf3, e3, Bd3: White prepares a central e4 break.',
    white: 'Play e4 and attack the king.',
    black: 'Develop the light-squared bishop actively before …e6.'
  },
  "Rat Defense": {
    tags: ['Offbeat'],
    summary: 'An early …d6 for a flexible but passive setup.',
    white: 'Take the center.',
    black: 'Transpose to a Pirc or Philidor setup.'
  },
  "King's Pawn Game": {
    tags: ['Open'],
    summary: '1.e4 grabs the center and opens lines for the queen and the king’s bishop. Black now chooses the kind of game: symmetrical with …e5, or a fighting defense like the Sicilian, French or Caro-Kann.',
    white: 'Follow up with d4 or Nf3 and develop quickly.',
    black: 'Contest the center at once — with a pawn (…e5, …c5, …d5) or by preparing it (…e6, …c6).'
  },
  "King's Knight Opening": {
    tags: ['Open'],
    summary: '1.e4 e5 2.Nf3: attack e5 and develop at the same time.',
    white: 'Develop and castle.',
    black: 'Defend e5 with …Nc6 or counter with …Nf6.'
  },
  "King's Pawn Opening": {
    tags: ['Open'],
    summary: '1.e4: the most direct first move, opening lines for the queen and bishop.',
    white: 'Take the center and develop.',
    black: 'Choose: symmetry with …e5 or a semi-open defense.'
  },
  "Mieses Opening": {
    tags: ['Offbeat'],
    summary: '1.d3: a modest, reversed-Pirc start.',
    white: 'Develop flexibly.',
    black: 'Take the center.'
  }
};
