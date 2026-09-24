/* Chem Playground — reactive molecular-dynamics engine.
 *
 * Units: length Å, time fs, mass amu (g/mol), energy kJ/mol, charge e, temperature K.
 *
 * Potential (a compact bond-order reactive force field in the spirit of Tersoff/Brenner/ReaxFF):
 *   - Covalent pairs: Morse  E = T(r)·De·[e^(-2y) − 2·b·e^(-y)],  y = a(r − re)
 *       re, De, a come from measured bond lengths / dissociation energies / stretch force constants
 *       and are interpolated in the continuous bond order n (1…3) of the pair.
 *   - Saturation b = p_i·p_j,  p_i = 1 / (1 + c·ramp(Z_i − f_ij + 1 − V_i)²)
 *       Z_i is the smooth coordination number of atom i, V_i its valence. An atom that already
 *       uses all its valence cannot attract a new partner, and a partner can only bind by weakening
 *       an existing bond. This produces reaction barriers, radical reactivity and bond exchange.
 *   - Bond order n: spare valence (V − Z) is shared between neighbours that also have spare
 *       valence → C=C, C≡C, O=O, N≡N, CO₂, aromatic 1.5 bonds emerge from connectivity.
 *   - VSEPR angles: harmonic in cos θ, θ0 from steric number (coordination + lone pairs).
 *   - Non-bonded: shielded Lennard-Jones (UFF) — the Pauli part is switched off (1 − b) between
 *       atoms that can bond — plus shielded, shifted-force Coulomb between bond-polarisation
 *       charges (electronegativity) and formal charges. 1-2 and 1-3 pairs are excluded smoothly.
 *   - Velocity-Verlet, fixed dt = 1 fs, local wall heat exchange (CSVR for conditioning), soft walls whose
 *       normal force gives the measured pressure.
 */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.ChemEngine = factory();
})(typeof self !== 'undefined' ? self : this, function () {
'use strict';

const KB = 0.0083144626;            // kJ/(mol·K)
const ACC = 1e-4;                   // (kJ/mol/Å)/amu → Å/fs²
const KEU = 1e4;                    // amu·(Å/fs)² → kJ/mol
const clampNum = (x, lo, hi) => Math.max(lo, Math.min(hi, x));
const COUL = 1389.35458;            // kJ·Å/(mol·e²)
const BAR = 16605.39;               // (kJ/mol)/Å³ → bar
const Q_MAX = 1.1;                  // saturation of bond-polarisation charge (e)
// Reactive-model constants, fitted to barriers of H + H₂, H + CH₄, H + Cl₂, F + H₂, O + H₂, OH + H₂ and
// H + O₂ → HO₂ (playground/tests/reactions.cjs) while keeping closed-shell dimers non-sticky:
//   saturation p(x) = 1/(1 + c1·x + c4·x^k);  valence-use decay β = kb·a (a = Morse exponent);
//   wpi  = how much a π bond blocks a newcomer (weaker and more polarizable than a σ bond);
//   wpiOO = the same for O=O, lower still (triplet O₂ is a diradical);
//   oo3e = extra O–O bond order when only one oxygen is unpaired (three-electron bond, HO₂·);
//   mu   = strength weighting of excess valence (Evans–Polanyi: exothermic transfers get lower barriers).
const TUNE = { kappa: 0.38, c1: 0.6, c4: 10.82, k: 6, kb: 1.3, yoff: 3.0, wpi: 0.85, wpiOO: 0.78, oo3e: 0.45, mu: 0.5, kbo: 0.4, tsStab: 0, pauliOpen: 0, share: 0.5, insert: 1, vstate: 1, cap: 0.10, gel: 1 };
const RAMP_W = 0.15;
const Y_ON = 3.6, Y_OFF = 5.6;      // Morse taper window (in units of a·(r − re))
const COORD_R1 = 1.22, COORD_R2 = 1.50; // structural coordination switch (× single-bond length)
const BO_CAP = 0.10, BO_W = 0.08, BO_BETA = 3.6; // saturation bond order: 1 up to r₁+cap, then exp(−β·Δr)
const BO_ON = 0.50, BO_OFF = 0.85;  // …tapered to zero between r₁+on and r₁+off (Å)
const COUL_D2 = 0.3;                // Coulomb short-range shielding (Å²)
const LJ_SHIELD = 0.55;             // LJ shielding radius as a fraction of r_min
const OPEN_W = 0.25;                // smoothing width of the free-valence ramp
const SUB_DX = 0.1;                 // Å an atom may travel in one sub-step
const SUB_HOLD = 200;               // fs a finer sub-step count is kept before it may relax
/* At most 64 sub-steps; redo a step whose energy error exceeds 1.5 kJ/mol plus 0.2% of the kinetic
   energy. It was 6, which let two hydrogens plunging into their 436 kJ/mol well lose 15 to 36 kJ/mol
   over a few steps each just under the line — enough to stay bound with nothing to take the energy,
   which an isolated pair cannot do. A strict guard only became affordable once a redone step stayed
   finer instead of dropping straight back, because the dropping back was itself a source of drift. */
const SUB_MAX = 64, SUB_ETOL = 1.5;
const SOLID_W = 0.05;               // Å over which a solid face's force eases in
const WALL_FOLD = 1.0;              // Å past a solid face before an atom is carried back rather than pushed
const CORE_A = 800;                 // kJ/mol, strength of the nuclear hard core
const VMAX = 0.6;                   // Å/fs safety speed limit (60 km/s)

/* ---------- element data ----------
 sym, Z, name, mass, χ(Pauling), rcov[single,double,triple] (Pyykkö), r_vdw (Bondi/Alvarez),
 UFF x (Å), UFF D (kcal/mol), valences, valence electrons, max bond order, angle K (kJ/mol/rad²), color */
const ELEMENT_ROWS = [
  ['H', 1, 'Hydrogen', 1.008, 2.20, [0.32], 1.20, 2.886, 0.044, [1], 1, 1, 0, '#f3f3ef'],
  ['He', 2, 'Helium', 4.0026, 0, [0.46], 1.40, 2.362, 0.056, [0], 2, 0, 0, '#d9ffff'],
  ['Li', 3, 'Lithium', 6.94, 0.98, [1.33, 1.24], 1.82, 2.451, 0.025, [1], 1, 1, 60, '#cc80ff'],
  ['Be', 4, 'Beryllium', 9.0122, 1.57, [1.02, 0.90, 0.85], 1.53, 2.745, 0.085, [2], 2, 2, 120, '#c2ff00'],
  ['B', 5, 'Boron', 10.81, 2.04, [0.85, 0.78, 0.73], 1.92, 4.083, 0.180, [3], 3, 2, 200, '#ffb5b5'],
  ['C', 6, 'Carbon', 12.011, 2.55, [0.75, 0.67, 0.60], 1.70, 3.851, 0.105, [4], 4, 3, 260, '#a9adb5'],
  ['N', 7, 'Nitrogen', 14.007, 3.04, [0.71, 0.60, 0.54], 1.55, 3.660, 0.069, [3], 5, 3, 250, '#5b8cff'],
  ['O', 8, 'Oxygen', 15.999, 3.44, [0.63, 0.57, 0.53], 1.52, 3.500, 0.060, [2], 6, 2, 200, '#ff4d3d'],
  ['F', 9, 'Fluorine', 18.998, 3.98, [0.64, 0.59, 0.53], 1.47, 3.364, 0.050, [1], 7, 1, 0, '#90e050'],
  ['Ne', 10, 'Neon', 20.180, 0, [0.67], 1.54, 3.243, 0.042, [0], 8, 0, 0, '#b3e3f5'],
  ['Na', 11, 'Sodium', 22.990, 0.93, [1.55, 1.60], 2.27, 2.983, 0.030, [1], 1, 1, 60, '#ab5cf2'],
  ['Mg', 12, 'Magnesium', 24.305, 1.31, [1.39, 1.32, 1.27], 1.73, 3.021, 0.111, [2], 2, 2, 100, '#8aff00'],
  ['Al', 13, 'Aluminium', 26.982, 1.61, [1.26, 1.13, 1.11], 1.84, 4.499, 0.505, [3], 3, 2, 150, '#bfa6a6'],
  ['Si', 14, 'Silicon', 28.085, 1.90, [1.16, 1.07, 1.02], 2.10, 4.295, 0.402, [4], 4, 2, 180, '#e0c09a'],
  ['P', 15, 'Phosphorus', 30.974, 2.19, [1.11, 1.02, 0.94], 1.80, 4.147, 0.305, [3, 5], 5, 3, 180, '#ff8a1f'],
  ['S', 16, 'Sulfur', 32.06, 2.58, [1.03, 0.94, 0.95], 1.80, 4.035, 0.274, [2, 4, 6], 6, 2, 180, '#ffd23f'],
  ['Cl', 17, 'Chlorine', 35.45, 3.16, [0.99, 0.95, 0.93], 1.75, 3.947, 0.227, [1, 3, 5, 7], 7, 1, 150, '#3fd06a'],
  ['Ar', 18, 'Argon', 39.948, 0, [0.96], 1.88, 3.868, 0.185, [0], 8, 0, 0, '#80d1e3'],
  ['K', 19, 'Potassium', 39.098, 0.82, [1.96, 1.93], 2.75, 3.812, 0.035, [1], 1, 1, 60, '#8f40d4'],
  ['Ca', 20, 'Calcium', 40.078, 1.00, [1.71, 1.47, 1.33], 2.31, 3.399, 0.238, [2], 2, 2, 100, '#3dff00'],
  ['Fe', 26, 'Iron', 55.845, 1.83, [1.16, 1.09, 1.02], 2.04, 2.912, 0.013, [2, 3], 2, 2, 100, '#e06633'],
  ['Cu', 29, 'Copper', 63.546, 1.90, [1.12, 1.15, 1.20], 1.40, 3.495, 0.005, [1, 2], 1, 1, 100, '#c88033'],
  ['Zn', 30, 'Zinc', 65.38, 1.65, [1.18, 1.20], 1.39, 2.763, 0.124, [2], 2, 1, 100, '#7d80b0'],
  ['Br', 35, 'Bromine', 79.904, 2.96, [1.14, 1.09, 1.10], 1.85, 4.189, 0.251, [1, 3, 5], 7, 1, 150, '#c9563a'],
  ['Kr', 36, 'Krypton', 83.798, 3.00, [1.17, 1.21, 1.08], 2.02, 4.141, 0.220, [0, 2], 8, 1, 0, '#5cb8d1'],
  ['I', 53, 'Iodine', 126.90, 2.66, [1.33, 1.29, 1.25], 1.98, 4.500, 0.339, [1, 3, 5, 7], 7, 1, 150, '#a064d8'],
  ['Xe', 54, 'Xenon', 131.29, 2.60, [1.31, 1.35, 1.22], 2.16, 4.404, 0.332, [0, 2, 4, 6], 8, 1, 100, '#429eb0']
];
// Homonuclear single-bond dissociation energies for fallback (Pauling-type) estimates, kJ/mol
const D_HOMO = { H: 436, Li: 105, Be: 60, B: 293, C: 348, N: 163, O: 146, F: 155, Na: 72, Mg: 30, Al: 186, Si: 222, P: 201, S: 266, Cl: 242, K: 55, Ca: 40, Fe: 100, Cu: 190, Zn: 30, Br: 193, Kr: 0, I: 151, Xe: 0 };

const ELEMENTS = ELEMENT_ROWS.map((r, t) => ({
  t, sym: r[0], Z: r[1], name: r[2], mass: r[3], chi: r[4], rcov: r[5], rvdw: r[6],
  ljX: r[7], ljD: r[8] * 4.184, valences: r[9], ve: r[10], maxOrder: r[11], kAngle: r[12], color: r[13],
  bonds: r[9].some(v => v > 0)
}));
const BY_SYM = Object.create(null);
for (const e of ELEMENTS) BY_SYM[e.sym] = e;
const NT = ELEMENTS.length;
const CHI = Float64Array.from(ELEMENTS, e => e.chi);
// chalcogens and halogens: with one of these in the formula, hydrogen is written first
const INS_MAX = 0.6;                 // how much of a sigma bond's competition a carbene can take over
const INS_ON = 0.8, INS_OFF = 1.8;   // Å past the bond length: fully engaged, and first felt
const INS_T1 = 0.8, INS_T2 = 1.15;   // an end's new bond + old bond order: lends fully, lends nothing
const INS_Q0 = 0.05, INS_Q1 = 0.15;  // the old bond's reach: not a bond, a bond
const INS_QC = 0.3;                  // Å: Pauling length of the old bond's order
const INS_U0 = 0.9, INS_U1 = 0.99;   // a new bond's reach: still forming, formed
const INS_REC = 35;                  // slots per insertion term
const INSERTS = new Uint8Array(ELEMENT_ROWS.length);
for (const sym of ['C', 'Si']) { const k = ELEMENT_ROWS.findIndex(r => r[0] === sym); if (k >= 0) INSERTS[k] = 1; }
const H_LEADS = new Set(['O', 'S', 'Se', 'Te', 'F', 'Cl', 'Br', 'I']);
const TETRA = [-1 / 3, -0.2924, -0.2504, -1 / 3]; // cos θ0 for steric number 4 with 0, 1, 2, 3 lone pairs (109.5°, 107°, 104.5°)

/* ---------- measured bonds: order → [re Å, De kJ/mol, k N/m] ---------- */
const BOND_DATA = {
  'H|H': { 1: [0.741, 436, 575] }, 'C|H': { 1: [1.09, 413, 500] }, 'H|N': { 1: [1.01, 391, 640] },
  'H|O': { 1: [0.96, 463, 780] }, 'F|H': { 1: [0.92, 567, 966] }, 'Cl|H': { 1: [1.27, 431, 516] },
  'Br|H': { 1: [1.41, 366, 412] }, 'H|I': { 1: [1.61, 299, 314] }, 'H|S': { 1: [1.34, 363, 420] },
  'H|P': { 1: [1.42, 322, 330] }, 'H|Si': { 1: [1.48, 318, 290] }, 'B|H': { 1: [1.19, 389, 350] },
  'H|Li': { 1: [1.60, 243, 103] }, 'H|Na': { 1: [1.89, 186, 78] }, 'H|K': { 1: [2.24, 175, 56] },
  'C|C': { 1: [1.54, 348, 450], 2: [1.34, 614, 950], 3: [1.20, 839, 1600] },
  'C|N': { 1: [1.47, 293, 490], 2: [1.28, 615, 1000], 3: [1.16, 891, 1790] },
  'C|O': { 1: [1.43, 358, 540], 2: [1.19, 780, 1250], 3: [1.13, 1072, 1900] }, // C=O: between ketone (745) and CO₂ (804)
  'C|F': { 1: [1.35, 485, 590] }, 'C|Cl': { 1: [1.77, 328, 340] }, 'Br|C': { 1: [1.94, 276, 290] },
  'C|I': { 1: [2.14, 240, 230] }, 'C|S': { 1: [1.82, 272, 300], 2: [1.60, 573, 700] },
  'C|Si': { 1: [1.87, 318, 300] }, 'C|P': { 1: [1.84, 264, 300], 2: [1.67, 513, 600] }, 'B|C': { 1: [1.56, 356, 380] },
  'N|N': { 1: [1.45, 163, 350], 2: [1.25, 418, 1000], 3: [1.10, 945, 2295] },
  'N|O': { 1: [1.40, 201, 400], 2: [1.21, 607, 1550] }, 'O|O': { 1: [1.48, 146, 400], 2: [1.21, 498, 1177] },
  'F|F': { 1: [1.42, 155, 470] }, 'Cl|Cl': { 1: [1.99, 242, 323] }, 'Br|Br': { 1: [2.28, 193, 246] },
  'I|I': { 1: [2.67, 151, 172] }, 'S|S': { 1: [2.05, 266, 250], 2: [1.89, 425, 500] },
  'O|S': { 1: [1.58, 265, 500], 2: [1.43, 522, 1000] }, 'O|P': { 1: [1.63, 335, 450], 2: [1.48, 544, 900] },
  'O|Si': { 1: [1.63, 452, 500], 2: [1.51, 640, 900] }, 'Si|Si': { 1: [2.33, 222, 170] },
  'P|P': { 1: [2.21, 201, 200], 2: [2.00, 351, 400], 3: [1.89, 489, 556] },
  'F|N': { 1: [1.36, 272, 450] }, 'Cl|N': { 1: [1.75, 200, 300] }, 'F|O': { 1: [1.42, 190, 400] },
  'Cl|O': { 1: [1.70, 203, 350] }, 'F|S': { 1: [1.56, 327, 450] }, 'Cl|S': { 1: [2.07, 271, 250] },
  'Cl|P': { 1: [2.04, 326, 250] }, 'F|P': { 1: [1.54, 490, 480] }, 'F|Si': { 1: [1.57, 565, 500] },
  'Cl|Si': { 1: [2.02, 381, 290] }, 'B|F': { 1: [1.31, 613, 600] }, 'B|O': { 1: [1.36, 536, 600], 2: [1.20, 800, 1100] },
  'Cl|F': { 1: [1.63, 253, 440] }, 'Br|Cl': { 1: [2.14, 218, 280] }, 'H|Xe': { 1: [1.7, 30, 60] },
  'Cl|Na': { 1: [2.36, 410, 110] }, 'F|Na': { 1: [1.93, 477, 180] }, 'Cl|K': { 1: [2.67, 433, 86] },
  'F|K': { 1: [2.17, 498, 140] }, 'F|Li': { 1: [1.56, 577, 250] }, 'Cl|Li': { 1: [2.02, 469, 143] },
  'Br|Na': { 1: [2.50, 363, 95] }, 'I|Na': { 1: [2.71, 304, 76] }, 'Li|O': { 1: [1.69, 341, 250] },
  'Na|O': { 1: [2.05, 270, 180] }, 'K|O': { 1: [2.17, 278, 150] },
  'Mg|O': { 1: [1.75, 363, 350], 2: [1.75, 363, 350] }, 'Ca|O': { 1: [1.82, 383, 300], 2: [1.82, 383, 300] },
  'Cl|Mg': { 1: [2.18, 312, 200] }, 'Ca|Cl': { 1: [2.44, 409, 180] }, 'F|Mg': { 1: [1.75, 463, 300] },
  'Al|O': { 1: [1.62, 502, 450], 2: [1.62, 502, 450] }, 'Al|Cl': { 1: [2.13, 502, 300] },
  'Fe|O': { 1: [1.62, 407, 480], 2: [1.62, 407, 480] }, 'Cu|O': { 1: [1.72, 287, 260] }, 'Zn|O': { 1: [1.70, 250, 300] },
  'Cl|Fe': { 1: [2.10, 330, 250] }, 'Cl|Cu': { 1: [2.05, 378, 230] }, 'Cl|Zn': { 1: [2.07, 229, 220] }
};

function pairKey(a, b) { return a < b ? a + '|' + b : b + '|' + a; }

function buildPairParams(ei, ej) {
  const canBond = ei.bonds && ej.bonds;
  const pp = {
    bond: canBond, maxOrder: Math.max(1, Math.min(ei.maxOrder, ej.maxOrder)),
    re: [0, 0, 0, 0], De: [0, 0, 0, 0], a: [0, 0, 0, 0], r1: 0, r2: 0,
    ljX: Math.sqrt(ei.ljX * ej.ljX), ljD: Math.sqrt(ei.ljD * ej.ljD),
    core: 0.6 * (ei.rcov[0] + ej.rcov[0]),
    oo: ei.sym === 'O' && ej.sym === 'O'
  };
  if (canBond) {
    const data = BOND_DATA[pairKey(ei.sym, ej.sym)] || {};
    const dchi = Math.abs(ei.chi - ej.chi);
    const d1 = data[1] ? data[1][1] : 0.5 * ((D_HOMO[ei.sym] || 150) + (D_HOMO[ej.sym] || 150)) + 60 * dchi * dchi;
    for (let n = 1; n <= 3; n++) {
      const rec = data[n];
      let re, De, k;
      if (rec) { re = rec[0]; De = rec[1]; k = rec[2]; }
      else {
        const ri = ei.rcov[n - 1] ?? ei.rcov[ei.rcov.length - 1] - 0.08 * (n - ei.rcov.length);
        const rj = ej.rcov[n - 1] ?? ej.rcov[ej.rcov.length - 1] - 0.08 * (n - ej.rcov.length);
        re = ri + rj - 0.05 * dchi;
        De = d1 * [1, 1.75, 2.4][n - 1];
        k = null;
      }
      const a = k ? Math.sqrt(k * 6.02214 / (2 * De)) : 1.95 + 0.3 * (n - 1);
      pp.re[n] = re; pp.De[n] = De; pp.a[n] = a;
    }
    pp.r1 = COORD_R1 * pp.re[1];
    pp.r2 = COORD_R2 * pp.re[1];
    pp.s1 = pp.re[1];
  }
  return pp;
}
const PAIR = new Array(NT * NT);
for (let i = 0; i < NT; i++) for (let j = 0; j < NT; j++) PAIR[i * NT + j] = buildPairParams(ELEMENTS[i], ELEMENTS[j]);
/* Saturation bond order decays like the pair's own Morse attraction, β = kb·a, so a stretched bond
   counts as much used valence as the attraction it still provides; window ends where s ≈ e^(−yoff). */
function refreshSaturation() {
  for (const pp of PAIR) {
    if (!pp.bond) continue;
    pp.sBeta = TUNE.kb * pp.a[1];
    pp.sOff = TUNE.cap + TUNE.yoff / pp.sBeta; pp.sOn = TUNE.cap + 0.6 * TUNE.yoff / pp.sBeta;
    pp.s2 = pp.s1 + pp.sOff;
  }
}
refreshSaturation();

/* ---------- helpers ---------- */
function mulberry(state) { // returns [value, newState]
  let t = (state + 0x6D2B79F5) >>> 0;
  let x = Math.imul(t ^ (t >>> 15), t | 1);
  x ^= x + Math.imul(x ^ (x >>> 7), x | 61);
  return [((x ^ (x >>> 14)) >>> 0) / 4294967296, t];
}
// Hot helpers return through module scratch variables (no allocation in the force loop).
let SW = 0, SWD = 0;          // smoothSwitch value, derivative
function smoothSwitch(x, x1, x2) { // 1 → 0 between x1 and x2
  if (x <= x1) { SW = 1; SWD = 0; return; }
  if (x >= x2) { SW = 0; SWD = 0; return; }
  const w = x2 - x1, t = (x - x1) / w;
  SW = 0.5 * (1 + Math.cos(Math.PI * t)); SWD = -0.5 * Math.PI * Math.sin(Math.PI * t) / w;
}
/* Pauling-type bond order used for valence saturation: 1 near the bond length, exp(−β·Δr)
   when stretched, exactly zero BO_OFF Å past it. Longer-ranged than the structural switch so a
   partner that is half-way into a bond already competes for valence. Output in SF, SFD. */
let SF = 0, SFD = 0;
function satF(r, pp) {
  const d = r - pp.s1;
  const cap = TUNE.cap;
  if (d <= cap) { SF = 1; SFD = 0; return; }
  if (d >= pp.sOff) { SF = 0; SFD = 0; return; }
  const x = d - cap, beta = pp.sBeta;
  let s, ds;
  if (x < BO_W) { s = x * x / (2 * BO_W); ds = x / BO_W; } else { s = x - BO_W / 2; ds = 1; }
  const e = Math.exp(-beta * s), de = -beta * ds * e;
  smoothSwitch(d, pp.sOn, pp.sOff);
  SF = e * SW; SFD = de * SW + e * SWD;
}
/* The multiplicity a pair counts with in coordination: sigma plus a partial pi, as Pass A does. */
function nEffOf(ti, tj, n) { return 1 + (PAIR[ti * NT + tj].oo ? TUNE.wpiOO : TUNE.wpi) * (n - 1); }
/* How much of a carbene an atom is: 1 with two spare valences or more, fading to 0 by 1.5, so a
   radical with one (CH3) is out while CH2, CH and a bare carbon atom are in. Smooth, so the forces
   are too. */
let HW = 0, HWD = 0;
function carbeneWeight(spare) {
  // two free valences or more once the two new bonds are set aside: CH2, CH, a bare C
  const lo = (spare - 1.5) / 0.5;
  if (lo <= 0) { HW = 0; HWD = 0; return 0; }
  if (lo >= 1) { HW = 1; HWD = 0; return 1; }
  HW = lo * lo * (3 - 2 * lo); HWD = 6 * lo * (1 - lo) / 0.5; return HW;
}
/* ...and carries no multiple bond. CO relaxes to C=O here, which also leaves carbon two spare
   valences, but beside a double bond — and CO does not insert into H2 (the barrier is some
   330 kJ/mol). CH2, CH and a bare carbon atom carry none, and all three insert. */
let GW = 0, GWD = 0;
function piWindow(pi) {
  // no multiple bond: full up to 0.2 of a π bond, gone by 0.5
  const t = (0.5 - pi) / 0.3;
  if (t >= 1) { GW = 1; GWD = 0; return 1; }
  if (t <= 0) { GW = 0; GWD = 0; return 0; }
  GW = t * t * (3 - 2 * t); GWD = -6 * t * (1 - t) / 0.3; return GW;
}
/* How much room a radical's unpaired electrons take around it, as extra steric number.
   Counting them as nothing — which is what the angles did — is right for CH3 and wrong for almost
   everything else. CH3's one electron sits in the p orbital above its plane of three bonds, so it
   stays flat. But an atom that already has a lone pair shares its nonbonding space with the
   unpaired electron just as water shares it between two lone pairs: NH2 is bent to 103°, like
   water, not 120°. And an atom with no lone pair and fewer than three bonds puts an unpaired
   electron into the open place in its plane, where it pushes about two thirds as hard as a bond:
   CH2 is bent to about 134° and BH2 to 129°, not straight. So:
     with a lone pair:    each unpaired electron is one more nonbonding domain, and the tetrahedral
                          angle follows the nonbonding count (NH2 -> the water angle);
     without one:         0.85 of a domain for each unpaired electron, up to the room left in the
                          plane (3 − bonds) — so CH3, with three bonds, gets none and stays flat.
   A closed shell has no free valence and is exactly as before. Smooth throughout, so the forces
   are the exact gradient. RD_E is the extra steric number, RD_EX its derivative in the free
   valence x = V − Zs, RD_EZ in the structural coordination Z; RD_C4 the tetrahedral-limit cos θ,
   RD_C4X its derivative in x. */
let RD_E = 0, RD_EX = 0, RD_EZ = 0, RD_C4 = 0, RD_C4X = 0;
const RD_W = 0.2, RD_SIGMA = 0.85, RD_P = 4;
function rdRamp(x) { if (x <= 0) return [0, 0]; if (x < RD_W) return [x * x / (2 * RD_W), x / RD_W]; return [x - RD_W / 2, 1]; }
function radicalDomains(x, Z, lp) {
  const [u, du] = rdRamp(x);
  if (lp >= 1) {
    RD_E = u; RD_EX = du; RD_EZ = 0;
    const L = lp + u;
    if (L >= 3) { RD_C4 = TETRA[3]; RD_C4X = 0; }
    else { const k = Math.floor(L), f = L - k, d = TETRA[k + 1] - TETRA[k]; RD_C4 = TETRA[k] + d * f; RD_C4X = d * du; }
    return;
  }
  RD_C4 = TETRA[0]; RD_C4X = 0;
  const [m, dm] = rdRamp(3 - Z);                     // the room left in the plane
  if (u <= 0 || m <= 0) { RD_E = 0; RD_EX = 0; RD_EZ = 0; return; }
  // a smooth min(u, m)
  const g = Math.pow(Math.pow(u, -RD_P) + Math.pow(m, -RD_P), -1 / RD_P);
  const gu = Math.pow(g / u, RD_P + 1), gm = Math.pow(g / m, RD_P + 1);
  RD_E = RD_SIGMA * g; RD_EX = RD_SIGMA * gu * du; RD_EZ = -RD_SIGMA * gm * dm;
}
/* Valence-state energy. A bond's energy is not the same whichever bond of the atom it is: carbon's
   four C–H bonds cost 439, 462, 424 and 338 kJ/mol to break in turn (CH4 → CH3 → CH2 → CH → C),
   nitrogen's 453, 391, 332, oxygen's 497 then 430. The atom rehybridises as it goes — promotion,
   spin pairing, lone pairs relaxing — and a partly bonded atom pays for being half-way. A pair
   potential gives every bond the average, so every radical came out too stable: OH by 33 kJ/mol,
   NH2 by 63, CH2 by 61, CH by 75. This adds, per atom, the difference between the average and the
   stepwise ladder, as a smooth function of its bonding Z (full multiplicity: N≡N counts three).
   Zero at no bonds and at a full shell, so atoms, every closed-shell molecule and every fitted
   bond energy are untouched. Tabulated at whole Z from the hydride ladders (ΔfH, NIST-JANAF),
   joined by a cubic that is flat at both ends. */
const VSTATE = { C: [0, 73, 59, 18, 0], N: [0, 51, 55, 0], O: [0, 32, 0], S: [0, 17, 0] };
const VS_TAB = new Array(ELEMENT_ROWS.length).fill(null);
for (const [sym, tab] of Object.entries(VSTATE)) { const k = ELEMENT_ROWS.findIndex(r => r[0] === sym); if (k >= 0) VS_TAB[k] = tab; }
let VS = 0, VSD = 0;
function valenceState(tab, z) {
  const V = tab.length - 1;
  if (z <= 0 || z >= V) { VS = 0; VSD = 0; return; }
  const k = Math.floor(z), t = z - k, p0 = tab[k], p1 = tab[k + 1];
  const m0 = k === 0 ? 0 : (tab[k + 1] - tab[k - 1]) / 2, m1 = k + 1 === V ? 0 : (tab[k + 2] - tab[k]) / 2;
  const t2 = t * t, t3 = t2 * t;
  VS = (2 * t3 - 3 * t2 + 1) * p0 + (t3 - 2 * t2 + t) * m0 + (-2 * t3 + 3 * t2) * p1 + (t3 - t2) * m1;
  VSD = (6 * t2 - 6 * t) * p0 + (3 * t2 - 4 * t + 1) * m0 + (-6 * t2 + 6 * t) * p1 + (3 * t2 - 2 * t) * m1;
}
/* Electronic degeneracy. A free atom has several electronic states at (nearly) the same energy — a
   halogen atom 4 (²P3/2, and 2 more ²P1/2 when hot), oxygen 9 (³P), hydrogen 2 — while the molecule
   it forms has one. Those states are entropy, R·ln g, and they are why real F2, Cl2, O2, I2 and H2
   come apart at flame temperatures 20–30 times more readily than a single potential surface says
   (measured against textbook statistical mechanics with spectroscopic constants: the whole gap was
   this factor). The equilibrium of the electronically averaged system is sampled exactly on the
   free-energy surface U − kT·ln g, with g following the atom's bonding Z smoothly between the free
   atom and its radicals (OH ²Π 4, CH ²Π 4, CH2 ³B1 3, NH ³Σ 3, CH3 and NH2 doublets 2) and 1 for a
   full shell. It vanishes at 0 K and for every closed-shell molecule. Levels: NIST ASD. */
const GEL = {
  H: T => [Math.log(2), 0],
  C: T => [Math.log(1 + 3 * Math.exp(-23.6 / T) + 5 * Math.exp(-62.4 / T)), Math.log(4), Math.log(3), Math.log(2), 0],
  N: T => [Math.log(4), Math.log(3), Math.log(2), 0],
  O: T => [Math.log(5 + 3 * Math.exp(-227.7 / T) + Math.exp(-326.6 / T)), Math.log(4), 0],
  F: T => [Math.log(4 + 2 * Math.exp(-581.0 / T)), 0],
  S: T => [Math.log(5 + 3 * Math.exp(-570.0 / T) + Math.exp(-825.0 / T)), Math.log(4), 0],
  Cl: T => [Math.log(4 + 2 * Math.exp(-1269.6 / T)), 0],
  Br: T => [Math.log(4 + 2 * Math.exp(-5302 / T)), 0],
  I: T => [Math.log(4), 0],
};
// 1 for an oxygen whose bonding is O2's own (Z = 2), fading to 0 by Z = 1 (coming apart) and Z = 2.6 (a third partner)
let OB = 0, OBD = 0;
function o2Bump(z) {
  if (z <= 1 || z >= 2.6) { OB = 0; OBD = 0; return; }
  if (z <= 2) { const t = z - 1; OB = t * t * (3 - 2 * t); OBD = 6 * t * (1 - t); return; }
  const t = (2.6 - z) / 0.6; OB = t * t * (3 - 2 * t); OBD = -6 * t * (1 - t) / 0.6;
}
const GEL_TAB = new Array(ELEMENT_ROWS.length).fill(null);
for (const [sym, f] of Object.entries(GEL)) { const k = ELEMENT_ROWS.findIndex(r => r[0] === sym); if (k >= 0) GEL_TAB[k] = f; }
/* An end with nothing to spare counts fully; one with 0.75 of a valence free is a ring partner. */
/* ...and an end lends only the valence it has. A hydrogen has one bond to give: while its reach
   to the carbene and to its partner sum to about one, it is on its way from one to the other and
   the carbene may take it over; past that it would be holding two bonds at once. Without this the
   screen let a hydrogen keep its H–H bond at full strength *and* bond to carbon, and CH2 + H2
   stuck as a CH2·H2 complex 550 kJ/mol more stable than it has any right to be. */
let CW = 0, CWD = 0;
function conserveWindow(T) {
  const t = (INS_T2 - T) / (INS_T2 - INS_T1);
  if (t >= 1) { CW = 1; CWD = 0; return; }
  if (t <= 0) { CW = 0; CWD = 0; return; }
  CW = t * t * (3 - 2 * t); CWD = -6 * t * (1 - t) / (INS_T2 - INS_T1);
}
/* Whether the two ends share a bond at all, by its saturation reach rather than its structural
   switch: the switch is gone by 1.1 Å of H–H while the bond still competes out to 1.6, and a screen
   that quit with the switch handed the competition back mid-stretch. Starting at a reach of 0.05
   keeps the 1-3 neighbours of a settled molecule out (geminal H in methane: 0.045). */
/* ...and only while the new bonds are still forming. Once both are whole the job is done, and a
   settled carbon — whose two bonds, set aside, would leave it two spare valences like any carbene —
   is exactly untouched however hot its neighbours get. */
let UW = 0, UWD = 0;
function formingWindow(x) {
  const t = (INS_U1 - x) / (INS_U1 - INS_U0);
  if (t >= 1) { UW = 1; UWD = 0; return; }
  if (t <= 0) { UW = 0; UWD = 0; return; }
  UW = t * t * (3 - 2 * t); UWD = -6 * t * (1 - t) / (INS_U1 - INS_U0);
}
let QW = 0, QWD = 0;
function bondGate(x) {
  const t = (x - INS_Q0) / (INS_Q1 - INS_Q0);
  if (t >= 1) { QW = 1; QWD = 0; return; }
  if (t <= 0) { QW = 0; QWD = 0; return; }
  QW = t * t * (3 - 2 * t); QWD = 6 * t * (1 - t) / (INS_Q1 - INS_Q0);
}
let EW = 0, EWD = 0;
function endWindow(se) {
  const t = (0.75 - se) / 0.4;
  if (t >= 1) { EW = 1; EWD = 0; return; }
  if (t <= 0) { EW = 0; EWD = 0; return; }
  EW = t * t * (3 - 2 * t); EWD = -6 * t * (1 - t) / 0.4;
}
let SP = 0, SPD = 0;          // saturation p(x), dp/dx
function sat(x) { // p(x) = 1 / (1 + c1·r + c4·r⁴), r = smooth ramp(x)
  if (x <= 0) { SP = 1; SPD = 0; return; }
  let r, dr;
  if (x < RAMP_W) { r = x * x / (2 * RAMP_W); dr = x / RAMP_W; } else { r = x - RAMP_W / 2; dr = 1; }
  const c1 = TUNE.c1, c4 = TUNE.c4, k = TUNE.k, r2 = r * r, rk1 = k === 6 ? r2 * r2 * r : Math.pow(r, k - 1), den = 1 + c1 * r + c4 * rk1 * r;
  SP = 1 / den; SPD = -(c1 + k * c4 * rk1) * dr / (den * den);
}
/* Bond-strength weight of a pair at bond order n: (De(n)/400)^μ. */
function bondWeight(pp, n, mu) {
  if (mu === 0) return 1;
  const lo = Math.min(2, Math.floor(n)), t = n - lo, De = pp.De[lo] + (pp.De[Math.min(lo + 1, 3)] - pp.De[lo]) * t;
  return mu === 0.5 ? Math.sqrt(De / 400) : mu === 1 ? De / 400 : Math.pow(De / 400, mu);
}
/* Strength-weighted excess valence: x = (C + 1 − V)·ratio, ratio = (W + εD)/((C + ε)·D), where C and W
   are the plain and strength-weighted valence used by the *other* bonds. XA = ∂x/∂u_k − D_k·XB, XB = ∂x/∂(D_k·u_k). */
let XS = 0, XA = 1, XB = 0;
const EX_EPS = 0.05;
function excess(C, W, V, D, mu) {
  const base = C + 1 - V;
  if (mu === 0) { XS = base; XA = 1; XB = 0; return; }
  const den = (C + EX_EPS) * D, ratio = (W + EX_EPS * D) / den;
  XS = base * ratio; XA = ratio - base * ratio / (C + EX_EPS); XB = base / den;
}
/* Shielded Lennard-Jones split into Pauli (LR) and shifted-force attraction (LA), each with d/dr.
   Shielding: r_s⁶ = r⁶ + s⁶, so the Pauli wall stays finite at contact. x2 = r_min², eps = depth. */
let LR = 0, LDR = 0, LA = 0, LDA = 0, LXR = 0, LXA = 0;
const S6 = Math.pow(LJ_SHIELD, 6);
function lj(r, x2, eps, rc) {
  const x6 = x2 * x2 * x2, r2 = r * r, r6 = r2 * r2 * r2, d6 = r6 + S6 * x6;
  const u6 = x6 / d6, du6 = -6 * u6 * r2 * r2 * r / d6; // d(u⁶)/dr
  const rc2 = rc * rc, rc6 = rc2 * rc2 * rc2, uc6 = x6 / (rc6 + S6 * x6);
  const ljc = eps * (uc6 * uc6 - 2 * uc6), dljc = eps * (2 * uc6 - 2) * (-6 * uc6 * rc2 * rc2 * rc / (rc6 + S6 * x6));
  const v = eps * (u6 * u6 - 2 * u6), dv = eps * (2 * u6 - 2) * du6;
  // Derivatives with respect to x2 are needed when a polar H radius changes with charge.
  const dx6 = 3 * x2 * x2, dc6 = rc6 + S6 * x6;
  const ux = dx6 * r6 / (d6 * d6), ucx = dx6 * rc6 / (dc6 * dc6);
  const ucr = -6 * rc2 * rc2 * rc * x6 / (dc6 * dc6);
  const ucrx = -6 * rc2 * rc2 * rc * dx6 * (rc6 - S6 * x6) / (dc6 * dc6 * dc6);
  const vcx = eps * (2 * uc6 - 2) * ucx;
  const dvcx = eps * (2 * ucx * ucr + (2 * uc6 - 2) * ucrx);
  const vx = eps * (2 * u6 - 2) * ux;
  if (u6 > 1) { // inside r_min: Pauli wall, attraction held at its (shifted) value at r_min
    LR = v + eps; LDR = dv;
    const rm = Math.sqrt(Math.cbrt(x6 - S6 * x6));
    LA = -eps - ljc - dljc * (rm - rc); LDA = 0;
    LXR = vx; LXA = -vcx - dvcx * (rm - rc) - dljc * rm / (2 * x2);
  } else { LR = 0; LDR = 0; LA = v - ljc - dljc * (r - rc); LDA = dv - dljc; LXR = 0; LXA = vx - vcx - dvcx * (r - rc); }
}

/* Effective valence and lone pairs for an element carrying a formal charge. */
function valenceFor(el, charge, want) {
  let V = want != null ? want : el.valences[0];
  if (!el.bonds) return 0;
  if (charge > 0) {
    if (el.ve >= 5 && el.ve <= 6) V += charge;          // N⁺, O⁺, S⁺, P⁺: extra bond
    else V = Math.max(0, V - charge);                   // metals, H⁺, C⁺: fewer bonds
  } else if (charge < 0) {
    if (el.ve === 3) V -= charge;                       // B⁻ (borate): four bonds
    else V = Math.max(0, V + charge);                   // O⁻, N⁻, halide: fewer bonds
  }
  return V;
}
function lonePairs(el, charge, V) {
  if (el.ve >= 8 || el.Z <= 2) return 0;
  return Math.max(0, Math.floor((el.ve - charge - V) / 2));
}

class Engine {
  constructor(opts = {}) {
    this.dt = 1.0;
    this.box = { x0: 0, x1: opts.width || 50, y0: 0, y1: opts.height || 32, z0: -(opts.depth || 10) / 2, z1: (opts.depth || 10) / 2 };
    this.sphere = null;             // {x,y,z,R} spherical container (used by the conditioning chamber)
    this.wallK = 60;                // kJ/mol/Å²
    // Boundary character. 'solid' is a stiff wall that starts exactly at the chamber face;
    // 'forcefield' is a soft cushion that begins fieldRange inside it and keeps growing outward,
    // so an atom is turned around gradually and may briefly cross the face.
    this.boundsMode = opts.boundsMode || 'solid';
    // Soft enough, and reaching in little enough, that ordinary thermal atoms still lean past
    // the face — otherwise a forcefield would hide the chamber's outside from every void wall.
    this.fieldK = opts.fieldK ?? 5;           // kJ/mol/Å², forcefield stiffness
    this.solidK = opts.solidK ?? 4000;        // kJ/mol/Å², a solid face: about as stiff as a C–H bond
    this.fieldRange = opts.fieldRange ?? 1.2; // Å the forcefield reaches inward
    // Void wall: what the region beyond the chamber face does to what reaches it.
    /* Void walls: what the boundary refuses to hand back. Each channel deletes one quantity the
       wall would otherwise return, so the chamber is open in that respect and nothing accumulates.
       They work the same for solid and soft bounds, because the contact shell — not the face —
       is where a wall takes its bite. */
    this.voidTemperature = opts.voidTemperature ?? false; // agitation relative to the atom's own molecule
    this.voidPressure = opts.voidPressure ?? false;       // the normal momentum, and the impulse it would report
    this.voidVelocity = opts.voidVelocity ?? false;       // the atom's motion entirely: it stops at the wall
    // 300 fs: firm enough that a chamber left alone settles instead of heating itself,
    // gentle enough that the wall heater can still hold a temperature against it.
    this.voidTau = opts.voidTau ?? 300;     // fs, absorption time at full contact
    this.voidSkin = opts.voidSkin ?? 0.8;   // Å of contact shell measured inward from the face
    this.voidHeat = 0; this.voidForce = 0;  // removed energy; voidForce retained only for legacy snapshots
    /* Optional pressure control. Passive measurement is always reported; when this is on, the
       chamber breathes toward the target the way a Berendsen barostat does, moving whole
       molecules rather than scaling the bonds inside them. */
    this.pressureControl = opts.pressureControl ?? false;
    this.pressureTarget = opts.pressureTarget ?? 1.0;   // bar
    this.pressureTau = opts.pressureTau ?? 4000;        // fs
    this.baroEvery = 100;                               // steps between adjustments
    this.T = opts.T ?? 298.15;
    this.tau = opts.tau ?? 250;     // thermostat coupling time, fs
    this.thermostat = opts.thermostat ?? true;
    // The sample exchanges heat only in a thin boundary layer. CSVR remains available
    // for preparing isolated molecules, not as the laboratory's thermostat.
    /* 'wall'   — a heater at the boundary; the interior warms through collisions.
       'kelvin' — surroundings at the setpoint touching every atom (a Langevin bath).
       'csvr'   — canonical sampling, kept for preparing isolated molecules.
       Kelvin is the default because it is the one that lets a molecule be built. A bond forming
       in open space releases its binding energy on the spot, and a boundary heater has no way to
       take that away from the middle of the chamber, so the molecule you just made blows itself
       apart — which is correct physics (recombination needs a third body) and useless as a
       default. The Kelvin stat is that third body, everywhere. */
    this.thermostatMode = opts.thermostatMode || 'kelvin';
    this.kelvinWork = 0;            // kJ/mol the Kelvin stat has put in (or taken out)
    this.wallT = opts.wallT ?? this.T;
    this.wallTarget = opts.wallTarget ?? Math.max(288.15, Math.min(623.15, this.T));
    this.wallTau = opts.wallTau ?? 10000; // fs, C_wall / conductance of heater
    this.wallCapacity = opts.wallCapacity ?? 1000 * KB; // kJ/(mol K), effective reservoir
    this.wallSkin = opts.wallSkin ?? 2; // Å; zero direct coupling in the interior
    this.wallCoupling = opts.wallCoupling ?? 100; // fs at the surface
    this.heatToSample = 0; this.heaterWork = 0;
    /* Ignition window. A spark is local and fast; coupling to the surroundings is neither. In a
       chamber of a dozen atoms a thermostat is instantaneous and global, so it would erase a
       spark on the step it landed and nothing could ever be lit. While this window is open the
       stat and the temperature void stand back and let the spark do its work. */
    this.sparkHold = 0;
    this.servoWork = 0; this.servoWorkTotal = 0;     // work the pointer has done on the sample
    this.wallMeasured = this.T; this.wallContact = 0; // what the fluid against the wall actually is
    /* Interaction cutoff. The electrostatics here are damped and already short-ranged, so 6.5 Å
       reaches everything that is worth more than a few kJ/mol: measured against the whole
       validation table, moving down from 8 Å changes no bond length by more than 0.02 Å and no
       reaction energy by more than 0.4 kJ/mol, brings the methane dimer closer to experiment, and
       leaves a third of the pairs — and a third of every step — behind. */
    this.rc = opts.rc || 6.5;
    this.skin = 1.0;
    this.rngState = (opts.seed ?? 12345) >>> 0;
    this.N = 0; this.cap = 0;
    this.time = 0; this.stepCount = 0;
    this.nextId = 1;
    this.Epot = 0; this.Ewall = 0; this.wallForce = 0; this.pressureBar = 0; this.pressureEMA = 0;
    this.tweezer = null;            // {i, x, y, k}
    this.nPairs = 0; this.pairCap = 0;
    this.checkpoints = []; this.ckEvery = 64; this.ckBudget = opts.historyBytes || 48e6; this.recording = true;
    this.needRebuild = true;
    this.clamped = 0; this.redone = 0;
    this._alloc(64);
    this._allocPairs(1024);
  }

  static get ELEMENTS() { return ELEMENTS; }
  static element(sym) { return BY_SYM[sym]; }
  static get UNITS() { return { KB, ACC, KEU, BAR, COUL }; }
  static pairParams(a, b) { return PAIR[BY_SYM[a].t * NT + BY_SYM[b].t]; }

  _alloc(cap) {
    const old = this.cap ? this : null;
    const f64 = n => new Float64Array(n);
    const arrays = {
      pos: f64(3 * cap), vel: f64(3 * cap), frc: f64(3 * cap), prev: f64(3 * cap), built: f64(3 * cap),
      mass: f64(cap), phi: f64(cap), Zs: f64(cap), openVal: f64(cap), Zw: f64(cap), Gb: f64(cap), kSh: f64(cap), kShP: f64(cap), qg: f64(cap), Gz: f64(cap), formal: f64(cap), q: f64(cap), Z: f64(cap), cos0: f64(cap), G: f64(cap), ljx: f64(cap), ljq: f64(cap), lje: f64(cap),
      type: new Int16Array(cap), val: new Int8Array(cap), lp: new Int8Array(cap), pinned: new Uint8Array(cap),
      ids: new Uint32Array(cap), cStart: new Int32Array(cap + 1), cCount: new Int32Array(cap)
    };
    for (const k in arrays) { if (old) arrays[k].set(this[k].subarray(0, Math.min(this[k].length, arrays[k].length))); this[k] = arrays[k]; }
    this.cap = cap;
  }
  _allocPairs(cap) {
    const f64 = n => new Float64Array(n);
    const keep = this.pairCap ? { pI: this.pI, pJ: this.pJ, pN: this.pN, pF: this.pF, pSraw: this.pSraw, pB: this.pB } : null;
    this.pI = new Int32Array(cap); this.pJ = new Int32Array(cap); this.pN = f64(cap);
    this.pR = f64(cap); this.pDx = f64(cap); this.pDy = f64(cap); this.pDz = f64(cap);
    this.pF = f64(cap); this.pFp = f64(cap); this.pS = f64(cap); this.pSp = f64(cap); this.pScr = f64(cap); this.pSraw = f64(cap); this.pSpRaw = f64(cap); this.pSig = f64(cap); this.pD = f64(cap); this.gAb = f64(cap); this.gBb = f64(cap);
    if (!this.tri) { this.tri = new Int32Array(4096); this.nTri = 0; } this.pB = f64(cap); this.gA = f64(cap); this.gB = f64(cap); this.gAs = f64(cap); this.gBs = f64(cap); this.gAbs = f64(cap); this.gBbs = f64(cap);
    this.cList = new Int32Array(2 * cap);
    if (keep) for (const k of ['pI', 'pJ', 'pN', 'pF', 'pSraw', 'pB']) this[k].set(keep[k].subarray(0, this.nPairs));
    this.pairCap = cap;
  }

  /* ---------- random numbers (deterministic, part of the saved state) ---------- */
  rand() { const r = mulberry(this.rngState); this.rngState = r[1]; return r[0]; }
  gauss() { let u = 0; while (u === 0) u = this.rand(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * this.rand()); }
  gamma(shape) { // Marsaglia–Tsang
    if (shape < 1) return this.gamma(shape + 1) * Math.pow(this.rand() || 1e-12, 1 / shape);
    const d = shape - 1 / 3, c = 1 / Math.sqrt(9 * d);
    for (;;) {
      let x, v;
      do { x = this.gauss(); v = 1 + c * x; } while (v <= 0);
      v = v * v * v; const u = this.rand();
      if (u < 1 - 0.0331 * x * x * x * x) return d * v;
      if (Math.log(u) < 0.5 * x * x + d * (1 - v + Math.log(v))) return d * v;
    }
  }

  /* ---------- atoms ---------- */
  addAtom(sym, x, y, z, o = {}) {
    const el = BY_SYM[sym]; if (!el) throw new Error('Unsupported element ' + sym);
    if (this.N >= this.cap) this._alloc(this.cap * 2);
    const i = this.N++, i3 = 3 * i;
    this.type[i] = el.t; this.mass[i] = el.mass;
    const charge = o.charge || 0;
    this.formal[i] = charge;
    const V = o.V != null ? o.V : valenceFor(el, charge, o.valence);
    this.val[i] = V; this.lp[i] = lonePairs(el, charge, V);
    this.pos[i3] = x; this.pos[i3 + 1] = y; this.pos[i3 + 2] = z || 0;
    this.prev[i3] = x; this.prev[i3 + 1] = y; this.prev[i3 + 2] = z || 0;
    if (o.v) { this.vel[i3] = o.v[0]; this.vel[i3 + 1] = o.v[1]; this.vel[i3 + 2] = o.v[2]; }
    else if (o.thermal !== false && this.T > 0) {
      const s = Math.sqrt(KB * (o.T ?? this.T) / (el.mass * KEU));
      this.vel[i3] = s * this.gauss(); this.vel[i3 + 1] = s * this.gauss(); this.vel[i3 + 2] = s * this.gauss();
    } else { this.vel[i3] = this.vel[i3 + 1] = this.vel[i3 + 2] = 0; }
    this.frc[i3] = this.frc[i3 + 1] = this.frc[i3 + 2] = 0;
    this.cos0[i] = -1 / 3; this.pinned[i] = 0; this.ids[i] = this.nextId++;
    this.ljx[i] = el.ljX; this.lje[i] = Math.sqrt(el.ljD);
    this.needRebuild = true; this.needForces = true;
    return i;
  }
  /* Remove atoms (array of indices); remaining atoms keep their ids. */
  removeAtoms(list) {
    if (!list.length) return;
    const kill = new Uint8Array(this.N); for (const i of list) if (i >= 0 && i < this.N) kill[i] = 1;
    const map = new Int32Array(this.N);
    let w = 0;
    const per1 = ['mass', 'formal', 'q', 'Z', 'cos0', 'type', 'val', 'lp', 'pinned', 'ids', 'ljx', 'lje'];
    const per3 = ['pos', 'vel', 'frc', 'prev', 'built'];
    for (let i = 0; i < this.N; i++) {
      if (kill[i]) { map[i] = -1; continue; }
      map[i] = w;
      if (w !== i) {
        for (const k of per1) this[k][w] = this[k][i];
        for (const k of per3) { this[k][3 * w] = this[k][3 * i]; this[k][3 * w + 1] = this[k][3 * i + 1]; this[k][3 * w + 2] = this[k][3 * i + 2]; }
      }
      w++;
    }
    // keep bond orders of surviving pairs
    let np = 0;
    for (let p = 0; p < this.nPairs; p++) {
      const a = map[this.pI[p]], b = map[this.pJ[p]];
      if (a < 0 || b < 0) continue;
      this.pI[np] = a; this.pJ[np] = b;
      this.pN[np] = this.pN[p]; this.pF[np] = this.pF[p]; this.pSraw[np] = this.pSraw[p]; this.pB[np] = this.pB[p];
      np++;
    }
    this.nPairs = np;
    this.N = w;
    if (this.tweezer) { const t = map[this.tweezer.i]; if (t == null || t < 0) this.tweezer = null; else this.tweezer.i = t; }
    this.needRebuild = true; this.needForces = true;
    return map;
  }
  /* Seed the continuous bond order of an existing pair (e.g. a drawn double bond). */
  setBondOrder(i, j, n) {
    this._checkRebuild();
    const a = Math.min(i, j), b = Math.max(i, j), p = this.pairMap.get(a * 1048576 + b);
    if (p !== undefined) this.pN[p] = Math.max(1, Math.min(3, n));
  }
  /* Smallest allowed valence that accommodates `bondSum` for an element with a formal charge. */
  static valenceForBonds(sym, charge, bondSum) {
    const el = BY_SYM[sym]; if (!el) return 0;
    for (const v of el.valences) { const V = valenceFor(el, charge, v); if (V >= bondSum - 1e-9) return V; }
    return Math.max(bondSum, valenceFor(el, charge, el.valences[el.valences.length - 1]));
  }
  clear() { this.N = 0; this.nPairs = 0; this.tweezer = null; this.needRebuild = true; this.checkpoints.length = 0; }
  indexOfId(id) { for (let i = 0; i < this.N; i++) if (this.ids[i] === id) return i; return -1; }

  /* ---------- neighbour list (cell grid, Verlet skin) ---------- */
  _rebuild() {
    const N = this.N, pos = this.pos, rl = this.rc + this.skin, rl2 = rl * rl;
    /* Remember each surviving pair by atom-id, and carry across everything anyone may read
       before the next force pass: the bond order and the three numbers bondStrength is made of.
       Leaving those behind let a rebuild — or an erase — hand the inventory, the feed and the
       renderer one pair's strength under another pair's name. */
    const keep = this._keepMap || (this._keepMap = new Map()); keep.clear();
    const kept = this._kept && this._kept.length >= 4 * this.nPairs ? this._kept : (this._kept = new Float64Array(4 * this.pairCap + 64));
    for (let p = 0; p < this.nPairs; p++) {
      keep.set(this._key(this.pI[p], this.pJ[p]), 4 * p);
      kept[4 * p] = this.pN[p]; kept[4 * p + 1] = this.pF[p]; kept[4 * p + 2] = this.pSraw[p]; kept[4 * p + 3] = this.pB[p];
    }
    let minx = Infinity, miny = Infinity, minz = Infinity, maxx = -Infinity, maxy = -Infinity, maxz = -Infinity;
    for (let i = 0; i < N; i++) {
      const x = pos[3 * i], y = pos[3 * i + 1], z = pos[3 * i + 2];
      if (x < minx) minx = x; if (x > maxx) maxx = x; if (y < miny) miny = y; if (y > maxy) maxy = y; if (z < minz) minz = z; if (z > maxz) maxz = z;
    }
    const nx = Math.max(1, Math.min(64, Math.floor((maxx - minx) / rl) + 1));
    const ny = Math.max(1, Math.min(64, Math.floor((maxy - miny) / rl) + 1));
    const nz = Math.max(1, Math.min(16, Math.floor((maxz - minz) / rl) + 1));
    const cx = Math.max(rl, (maxx - minx) / nx + 1e-9), cy = Math.max(rl, (maxy - miny) / ny + 1e-9), cz = Math.max(rl, (maxz - minz) / nz + 1e-9);
    const nc = nx * ny * nz;
    let head = this._clHead, next = this._clNext, cell = this._clCell;
    if (!head || head.length < nc) head = this._clHead = new Int32Array(Math.max(nc, 1024));
    if (!next || next.length < N) { next = this._clNext = new Int32Array(N + 128); cell = this._clCell = new Int32Array(N + 128); }
    head.fill(-1, 0, nc);
    for (let i = N - 1; i >= 0; i--) {
      const ix = Math.min(nx - 1, Math.floor((pos[3 * i] - minx) / cx)), iy = Math.min(ny - 1, Math.floor((pos[3 * i + 1] - miny) / cy)), iz = Math.min(nz - 1, Math.floor((pos[3 * i + 2] - minz) / cz));
      const c = (iz * ny + iy) * nx + ix; cell[i] = c; next[i] = head[c]; head[c] = i;
    }
    let np = 0;
    for (let i = 0; i < N; i++) {
      const c = cell[i], ix = c % nx, iy = Math.floor(c / nx) % ny, iz = Math.floor(c / (nx * ny));
      const xi = pos[3 * i], yi = pos[3 * i + 1], zi = pos[3 * i + 2];
      for (let dz = -1; dz <= 1; dz++) {
        const jz = iz + dz; if (jz < 0 || jz >= nz) continue;
        for (let dy = -1; dy <= 1; dy++) {
          const jy = iy + dy; if (jy < 0 || jy >= ny) continue;
          for (let dx = -1; dx <= 1; dx++) {
            const jx = ix + dx; if (jx < 0 || jx >= nx) continue;
            for (let j = head[(jz * ny + jy) * nx + jx]; j >= 0; j = next[j]) {
              if (j <= i) continue;
              const ddx = pos[3 * j] - xi, ddy = pos[3 * j + 1] - yi, ddz = pos[3 * j + 2] - zi;
              if (ddx * ddx + ddy * ddy + ddz * ddz > rl2) continue;
              if (np >= this.pairCap) { this.nPairs = np; this._allocPairs(this.pairCap * 2); }
              this.pI[np] = i; this.pJ[np] = j;
              const k = keep.get(this._key(i, j));
              if (k === undefined) { this.pN[np] = 1; this.pF[np] = this.pSraw[np] = this.pB[np] = 0; }
              else { this.pN[np] = kept[k]; this.pF[np] = kept[k + 1]; this.pSraw[np] = kept[k + 2]; this.pB[np] = kept[k + 3]; }
              np++;
            }
          }
        }
      }
    }
    this.nPairs = np;
    const pm = this.pairMap = new Map();
    for (let p = 0; p < np; p++) pm.set(this.pI[p] * 1048576 + this.pJ[p], p);
    this.built.set(pos.subarray(0, 3 * N));
    this.needRebuild = false;
  }
  _key(i, j) { const a = this.ids[i], b = this.ids[j]; return a < b ? a * 67108864 + b : b * 67108864 + a; }
  _checkRebuild() {
    if (this.needRebuild) return this._rebuild();
    const lim = 0.25 * this.skin * this.skin, pos = this.pos, bl = this.built;
    for (let k = 0; k < 3 * this.N; k += 3) {
      const dx = pos[k] - bl[k], dy = pos[k + 1] - bl[k + 1], dz = pos[k + 2] - bl[k + 2];
      if (dx * dx + dy * dy + dz * dz > lim) return this._rebuild();
    }
  }

  /* ---------- forces ---------- */
  computeForces(relaxDt = 0) {
    const N = this.N, P = this.nPairs, pos = this.pos, F = this.frc, Z = this.Z, q = this.q, type = this.type, val = this.val;
    const pI = this.pI, pJ = this.pJ, pR = this.pR, pDx = this.pDx, pDy = this.pDy, pDz = this.pDz, pF = this.pF, pFp = this.pFp;
    const rc = this.rc, rc2 = rc * rc;
    F.fill(0, 0, 3 * N);
    const phi = this.phi, Zs = this.Zs, pS = this.pS, pSp = this.pSp, qg = this.qg;
    const Zw = this.Zw, Gb = this.Gb, pD = this.pD, mu = TUNE.mu;
    let Zf = this._Zf, Gf = this._Gf;
    if (!Zf || Zf.length < N) { Zf = this._Zf = new Float64Array(this.cap + 64); Gf = this._Gf = new Float64Array(this.cap + 64); }
    Zf.fill(0, 0, N); Gf.fill(0, 0, N);
    for (let i = 0; i < N; i++) { Z[i] = 0; Zs[i] = 0; Zw[i] = 0; Gb[i] = 0; q[i] = 0; this.G[i] = 0; this.cCount[i] = 0; phi[i] = 0; }
    let E = 0;

    // Pass A — geometry, coordination numbers, bond-polarisation charges
    // Clearing a whole array at once is a memset; clearing nine of them one element at a time,
    // inside a loop that is mostly skipped anyway, is nine scattered writes per pair.
    pF.fill(0, 0, P); pFp.fill(0, 0, P); pS.fill(0, 0, P); pSp.fill(0, 0, P); this.pSpRaw.fill(0, 0, P);
    for (let p = 0; p < P; p++) {
      const i = pI[p], j = pJ[p];
      const dx = pos[3 * j] - pos[3 * i], dy = pos[3 * j + 1] - pos[3 * i + 1], dz = pos[3 * j + 2] - pos[3 * i + 2];
      const r2 = dx * dx + dy * dy + dz * dz;
      if (r2 > rc2) { pR[p] = -1; continue; }
      const r = Math.sqrt(r2) || 1e-6;
      pR[p] = r; pDx[p] = dx; pDy[p] = dy; pDz[p] = dz;
      const pp = PAIR[type[i] * NT + type[j]];
      if (pp.bond && r < pp.s2) {
        satF(r, pp);
        if (SF > 0) { pS[p] = SF; pSp[p] = SFD; this.pSpRaw[p] = SFD; }
      }
      if (pp.bond && r < pp.r2) {
        smoothSwitch(r, pp.r1, pp.r2);
        if (SW > 0) {
          pF[p] = SW; pFp[p] = SWD;
          Z[i] += SW; Z[j] += SW;
          this.cCount[i]++; this.cCount[j]++;
          const raw = TUNE.kappa * SW * (CHI[type[j]] - CHI[type[i]]);
          const dq = raw / Math.pow(1 + (raw / Q_MAX) ** 4, 0.25);
          q[i] += dq; q[j] -= dq;
        }
      }
    }
    // Saturate each bond's transfer, then add equal and opposite charges.
    // Saturating net atomic charges separately used to create spurious total charge.
    for (let i = 0; i < N; i++) q[i] += this.formal[i];
    // coordination adjacency (CSR of pair indices)
    const cStart = this.cStart, cList = this.cList;
    cStart[0] = 0;
    for (let i = 0; i < N; i++) cStart[i + 1] = cStart[i] + this.cCount[i];
    for (let i = 0; i < N; i++) this.cCount[i] = 0;
    for (let p = 0; p < P; p++) if (pF[p] > 0) {
      const i = pI[p], j = pJ[p];
      cList[cStart[i] + this.cCount[i]++] = p; cList[cStart[j] + this.cCount[j]++] = p;
    }
    // Screening: two atoms bonded to a common neighbour i (a 1-3 pair, e.g. the H···H of water) do not
    // compete for each other's valence: S_jk = Π_i (1 − f_ij·f_ik). Its gradient is applied in
    // _screenForces once the saturation coefficients G are known.
    const scr = this.pScr, pSraw = this.pSraw;
    scr.fill(1, 0, P); pSraw.set(pS.subarray(0, P));
    let nt = 0;
    for (let i = 0; i < N; i++) {
      for (let a = cStart[i]; a < cStart[i + 1]; a++) {
        const pa = cList[a], j = pI[pa] === i ? pJ[pa] : pI[pa];
        for (let b = a + 1; b < cStart[i + 1]; b++) {
          const pb = cList[b], k = pI[pb] === i ? pJ[pb] : pI[pb];
          const q2 = this.pairMap.get(j < k ? j * 1048576 + k : k * 1048576 + j);
          if (q2 === undefined || pS[q2] === 0) continue;
          scr[q2] *= 1 - pF[pa] * pF[pb] * (1 - pF[q2]); // a pair that is itself bonded (3-ring) is not screened
          if (nt * 4 >= this.tri.length) { const t = new Int32Array(this.tri.length * 2); t.set(this.tri); this.tri = t; }
          this.tri[4 * nt] = q2; this.tri[4 * nt + 1] = pa; this.tri[4 * nt + 2] = pb; this.tri[4 * nt + 3] = i; nt++;
        }
      }
    }
    this.nTri = nt;
    this._insertionScreen();
    // Only integration / minimisation advances this internal model variable.
    // A force query or display refresh must not change the physical state.
    if (relaxDt > 0) this._updateBondOrders(relaxDt);
    for (let p = 0; p < P; p++) {
      if (pS[p] === 0) continue;
      // valence used = screened bond order × multiplicity; the π part of O=O counts only partly,
      // because ground-state O₂ is a triplet diradical that radicals add to without a barrier
      const n = this.pN[p], S = scr[p] * (1 + (PAIR[type[pI[p]] * NT + type[pJ[p]]].oo ? TUNE.wpiOO : TUNE.wpi) * (n - 1));
      pS[p] *= S; pSp[p] *= S;
      Zs[pI[p]] += pS[p]; Zs[pJ[p]] += pS[p];
      const zf = pSraw[p] * scr[p] * n; Zf[pI[p]] += zf; Zf[pJ[p]] += zf;
      const Dp = pD[p] = bondWeight(PAIR[type[pI[p]] * NT + type[pJ[p]]], this.pN[p], mu);
      Zw[pI[p]] += pS[p] * Dp; Zw[pJ[p]] += pS[p] * Dp;
    }
    /* Valence sharing. A partner's claim on an atom is only meaningful relative to what the atom
       has to give: an atom cannot contribute a whole bond to two neighbours at once. Where the
       raw claims exceed the valence they are scaled back in proportion, k = 1 − share·(1 − V/Z),
       so the competing coordination each bond sees is what its rivals can actually hold rather
       than what they would claim alone. k is exactly 1 for every atom at or under its valence,
       which is every equilibrium structure — so bond lengths, bond energies and thermochemistry
       are untouched by construction, and only the half-made crossing of a reaction changes. */
    const kSh = this.kSh, kShP = this.kShP, share = TUNE.share;
    for (let i = 0; i < N; i++) {
      const Z = Zs[i], V = val[i];
      if (share === 0 || V <= 0 || Z <= V) { kSh[i] = 1; kShP[i] = 0; continue; }
      kSh[i] = 1 - share + share * V / Z;
      kShP[i] = -share * V / (Z * Z);
    }
    if (TUNE.vstate) for (let i = 0; i < N; i++) {
      const tab = VS_TAB[type[i]];
      if (!tab || val[i] !== tab.length - 1) continue;      // a charged atom has another ladder
      valenceState(tab, Zf[i]);
      if (VS === 0 && VSD === 0) continue;
      E += TUNE.vstate * VS; Gf[i] = TUNE.vstate * VSD;
    }
    if (TUNE.gel && this.T > 0) {
      const kT = KB * this.T;
      let tabs = this._gelTabs;
      if (!tabs || this._gelT !== this.T) { tabs = this._gelTabs = GEL_TAB.map(f => f && f(this.T)); this._gelT = this.T; }
      for (let i = 0; i < N; i++) {
        const tab = tabs[type[i]];
        if (!tab || val[i] !== tab.length - 1) continue;
        valenceState(tab, Zf[i]);
        let L = VS, dL = VSD;
        if (Zf[i] <= 0) { L = tab[0]; dL = 0; }                 // the free atom itself
        if (L === 0 && dL === 0) continue;
        E -= TUNE.gel * kT * L; Gf[i] -= TUNE.gel * kT * dL;
      }
      /* Ground-state O2 is itself a triplet (³Σg⁻, g = 3), which no other common diatomic is: an
         O=O double bond between two oxygens with no other partner carries ln 3, fading as either
         oxygen takes on more (ozone, peroxides) or the bond comes apart. */
      const tO = BY_SYM.O.t;
      for (let p = 0; p < P; p++) {
        const i = pI[p], j = pJ[p];
        if (type[i] !== tO || type[j] !== tO || !(this.pN[p] > 1) || pSraw[p] <= 0) continue;
        o2Bump(Zf[i]); const bi = OB, dbi = OBD; if (bi === 0) continue;
        o2Bump(Zf[j]); const bj = OB, dbj = OBD; if (bj === 0) continue;
        const c = TUNE.gel * kT * Math.log(3) * Math.min(1, this.pN[p] - 1);
        E -= c * bi * bj; Gf[i] -= c * dbi * bj; Gf[j] -= c * bi * dbj;
      }
    }
    // Free (unpaired) valence per atom, smoothed so forces stay continuous. A half-filled orbital
    // feels much less Pauli repulsion than a closed shell, which is why radicals add without a barrier.
    const open = this.openVal;
    for (let i = 0; i < N; i++) {
      const x = val[i] - Zs[i];
      open[i] = x <= 0 ? 0 : x < OPEN_W ? x * x / (2 * OPEN_W) : Math.min(1, x - OPEN_W / 2);
    }
    // Partly-formed valence is stabilised (the delocalisation a half-made bond enjoys, ReaxFF's
    // under-coordination term). It vanishes at whole valence, so reactants, products and every fitted
    // bond energy are untouched; it only lowers the cost of being half-way through a reaction.
    if (TUNE.tsStab > 0) {
      const lam = TUNE.tsStab;
      for (let i = 0; i < N; i++) {
        const spare = val[i] - Zs[i];
        if (spare <= 0 || spare >= 1) continue;
        E -= lam * 4 * spare * (1 - spare);
        this.G[i] += 4 * lam * (1 - 2 * spare);
      }
    }
    // polar hydrogens get a small Pauli radius (hydrogen bonding)
    const ljx = this.ljx;
    for (let i = 0; i < N; i++) {
      const el = ELEMENTS[type[i]];
      this.ljq[i] = 0;
      if (el.Z === 1) {
        const t = (q[i] - this.formal[i] - 0.15) / 0.15;
        ljx[i] = el.ljX + (1.7 - el.ljX) * Math.min(1, Math.max(0, t));
        if (t > 0 && t < 1) this.ljq[i] = (1.7 - el.ljX) / 0.15;
      }
    }

    // Pass C — pair energies
    const pB = this.pB, gA = this.gA, gB = this.gB, gAs = this.gAs, gBs = this.gBs, gAb = this.gAb, gBb = this.gBb, pN = this.pN, G = this.G, lje = this.lje;
    const gc = 1 / Math.sqrt(rc * rc + COUL_D2), dgc = -rc * gc * gc * gc;
    gA.fill(0, 0, P); gB.fill(0, 0, P); gAs.fill(0, 0, P); gBs.fill(0, 0, P);
    this.gAbs.fill(0, 0, P); this.gBbs.fill(0, 0, P); gAb.fill(0, 0, P); gBb.fill(0, 0, P); pB.fill(0, 0, P);
    for (let p = 0; p < P; p++) {
      const r = pR[p];
      if (r < 0) continue;
      const i = pI[p], j = pJ[p], ti = type[i], tj = type[j];
      const pp = PAIR[ti * NT + tj], f = pF[p], fp = pFp[p], fs = pS[p];
      let e = 0, dEdr = 0, lam = 0, pi = 0, pj = 0, dpi = 0, dpj = 0, b = 0, ai = 1, bi = 0, aj = 1, bj = 0, asi = 1, bsi = 0, asj = 1, bsj = 0;
      // Lennard-Jones (shielded) first: its Pauli part decides whether b matters for a distant pair
      lj(r, ljx[i] * ljx[j], lje[i] * lje[j], rc);
      if (pp.bond) {
        // Morse with bond-order-interpolated parameters
        const n = pN[p], lo = Math.min(2, Math.floor(n)), t = n - lo, hi = lo + 1;
        const tl = t === 0 ? 0 : 1 - Math.pow(1 - t, 1.6); // bond length contracts fastest at low fractional order (Pauling/resonance)
        const re = pp.re[lo] + (pp.re[Math.min(hi, 3)] - pp.re[lo]) * tl;
        const De = pp.De[lo] + (pp.De[Math.min(hi, 3)] - pp.De[lo]) * t;
        const a = pp.a[lo] + (pp.a[Math.min(hi, 3)] - pp.a[lo]) * t;
        const y = a * (r - re);
        if (y < Y_OFF || LR !== 0) { // b only matters inside Morse range or the Pauli wall
          // Excess valence x = (used − this pair + 1 − V), scaled when over-crowded by the strength of the
          // competing bonds relative to this one (weights D = (De/400)^μ): a stronger incoming bond
          // displaces a weaker one more easily (Evans–Polanyi); an atom at its normal valence is unaffected.
          const Dj = mu === 0 ? 1 : mu === 0.5 ? Math.sqrt(De / 400) : bondWeight(pp, n, mu);
          // C = (Z − this pair)·k(Z): ∂C/∂Z = k + (Z − f)·k′, while the pair's own explicit
          // appearance contributes −k. The two differ once k varies, so they are carried apart.
          /* C = c·k(Z) and W = w·k(Z) both carry k, and k varies with Z, so the Zs channel picks
             up ∂/∂Z through both: XA·(k + c·k′) + XB·w·k′. The Zw channel sees only ∂W/∂Zw = k.
             The pair's own explicit appearance contributes −k to each. */
          if (val[i] > 0) {
            const ki = kSh[i], kp = kShP[i], ci = Zs[i] - fs, wi = Zw[i] - fs * Dj;
            excess(ci * ki, wi * ki, val[i], Dj, mu); sat(XS); pi = SP; dpi = SPD;
            ai = XA * (ki + ci * kp) + XB * wi * kp; bi = XB * ki; asi = XA * ki; bsi = XB * ki;
          }
          if (val[j] > 0) {
            const kj = kSh[j], kp = kShP[j], cj = Zs[j] - fs, wj = Zw[j] - fs * Dj;
            excess(cj * kj, wj * kj, val[j], Dj, mu); sat(XS); pj = SP; dpj = SPD;
            aj = XA * (kj + cj * kp) + XB * wj * kp; bj = XB * kj; asj = XA * kj; bsj = XB * kj;
          }
          b = pi * pj;
        }
        if (y < Y_OFF) {
          smoothSwitch(y, Y_ON, Y_OFF);
          const ey = Math.exp(-y), VR = De * ey * ey, VA = 2 * De * ey;
          e += SW * (VR - b * VA);
          dEdr += SWD * a * (VR - b * VA) + SW * (-2 * a * VR + b * a * VA);
          lam -= SW * VA;
        }
      }
      pB[p] = b;
      // hard core: nuclei never overlap, E = A·(rc/r − 1)² inside rc = 0.6·r(single bond)
      if (pp.core > 0 && r < pp.core) { const u = pp.core / r - 1; e += CORE_A * u * u; dEdr += -2 * CORE_A * u * pp.core / (r * r); }
      // LJ Pauli wall acts at contact distances; inside the bonding window (s > 0) the Morse term
      // already carries the repulsion, so the wall fades there: weight (1 − b)·(1 − s)
      const sr = this.pSraw[p], srp = this.pSpRaw[p];
      // open-shell partners screen each other less: the wall fades with their free valence
      let openF = 1, dOpenI = 0, dOpenJ = 0;
      if (TUNE.pauliOpen > 0) {
        const o = open[i] + open[j];
        if (o > 0) {
          const t = o < 1 ? o * o * (3 - 2 * o) : 1, dt = o < 1 ? 6 * o * (1 - o) : 0;
          openF = 1 - TUNE.pauliOpen * t;
          const d = -TUNE.pauliOpen * dt;
          dOpenI = open[i] > 0 && open[i] < 1 ? d : 0; dOpenJ = open[j] > 0 && open[j] < 1 ? d : 0;
        }
      }
      const wP = (1 - b) * (1 - sr) * openF;
      e += wP * LR; dEdr += wP * LDR - (1 - b) * openF * LR * srp; lam -= (1 - sr) * openF * LR;
      // dE/dZs through the open-shell factor (dopen/dZs = −1 inside the ramp)
      if (dOpenI !== 0) G[i] -= dOpenI * (1 - b) * (1 - sr) * LR;
      if (dOpenJ !== 0) G[j] -= dOpenJ * (1 - b) * (1 - sr) * LR;
      // excluded (1-2) dispersion + electrostatics
      const qq = q[i] * q[j];
      let h = 0, dh = 0;
      if (qq !== 0 || q[i] !== 0 || q[j] !== 0) { const g = 1 / Math.sqrt(r * r + COUL_D2); h = g - gc - dgc * (r - rc); dh = -r * g * g * g - dgc; }
      const ex = LA + COUL * qq * h, dex = LDA + COUL * qq * dh;
      e += (1 - f) * ex; dEdr += (1 - f) * dex - fp * ex;
      const radiusGrad = wP * LXR + (1 - f) * LXA;
      phi[i] += radiusGrad * ljx[j] * this.ljq[i];
      phi[j] += radiusGrad * ljx[i] * this.ljq[j];
      if (h !== 0) { const c = (1 - f) * COUL * h; phi[i] += c * q[j]; phi[j] += c * q[i]; }
      if (lam !== 0 && pp.bond) {
        const ga = lam * pj * dpi, gb = lam * pi * dpj; // dE/dx on each side
        gA[p] = ga * ai; gAb[p] = ga * bi; gB[p] = gb * aj; gBb[p] = gb * bj;
        gAs[p] = ga * asi; gBs[p] = gb * asj;              // the pair's own explicit −k term
        this.gAbs[p] = ga * bsi; this.gBbs[p] = gb * bsj;
        G[i] += gA[p]; Gb[i] += gAb[p]; G[j] += gB[p]; Gb[j] += gBb[p];
      }
      E += e;
      if (dEdr !== 0) {
        const s = dEdr / r, fx = s * pDx[p], fy = s * pDy[p], fz = s * pDz[p];
        F[3 * i] += fx; F[3 * i + 1] += fy; F[3 * i + 2] += fz;
        F[3 * j] -= fx; F[3 * j + 1] -= fy; F[3 * j + 2] -= fz;
      }
    }
    // Pass E — angles (VSEPR) and 1-3 exclusions (the exclusions also add to G)
    E += this._angles();
    // Pass S — forces from the screening factor
    this._screenForces();
    this._insertionForces();
    // Pass D — many-body forces from the coordination dependence of b
    for (let p = 0; p < P; p++) {
      const fp = pSp[p]; if (fp === 0) continue;
      const i = pI[p], j = pJ[p];
      const D = pD[p];
      let coef = (G[i] + Gb[i] * D - gAs[p] - this.gAbs[p] * D) + (G[j] + Gb[j] * D - gBs[p] - this.gBbs[p] * D);
      if (Gf[i] !== 0 || Gf[j] !== 0) coef += (Gf[i] + Gf[j]) * this.pN[p] / (1 + (PAIR[type[i] * NT + type[j]].oo ? TUNE.wpiOO : TUNE.wpi) * (this.pN[p] - 1));
      if (coef === 0) continue;
      const s = coef * fp / pR[p], fx = s * pDx[p], fy = s * pDy[p], fz = s * pDz[p];
      F[3 * i] += fx; F[3 * i + 1] += fy; F[3 * i + 2] += fz;
      F[3 * j] -= fx; F[3 * j + 1] -= fy; F[3 * j + 2] -= fz;
    }
    // Pass E' — charges and VSEPR targets depend on the structural bond orders f:
    // dE/dr = κ·f'·Δχ·(φ_i·q'_i − φ_j·q'_j) + (dE/dZ_i + dE/dZ_j)·f'
    const Gz = this.Gz;
    for (let p = 0; p < P; p++) {
      const fp = pFp[p]; if (fp === 0) continue;
      const i = pI[p], j = pJ[p];
      const dchi = CHI[type[j]] - CHI[type[i]];
      const raw = TUNE.kappa * pF[p] * dchi;
      const transferGrad = Math.pow(1 + (raw / Q_MAX) ** 4, -1.25);
      const dEdr = TUNE.kappa * fp * dchi * transferGrad * (phi[i] - phi[j]) + (Gz[i] + Gz[j]) * fp;
      if (dEdr === 0) continue;
      const s = dEdr / pR[p], fx = s * pDx[p], fy = s * pDy[p], fz = s * pDz[p];
      F[3 * i] += fx; F[3 * i + 1] += fy; F[3 * i + 2] += fz;
      F[3 * j] -= fx; F[3 * j + 1] -= fy; F[3 * j + 2] -= fz;
    }
    // Pass F — container walls and tweezers
    E += this._walls();
    /* Dragging moves a molecule; it does not whip one atom until the molecule comes off it.
       A stiff spring on a single atom reaches the force cap the moment the pointer is a few
       ångström ahead, tears the bond, and flings the pieces at kilometres a second — which then
       counts as the sample's heat and makes a stat stop everything else. So the pointer drives
       the whole dragged cluster as one: a damped servo on its centre of mass, distributed
       mass-weighted so every atom takes the same acceleration and the cluster feels no internal
       stress at all, with a speed limit that keeps a drag from turning into a projectile. */
    /* The servo does real work on the sample, so the step's energy check has to know how much:
       without it the check reads every drag as a failure, which is why it used to be switched off
       while dragging — and switched off is exactly when a reaction needs it most. */
    this._servoActive = false;
    if (this.tweezer) {
      const tw = this.tweezer, i = tw.i;
      if (i < N) {
        this._markDriven();
        const m = this.drivenMask;
        let M = 0, cx = 0, cy = 0, cz = 0, vx = 0, vy = 0, vz = 0;
        for (let a2 = 0; a2 < N; a2++) {
          if (!m[a2]) continue;
          const w = this.mass[a2]; M += w;
          cx += w * pos[3 * a2]; cy += w * pos[3 * a2 + 1]; cz += w * pos[3 * a2 + 2];
          vx += w * this.vel[3 * a2]; vy += w * this.vel[3 * a2 + 1]; vz += w * this.vel[3 * a2 + 2];
        }
        if (M > 0) {
          cx /= M; cy /= M; cz /= M; vx /= M; vy /= M; vz /= M;
          const tau = tw.tau || 120, vmax = tw.vmax || 0.03;   // fs to close the gap, Å/fs ceiling
          let wx = (tw.x - cx) / tau, wy = (tw.y - cy) / tau, wz = -cz / tau * 0.2;
          const ws = Math.hypot(wx, wy, wz);
          if (ws > vmax) { const f = vmax / ws; wx *= f; wy *= f; wz *= f; }
          const g = M / (tw.resp || 60) * (ACC > 0 ? 1 / ACC : 1);   // force to reach that speed
          const fx = g * (wx - vx), fy = g * (wy - vy), fz = g * (wz - vz);
          let tf = this._twF;
          if (!tf || tf.length < 3 * N) tf = this._twF = new Float64Array(3 * (N + 64));
          tf.fill(0, 0, 3 * N);
          for (let a2 = 0; a2 < N; a2++) {
            if (!m[a2]) continue;
            const w = this.mass[a2] / M;
            const ax = fx * w, ay = fy * w, az = fz * w;
            F[3 * a2] += ax; F[3 * a2 + 1] += ay; F[3 * a2 + 2] += az;
            tf[3 * a2] = ax; tf[3 * a2 + 1] = ay; tf[3 * a2 + 2] = az;
          }
          this._servoActive = true;
        }
      }
    }
    this.Epot = E;
    return E;
  }
  _coulH(r) {
    const rc = this.rc, g = 1 / Math.sqrt(r * r + COUL_D2), gc = 1 / Math.sqrt(rc * rc + COUL_D2), dgc = -rc * gc * gc * gc;
    return g - gc - dgc * (r - rc);
  }
  _coul(r, qq) {
    if (qq === 0) return 0;
    const rc = this.rc, g = 1 / Math.sqrt(r * r + COUL_D2), gc = 1 / Math.sqrt(rc * rc + COUL_D2), dgc = -rc * gc * gc * gc;
    return COUL * qq * (g - gc - dgc * (r - rc));
  }
  _dcoul(r, qq) {
    if (qq === 0) return 0;
    const rc = this.rc, g = 1 / Math.sqrt(r * r + COUL_D2), gc = 1 / Math.sqrt(rc * rc + COUL_D2);
    return COUL * qq * (-r * g * g * g + rc * gc * gc * gc);
  }
  /* Multiple bonds follow the sigma bonding each atom already has, measured with the same
     long-ranged bond order the saturation uses (screened, so 1-3 neighbours do not count). A radical
     approaching a double bond therefore weakens the pi bond while its own bond forms — bond order is
     conserved along the path instead of the pi bond having to break first. */
  _updateBondOrders(dt) {
    const N = this.N, val = this.val, cStart = this.cStart, cList = this.cList, pI = this.pI, pJ = this.pJ, pN = this.pN, type = this.type;
    const sig = this.pSig, scr = this.pScr, praw = this.pSraw;
    const spare = this._spare && this._spare.length >= N ? this._spare : (this._spare = new Float64Array(this.cap));
    const D = this._D && this._D.length >= N ? this._D : (this._D = new Float64Array(this.cap));
    const Zs = this._Zsig && this._Zsig.length >= N ? this._Zsig : (this._Zsig = new Float64Array(this.cap));
    for (let p = 0; p < this.nPairs; p++) {
      const raw = scr[p] * praw[p];
      // kbo < 1 starts the response while the partner is still approaching, so the pi bond gives way
      // as the new bond forms instead of afterwards
      sig[p] = raw <= 0 || raw >= 1 ? raw : Math.pow(raw, TUNE.kbo);
    }
    // Two passes: a partner consumes valence only to the extent that it can bond at all, so a radical
    // approaching a double bond frees the pi bond while a saturated molecule drifting past does not.
    // An already-formed bond (sigma ≈ 1) always counts in full.
    for (let pass = 0; pass < 2; pass++) {
      Zs.fill(0, 0, N);
      for (let p = 0; p < this.nPairs; p++) {
        const v = sig[p]; if (v <= 0) continue;
        const i = pI[p], j = pJ[p];
        const wi = pass === 0 ? 1 : Math.min(1, spare[j] + v), wj = pass === 0 ? 1 : Math.min(1, spare[i] + v);
        Zs[i] += v * wi; Zs[j] += v * wj;
      }
      for (let i = 0; i < N; i++) spare[i] = Math.max(0, val[i] - Zs[i]);
    }
    for (let i = 0; i < N; i++) {
      let s = 0;
      for (let c = cStart[i]; c < cStart[i + 1]; c++) { const p = cList[c], k = pI[p] === i ? pJ[p] : pI[p]; s += sig[p] * Math.min(2, spare[k]); }
      D[i] = s;
    }
    const rate = -Math.expm1(Math.log(0.92) * dt); // identical decay over any subdivision
    for (let p = 0; p < this.nPairs; p++) {
      const f = sig[p];
      let target = 1;
      if (f > 0) {
        const i = pI[p], j = pJ[p];
        const pp = PAIR[type[i] * NT + type[j]];
        if (pp.maxOrder > 1 && D[i] > 1e-9 && D[j] > 1e-9) {
          const si = spare[i] * f * Math.min(2, spare[j]) / D[i];
          const sj = spare[j] * f * Math.min(2, spare[i]) / D[j];
          target = 1 + Math.min(si, sj, pp.maxOrder - 1);
        }
        // three-electron O–O bond (HO₂·, RO₂·): one oxygen keeps an unpaired valence, the other does not
        if (pp.oo) target = Math.min(pp.maxOrder, target + TUNE.oo3e * f * Math.min(1, Math.abs(spare[i] - spare[j])));
      }
      if (pN[p] !== target) { const d = target - pN[p]; pN[p] = Math.abs(d) < 1e-4 ? target : pN[p] + d * rate; }
    }
  }
  _angles() {
    const N = this.N, pos = this.pos, F = this.frc, Z = this.Z, cStart = this.cStart, cList = this.cList, pI = this.pI, pJ = this.pJ, pF = this.pF, pFp = this.pFp, type = this.type, cos0 = this.cos0, Gz = this.Gz;
    let E = 0;
    for (let i = 0; i < N; i++) Gz[i] = 0;
    for (let i = 0; i < N; i++) {
      const c0 = cStart[i], c1 = cStart[i + 1], deg = c1 - c0;
      const el = ELEMENTS[type[i]];
      if (deg < 2 || el.kAngle === 0) continue;
      // VSEPR: ideal cos θ0 is a smooth function of the continuous steric number SN = Z + lone
      // pairs + whatever room the atom's unpaired electrons take (radicalDomains, below)
      const lp = this.lp[i];
      radicalDomains(this.val[i] - this.Zs[i], Z[i], lp);
      const sn = Z[i] + lp + RD_E, c4 = RD_C4;
      let c0v, dc0, dc0c4 = 0;
      if (sn <= 2) { c0v = -1; dc0 = 0; }
      else if (sn < 3) { const t = sn - 2; c0v = -1 + 0.5 * (3 * t * t - 2 * t * t * t); dc0 = 0.5 * (6 * t - 6 * t * t); }
      else if (sn < 4) { const t = sn - 3, st = 3 * t * t - 2 * t * t * t; c0v = -0.5 + (c4 + 0.5) * st; dc0 = (c4 + 0.5) * (6 * t - 6 * t * t); dc0c4 = st; }
      else { c0v = c4; dc0 = 0; dc0c4 = 1; }
      cos0[i] = c0v;
      // E = Kc·(c − c0)² + Kl·ℓ(c0)·(1 + c); ℓ switches on a harmonic pull toward 180° for sp centres
      const K = el.kAngle, Kc = 1.2 * K, Kl = 2 * K;
      let ell = 0, dell = 0;
      if (c0v < -0.5) { const u = (-c0v - 0.5) / 0.5; ell = u * u; dell = -4 * u; }
      let dEdc0sum = 0;
      const xi = pos[3 * i], yi = pos[3 * i + 1], zi = pos[3 * i + 2];
      for (let a = c0; a < c1; a++) {
        const pa = cList[a], j = pI[pa] === i ? pJ[pa] : pI[pa];
        const ux = pos[3 * j] - xi, uy = pos[3 * j + 1] - yi, uz = pos[3 * j + 2] - zi;
        const ru = Math.sqrt(ux * ux + uy * uy + uz * uz) || 1e-6;
        const fj = pF[pa], fpj = pFp[pa];
        for (let b = a + 1; b < c1; b++) {
          const pb = cList[b], k = pI[pb] === i ? pJ[pb] : pI[pb];
          const vx = pos[3 * k] - xi, vy = pos[3 * k + 1] - yi, vz = pos[3 * k + 2] - zi;
          const rv = Math.sqrt(vx * vx + vy * vy + vz * vz) || 1e-6;
          const fk = pF[pb], fpk = pFp[pb];
          const w = fj * fk;
          let c = (ux * vx + uy * vy + uz * vz) / (ru * rv);
          if (c > 1) c = 1; else if (c < -1) c = -1;
          const d = c - c0v, eu = Kc * d * d + Kl * ell * (1 + c), dEdc = 2 * Kc * d + Kl * ell;
          E += w * eu;
          dEdc0sum += w * (-2 * Kc * d + Kl * dell * (1 + c));
          const g = w * dEdc, iuv = 1 / (ru * rv);
          // dc/du, dc/dv
          let Fjx = -g * (vx * iuv - c * ux / (ru * ru)), Fjy = -g * (vy * iuv - c * uy / (ru * ru)), Fjz = -g * (vz * iuv - c * uz / (ru * ru));
          let Fkx = -g * (ux * iuv - c * vx / (rv * rv)), Fky = -g * (uy * iuv - c * vy / (rv * rv)), Fkz = -g * (uz * iuv - c * vz / (rv * rv));
          // weight derivatives
          const sj = -eu * fpj * fk / ru, sk = -eu * fj * fpk / rv;
          Fjx += sj * ux; Fjy += sj * uy; Fjz += sj * uz;
          Fkx += sk * vx; Fky += sk * vy; Fkz += sk * vz;
          F[3 * j] += Fjx; F[3 * j + 1] += Fjy; F[3 * j + 2] += Fjz;
          F[3 * k] += Fkx; F[3 * k + 1] += Fky; F[3 * k + 2] += Fkz;
          F[3 * i] -= Fjx + Fkx; F[3 * i + 1] -= Fjy + Fky; F[3 * i + 2] -= Fjz + Fkz;
          // 1-3 exclusion of the non-bonded pair (j,k)
          E += this._exclude13(j, k, w, fpj * fk, fj * fpk, ux, uy, uz, ru, vx, vy, vz, rv, i);
        }
      }
      Gz[i] = dEdc0sum * dc0 * (1 + RD_EZ);                  // dE/dZ_i through θ0
      // and through the free valence x = V − Zs, which sets how much room the unpaired electrons take
      const dEdx = dEdc0sum * (dc0 * RD_EX + dc0c4 * RD_C4X);
      if (dEdx !== 0) this.G[i] -= dEdx;
    }
    return E;
  }
  /* Insertion. A carbene - a carbon or silicon with two or more free valences, CH2, CH, C, SiH2 -
     can do what no radical can: take a sigma bond's own electron pair into its empty orbital, so it
     bonds to both ends at once while the bond between them lets go. That is why singlet methylene
     inserts into H2 with no barrier at all and makes methane, while an H atom, an OH or a CH3 has to
     pull one end off first and pays for it. In this model the pair shows up as the sigma bond's
     competition at its two ends: an H already bonded to its partner offers almost nothing to a
     newcomer. So when a carbene reaches both ends of one bond, that bond competes less at those
     ends, in proportion to how far the carbene has reached. Its reach is its own - first felt
     1.8 A past the bond length, fully engaged by 0.8 A past it, so 2.9 to 1.9 A for C-H - because
     the pull of an empty orbital on a bond's electrons starts before the ends' repulsion does. Built
     on the bond switch (1.6 A) or the saturation value (a tenth by 2 A) it arrived too late: the
     carbon was stopped at 2.2 A, 140 kJ/mol up the ends' repulsion, before it could engage.
       Only an atom with two or more spare valences once these two new bonds are set aside, and no
     multiple bond, qualifies — CH2, CH, a bare C — so a radical with one (H, OH, CH3) cannot, and a
     settled molecule - every atom at its valence - is untouched: competition only matters to an
     atom over its valence.
       And a hydrogen has one bond to give, not two. The crossing is three-centre: the H–H bond
     hands its order over to the two C–H bonds, the sum staying about one on each hydrogen. A screen
     that ignored that (full strength, whatever either bond was doing) let an H keep its H–H bond
     whole and bond to carbon besides - 1.4 bonds on one hydrogen - and CH2 + H2 settled as a
     CH2·H2 complex held together by 550 kJ/mol of double counting, never becoming methane. So an
     end lends only while its new bond and its old one's order sum to about one (conserveWindow),
     only while the new bonds are still forming (formingWindow), and the old bond is measured by
     its reach, not its switch (bondGate), so the screen never quits halfway through the stretch. */
  _insertionScreen() {
    this.nIns = 0;
    if (!TUNE.insert) return;
    const N = this.N, P = this.nPairs, type = this.type, val = this.val;
    const pI = this.pI, pJ = this.pJ, pF = this.pF, sr = this.pSraw, pN = this.pN, scr = this.pScr;
    let zr = this._insZr, zs = this._insZs, head = this._insHead;
    if (!zr || zr.length < N) { zr = this._insZr = new Float64Array(this.cap + 64); zs = this._insZs = new Float64Array(this.cap + 64); head = this._insHead = new Int32Array(this.cap + 64); }
    zr.fill(0, 0, N); zs.fill(0, 0, N); head.fill(-1, 0, N);
    for (let p = 0; p < P; p++) {
      const v = sr[p]; if (v <= 0) continue;
      const w = v * nEffOf(type[pI[p]], type[pJ[p]], pN[p]);
      zr[pI[p]] += w; zr[pJ[p]] += w; zs[pI[p]] += v; zs[pJ[p]] += v;
    }
    let cand = false;
    for (let i = 0; i < N; i++) if (INSERTS[type[i]] && val[i] - zr[i] > -1) { cand = true; break; }
    if (!cand) return;
    // pairs within saturation reach of each candidate centre, as linked lists
    let nxt = this._insNext, lst = this._insList;
    if (!nxt || nxt.length < 2 * P) { nxt = this._insNext = new Int32Array(2 * this.pairCap + 64); lst = this._insList = new Int32Array(2 * this.pairCap + 64); }
    let m = 0;
    const pR = this.pR;
    for (let p = 0; p < P; p++) {
      const r = pR[p]; if (!(r > 0)) continue;
      const pp = PAIR[type[pI[p]] * NT + type[pJ[p]]]; if (!pp.bond || r >= pp.re[1] + INS_OFF) continue;
      for (let e = 0; e < 2; e++) {
        const c = e ? pJ[p] : pI[p];
        if (!INSERTS[type[c]] || val[c] - zr[c] <= -1) continue;
        lst[m] = p; nxt[m] = head[c]; head[c] = m; m++;
      }
    }
    let rec = this._ins;
    for (let c = 0; c < N; c++) {
      if (head[c] < 0) continue;
      for (let x = head[c]; x >= 0; x = nxt[x]) {
        const pa = lst[x], a = pI[pa] === c ? pJ[pa] : pI[pa];
        for (let y = nxt[x]; y >= 0; y = nxt[y]) {
          const pb = lst[y], b = pI[pb] === c ? pJ[pb] : pI[pb];
          if (a === b) continue;
          const q = this.pairMap.get(a < b ? a * 1048576 + b : b * 1048576 + a);
          if (q === undefined || sr[q] <= INS_Q0) continue;      // the two ends must share a bond
          bondGate(sr[q]); const Q = QW, dQ = QWD;
          // Only the σ part of each arm is a new bond. Setting an arm's π part aside too made an
          // ethylene carbon — two H and its own C=C — look like it had two free valences, and it
          // was treated as a carbene inserting into its neighbour's C–H: pyramidal, C=C at 1.24 Å.
          const na = 1, nb = 1;
          const spare = val[c] - (zr[c] - sr[pa] - sr[pb]);
          carbeneWeight(spare); if (HW <= 0) continue;
          // and no multiple bond on the centre: CH2, CH and a bare C atom insert, C=O and C=C do not.
          // pi = the π part of the centre's other bonds, V − Zs counted with and without multiplicity
          const nfa = nEffOf(type[c], type[a], pN[pa]), nfb = nEffOf(type[c], type[b], pN[pb]);
          const pi = (zr[c] - zs[c]) - sr[pa] * (nfa - 1) - sr[pb] * (nfb - 1);
          piWindow(pi); if (GW <= 0) continue;
          const hw = HW * GW, hwd = HWD * GW, hgd = HW * GWD;
          // how far the carbene has reached each end, on its own reach rather than the bond's
          const ppa = PAIR[type[c] * NT + type[a]], ppb = PAIR[type[c] * NT + type[b]];
          smoothSwitch(pR[pa], ppa.re[1] + INS_ON, ppa.re[1] + INS_OFF); const Sa = SW, dSa = SWD;
          smoothSwitch(pR[pb], ppb.re[1] + INS_ON, ppb.re[1] + INS_OFF); const Sb = SW, dSb = SWD;
          const u = Sa * Sb, om = 1 - u, A = 1 - om * om * om;
          // and the bond's two ends have nothing to spare without it: H2's hydrogens, a C-H.
          // A ring bond seen from the third atom of a three-membered ring has ends with a free
          // valence each (the third atom is already their partner) - that is cyclopropane, not an
          // insertion, and treating it as one made cyclopropane 180 kJ/mol too stable.
          // Counted with the bond at full strength, so the window stays put while the bond stretches.
          // Counting its current strength switched the screen off halfway through the stretch, while
          // the bond still competed: the competition the carbene saw went back up, and the reaction
          // had to climb out of a pit the screen itself had dug.
          const nq = nEffOf(type[a], type[b], pN[q]), gap = nq * (1 - sr[q]);
          const sea = val[a] - (zr[a] - sr[pa]) - gap, seb = val[b] - (zr[b] - sr[pb]) - gap;
          endWindow(sea); const Ea = EW, dEa = EWD; if (Ea <= 0) continue;
          endWindow(seb); const Eb = EW, dEb = EWD; if (Eb <= 0) continue;
          // the bond's order, Pauling-style, starting to fall the moment it stretches: the saturation
          // reach sits flat at 1 until 0.82 Å of H–H, and the carbene had nothing to pull on until then
          const ppq = PAIR[type[a] * NT + type[b]], xq = pR[q] - (ppq.re[1] - 0.1);
          let rq = 0, drq = 0;
          if (xq > 0.2) { rq = xq - 0.1; drq = 1; } else if (xq > 0) { rq = xq * xq / 0.4; drq = xq / 0.2; }
          const mq = Math.exp(-rq / INS_QC), dmq = -mq * drq / INS_QC;
          formingWindow(sr[pa]); const Ua = UW, dUa = UWD; if (Ua <= 0) continue;
          formingWindow(sr[pb]); const Ub = UW, dUb = UWD; if (Ub <= 0) continue;
          conserveWindow(sr[pa] + mq); const Ca = CW * Ua, dCa = CWD * Ua, dUCa = CW * dUa; if (Ca <= 0) continue;
          conserveWindow(sr[pb] + mq); const Cb = CW * Ub, dCb = CWD * Ub, dUCb = CW * dUb; if (Cb <= 0) continue;
          const w = INS_MAX * hw * A * Q * Ea * Eb * Ca * Cb;
          if (w <= 0) continue;
          scr[q] *= 1 - w;
          if (!rec || rec.length < INS_REC * (this.nIns + 1)) { const t = new Float64Array(INS_REC * Math.max(8, 2 * (this.nIns + 1))); if (rec) t.set(rec); rec = this._ins = t; }
          const o = INS_REC * this.nIns++;
          rec[o] = q; rec[o + 1] = pa; rec[o + 2] = pb; rec[o + 3] = c;
          rec[o + 4] = w; rec[o + 5] = hw; rec[o + 6] = hwd; rec[o + 7] = A; rec[o + 8] = 3 * om * om;
          rec[o + 9] = na; rec[o + 10] = nb; rec[o + 11] = Sa; rec[o + 12] = dSa; rec[o + 13] = Sb; rec[o + 14] = dSb; rec[o + 15] = hgd;
          rec[o + 16] = a; rec[o + 17] = b; rec[o + 18] = Ea * Ca; rec[o + 19] = dEa * Ca; rec[o + 20] = Eb * Cb; rec[o + 21] = dEb * Cb;
          rec[o + 22] = nfa; rec[o + 23] = nfb;
          rec[o + 24] = Ca; rec[o + 25] = dCa; rec[o + 26] = Cb; rec[o + 27] = dCb; rec[o + 28] = Ea * Eb; rec[o + 29] = nq;
          rec[o + 30] = Q; rec[o + 31] = dQ; rec[o + 32] = dmq; rec[o + 33] = dUCa; rec[o + 34] = dUCb;
        }
      }
    }
  }
  _insertionForces() {
    if (!this.nIns) return;
    const rec = this._ins, sr = this.pSraw, spr = this.pSpRaw, pF = this.pF, pFp = this.pFp, scr = this.pScr;
    const pI = this.pI, pJ = this.pJ, pN = this.pN, type = this.type, G = this.G, Gb = this.Gb;
    const pR = this.pR, pDx = this.pDx, pDy = this.pDy, pDz = this.pDz, F = this.frc, N = this.N, P = this.nPairs;
    let Hc = this._insHc, Hg = this._insHg;
    if (!Hc || Hc.length < N) { Hc = this._insHc = new Float64Array(this.cap + 64); Hg = this._insHg = new Float64Array(this.cap + 64); }
    Hc.fill(0, 0, N); Hg.fill(0, 0, N);
    const push = (p, dEdr) => {
      if (dEdr === 0 || !(pR[p] > 0)) return;
      const s = dEdr / pR[p], i = pI[p], j = pJ[p], fx = s * pDx[p], fy = s * pDy[p], fz = s * pDz[p];
      F[3 * i] += fx; F[3 * i + 1] += fy; F[3 * i + 2] += fz;
      F[3 * j] -= fx; F[3 * j + 1] -= fy; F[3 * j + 2] -= fz;
    };
    let anyH = false;
    for (let t = 0; t < this.nIns; t++) {
      const o = INS_REC * t, q = rec[o], pa = rec[o + 1], pb = rec[o + 2], c = rec[o + 3];
      const ea = rec[o + 16], eb = rec[o + 17], Ea = rec[o + 18], dEa = rec[o + 19], Eb = rec[o + 20], dEb = rec[o + 21], EE = Ea * Eb;
      const w = rec[o + 4], h = rec[o + 5], hd = rec[o + 6], A = rec[o + 7], Ad = rec[o + 8], na = rec[o + 9], nb = rec[o + 10];
      const Sa = rec[o + 11], dSa = rec[o + 12], Sb = rec[o + 13], dSb = rec[o + 14], hg = rec[o + 15];
      const j = pI[q], k = pJ[q];
      const nEff = nEffOf(type[j], type[k], pN[q]);
      const D = this.pD[q];
      // dE/d(scr[q]), exactly as the ordinary screening forces take it
      const dEdS = sr[q] * nEff * ((G[j] + Gb[j] * D - this.gAs[q] - this.gAbs[q] * D) + (G[k] + Gb[k] * D - this.gBs[q] - this.gBbs[q] * D))
        + sr[q] * pN[q] * (this._Gf[j] + this._Gf[k]);
      if (dEdS === 0) continue;
      const dEdw = -dEdS * scr[q] / (1 - w);
      const fq = rec[o + 30], nq = rec[o + 29];
      // w = INS_MAX * h(spare, pi) * A(Sa*Sb) * Q(sr_q) * Ea(spare of a) * Eb(spare of b) * Ca * Cb
      push(pa, dEdw * INS_MAX * h * Ad * Sb * fq * EE * dSa);
      push(pb, dEdw * INS_MAX * h * Ad * Sa * fq * EE * dSb);
      push(q, dEdw * INS_MAX * h * A * EE * rec[o + 31] * spr[q]);
      // each end: Ca = C(sr_a + m_q)·U(sr_a), and likewise b
      const base0 = dEdw * INS_MAX * h * A * fq * rec[o + 28];
      const gTa = base0 * rec[o + 25] * rec[o + 26], gTb = base0 * rec[o + 24] * rec[o + 27];
      const gUa = base0 * rec[o + 33] * rec[o + 26], gUb = base0 * rec[o + 24] * rec[o + 34];
      if (gTa !== 0 || gTb !== 0 || gUa !== 0 || gUb !== 0) {
        push(pa, (gTa + gUa) * spr[pa]); push(pb, (gTb + gUb) * spr[pb]); push(q, (gTa + gTb) * rec[o + 32]);
      }
      // each end's spare = V - (its coordination - its arm): all its pairs lower it, its arm not
      for (let e2 = 0; e2 < 2; e2++) {
        const at = e2 ? eb : ea, arm = e2 ? pb : pa;
        const dEdse = dEdw * INS_MAX * h * A * fq * (e2 ? Ea * dEb : dEa * Eb);
        if (dEdse === 0) continue;
        Hc[at] -= dEdse; anyH = true;
        push(arm, dEdse * spr[arm]);
        push(q, dEdse * nq * spr[q]);               // the bond counted at full strength
      }
      // spare = V - (all of the centre's coordination - the two arms): every other pair of the
      // centre lowers it one-for-one, the two arms not at all
      const dEdspare = dEdw * INS_MAX * hd * A * fq * EE;
      if (dEdspare !== 0) {
        Hc[c] -= dEdspare; anyH = true;
        push(pa, dEdspare * na * spr[pa]);
        push(pb, dEdspare * nb * spr[pb]);
      }
      // pi = Σ over the centre's other pairs of S·(n − 1): up with the weighted count, down with the
      // plain one, and the two arms' own π parts set aside
      const dEdpi = dEdw * INS_MAX * hg * A * fq * EE;
      if (dEdpi !== 0) {
        Hc[c] += dEdpi; Hg[c] -= dEdpi; anyH = true;
        push(pa, -dEdpi * (rec[o + 22] - 1) * spr[pa]);
        push(pb, -dEdpi * (rec[o + 23] - 1) * spr[pb]);
      }
    }
    if (!anyH) return;
    for (let p = 0; p < P; p++) {
      if (spr[p] === 0) continue;
      const hsum = Hc[pI[p]] + Hc[pJ[p]], gsum = Hg[pI[p]] + Hg[pJ[p]];
      if (hsum === 0 && gsum === 0) continue;
      push(p, (hsum * nEffOf(type[pI[p]], type[pJ[p]], pN[p]) + gsum) * spr[p]);
    }
  }
  _screenForces() {
    const tri = this.tri, pS = this.pSraw, pN = this.pN, scr = this.pScr, G = this.G, gA = this.gA, gB = this.gB;
    const pI = this.pI, pJ = this.pJ, pF = this.pF, pFp = this.pFp, pR = this.pR, pDx = this.pDx, pDy = this.pDy, pDz = this.pDz, F = this.frc;
    for (let t = 0; t < this.nTri; t++) {
      const q = tri[4 * t], pa = tri[4 * t + 1], pb = tri[4 * t + 2];
      const j = pI[q], k = pJ[q];
      const nEff = 1 + (PAIR[this.type[j] * NT + this.type[k]].oo ? TUNE.wpiOO : TUNE.wpi) * (pN[q] - 1);
      const D = this.pD[q], Gb = this.Gb;
      const dEdS = pS[q] * nEff * ((G[j] + Gb[j] * D - this.gAs[q] - this.gAbs[q] * D) + (G[k] + Gb[k] * D - this.gBs[q] - this.gBbs[q] * D))
        + pS[q] * pN[q] * (this._Gf[j] + this._Gf[k]);
      if (dEdS === 0) continue;
      const m = 1 - pF[q], w = pF[pa] * pF[pb] * m, one = 1 - w;
      let others; // Π over the other shared neighbours
      if (one > 1e-9) others = scr[q] / one;
      else { others = 1; for (let u = 0; u < this.nTri; u++) if (u !== t && tri[4 * u] === q) others *= 1 - pF[tri[4 * u + 1]] * pF[tri[4 * u + 2]] * m; }
      const dEdw = -dEdS * others;
      // w = f(r_a)·f(r_b)·(1 − f(r_jk)): push along each arm of the triple and along j–k
      for (let arm = 0; arm < 3; arm++) {
        const p = arm === 0 ? pa : arm === 1 ? pb : q;
        const dw = arm === 0 ? pFp[pa] * pF[pb] * m : arm === 1 ? pF[pa] * pFp[pb] * m : -pF[pa] * pF[pb] * pFp[q];
        if (dw === 0) continue;
        const s = dEdw * dw / pR[p], a = pI[p], b = pJ[p];
        const fx = s * pDx[p], fy = s * pDy[p], fz = s * pDz[p];
        F[3 * a] += fx; F[3 * a + 1] += fy; F[3 * a + 2] += fz;
        F[3 * b] -= fx; F[3 * b + 1] -= fy; F[3 * b + 2] -= fz;
      }
    }
  }
  _exclude13(j, k, w, dwdru, dwdrv, ux, uy, uz, ru, vx, vy, vz, rv, i) {
    const pos = this.pos, F = this.frc, type = this.type, tj = type[j], tk = type[k];
    const dx = pos[3 * k] - pos[3 * j], dy = pos[3 * k + 1] - pos[3 * j + 1], dz = pos[3 * k + 2] - pos[3 * j + 2];
    const r = Math.sqrt(dx * dx + dy * dy + dz * dz) || 1e-6;
    if (r > this.rc) return 0;
    const pp = PAIR[tj * NT + tk];
    // b for (j,k) with current coordination (held fixed for the exclusion)
    let b = 0, fjk = 0, fpjk = 0, X = 0, dX = 0, dpj = 0, dpk = 0, bpj = 0, bpk = 0, sr = 0, srp = 0;
    if (pp.bond) {
      if (r < pp.s2) { satF(r, pp); sr = SF; srp = SFD; }
      smoothSwitch(r, pp.r1, pp.r2); fjk = SW; fpjk = SWD;
      const sjk = 0; // a 1-3 pair is screened by the shared neighbour
      let pj = 0, pk = 0;
      if (this.val[j] > 0) { sat(this.Zs[j] - sjk + 1 - this.val[j]); pj = SP; dpj = SPD; }
      if (this.val[k] > 0) { sat(this.Zs[k] - sjk + 1 - this.val[k]); pk = SP; dpk = SPD; }
      b = pj * pk; bpj = pj; bpk = pk;
      const re = pp.re[1], De = pp.De[1], a = pp.a[1], y = a * (r - re);
      if (y < Y_OFF) { smoothSwitch(y, Y_ON, Y_OFF); const ey = Math.exp(-y), VR = De * ey * ey; X += SW * VR; dX += SWD * a * VR - 2 * a * SW * VR; }
    }
    lj(r, this.ljx[j] * this.ljx[k], this.lje[j] * this.lje[k], this.rc);
    const qq = this.q[j] * this.q[k], h = this._coulH(r);
    const wP = (1 - b) * (1 - sr);
    X += wP * LR + LA + COUL * qq * h;
    dX += wP * LDR - (1 - b) * LR * srp + LDA + this._dcoul(r, qq);
    const m = 1 - fjk; // do not double-exclude a pair that is itself bonded (3-rings)
    const radiusGrad = -w * m * (wP * LXR + LXA);
    this.phi[j] += radiusGrad * this.ljx[k] * this.ljq[j];
    this.phi[k] += radiusGrad * this.ljx[j] * this.ljq[k];
    const Esub = -w * m * X;
    // Esub depends on b = p_j·p_k through the Pauli term: dEsub/db = w·m·LR → saturation coefficients
    if (LR !== 0 && (dpj !== 0 || dpk !== 0)) { const c = w * m * (1 - sr) * LR; this.G[j] += c * bpk * dpj; this.G[k] += c * bpj * dpk; }
    if (h !== 0) { const c = w * m * COUL * h; this.phi[j] -= c * this.q[k]; this.phi[k] -= c * this.q[j]; }
    // gradient wrt r_jk
    const dEdr = -w * (m * dX - fpjk * X);
    const s = dEdr / r;
    F[3 * j] += s * dx; F[3 * j + 1] += s * dy; F[3 * j + 2] += s * dz;
    F[3 * k] -= s * dx; F[3 * k + 1] -= s * dy; F[3 * k + 2] -= s * dz;
    // gradient wrt weights w = f_ij·f_ik
    const cu = m * X * dwdru / ru, cv = m * X * dwdrv / rv; // −dE/dr_ij = +m X dw/dr
    F[3 * j] += cu * ux; F[3 * j + 1] += cu * uy; F[3 * j + 2] += cu * uz;
    F[3 * k] += cv * vx; F[3 * k + 1] += cv * vy; F[3 * k + 2] += cv * vz;
    F[3 * i] -= cu * ux + cv * vx; F[3 * i + 1] -= cu * uy + cv * vy; F[3 * i + 2] -= cu * uz + cv * vz;
    return Esub;
  }
  /* Rectangular solid walls reflect the drift at the face. Soft fields and the
     spherical conditioning vessel retain conservative harmonic potentials. */
  /* Solid faces reflect what reaches them. Only a velocity void changes that: an atom that
     reaches the face stops there instead of bouncing. A pressure void changes nothing about how
     atoms move — it voids the pressure, not the motion — so the wall simply records nothing. */
  /* Which molecule each atom belongs to, cached for the step. Bonds do not come and go inside one
     step, so the sub-steps of a step can share one answer. */
  _fragRoots() {
    const N = this.N;
    if (this._fragStamp === this.stepCount && this._fragN === N && this._fragRoot && this._fragRoot.length >= N) return this._fragRoot;
    let root = this._fragRoot;
    if (!root || root.length < N) root = this._fragRoot = new Int32Array(N + 128);
    for (let i = 0; i < N; i++) root[i] = i;
    const find = i => { while (root[i] !== i) { root[i] = root[root[i]]; i = root[i]; } return i; };
    for (let p = 0; p < this.nPairs; p++) {
      if (this.bondStrength(p) <= 0.25) continue;
      const a = find(this.pI[p]), b = find(this.pJ[p]); if (a !== b) root[a] = b;
    }
    for (let i = 0; i < N; i++) root[i] = find(i);
    this._fragStamp = this.stepCount; this._fragN = N;
    return root;
  }
  /* The safety net behind a solid face. Ordinary contact is the face's own force, in _walls. What
     is left for this is an atom the face never had the chance to push — one put far outside by
     an edit, a smaller chamber or the barostat. That is carried back inside with the molecule it
     belongs to as one rigid piece, so no bond is stretched on the way, and its travel into the
     face, if it is still travelling that way, is reversed. The velocity void also lives here:
     an atom that touches a face while it is on stops where it stands. */
  _reflectWalls() {
    if (this.sphere || this.boundsMode !== 'solid') return;
    const N = this.N, pos = this.pos, vel = this.vel, b = this.box;
    const lo = this._wLo || (this._wLo = [0, 0, 0]), hi = this._wHi || (this._wHi = [0, 0, 0]);
    lo[0] = b.x0; lo[1] = b.y0; lo[2] = b.z0; hi[0] = b.x1; hi[1] = b.y1; hi[2] = b.z1;
    // the face's force handles ordinary contact; this only acts on contact when the velocity void
    // is on, and otherwise only on an atom left far outside by an edit, a resize or the barostat
    const margin = this.voidVelocity ? 0 : WALL_FOLD;
    let any = false;
    for (let i = 0; i < N && !any; i++) {
      if (this.pinned[i]) continue;
      const k = 3 * i;
      for (let d = 0; d < 3; d++) if (pos[k + d] < lo[d] - margin || pos[k + d] > hi[d] + margin) { any = true; break; }
    }
    if (!any) return;
    const report = !this.voidPressure;
    if (this.voidVelocity) {          // the velocity void: whatever touches the wall stops there
      for (let i = 0; i < N; i++) {
        if (this.pinned[i]) continue;
        const k = 3 * i;
        let touched = false;
        for (let d = 0; d < 3; d++) {
          const j = k + d;
          if (pos[j] < lo[d] || pos[j] > hi[d]) {
            touched = true;
            pos[j] = clampNum(pos[j], lo[d], hi[d]);
            if (report) this._wallImpulse += this.mass[i] * Math.abs(vel[j]) * KEU;  // absorbed: m*v
          }
        }
        if (touched) for (let d = 0; d < 3; d++) {
          const j = k + d;
          this.voidHeat += 0.5 * KEU * this.mass[i] * vel[j] * vel[j];
          vel[j] = 0;
        }
      }
      return;
    }
    const root = this._fragRoots();
    let acc = this._wAcc;
    if (!acc || acc.length < 7 * N) acc = this._wAcc = new Float64Array(7 * (N + 128));
    acc.fill(0, 0, 7 * N);            // per molecule: mass, momentum(3), deepest overshoot(3)
    for (let i = 0; i < N; i++) {
      if (this.pinned[i]) continue;
      const r = 7 * root[i], k = 3 * i, m = this.mass[i];
      acc[r] += m;
      for (let d = 0; d < 3; d++) {
        acc[r + 1 + d] += m * vel[k + d];
        const x = pos[k + d];
        const o = x < lo[d] - WALL_FOLD ? x - lo[d] : x > hi[d] + WALL_FOLD ? x - hi[d] : 0;
        // remember the shift that would carry the atom furthest out back where it belongs, folded
        // the same way however many chamber-lengths it has travelled
        if (Math.abs(o) > Math.abs(acc[r + 4 + d])) {
          const L = hi[d] - lo[d], u = (x - lo[d]) / L, folded = ((u % 2) + 2) % 2;
          acc[r + 4 + d] = (lo[d] + L * (folded <= 1 ? folded : 2 - folded)) - x;
        }
      }
    }
    // a molecule still travelling into the wall has that travel reversed; one that has already
    // turned around, or merely straddles a face, is only carried back inside
    for (let i = 0; i < N; i++) {
      const r = 7 * i;
      if (root[i] !== i || acc[r] <= 0) continue;
      const M = acc[r];
      for (let d = 0; d < 3; d++) {
        const V = acc[r + 1 + d] / M;
        const shift = acc[r + 4 + d];
        acc[r + 1 + d] = shift !== 0 && (shift > 0) !== (V > 0) && V !== 0 ? 2 * V : 0;
        if (report && acc[r + 1 + d] !== 0) this._wallImpulse += Math.abs(acc[r + 1 + d]) * M * KEU;
      }
    }
    for (let i = 0; i < N; i++) {
      if (this.pinned[i]) continue;
      const r = 7 * root[i], k = 3 * i;
      for (let d = 0; d < 3; d++) {
        pos[k + d] += acc[r + 4 + d];
        vel[k + d] -= acc[r + 1 + d];
      }
    }
  }
  _walls() {
    const N = this.N, pos = this.pos, F = this.frc;
    const field = this.boundsMode === 'forcefield' && !this.sphere;
    const K = field ? this.fieldK : this.wallK, reach = field ? this.fieldRange : 0;
    let E = 0, Fsum = 0;
    if (this.sphere) {
      const s = this.sphere;
      for (let i = 0; i < N; i++) {
        const dx = pos[3 * i] - s.x, dy = pos[3 * i + 1] - s.y, dz = pos[3 * i + 2] - s.z, r = Math.sqrt(dx * dx + dy * dy + dz * dz);
        if (r > s.R) { const d = r - s.R, f = K * d / r; E += 0.5 * K * d * d; Fsum += K * d; F[3 * i] -= f * dx; F[3 * i + 1] -= f * dy; F[3 * i + 2] -= f * dz; }
      }
      this.wallArea = 4 * Math.PI * s.R * s.R;
    } else {
      const b = this.box, lo = [b.x0, b.y0, b.z0], hi = [b.x1, b.y1, b.z1];
      if (field) for (let i = 0; i < N; i++) for (let d = 0; d < 3; d++) {
        const x = pos[3 * i + d];
        let e = 0, dir = 0;
        if (x < lo[d] + reach) { e = lo[d] + reach - x; dir = 1; }
        else if (x > hi[d] - reach) { e = x - (hi[d] - reach); dir = -1; }
        if (e <= 0) continue;
        E += 0.5 * K * e * e; F[3 * i + d] += dir * K * e;
        Fsum += K * e;
      }
      /* A solid face is a force too, just a stiff one that starts exactly at the face. It used to
         be a rule instead — anything past the face was put back and had its velocity turned
         round — and a rule applied to positions is energy nobody pays for. Mirroring one atom of
         a molecule stretched its bonds; carrying the whole molecule back pushed it into its
         neighbours; either way a crowded chamber against a solid wall heated itself, by 170 kJ/mol
         in 20 ps for fifteen waters. As a force it is integrated like every other force: it
         conserves energy, it pushes on the one atom that touches so a spinning molecule is turned
         by the wall rather than shoved off it, and the pressure is simply what it pushes with.
         At this stiffness a room-temperature atom goes about 0.04 Å past the face and comes back. */
      if (!field && this.boundsMode === 'solid') {
        // the force grows from zero with zero slope over the first SOLID_W, then linearly: a force
        // that switched on with a kink cost a fast atom a quarter of a percent of its energy per bounce
        const Ks = this.solidK, w = SOLID_W;
        for (let i = 0; i < N; i++) for (let d = 0; d < 3; d++) {
          const x = pos[3 * i + d];
          let e, dir;
          if (x < lo[d]) { e = lo[d] - x; dir = 1; } else if (x > hi[d]) { e = x - hi[d]; dir = -1; } else continue;
          let f;
          if (e < w) { f = Ks * e * e / (2 * w); E += Ks * e * e * e / (6 * w); }
          else { f = Ks * (e - w / 2); E += Ks * (e * e / 2 - w * e / 2 + w * w / 6); }
          F[3 * i + d] += dir * f;
          Fsum += f;
        }
      }
      const Lx = b.x1 - b.x0, Ly = b.y1 - b.y0, Lz = b.z1 - b.z0;
      this.wallArea = 2 * (Lx * Ly + Lx * Lz + Ly * Lz);
    }
    this.Ewall = E; this.wallForce = this.voidPressure ? 0 : Fsum; this.voidForce = 0; // legacy snapshot field; reflecting walls always transfer momentum
    return E;
  }

  /* Void walls. Each channel deletes one thing the chamber would otherwise keep. They are not
     gentle absorbers: what is voided is gone the moment it appears, which is what "void" means.

       velocity    — an atom that reaches the wall stops dead and stays stopped. Not damped, not
                     reflected: every component is zeroed for as long as it is in contact, so it
                     rests there until a neighbour pushes it off.
       pressure    — only the component along the face it touched is absorbed, so the atom can
                     still slide along the wall. The wall registers no impulse, so the pressure
                     it would have reported radiates away without the chamber changing size.
       temperature — handled separately, in _voidHeat: heat is radiated off the instant it
                     appears, everywhere, not only where an atom happens to touch a wall. */
  /* Velocity void. An atom that reaches a face stops there: every component zeroed, once, at the
     moment of contact, and then left alone — free to be moved by whatever else acts on it. Only
     motion driving into a face it has actually reached counts, so an atom gliding past a wall or
     still approaching one is untouched, and an atom held at the wall is never clamped against its
     own bonded neighbours. Pressure is not handled here: voiding a pressure reading must not
     change how anything moves. */
  _voidWalls() {
    if (this.sphere || !this.voidVelocity) return;
    const b = this.box, lo = [b.x0, b.y0, b.z0], hi = [b.x1, b.y1, b.z1];
    const N = this.N, pos = this.pos, vel = this.vel;
    let lost = 0;
    for (let i = 0; i < N; i++) {
      if (this.pinned[i]) continue;
      const k = 3 * i;
      let driving = false;
      for (let d = 0; d < 3 && !driving; d++) {
        const p = pos[k + d], v = vel[k + d];
        if ((p <= lo[d] && v < 0) || (p >= hi[d] && v > 0)) driving = true;
      }
      if (!driving) continue;
      const m = this.mass[i];
      for (let d = 0; d < 3; d++) {
        const j = k + d;
        lost += 0.5 * KEU * m * vel[j] * vel[j];
        vel[j] = 0;
      }
    }
    this.voidHeat += lost;
  }
  /* Void temperature. A sample that heats itself — friction, a reaction, work done on it —
     radiates that heat away as fast as it is made, so the chamber never runs hotter than the
     temperature it was set to. One-sided: it only ever removes, so a cold chamber stays cold
     and nothing here can drive the sample. */
  _voidHeat() {
    if (!this.voidTemperature) return;
    const Nf = this.dof(); if (!Nf) return;
    const K = this.thermalKinetic(), Kt = 0.5 * Nf * KB * Math.max(0, this.T);
    if (!(K > Kt)) return;
    const sc = Math.sqrt(Kt / K), v = this.vel;
    const vd = this._vdTmp || (this._vdTmp = [0, 0, 0]), dragging = this.dragVelocity(vd);
    for (let i = 0; i < this.N; i++) {
      if (this.pinned[i]) { v.fill(0, 3 * i, 3 * i + 3); continue; }
      const k = 3 * i;
      if (dragging && this.driven(i)) {
        for (let d = 0; d < 3; d++) v[k + d] = vd[d] + (v[k + d] - vd[d]) * sc;
      } else for (let d = 0; d < 3; d++) v[k + d] *= sc;
    }
    this.voidHeat += K - Kt;
  }
  /* What the wall actually is, rather than what a heater was told to make it: the kinetic
     temperature of the fluid lying against it. With heat radiating away this is the only
     honest reading, and it is what the gauge shows. */
  _measureWall() {
    const b = this.box, p = this.pos, v = this.vel, skin = this.wallSkin;
    let K = 0, n = 0;
    for (let i = 0; i < this.N; i++) {
      if (this.pinned[i] || this.driven(i)) continue;
      const k = 3 * i;
      const distance = this.sphere ? this.sphere.R - Math.hypot(p[k] - this.sphere.x, p[k + 1] - this.sphere.y, p[k + 2] - this.sphere.z) :
        Math.min(p[k] - b.x0, b.x1 - p[k], p[k + 1] - b.y0, b.y1 - p[k + 1], p[k + 2] - b.z0, b.z1 - p[k + 2]);
      if (distance >= skin) continue;
      K += 0.5 * KEU * this.mass[i] * (v[k] * v[k] + v[k + 1] * v[k + 1] + v[k + 2] * v[k + 2]);
      n += 3;
    }
    this.wallMeasured = n ? 2 * K / (n * KB) : 0;
    this.wallContact = n / 3;
  }

  /* Pressure control: the chamber breathes toward the target instead of being held by hand.
     Whole molecules move with the walls; bond lengths are never scaled. */
  _barostat() {
    const P = this.pressureEMA, target = this.pressureTarget;
    const span = Math.max(1, Math.abs(target)) + 50;
    const gain = this.baroEvery * this.dt / Math.max(1, this.pressureTau);
    const mu = 1 + clampNum(gain * (P - target) / span, -0.004, 0.004);
    const b = this.box, cx = (b.x0 + b.x1) / 2, cy = (b.y0 + b.y1) / 2;
    const { list } = this.fragments();
    let minW=10,minH=10;
    // Pinned fragments stay fixed; a moving wall must not pass through them.
    for(const g of list) if(g.some(i=>this.pinned[i])) for(const i of g){
      minW=Math.max(minW,2*Math.abs(this.pos[3*i]-cx));
      minH=Math.max(minH,2*Math.abs(this.pos[3*i+1]-cy));
    }
    const wNew = clampNum((b.x1 - b.x0) * mu, Math.min(500,minW), 500), hNew = clampNum((b.y1 - b.y0) * mu, Math.min(500,minH), 500);
    const sx = wNew / (b.x1 - b.x0), sy = hNew / (b.y1 - b.y0);
    if (sx === 1 && sy === 1) return;
    b.x0 = cx - wNew / 2; b.x1 = cx + wNew / 2;
    b.y0 = cy - hNew / 2; b.y1 = cy + hNew / 2;
    for (const g of list) {
      if(g.some(i=>this.pinned[i]))continue;
      let M = 0, gx = 0, gy = 0;
      for (const i of g) { const m = this.mass[i]; M += m; gx += m * this.pos[3 * i]; gy += m * this.pos[3 * i + 1]; }
      gx /= M; gy /= M;
      const dx = cx + (gx - cx) * sx - gx, dy = cy + (gy - cy) * sy - gy;
      for (const i of g) { this.pos[3 * i] += dx; this.pos[3 * i + 1] += dy; this.prev[3 * i] += dx; this.prev[3 * i + 1] += dy; }
    }
    this.needRebuild = true; this.needForces = true;
  }

  /* ---------- integration ---------- */
  /* One simulation step = exactly dt (1 fs), velocity Verlet. A step that meets extreme curvature
     (a hydrogen squeezed between competing bonds, a hard collision) is split into sub-steps: the
     sub-step count keeps ω·dt ≤ 1 for the curvature each atom felt in the previous step, and a step whose
     energy error jumps is redone with finer sub-steps. Deterministic, so replay stays exact. */
  step() {
    if (this.recording && this.stepCount % this.ckEvery === 0) this._checkpoint();
    const N = this.N, n3 = 3 * N, pos = this.pos, vel = this.vel, F = this.frc, dt = this.dt;
    this.prev.set(pos.subarray(0, n3));
    if (this.needForces) { this._checkRebuild(); this.computeForces(); this.needForces = false; }
    const forceRebuild = this.stepCount % this.ckEvery === 0; // makes replay bit-identical
    let nsub = Math.max(1, Math.min(SUB_MAX, this.nextSub || 1));
    const canCheck = N > 0;
    let E0 = 0, K0 = 0;
    this._save();
    if (canCheck) { K0 = this.kinetic(); E0 = this.Epot + K0; }
    for (;;) {
      this._wallImpulse = 0;
      this._integrate(nsub, forceRebuild);
      if (!canCheck || nsub >= SUB_MAX) break;
      const dE = Math.abs(this.Epot + this.kinetic() + this.voidHeat - this._sv.voidHeat - this.servoWork - E0);
      if (dE <= SUB_ETOL + 0.002 * K0) break;
      this._load(); nsub = Math.min(SUB_MAX, nsub * 2); this.redone++; // redo this step at half the sub-step: a gentler jump shifts the conserved energy less
    }
    this.lastSub = nsub;
    this.servoWorkTotal += this.servoWork;
    // a step that had to be redone finer stays finer for a while, rather than dropping straight back
    if (nsub > (this.nextSub || 1)) { this.nextSub = nsub; this.subHold = SUB_HOLD; }
    // curvature felt by each atom this step → sub-steps for the next one (ω·dt ≤ 0.6)
    let w2max = 0, vmax2 = 0;
    for (let i = 0; i < N; i++) {
      if (this.pinned[i]) continue;
      const v2 = vel[3 * i] ** 2 + vel[3 * i + 1] ** 2 + vel[3 * i + 2] ** 2;
      if (v2 > vmax2) vmax2 = v2;
    }
    /* How finely to step comes from the things that are actually stiff, read from where they are
       now. Each formed bond sets a floor from its own curvature at its current length — a hot
       bond squeezed to 0.8 Å sits on a Morse wall nearly three times stiffer than at rest — over
       its reduced mass; so does a solid face an atom is about to meet; and anything moving fast
       enough to cross SUB_DX in a sub-step sets one too, which is what catches a spark or a
       violent impact driving atoms into each other's cores. What this replaced was an estimate
       read off how far each atom had just moved. A molecule flying across the chamber moves a long
       way whatever its bonds are doing, so it looked soft; and in a big chamber some atom always
       spiked it, so it flickered — and flicker is itself a source of drift. */
    { const need = Math.sqrt(vmax2) * dt / SUB_DX, w2 = need > 1 ? (0.5 * need / dt) ** 2 : 0; if (w2 > w2max) w2max = w2; }
    const type = this.type, mass = this.mass;
    for (let p = 0; p < this.nPairs; p++) {
      if (this.pF[p] < 0.3 || this.pR[p] <= 0) continue;
      const i = this.pI[p], j = this.pJ[p];
      const pp = PAIR[type[i] * NT + type[j]]; if (!pp.bond) continue;
      const lo = Math.min(2, Math.max(1, Math.floor(this.pN[p]))), a = pp.a[lo];
      const ey = Math.exp(-a * (this.pR[p] - pp.re[lo]));
      const curv = a * a * pp.De[lo] * Math.max(2, 4 * ey * ey - 2 * this.pB[p] * ey);
      const w2 = curv * (1 / mass[i] + 1 / mass[j]) * ACC;
      if (w2 > w2max) w2max = w2;
    }
    // an atom that could reach a solid face during the next step is stepped as finely as the face
    // needs from its first femtosecond of contact, not from the step after it has already arrived
    if (!this.sphere && this.boundsMode === 'solid') {
      const b = this.box, lo = [b.x0, b.y0, b.z0], hi = [b.x1, b.y1, b.z1];
      for (let i = 0; i < N; i++) {
        if (this.pinned[i]) continue;
        const k = 3 * i;
        let near = false;
        for (let d = 0; d < 3 && !near; d++) {
          const reach = Math.abs(vel[k + d]) * dt * 1.5 + SOLID_W;
          if (pos[k + d] < lo[d] + reach || pos[k + d] > hi[d] - reach) near = true;
        }
        if (near) { const w2 = this.solidK * ACC / mass[i]; if (w2 > w2max) w2max = w2; }
      }
    }
    /* ω·dt ≤ 0.5, about a dozen steps across the fastest vibration. At ω·dt ≤ 1 a cluster of
       fifteen water molecules gained 180 kJ/mol in 20 ps with no thermostat — the O–H stretch
       was being integrated with nine points per period. At 0.5 the same cluster keeps its energy
       to half a kJ/mol, and because an O–H or H–H bond is always present the count settles at a
       steady two instead of flickering between one and two, which is its own source of drift. */
    const need = Math.sqrt(w2max) * dt / 0.5;
    /* Go finer at once; come back down only after SUB_HOLD fs in which nothing asked for more.
       Velocity Verlet keeps energy only while its step length holds still: a bouncing water
       molecule stepped at a steady two sub-steps kept its energy to 1 kJ/mol over 20 ps, and the
       same molecule stepped at two-and-occasionally-three gained 128. The count has to change
       sometimes, but it should change rarely, not flicker with every compression. */
    const want = need > 1 ? Math.min(SUB_MAX, Math.ceil(need)) : 1, cur = this.nextSub || 1;
    if (want >= cur) { this.nextSub = want; this.subHold = SUB_HOLD; }
    else if ((this.subHold -= dt) <= 0) { this.nextSub = want; this.subHold = SUB_HOLD; }
    // Safety net: a numerically broken state is rolled back and paused; a merely very fast atom
    // (rare after sub-stepping) is capped at 60 km/s and counted, so hot gases keep running.
    let clamped = 0;
    for (let i = 0; i < N; i++) {
      if (this.pinned[i]) continue;
      const v2 = vel[3 * i] ** 2 + vel[3 * i + 1] ** 2 + vel[3 * i + 2] ** 2;
      if (!Number.isFinite(v2) || !Number.isFinite(pos[3 * i] + pos[3 * i + 1] + pos[3 * i + 2])) {
        this._load(); this.needForces = true;
        throw new Error('Simulation paused: numerical failure. Separate overlapping atoms or lower the temperature.');
      }
      if (v2 >= VMAX * VMAX) { clamped++; const s = VMAX / Math.sqrt(v2); vel[3 * i] *= s; vel[3 * i + 1] *= s; vel[3 * i + 2] *= s; }
    }
    this.clamped += clamped;
    this._voidWalls();
    const igniting = this.time < this.sparkHold;
    if (this.thermostat && !igniting) {
      if (this.thermostatMode === 'csvr') this._csvr();
      else if (this.thermostatMode === 'kelvin') this._kelvin();
      // A heater cannot warm a wall that radiates everything away, so it stands down and the
      // wall simply reports the fluid against it.
      else if (!this.voidTemperature) this._wallBath(dt);
    }
    // Radiating last gives it the final say: heat made during this step never survives it.
    if (!igniting) this._voidHeat();
    this._measureWall();
    if (this.voidTemperature) this.wallT = this.wallMeasured;
    this.time += dt; this.stepCount++;
    if (!this.sphere && this.boundsMode === 'solid') this.wallForce += this._wallImpulse / dt;
    // Normal momentum flux, exponentially averaged over ~1 ps
    const Pinst = this.wallArea > 0 ? this.wallForce / this.wallArea * BAR : 0;
    this.pressureBar = Pinst;
    this.pressureEMA += (Pinst - this.pressureEMA) * 0.001;
    if (this.pressureControl && this.stepCount % this.baroEvery === 0) this._barostat();
  }
  _integrate(nsub, forceRebuild) {
    const N = this.N, pos = this.pos, vel = this.vel, F = this.frc, h = this.dt / nsub;
    this.servoWork = 0;
    for (let sub = 0; sub < nsub; sub++) {
      const tf = this._servoActive ? this._twF : null;   // the servo force this sub-step starts with
      let W = 0;
      for (let i = 0; i < N; i++) {
        if (this.pinned[i]) { vel[3 * i] = vel[3 * i + 1] = vel[3 * i + 2] = 0; continue; }
        const hk = 0.5 * h * ACC / this.mass[i];
        for (let d = 3 * i; d < 3 * i + 3; d++) { vel[d] += hk * F[d]; const dx = vel[d] * h; pos[d] += dx; if (tf) W += tf[d] * dx; }
      }
      this.servoWork += W;
      this._reflectWalls();
      if (sub === 0 && forceRebuild) this.needRebuild = true;
      this._checkRebuild();
      this.computeForces(h);
      for (let i = 0; i < N; i++) {
        if (this.pinned[i]) continue;
        const hk = 0.5 * h * ACC / this.mass[i];
        for (let d = 3 * i; d < 3 * i + 3; d++) vel[d] += hk * F[d];
      }
    }
  }
  _copy(dst, src, n) { if (!dst || dst.length < n) dst = new Float64Array(Math.max(n, 64)); dst.set(src.subarray(0, n)); return dst; }
  _save() { // lightweight state for redoing one step
    const n3 = 3 * this.N, P = this.nPairs, sv = this._sv || (this._sv = {});
    sv.pos = this._copy(sv.pos, this.pos, n3); sv.vel = this._copy(sv.vel, this.vel, n3); sv.frc = this._copy(sv.frc, this.frc, n3); sv.built = this._copy(sv.built, this.built, n3);
    sv.pN = this._copy(sv.pN, this.pN, P);
    if (!sv.pI || sv.pI.length < P) { sv.pI = new Int32Array(Math.max(P, 64)); sv.pJ = new Int32Array(Math.max(P, 64)); }
    sv.pI.set(this.pI.subarray(0, P)); sv.pJ.set(this.pJ.subarray(0, P));
    sv.voidHeat = this.voidHeat;
    sv.P = P; sv.Epot = this.Epot; sv.map = this.pairMap; sv.needRebuild = this.needRebuild;
    sv.servoActive = this._servoActive;
    if (this._servoActive) sv.twF = this._copy(sv.twF, this._twF, n3);
    this.Fold = this._copy(this.Fold, this.frc, n3);
  }
  _load() {
    const sv = this._sv, n3 = 3 * this.N, P = sv.P;
    this.pos.set(sv.pos.subarray(0, n3)); this.vel.set(sv.vel.subarray(0, n3)); this.frc.set(sv.frc.subarray(0, n3)); this.built.set(sv.built.subarray(0, n3));
    if (P > this.pairCap) this._allocPairs(P);
    this.pI.set(sv.pI.subarray(0, P)); this.pJ.set(sv.pJ.subarray(0, P)); this.pN.set(sv.pN.subarray(0, P)); this.nPairs = P;
    this.voidHeat = sv.voidHeat;
    this.pairMap = sv.map; this.Epot = sv.Epot; this.needRebuild = sv.needRebuild;
    this._servoActive = sv.servoActive;
    if (sv.servoActive) { if (!this._twF || this._twF.length < n3) this._twF = new Float64Array(n3 + 192); this._twF.set(sv.twF.subarray(0, n3)); }
  }
  kinetic() {
    let K = 0; const v = this.vel;
    for (let i = 0; i < this.N; i++) { if (this.pinned[i]) continue; K += this.mass[i] * (v[3 * i] ** 2 + v[3 * i + 1] ** 2 + v[3 * i + 2] ** 2); }
    return 0.5 * KEU * K;
  }
  /* An atom held by the pointer is being driven from outside, so its motion is not thermal and
     must not be read as temperature. Counting it was what made dragging one molecule stop every
     other: a stat holding the total kinetic energy saw the drag as a huge excess and scaled the
     whole chamber down to compensate, and the dragged atom — re-accelerated by the tweezer every
     step — was the only thing left moving. */
  /* Everything the pointer is dragging, not just the atom it holds. Grab one atom of a molecule
     and its bonded partners are hauled along by their bonds; that motion is the drag's, not the
     sample's, and counting it as temperature let a stat scale the rest of the chamber to a halt. */
  _markDriven() {
    const N = this.N;
    if (!this.drivenMask || this.drivenMask.length < N) this.drivenMask = new Uint8Array(N + 100);
    const m = this.drivenMask;
    m.fill(0, 0, N);
    this._drivenFor = this.tweezer ? this.tweezer.i : -1;
    this._drivenAt = this.stepCount;
    if (!this.tweezer || this.tweezer.i >= N) return;
    m[this.tweezer.i] = 1;
    for (let pass = 0; pass < 8; pass++) {           // flood along bonds to the whole molecule
      let changed = false;
      for (let p = 0; p < this.nPairs; p++) {
        if (this.bondStrength(p) <= 0.25) continue;
        const i = this.pI[p], j = this.pJ[p];
        if (m[i] !== m[j]) { m[i] = m[j] = 1; changed = true; }
      }
      if (!changed) break;
    }
  }
  driven(i) {
    if (!this.tweezer) return false;
    if (this._drivenFor !== this.tweezer.i || this._drivenAt !== this.stepCount || !this.drivenMask) this._markDriven();
    return this.drivenMask[i] === 1;
  }
  /* Only the pointer's own contribution is excluded, which is the dragged cluster's centre-of-mass
     velocity — the coherent part the servo imposes. Everything else a dragged molecule does is
     still the sample's: it still vibrates, still rotates, and still has to give up the energy a
     new bond releases. Exempting those atoms wholesale meant a bond formed while dragging kept
     its own binding energy and tore itself straight back apart. */
  dragVelocity(out) {
    out[0] = out[1] = out[2] = 0;
    if (!this.tweezer) return false;
    let M = 0;
    for (let i = 0; i < this.N; i++) {
      if (this.pinned[i] || !this.driven(i)) continue;
      const w = this.mass[i]; M += w;
      for (let d = 0; d < 3; d++) out[d] += w * this.vel[3 * i + d];
    }
    if (M <= 0) return false;
    for (let d = 0; d < 3; d++) out[d] /= M;
    return true;
  }
  thermalKinetic() {
    const v = this.vel, vd = this._vdTmp || (this._vdTmp = [0, 0, 0]);
    const dragging = this.dragVelocity(vd);
    let K = 0;
    for (let i = 0; i < this.N; i++) {
      if (this.pinned[i]) continue;
      const k = 3 * i, d = dragging && this.driven(i);
      const ax = v[k] - (d ? vd[0] : 0), ay = v[k + 1] - (d ? vd[1] : 0), az = v[k + 2] - (d ? vd[2] : 0);
      K += this.mass[i] * (ax * ax + ay * ay + az * az);
    }
    return 0.5 * KEU * K;
  }
  // three degrees of freedom belong to the pointer while it holds a cluster; the rest are thermal
  dof() {
    let n = 0, held = false;
    for (let i = 0; i < this.N; i++) {
      if (this.pinned[i]) continue;
      n += 3;
      if (this.driven(i)) held = true;
    }
    return held ? Math.max(0, n - 3) : n;
  }
  temperature() { const nf = this.dof(); return nf ? 2 * this.thermalKinetic() / (nf * KB) : 0; }
  /* Finite-capacity wall reservoir with a first-order heater, then an exact
     Ornstein–Uhlenbeck velocity update at fixed position in the boundary layer.
     Noise and drag obey fluctuation–dissipation. Both act on the SAME atoms.
     Wall heat is accounted separately from mechanical work and reaction energy.
     All time constants use simulated fs, independent of rendering/playback. */
  _wallBath(dt) {
    const target = Math.max(288.15, Math.min(623.15, this.wallTarget));
    const dT = (target - this.wallT) * -Math.expm1(-dt / this.wallTau);
    this.wallT += dT;
    this.heaterWork += this.wallCapacity * dT;
    let heat = 0;
    this._bathBefore = this._copy(this._bathBefore, this.vel, 3 * this.N);
    const b = this.box, v = this.vel, p = this.pos;
    for (let i = 0; i < this.N; i++) {
      if (this.pinned[i] || this.driven(i)) continue;   // the pointer drives this one, not the bath
      const k = 3 * i;
      const distance = this.sphere ? this.sphere.R - Math.hypot(p[k] - this.sphere.x, p[k + 1] - this.sphere.y, p[k + 2] - this.sphere.z) :
        Math.min(p[k] - b.x0, b.x1 - p[k], p[k + 1] - b.y0, b.y1 - p[k + 1], p[k + 2] - b.z0, b.z1 - p[k + 2]);
      if (distance >= this.wallSkin) continue;
      const u = Math.max(0, distance / this.wallSkin);
      const weight = (1 - u) ** 2;
      const c = Math.exp(-dt * weight / this.wallCoupling);
      const variance = -Math.expm1(-2 * dt * weight / this.wallCoupling) * KB * this.wallT / (this.mass[i] * KEU);
      const sigma = Math.sqrt(Math.max(0, variance));
      for (let d = k; d < k + 3; d++) {
        const before = v[d];
        v[d] = c * before + sigma * this.gauss();
        heat += 0.5 * KEU * this.mass[i] * (v[d] ** 2 - before ** 2);
      }
    }
    // Positive heat enters the sample and leaves the wall. A finite reservoir
    // must not provide more energy than it owns; callers choose a macroscopic C.
    if (heat > this.wallCapacity * this.wallT) {
      // Reject an impossible reservoir exchange without corrupting the state.
      this.vel.set(this._bathBefore.subarray(0, 3 * this.N));
      return;
    }
    this.wallT -= heat / this.wallCapacity;
    this.heatToSample += heat;
  }
  setTemperature(T) {
    if (!Number.isFinite(T) || T < 0) throw new RangeError('Temperature must be finite and non-negative');
    this.checkpoints.length = 0;
    if (this.thermostat && this.thermostatMode === 'wall' && !this.voidTemperature) {
      this.wallTarget = Math.max(288.15, Math.min(623.15, T));
      this.T = this.wallTarget;
      return; // heater setpoint changes; neither walls nor sample jump
    }
    this.T = T;
    if (this.thermostat && this.thermostatMode === 'kelvin') { this._kelvinSet(); return; } // instant, by definition
    if (!this.thermostat) {
      // This is an explicit user intervention, not thermostatted dynamics.
      const current = this.temperature();
      if (T === 0) this.vel.fill(0, 0, 3 * this.N);
      else if (this.dof()) {
        if (current <= 1e-20) this.thermalize(T);
        const scale = Math.sqrt(T / this.temperature());
        for (let i = 0; i < this.N; i++) {
          if (this.pinned[i]) this.vel.fill(0, 3 * i, 3 * i + 3);
          else for (let d = 3 * i; d < 3 * i + 3; d++) this.vel[d] *= scale;
        }
      }
    }
  }
  /* Kelvin bath. The sample sits in surroundings at the set temperature that touch every atom — a
     solvent, a dense buffer gas — and each atom feels them as friction plus random kicks in exact
     balance (Langevin; an exact Ornstein–Uhlenbeck update of each velocity, the same law the wall
     bath uses in its boundary layer, here applied everywhere). The friction is that of a light
     buffer gas: τ for hydrogen, proportionally longer for heavier atoms (Epstein drag, γ ∝ 1/m).
     The same τ for every atom made xenon wade through the chamber like a liquid, twelve times
     slower than a real gas. That is what gives a real sample
     its Boltzmann spread of energies. What this did before — rescale every atom toward exactly its
     share, then the total toward exactly the setpoint — held the reading perfect and left every
     atom with the same energy: no fast tail at all (not one atom in 10^5 above 5 kT, against 1.9%
     in any real gas), and with it no activated chemistry at the rate a real sample has.
       The reading now fluctuates about the setpoint the way a real handful of molecules does, and a
     hot new molecule hands its excess to the surroundings over about τ. */
  _kelvin() {
    const v = this.vel, target = Math.max(0, this.T);
    const initialK = this.thermalKinetic();
    const vd = this._vdTmp || (this._vdTmp = [0, 0, 0]), dragging = this.dragVelocity(vd);
    for (let i = 0; i < this.N; i++) {
      if (this.pinned[i]) { v.fill(0, 3 * i, 3 * i + 3); continue; }
      // friction from a light buffer gas scales as 1/mass (Epstein drag): τ is set for hydrogen and
      // a heavier atom relaxes proportionally slower, so it still flies rather than wading
      const tau = Math.max(this.dt, this.tau * this.mass[i] / 1.008);
      const c = Math.exp(-this.dt / tau), c2 = -Math.expm1(-2 * this.dt / tau);
      const k = 3 * i, sigma = Math.sqrt(c2 * KB * target / (this.mass[i] * KEU));
      // an atom the pointer drags keeps the drag's velocity; the bath acts on its motion relative to it
      const d0 = dragging && this.driven(i);
      for (let d = 0; d < 3; d++) {
        const base = d0 ? vd[d] : 0;
        v[k + d] = base + c * (v[k + d] - base) + sigma * this.gauss();
      }
    }
    this.kelvinWork += this.thermalKinetic() - initialK;
  }
  /* Setting a temperature under the bath takes effect at once, as a one-off rescale; the bath then
     keeps it there with its natural fluctuations. */
  _kelvinSet() {
    const Nf = this.dof(); if (!Nf) return;
    const target = Math.max(0, this.T), initialK = this.thermalKinetic();
    if (target === 0) { for (let i = 0; i < this.N; i++) if (!this.driven(i)) this.vel.fill(0, 3 * i, 3 * i + 3); this.kelvinWork -= initialK; return; }
    let K = initialK;
    if (K <= 1e-14) { this.thermalize(target); K = this.thermalKinetic(); if (K <= 1e-14) return; }
    const sc = Math.sqrt(0.5 * Nf * KB * target / K), v = this.vel;
    const vd = this._vdTmp || (this._vdTmp = [0, 0, 0]), dragging = this.dragVelocity(vd);
    for (let i = 0; i < this.N; i++) {
      if (this.pinned[i]) { v.fill(0, 3 * i, 3 * i + 3); continue; }
      const k = 3 * i, d0 = dragging && this.driven(i);
      for (let d = 0; d < 3; d++) { const base = d0 ? vd[d] : 0; v[k + d] = base + (v[k + d] - base) * sc; }
    }
    this.kelvinWork += this.thermalKinetic() - initialK;
  }
  _csvr() {
    const Nf = this.dof(); if (!Nf) return;
    const K = this.thermalKinetic(), Kt = 0.5 * Nf * KB * Math.max(0, this.T);
    const c = Math.exp(-this.dt / Math.max(this.dt, this.tau));
    if (K <= 1e-12) { if (Kt > 0) this.thermalize(this.T * 0.05); return; }
    const r1 = this.gauss(), sum = Nf > 1 ? 2 * this.gamma((Nf - 1) / 2) : 0;
    let Kn = K + (1 - c) * (Kt * (r1 * r1 + sum) / Nf - K) + 2 * r1 * Math.sqrt(c * (1 - c) * Kt * K / Nf);
    if (Kn < 0) Kn = 0;
    // CSVR includes the sign of alpha, required for the canonical transition kernel.
    const sign = Kt > 0 && r1 + Math.sqrt(c * K * Nf / ((1 - c) * Kt)) < 0 ? -1 : 1;
    const s = sign * Math.sqrt(Kn / K), v = this.vel;
    for (let i = 0; i < this.N; i++) {
      if (this.pinned[i]) { v.fill(0, 3 * i, 3 * i + 3); continue; }
      if (this.driven(i)) continue;
      for (let k = 3 * i; k < 3 * i + 3; k++) v[k] *= s;
    }
  }
  thermalize(T, list) {
    if (!Number.isFinite(T) || T < 0) throw new RangeError('Temperature must be finite and non-negative');
    const idx = list || Array.from({ length: this.N }, (_, i) => i);
    for (const i of idx) {
      if (this.pinned[i]) { this.vel.fill(0, 3 * i, 3 * i + 3); continue; }
      const s = Math.sqrt(KB * T / (this.mass[i] * KEU));
      this.vel[3 * i] = s * this.gauss(); this.vel[3 * i + 1] = s * this.gauss(); this.vel[3 * i + 2] = s * this.gauss();
    }
  }
  scaleVelocities(list, f) { for (const i of list) { this.vel[3 * i] *= f; this.vel[3 * i + 1] *= f; this.vel[3 * i + 2] *= f; } }
  zeroMomentum(list) {
    const idx = (list || Array.from({ length: this.N }, (_, i) => i)).filter(i => !this.pinned[i]);
    let px = 0, py = 0, pz = 0, M = 0;
    for (const i of idx) { const m = this.mass[i]; px += m * this.vel[3 * i]; py += m * this.vel[3 * i + 1]; pz += m * this.vel[3 * i + 2]; M += m; }
    if (!M) return;
    for (const i of idx) { this.vel[3 * i] -= px / M; this.vel[3 * i + 1] -= py / M; this.vel[3 * i + 2] -= pz / M; }
  }
  /* Invalidate cached forces after editing positions/atoms. */
  touch() { this.needRebuild = true; this.needForces = true; this.checkpoints.length = 0; }
  refresh() { this._checkRebuild(); this.computeForces(); this.needForces = false; this.prev.set(this.pos.subarray(0, 3 * this.N)); }

  /* FIRE energy minimisation (used for absolute-zero conditioning). */
  minimize(maxIter = 2000, fTol = 0.5) {
    const N = this.N, v = this.vel, F = this.frc, pos = this.pos;
    let dt = 0.2, alpha = 0.1, npos = 0;
    const dtMax = 1.2;
    v.fill(0, 0, 3 * N);
    this._checkRebuild(); this.computeForces(1);
    for (let it = 0; it < maxIter; it++) {
      let P = 0, vn = 0, fn = 0, fmax = 0;
      for (let k = 0; k < 3 * N; k++) { P += F[k] * v[k]; vn += v[k] * v[k]; fn += F[k] * F[k]; fmax = Math.max(fmax, Math.abs(F[k])); }
      if (fmax < fTol) return it;
      vn = Math.sqrt(vn); fn = Math.sqrt(fn) || 1;
      if (P > 0) {
        for (let k = 0; k < 3 * N; k++) v[k] = (1 - alpha) * v[k] + alpha * F[k] / fn * vn;
        if (++npos > 5) { dt = Math.min(dt * 1.1, dtMax); alpha *= 0.99; }
      } else { v.fill(0, 0, 3 * N); dt *= 0.5; alpha = 0.1; npos = 0; }
      for (let i = 0; i < N; i++) {
        if (this.pinned[i]) continue;
        const h = dt * ACC / this.mass[i];
        for (let d = 3 * i; d < 3 * i + 3; d++) { v[d] += h * F[d]; let s = v[d] * dt; if (s > 0.1) s = 0.1; else if (s < -0.1) s = -0.1; pos[d] += s; }
      }
      this._checkRebuild(); this.computeForces(1);
    }
    v.fill(0, 0, 3 * N);
    return maxIter;
  }

  /* ---------- analysis ---------- */
  /* Bond strength in 0…1: the pair's saturation b times how much of the bond is left. Uses the wider
     saturation bond order as well as the structural switch, so a stretched or soft bond (a hot diatomic,
     an ionic contact) still counts as a bond. */
  bondStrength(p) {
    const f = this.pF[p], sr = this.pSraw[p];
    return this.pB[p] * (f > sr ? f : sr);
  }
  /* Bonds as {i, j, order, strength}. */
  bonds(threshold = 0.2) {
    const out = [];
    for (let p = 0; p < this.nPairs; p++) {
      if (this.pF[p] <= 0 && this.pSraw[p] <= 0) continue;
      const s = this.bondStrength(p); if (s < threshold) continue;
      out.push({ i: this.pI[p], j: this.pJ[p], order: this.pN[p], strength: s });
    }
    return out;
  }
  /* Connected fragments using bonds with strength > 0.25. Returns {comp: Int32Array, list: [[indices]]}. */
  fragments(threshold = 0.25) {
    const N = this.N, parent = new Int32Array(N);
    for (let i = 0; i < N; i++) parent[i] = i;
    const find = i => { while (parent[i] !== i) { parent[i] = parent[parent[i]]; i = parent[i]; } return i; };
    for (let p = 0; p < this.nPairs; p++) {
      if (this.bondStrength(p) > threshold) { const a = find(this.pI[p]), b = find(this.pJ[p]); if (a !== b) parent[a] = b; }
    }
    const comp = new Int32Array(N), groups = new Map();
    for (let i = 0; i < N; i++) { const r = find(i); if (!groups.has(r)) groups.set(r, []); groups.get(r).push(i); }
    const list = [...groups.values()];
    list.forEach((g, k) => { for (const i of g) comp[i] = k; });
    return { comp, list };
  }
  /* A written formula is not alphabetical, it is chemical. Two habits cover almost everything a
     chamber can make: the least electronegative element is written first, which gives NaCl, SO2
     and H2SO4; and hydrogen leads only in front of a chalcogen or a halogen, which is why water
     is H2O and ammonia is NH3, and why sodium hydroxide is NaOH rather than either. Carbon still
     comes first whenever there is any, as everyone writes it. */
  formulaOf(indices) {
    const cnt = {}; let charge = 0;
    for (const i of indices) { const s = ELEMENTS[this.type[i]].sym; cnt[s] = (cnt[s] || 0) + 1; charge += this.formal[i]; }
    const keys = Object.keys(cnt);
    const chi = k => BY_SYM[k].chi;
    const byChi = list => list.sort((a, b) => chi(a) - chi(b) || (a < b ? -1 : 1));
    let order;
    if (cnt.C) order = ['C', 'H', ...byChi(keys.filter(k => k !== 'C' && k !== 'H'))];
    else if (cnt.H && keys.length > 1) {
      const rest = byChi(keys.filter(k => k !== 'H')), hchi = chi('H');
      const leads = rest.some(k => H_LEADS.has(k)) && !rest.some(k => chi(k) < hchi);
      order = leads ? ['H', ...rest] : [...rest, 'H'];
    } else order = byChi(keys);
    let s = order.filter(k => cnt[k]).map(k => k + (cnt[k] > 1 ? cnt[k] : '')).join('');
    if (charge) s += (Math.abs(charge) > 1 ? Math.abs(Math.round(charge)) : '') + (charge > 0 ? '+' : '−');
    return s;
  }
  /* A fragment is a radical if any atom has unused valence: valence − Σ(bond order) ≥ 0.5. */
  isRadical(indices) {
    if (!this._boSum || this._boSum.length < this.N) this._boSum = new Float64Array(this.cap);
    if (this._boStamp !== this.stepCount + ':' + this.N + ':' + this.Epot) {
      const s = this._boSum; s.fill(0, 0, this.N);
      for (let p = 0; p < this.nPairs; p++) {
        if (this.bondStrength(p) <= 0.25) continue; // a bond counts in full, however stretched it is
        const v = this.pN[p]; s[this.pI[p]] += v; s[this.pJ[p]] += v;
      }
      this._boStamp = this.stepCount + ':' + this.N + ':' + this.Epot;
    }
    for (const i of indices) if (this.val[i] - this._boSum[i] > 0.5) return true;
    return false;
  }

  /* ---------- state, history ---------- */
  snapshot() {
    const N = this.N, n3 = 3 * N;
    const bo = [];
    for (let p = 0; p < this.nPairs; p++) if (this.pN[p] !== 1) bo.push(this._key(this.pI[p], this.pJ[p]), this.pN[p]);
    return {
      N, time: this.time, stepCount: this.stepCount, rng: this.rngState, nextId: this.nextId,
      box: { ...this.box }, sphere: this.sphere && { ...this.sphere }, T: this.T, tau: this.tau, thermostat: this.thermostat,
      thermostatMode: this.thermostatMode, kelvinWork: this.kelvinWork, sparkHold: this.sparkHold, wallT: this.wallT, wallTarget: this.wallTarget,
      wallMeasured: this.wallMeasured, wallContact: this.wallContact,
      boundsMode: this.boundsMode, fieldK: this.fieldK, fieldRange: this.fieldRange,
      voidTemperature: this.voidTemperature, voidPressure: this.voidPressure, voidVelocity: this.voidVelocity, voidTau: this.voidTau, voidSkin: this.voidSkin,
      pressureControl: this.pressureControl, pressureTarget: this.pressureTarget, pressureTau: this.pressureTau,
      voidHeat: this.voidHeat, voidForce: this.voidForce,
      servoWork: this.servoWork, servoWorkTotal: this.servoWorkTotal,
      wallTau: this.wallTau, wallCapacity: this.wallCapacity, wallSkin: this.wallSkin, wallCoupling: this.wallCoupling,
      heatToSample: this.heatToSample, heaterWork: this.heaterWork,
      nextSub: this.nextSub || 1, subHold: this.subHold || 0, lastSub: this.lastSub || 1, Epot: this.Epot, Ewall: this.Ewall,
      wallForce: this.wallForce, wallArea: this.wallArea, pressureBar: this.pressureBar, pressureEMA: this.pressureEMA,
      needForces: this.needForces, needRebuild: this.needRebuild, redone: this.redone, clamped: this.clamped,
      built: this.built.slice(0, n3), pI: this.pI.slice(0, this.nPairs), pJ: this.pJ.slice(0, this.nPairs), pN: this.pN.slice(0, this.nPairs),
      pos: this.pos.slice(0, n3), vel: this.vel.slice(0, n3), frc: this.frc.slice(0, n3),
      type: this.type.slice(0, N), formal: this.formal.slice(0, N), val: this.val.slice(0, N), lp: this.lp.slice(0, N),
      pinned: this.pinned.slice(0, N), ids: this.ids.slice(0, N), cos0: this.cos0.slice(0, N), bo: Float64Array.from(bo)
    };
  }
  restore(s) {
    if (s.N > this.cap) this._alloc(Math.max(s.N, 2 * this.cap));
    this.N = s.N; this.time = s.time; this.stepCount = s.stepCount; this.rngState = s.rng; this.nextId = s.nextId;
    if (s.box) this.box = { ...s.box };
    this.sphere = s.sphere ? { ...s.sphere } : null;
    for (const key of ['T', 'tau', 'thermostat', 'thermostatMode', 'kelvinWork', 'sparkHold', 'wallMeasured', 'wallContact', 'wallT', 'wallTarget', 'wallTau', 'wallCapacity', 'wallSkin', 'wallCoupling', 'boundsMode', 'fieldK', 'fieldRange', 'voidTemperature', 'voidPressure', 'voidVelocity', 'voidTau', 'voidSkin', 'voidHeat', 'voidForce', 'servoWork', 'servoWorkTotal', 'pressureControl', 'pressureTarget', 'pressureTau', 'heatToSample', 'heaterWork', 'nextSub', 'subHold', 'lastSub', 'redone', 'clamped']) if (s[key] !== undefined) this[key] = s[key];
    this.tweezer = null;
    this.pos.set(s.pos); this.vel.set(s.vel); this.frc.set(s.frc); this.prev.set(s.pos);
    this.type.set(s.type); this.formal.set(s.formal); this.val.set(s.val); this.lp.set(s.lp); this.pinned.set(s.pinned); this.ids.set(s.ids); this.cos0.set(s.cos0);
    for (let i = 0; i < this.N; i++) { const el = ELEMENTS[this.type[i]]; this.mass[i] = el.mass; this.ljx[i] = el.ljX; this.lje[i] = Math.sqrt(el.ljD); }
    // rebuild pair list with the saved bond orders
    this.nPairs = 0;
    this._rebuild();
    const bo = new Map(); for (let k = 0; k < s.bo.length; k += 2) bo.set(s.bo[k], s.bo[k + 1]);
    for (let p = 0; p < this.nPairs; p++) { const v = bo.get(this._key(this.pI[p], this.pJ[p])); this.pN[p] = v === undefined ? 1 : v; }
    this.needRebuild = true; this.needForces = false;
    // recompute derived per-pair data (Z, q, b) for display without touching forces
    const f = this.frc.slice(0, 3 * this.N), pn = this.pN.slice(0, this.nPairs), c0 = this.cos0.slice(0, this.N);
    if (s.pI) {
      if (s.pI.length > this.pairCap) this._allocPairs(s.pI.length);
      this.nPairs = s.pI.length; this.pI.set(s.pI); this.pJ.set(s.pJ); this.pN.set(s.pN); this.built.set(s.built);
      this.pairMap = new Map();
      for (let p = 0; p < this.nPairs; p++) this.pairMap.set(this.pI[p] * 1048576 + this.pJ[p], p);
    }
    this.computeForces(); this.frc.set(f); this.cos0.set(c0);
    if (!s.pI) this.pN.set(pn.subarray(0, Math.min(pn.length, this.nPairs)));
    for (const key of ['Epot', 'Ewall', 'wallForce', 'wallArea', 'pressureBar', 'pressureEMA']) if (s[key] !== undefined) this[key] = s[key];
    this.needRebuild = s.needRebuild ?? true; this.needForces = s.needForces ?? false;
  }
  _checkpoint() {
    const s = this.snapshot();
    const bytes = s.pos.byteLength * 4 + s.N * 40 + s.bo.byteLength + s.pN.byteLength * 2 + 400;
    this.checkpoints.push(s);
    const max = Math.max(8, Math.floor(this.ckBudget / bytes));
    if (this.checkpoints.length > max) this.checkpoints.splice(0, this.checkpoints.length - max);
  }
  /* Go back exactly one step: restore the nearest checkpoint and replay. */
  stepBack(n = 1) {
    const target = this.stepCount - n;
    if (target < 0) return false;
    let ck = null;
    while (this.checkpoints.length) {
      const c = this.checkpoints[this.checkpoints.length - 1];
      if (c.stepCount <= target) { ck = c; break; }
      this.checkpoints.pop();
    }
    if (!ck) return false;
    const last = this.snapshot();
    this.restore(ck);
    this.checkpoints.pop(); // step() re-records it
    const rec = this.recording;
    while (this.stepCount < target) this.step();
    this.recording = rec;
    // interpolate from the later state back to this one for display
    this.prev.set(last.pos.subarray(0, Math.min(last.pos.length, 3 * this.N)));
    return true;
  }
  canStepBack() { return this.checkpoints.length > 0 && this.checkpoints[0].stepCount <= this.stepCount - 1; }
  historySpan() { return this.checkpoints.length ? this.stepCount - this.checkpoints[0].stepCount : 0; }

  /* Serialisable scene (for save/load). */
  toJSON() {
    const atoms = [];
    for (let i = 0; i < this.N; i++) atoms.push([ELEMENTS[this.type[i]].sym, +this.pos[3 * i].toFixed(4), +this.pos[3 * i + 1].toFixed(4), +this.pos[3 * i + 2].toFixed(4), +this.vel[3 * i].toFixed(6), +this.vel[3 * i + 1].toFixed(6), +this.vel[3 * i + 2].toFixed(6), this.formal[i], this.val[i]]);
    return { format: 'chem-playground/scene@1', box: this.box, T: this.T, tau: this.tau, thermostat: this.thermostat, thermostatMode: this.thermostatMode, kelvinWork: this.kelvinWork, boundsMode: this.boundsMode, voidTemperature: this.voidTemperature, voidPressure: this.voidPressure, voidVelocity: this.voidVelocity, dampingVersion: 3, voidTau: this.voidTau, voidSkin: this.voidSkin, voidHeat: this.voidHeat, pressureControl: this.pressureControl, pressureTarget: this.pressureTarget, wallT: this.wallT, wallTarget: this.wallTarget, wallTau: this.wallTau, heatToSample: this.heatToSample, heaterWork: this.heaterWork, time: this.time, atoms };
  }
}

return { Engine, ELEMENTS, BY_SYM, PAIR, TUNE, refreshSaturation, KB, KEU, BAR, valenceFor, lonePairs };
});
