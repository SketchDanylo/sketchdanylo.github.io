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
 *   - Velocity-Verlet, fixed dt = 1 fs, Bussi CSVR thermostat, soft container walls whose
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
const COUL = 1389.35458;            // kJ·Å/(mol·e²)
const BAR = 16605.39;               // (kJ/mol)/Å³ → bar
const KAPPA = 0.32;                 // e per unit Pauling electronegativity difference per bond
const Q_MAX = 1.1;                  // saturation of bond-polarisation charge (e)
const TUNE = { c1: 1.3, c4: 20 };   // saturation p(x) = 1/(1 + c1·x + c4·x⁴)
const RAMP_W = 0.15;
const Y_ON = 3.6, Y_OFF = 5.6;      // Morse taper window (in units of a·(r − re))
const COORD_R1 = 1.22, COORD_R2 = 1.50; // structural coordination switch (× single-bond length)
const BO_CAP = 0.10, BO_W = 0.08, BO_BETA = 3.6; // saturation bond order: 1 up to r₁+cap, then exp(−β·Δr)
const BO_ON = 0.50, BO_OFF = 0.85;  // …tapered to zero between r₁+on and r₁+off (Å)
const COUL_D2 = 0.3;                // Coulomb short-range shielding (Å²)
const LJ_SHIELD = 0.55;             // LJ shielding radius as a fraction of r_min
const SUB_MAX = 16, SUB_ETOL = 6;  // at most 16 sub-steps; redo a step whose energy error exceeds 6 kJ/mol
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
  'C|O': { 1: [1.43, 358, 540], 2: [1.21, 745, 1200], 3: [1.13, 1072, 1900] },
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
    core: 0.6 * (ei.rcov[0] + ej.rcov[0])
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
    pp.s2 = pp.re[1] + BO_OFF;
  }
  return pp;
}
const PAIR = new Array(NT * NT);
for (let i = 0; i < NT; i++) for (let j = 0; j < NT; j++) PAIR[i * NT + j] = buildPairParams(ELEMENTS[i], ELEMENTS[j]);

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
  if (d <= BO_CAP) { SF = 1; SFD = 0; return; }
  if (d >= BO_OFF) { SF = 0; SFD = 0; return; }
  const x = d - BO_CAP;
  let s, ds;
  if (x < BO_W) { s = x * x / (2 * BO_W); ds = x / BO_W; } else { s = x - BO_W / 2; ds = 1; }
  const e = Math.exp(-BO_BETA * s), de = -BO_BETA * ds * e;
  smoothSwitch(d, BO_ON, BO_OFF);
  SF = e * SW; SFD = de * SW + e * SWD;
}
let SP = 0, SPD = 0;          // saturation p(x), dp/dx
function sat(x) { // p(x) = 1 / (1 + c1·r + c4·r⁴), r = smooth ramp(x)
  if (x <= 0) { SP = 1; SPD = 0; return; }
  let r, dr;
  if (x < RAMP_W) { r = x * x / (2 * RAMP_W); dr = x / RAMP_W; } else { r = x - RAMP_W / 2; dr = 1; }
  const c1 = TUNE.c1, c4 = TUNE.c4, r3 = r * r * r, den = 1 + c1 * r + c4 * r3 * r;
  SP = 1 / den; SPD = -(c1 + 4 * c4 * r3) * dr / (den * den);
}
/* Shielded Lennard-Jones split into Pauli (LR) and shifted-force attraction (LA), each with d/dr.
   Shielding: r_s⁶ = r⁶ + s⁶, so the Pauli wall stays finite at contact. x2 = r_min², eps = depth. */
let LR = 0, LDR = 0, LA = 0, LDA = 0;
const S6 = Math.pow(LJ_SHIELD, 6);
function lj(r, x2, eps, rc) {
  const x6 = x2 * x2 * x2, r2 = r * r, r6 = r2 * r2 * r2, d6 = r6 + S6 * x6;
  const u6 = x6 / d6, du6 = -6 * u6 * r2 * r2 * r / d6; // d(u⁶)/dr
  const rc2 = rc * rc, rc6 = rc2 * rc2 * rc2, uc6 = x6 / (rc6 + S6 * x6);
  const ljc = eps * (uc6 * uc6 - 2 * uc6), dljc = eps * (2 * uc6 - 2) * (-6 * uc6 * rc2 * rc2 * rc / (rc6 + S6 * x6));
  const v = eps * (u6 * u6 - 2 * u6), dv = eps * (2 * u6 - 2) * du6;
  if (u6 > 1) { // inside r_min: Pauli wall, attraction held at its (shifted) value at r_min
    LR = v + eps; LDR = dv;
    const rm = Math.sqrt(Math.cbrt(x6 - S6 * x6));
    LA = -eps - ljc - dljc * (rm - rc); LDA = 0;
  } else { LR = 0; LDR = 0; LA = v - ljc - dljc * (r - rc); LDA = dv - dljc; }
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
    this.T = opts.T ?? 298.15;
    this.tau = opts.tau ?? 250;     // thermostat coupling time, fs
    this.thermostat = opts.thermostat ?? true;
    this.rc = opts.rc || 8.0;
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
      mass: f64(cap), phi: f64(cap), Zs: f64(cap), qg: f64(cap), Gz: f64(cap), formal: f64(cap), q: f64(cap), Z: f64(cap), cos0: f64(cap), G: f64(cap), ljx: f64(cap), lje: f64(cap),
      type: new Int16Array(cap), val: new Int8Array(cap), lp: new Int8Array(cap), pinned: new Uint8Array(cap),
      ids: new Uint32Array(cap), cStart: new Int32Array(cap + 1), cCount: new Int32Array(cap)
    };
    for (const k in arrays) { if (old) arrays[k].set(this[k].subarray(0, Math.min(this[k].length, arrays[k].length))); this[k] = arrays[k]; }
    this.cap = cap;
  }
  _allocPairs(cap) {
    const f64 = n => new Float64Array(n);
    const keep = this.pairCap ? { pI: this.pI, pJ: this.pJ, pN: this.pN } : null;
    this.pI = new Int32Array(cap); this.pJ = new Int32Array(cap); this.pN = f64(cap);
    this.pR = f64(cap); this.pDx = f64(cap); this.pDy = f64(cap); this.pDz = f64(cap);
    this.pF = f64(cap); this.pFp = f64(cap); this.pS = f64(cap); this.pSp = f64(cap); this.pScr = f64(cap); this.pSraw = f64(cap);
    if (!this.tri) { this.tri = new Int32Array(4096); this.nTri = 0; } this.pB = f64(cap); this.gA = f64(cap); this.gB = f64(cap);
    this.cList = new Int32Array(2 * cap);
    if (keep) { this.pI.set(keep.pI.subarray(0, this.nPairs)); this.pJ.set(keep.pJ.subarray(0, this.nPairs)); this.pN.set(keep.pN.subarray(0, this.nPairs)); }
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
      this.pI[np] = a; this.pJ[np] = b; this.pN[np] = this.pN[p]; np++;
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
    // remember bond orders by atom-id pair
    const keep = new Map();
    for (let p = 0; p < this.nPairs; p++) if (this.pN[p] !== 1) keep.set(this._key(this.pI[p], this.pJ[p]), this.pN[p]);
    let minx = Infinity, miny = Infinity, minz = Infinity, maxx = -Infinity, maxy = -Infinity, maxz = -Infinity;
    for (let i = 0; i < N; i++) {
      const x = pos[3 * i], y = pos[3 * i + 1], z = pos[3 * i + 2];
      if (x < minx) minx = x; if (x > maxx) maxx = x; if (y < miny) miny = y; if (y > maxy) maxy = y; if (z < minz) minz = z; if (z > maxz) maxz = z;
    }
    const nx = Math.max(1, Math.min(64, Math.floor((maxx - minx) / rl) + 1));
    const ny = Math.max(1, Math.min(64, Math.floor((maxy - miny) / rl) + 1));
    const nz = Math.max(1, Math.min(16, Math.floor((maxz - minz) / rl) + 1));
    const cx = Math.max(rl, (maxx - minx) / nx + 1e-9), cy = Math.max(rl, (maxy - miny) / ny + 1e-9), cz = Math.max(rl, (maxz - minz) / nz + 1e-9);
    const nc = nx * ny * nz, head = new Int32Array(nc).fill(-1), next = new Int32Array(N), cell = new Int32Array(N);
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
              const n = keep.get(this._key(i, j)); this.pN[np] = n === undefined ? 1 : n;
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
  computeForces() {
    const N = this.N, P = this.nPairs, pos = this.pos, F = this.frc, Z = this.Z, q = this.q, type = this.type, val = this.val;
    const pI = this.pI, pJ = this.pJ, pR = this.pR, pDx = this.pDx, pDy = this.pDy, pDz = this.pDz, pF = this.pF, pFp = this.pFp;
    const rc = this.rc, rc2 = rc * rc;
    F.fill(0, 0, 3 * N);
    const phi = this.phi, Zs = this.Zs, pS = this.pS, pSp = this.pSp, qg = this.qg;
    for (let i = 0; i < N; i++) { Z[i] = 0; Zs[i] = 0; q[i] = 0; this.G[i] = 0; this.cCount[i] = 0; phi[i] = 0; }
    let E = 0;

    // Pass A — geometry, coordination numbers, bond-polarisation charges
    for (let p = 0; p < P; p++) {
      const i = pI[p], j = pJ[p];
      const dx = pos[3 * j] - pos[3 * i], dy = pos[3 * j + 1] - pos[3 * i + 1], dz = pos[3 * j + 2] - pos[3 * i + 2];
      const r2 = dx * dx + dy * dy + dz * dz;
      pF[p] = 0; pFp[p] = 0; pS[p] = 0; pSp[p] = 0;
      if (r2 > rc2) { pR[p] = -1; continue; }
      const r = Math.sqrt(r2) || 1e-6;
      pR[p] = r; pDx[p] = dx; pDy[p] = dy; pDz[p] = dz;
      const pp = PAIR[type[i] * NT + type[j]];
      if (pp.bond && r < pp.s2) {
        satF(r, pp);
        if (SF > 0) { pS[p] = SF; pSp[p] = SFD; }
      }
      if (pp.bond && r < pp.r2) {
        smoothSwitch(r, pp.r1, pp.r2);
        if (SW > 0) {
          pF[p] = SW; pFp[p] = SWD;
          Z[i] += SW; Z[j] += SW;
          this.cCount[i]++; this.cCount[j]++;
          const dq = KAPPA * SW * (CHI[type[j]] - CHI[type[i]]);
          q[i] += dq; q[j] -= dq;
        }
      }
    }
    // Bond-polarisation charge saturates smoothly (an atom cannot give away unlimited charge):
    // q = formal + raw / (1 + (raw/qmax)⁴)^¼,  dq/draw = (1 + (raw/qmax)⁴)^(−5/4)
    // (qg declared above)
    for (let i = 0; i < N; i++) {
      const raw = q[i], u = raw / Q_MAX, u4 = u * u * u * u, den = Math.pow(1 + u4, 0.25);
      q[i] = this.formal[i] + raw / den; qg[i] = 1 / (den * (1 + u4));
    }
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
    for (let p = 0; p < P; p++) { scr[p] = 1; pSraw[p] = pS[p]; }
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
    for (let p = 0; p < P; p++) {
      if (pS[p] === 0) continue;
      const S = scr[p] * this.pN[p]; // valence used = screened bond order × multiplicity
      pS[p] *= S; pSp[p] *= S;
      Zs[pI[p]] += pS[p]; Zs[pJ[p]] += pS[p];
    }
    // polar hydrogens get a small Pauli radius (hydrogen bonding)
    const ljx = this.ljx;
    for (let i = 0; i < N; i++) {
      const el = ELEMENTS[type[i]];
      if (el.Z === 1) { const t = Math.min(1, Math.max(0, (q[i] - this.formal[i] - 0.15) / 0.15)); ljx[i] = el.ljX + (1.7 - el.ljX) * t; }
    }

    // Pass B — relax continuous bond orders toward the spare-valence target
    this._updateBondOrders();

    // Pass C — pair energies
    const pB = this.pB, gA = this.gA, gB = this.gB, pN = this.pN, G = this.G, lje = this.lje;
    const gc = 1 / Math.sqrt(rc * rc + COUL_D2), dgc = -rc * gc * gc * gc;
    for (let p = 0; p < P; p++) {
      const r = pR[p];
      gA[p] = 0; gB[p] = 0; pB[p] = 0;
      if (r < 0) continue;
      const i = pI[p], j = pJ[p], ti = type[i], tj = type[j];
      const pp = PAIR[ti * NT + tj], f = pF[p], fp = pFp[p], fs = pS[p];
      let e = 0, dEdr = 0, lam = 0, pi = 0, pj = 0, dpi = 0, dpj = 0, b = 0;
      if (pp.bond) {
        if (val[i] > 0) { sat(Zs[i] - fs + 1 - val[i]); pi = SP; dpi = SPD; }
        if (val[j] > 0) { sat(Zs[j] - fs + 1 - val[j]); pj = SP; dpj = SPD; }
        b = pi * pj;
        // Morse with bond-order-interpolated parameters
        const n = pN[p], lo = Math.min(2, Math.floor(n)), t = n - lo, hi = lo + 1;
        const tl = 1 - Math.pow(1 - t, 1.6); // bond length contracts fastest at low fractional order (Pauling/resonance)
        const re = pp.re[lo] + (pp.re[Math.min(hi, 3)] - pp.re[lo]) * tl;
        const De = pp.De[lo] + (pp.De[Math.min(hi, 3)] - pp.De[lo]) * t;
        const a = pp.a[lo] + (pp.a[Math.min(hi, 3)] - pp.a[lo]) * t;
        const y = a * (r - re);
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
      // Lennard-Jones (shielded); Pauli part weighted by (1 − b)
      lj(r, ljx[i] * ljx[j], lje[i] * lje[j], rc);
      e += (1 - b) * LR; dEdr += (1 - b) * LDR; lam -= LR;
      // excluded (1-2) dispersion + electrostatics
      const qq = q[i] * q[j];
      let h = 0, dh = 0;
      if (qq !== 0 || q[i] !== 0 || q[j] !== 0) { const g = 1 / Math.sqrt(r * r + COUL_D2); h = g - gc - dgc * (r - rc); dh = -r * g * g * g - dgc; }
      const ex = LA + COUL * qq * h, dex = LDA + COUL * qq * dh;
      e += (1 - f) * ex; dEdr += (1 - f) * dex - fp * ex;
      if (h !== 0) { const c = (1 - f) * COUL * h; phi[i] += c * q[j]; phi[j] += c * q[i]; }
      if (lam !== 0 && pp.bond) {
        const ga = lam * pj * dpi, gb = lam * pi * dpj;
        gA[p] = ga; gB[p] = gb; G[i] += ga; G[j] += gb;
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
    // Pass D — many-body forces from the coordination dependence of b
    for (let p = 0; p < P; p++) {
      const fp = pSp[p]; if (fp === 0) continue;
      const i = pI[p], j = pJ[p];
      const coef = (G[i] - gA[p]) + (G[j] - gB[p]);
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
      const dEdr = KAPPA * fp * dchi * (phi[i] * qg[i] - phi[j] * qg[j]) + (Gz[i] + Gz[j]) * fp;
      if (dEdr === 0) continue;
      const s = dEdr / pR[p], fx = s * pDx[p], fy = s * pDy[p], fz = s * pDz[p];
      F[3 * i] += fx; F[3 * i + 1] += fy; F[3 * i + 2] += fz;
      F[3 * j] -= fx; F[3 * j + 1] -= fy; F[3 * j + 2] -= fz;
    }
    // Pass F — container walls and tweezers
    E += this._walls();
    if (this.tweezer) {
      const tw = this.tweezer, i = tw.i;
      if (i < N) {
        const dx = pos[3 * i] - tw.x, dy = pos[3 * i + 1] - tw.y, k = tw.k || 25;
        const d = Math.hypot(dx, dy), cap = tw.maxF || 400;
        let s = k; if (k * d > cap) s = cap / d;
        F[3 * i] -= s * dx; F[3 * i + 1] -= s * dy;
        F[3 * i + 2] -= 0.5 * k * pos[3 * i + 2] * 0.2;
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
  _updateBondOrders() {
    const N = this.N, Z = this.Z, val = this.val, cStart = this.cStart, cList = this.cList, pF = this.pF, pI = this.pI, pJ = this.pJ, pN = this.pN, type = this.type;
    const spare = this._spare && this._spare.length >= N ? this._spare : (this._spare = new Float64Array(this.cap));
    const D = this._D && this._D.length >= N ? this._D : (this._D = new Float64Array(this.cap));
    for (let i = 0; i < N; i++) spare[i] = Math.max(0, val[i] - Z[i]);
    for (let i = 0; i < N; i++) {
      let s = 0;
      for (let c = cStart[i]; c < cStart[i + 1]; c++) { const p = cList[c], k = pI[p] === i ? pJ[p] : pI[p]; s += pF[p] * Math.min(2, spare[k]); }
      D[i] = s;
    }
    const rate = 0.08 * (this._rate || 1); // relaxation per fs, independent of sub-stepping
    for (let p = 0; p < this.nPairs; p++) {
      const f = pF[p];
      let target = 1;
      if (f > 0) {
        const i = pI[p], j = pJ[p];
        const pp = PAIR[type[i] * NT + type[j]];
        if (pp.maxOrder > 1 && D[i] > 1e-9 && D[j] > 1e-9) {
          const si = spare[i] * f * Math.min(2, spare[j]) / D[i];
          const sj = spare[j] * f * Math.min(2, spare[i]) / D[j];
          target = 1 + Math.min(si, sj, pp.maxOrder - 1);
        }
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
      // VSEPR: ideal cos θ0 is a smooth function of the continuous steric number SN = Z + lone pairs
      const lp = this.lp[i], sn = Z[i] + lp, c4 = TETRA[Math.min(3, lp)];
      let c0v, dc0;
      if (sn <= 2) { c0v = -1; dc0 = 0; }
      else if (sn < 3) { const t = sn - 2; c0v = -1 + 0.5 * (3 * t * t - 2 * t * t * t); dc0 = 0.5 * (6 * t - 6 * t * t); }
      else if (sn < 4) { const t = sn - 3; c0v = -0.5 + (c4 + 0.5) * (3 * t * t - 2 * t * t * t); dc0 = (c4 + 0.5) * (6 * t - 6 * t * t); }
      else { c0v = c4; dc0 = 0; }
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
      Gz[i] = dEdc0sum * dc0; // dE/dZ_i through θ0
    }
    return E;
  }
  _screenForces() {
    const tri = this.tri, pS = this.pSraw, pN = this.pN, scr = this.pScr, G = this.G, gA = this.gA, gB = this.gB;
    const pI = this.pI, pJ = this.pJ, pF = this.pF, pFp = this.pFp, pR = this.pR, pDx = this.pDx, pDy = this.pDy, pDz = this.pDz, F = this.frc;
    for (let t = 0; t < this.nTri; t++) {
      const q = tri[4 * t], pa = tri[4 * t + 1], pb = tri[4 * t + 2];
      const j = pI[q], k = pJ[q];
      const dEdS = pS[q] * pN[q] * ((G[j] - gA[q]) + (G[k] - gB[q]));
      if (dEdS === 0) continue;
      const m = 1 - pF[q], w = pF[pa] * pF[pb] * m, one = 1 - w;
      let others; // Π over the other shared neighbours
      if (one > 1e-9) others = scr[q] / one;
      else { others = 1; for (let u = 0; u < this.nTri; u++) if (u !== t && tri[4 * u] === q) others *= 1 - pF[tri[4 * u + 1]] * pF[tri[4 * u + 2]] * m; }
      const dEdw = -dEdS * others;
      // w = f(r_a)·f(r_b)·(1 − f(r_jk)): push along each arm of the triple and along j–k
      for (const [p, dw] of [[pa, pFp[pa] * pF[pb] * m], [pb, pF[pa] * pFp[pb] * m], [q, -pF[pa] * pF[pb] * pFp[q]]]) {
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
    let b = 0, fjk = 0, fpjk = 0, X = 0, dX = 0, dpj = 0, dpk = 0, bpj = 0, bpk = 0;
    if (pp.bond) {
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
    X += (1 - b) * LR + LA + COUL * qq * h;
    dX += (1 - b) * LDR + LDA + this._dcoul(r, qq);
    const m = 1 - fjk; // do not double-exclude a pair that is itself bonded (3-rings)
    const Esub = -w * m * X;
    // Esub depends on b = p_j·p_k through the Pauli term: dEsub/db = w·m·LR → saturation coefficients
    if (LR !== 0 && (dpj !== 0 || dpk !== 0)) { const c = w * m * LR; this.G[j] += c * bpk * dpj; this.G[k] += c * bpj * dpk; }
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
  _walls() {
    const N = this.N, pos = this.pos, F = this.frc, K = this.wallK;
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
      for (let i = 0; i < N; i++) for (let d = 0; d < 3; d++) {
        const x = pos[3 * i + d];
        if (x < lo[d]) { const e = lo[d] - x; E += 0.5 * K * e * e; F[3 * i + d] += K * e; Fsum += K * e; }
        else if (x > hi[d]) { const e = x - hi[d]; E += 0.5 * K * e * e; F[3 * i + d] -= K * e; Fsum += K * e; }
      }
      const Lx = b.x1 - b.x0, Ly = b.y1 - b.y0, Lz = b.z1 - b.z0;
      this.wallArea = 2 * (Lx * Ly + Lx * Lz + Ly * Lz);
    }
    this.Ewall = E; this.wallForce = Fsum;
    return E;
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
    const canCheck = !this.tweezer && N > 0;
    let E0 = 0, K0 = 0;
    if (canCheck) { K0 = this.kinetic(); E0 = this.Epot + K0; this._save(); }
    else this.Fold = this._copy(this.Fold, F, n3);
    for (;;) {
      this._integrate(nsub, forceRebuild);
      if (!canCheck || nsub >= SUB_MAX) break;
      const dE = Math.abs(this.Epot + this.kinetic() - E0);
      if (dE <= SUB_ETOL + 0.005 * K0) break;
      this._load(); nsub = Math.min(SUB_MAX, nsub * 4); this.redone++; // redo this step more finely
    }
    this.lastSub = nsub;
    // curvature felt by each atom this step → sub-steps for the next one (ω·dt ≤ 0.6)
    let w2max = 0;
    const Fo = this.Fold, prev = this.prev;
    for (let i = 0; i < N; i++) {
      if (this.pinned[i]) continue;
      const dx = pos[3 * i] - prev[3 * i], dy = pos[3 * i + 1] - prev[3 * i + 1], dz = pos[3 * i + 2] - prev[3 * i + 2];
      const d2 = dx * dx + dy * dy + dz * dz; if (d2 < 1e-10) continue;
      const fx = F[3 * i] - Fo[3 * i], fy = F[3 * i + 1] - Fo[3 * i + 1], fz = F[3 * i + 2] - Fo[3 * i + 2];
      const w2 = Math.sqrt((fx * fx + fy * fy + fz * fz) / d2) * ACC / this.mass[i];
      if (w2 > w2max) w2max = w2;
    }
    const need = Math.sqrt(w2max) * dt / 1.0;
    this.nextSub = need > 1 ? Math.min(SUB_MAX, Math.ceil(need)) : 1;
    let clamped = 0;
    for (let i = 0; i < N; i++) {
      if (this.pinned[i]) continue;
      const v2 = vel[3 * i] ** 2 + vel[3 * i + 1] ** 2 + vel[3 * i + 2] ** 2;
      if (!(v2 < VMAX * VMAX)) { // also catches NaN
        clamped++;
        if (!Number.isFinite(v2)) { vel[3 * i] = vel[3 * i + 1] = vel[3 * i + 2] = 0; for (let d = 0; d < 3; d++) if (!Number.isFinite(pos[3 * i + d])) pos[3 * i + d] = this.prev[3 * i + d]; }
        else { const s = VMAX / Math.sqrt(v2); vel[3 * i] *= s; vel[3 * i + 1] *= s; vel[3 * i + 2] *= s; }
      }
    }
    this.clamped += clamped;
    if (this.thermostat) this._csvr();
    this.time += dt; this.stepCount++;
    // pressure from wall forces, exponentially averaged over ~1 ps
    const Pinst = this.wallArea > 0 ? this.wallForce / this.wallArea * BAR : 0;
    this.pressureBar = Pinst;
    this.pressureEMA += (Pinst - this.pressureEMA) * 0.001;
  }
  _integrate(nsub, forceRebuild) {
    const N = this.N, pos = this.pos, vel = this.vel, F = this.frc, h = this.dt / nsub;
    for (let sub = 0; sub < nsub; sub++) {
      for (let i = 0; i < N; i++) {
        if (this.pinned[i]) { vel[3 * i] = vel[3 * i + 1] = vel[3 * i + 2] = 0; continue; }
        const hk = 0.5 * h * ACC / this.mass[i];
        for (let d = 3 * i; d < 3 * i + 3; d++) { vel[d] += hk * F[d]; pos[d] += vel[d] * h; }
      }
      if (sub === 0 && forceRebuild) this.needRebuild = true;
      this._checkRebuild();
      this._rate = h;
      this.computeForces();
      for (let i = 0; i < N; i++) {
        if (this.pinned[i]) continue;
        const hk = 0.5 * h * ACC / this.mass[i];
        for (let d = 3 * i; d < 3 * i + 3; d++) vel[d] += hk * F[d];
      }
    }
    this._rate = 1;
  }
  _copy(dst, src, n) { if (!dst || dst.length < n) dst = new Float64Array(Math.max(n, 64)); dst.set(src.subarray(0, n)); return dst; }
  _save() { // lightweight state for redoing one step
    const n3 = 3 * this.N, P = this.nPairs, sv = this._sv || (this._sv = {});
    sv.pos = this._copy(sv.pos, this.pos, n3); sv.vel = this._copy(sv.vel, this.vel, n3); sv.frc = this._copy(sv.frc, this.frc, n3); sv.built = this._copy(sv.built, this.built, n3);
    sv.pN = this._copy(sv.pN, this.pN, P);
    if (!sv.pI || sv.pI.length < P) { sv.pI = new Int32Array(Math.max(P, 64)); sv.pJ = new Int32Array(Math.max(P, 64)); }
    sv.pI.set(this.pI.subarray(0, P)); sv.pJ.set(this.pJ.subarray(0, P));
    sv.P = P; sv.Epot = this.Epot; sv.map = this.pairMap; sv.needRebuild = this.needRebuild;
    this.Fold = this._copy(this.Fold, this.frc, n3);
  }
  _load() {
    const sv = this._sv, n3 = 3 * this.N, P = sv.P;
    this.pos.set(sv.pos.subarray(0, n3)); this.vel.set(sv.vel.subarray(0, n3)); this.frc.set(sv.frc.subarray(0, n3)); this.built.set(sv.built.subarray(0, n3));
    if (P > this.pairCap) this._allocPairs(P);
    this.pI.set(sv.pI.subarray(0, P)); this.pJ.set(sv.pJ.subarray(0, P)); this.pN.set(sv.pN.subarray(0, P)); this.nPairs = P;
    this.pairMap = sv.map; this.Epot = sv.Epot; this.needRebuild = sv.needRebuild;
  }
  kinetic() {
    let K = 0; const v = this.vel;
    for (let i = 0; i < this.N; i++) { if (this.pinned[i]) continue; K += this.mass[i] * (v[3 * i] ** 2 + v[3 * i + 1] ** 2 + v[3 * i + 2] ** 2); }
    return 0.5 * KEU * K;
  }
  dof() { let n = 0; for (let i = 0; i < this.N; i++) if (!this.pinned[i]) n += 3; return n; }
  temperature() { const nf = this.dof(); return nf ? 2 * this.kinetic() / (nf * KB) : 0; }
  _csvr() {
    const Nf = this.dof(); if (!Nf) return;
    const K = this.kinetic(), Kt = 0.5 * Nf * KB * Math.max(0, this.T);
    const c = Math.exp(-this.dt / Math.max(this.dt, this.tau));
    if (K <= 1e-12) { if (Kt > 0) this.thermalize(this.T * 0.05); return; }
    const r1 = this.gauss(), sum = Nf > 1 ? 2 * this.gamma((Nf - 1) / 2) : 0;
    let Kn = K + (1 - c) * (Kt * (r1 * r1 + sum) / Nf - K) + 2 * r1 * Math.sqrt(c * (1 - c) * Kt * K / Nf);
    if (Kn < 0) Kn = 0;
    const s = Math.sqrt(Kn / K), v = this.vel;
    for (let k = 0; k < 3 * this.N; k++) v[k] *= s;
  }
  thermalize(T, list) {
    const idx = list || Array.from({ length: this.N }, (_, i) => i);
    for (const i of idx) {
      const s = Math.sqrt(KB * T / (this.mass[i] * KEU));
      this.vel[3 * i] = s * this.gauss(); this.vel[3 * i + 1] = s * this.gauss(); this.vel[3 * i + 2] = s * this.gauss();
    }
  }
  scaleVelocities(list, f) { for (const i of list) { this.vel[3 * i] *= f; this.vel[3 * i + 1] *= f; this.vel[3 * i + 2] *= f; } }
  zeroMomentum(list) {
    const idx = list || Array.from({ length: this.N }, (_, i) => i);
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
    this._checkRebuild(); this.computeForces();
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
      this._checkRebuild(); this.computeForces();
    }
    v.fill(0, 0, 3 * N);
    return maxIter;
  }

  /* ---------- analysis ---------- */
  /* Bonds as {i, j, order, strength}: strength = f·b in 0…1, order = continuous bond order. */
  bonds(threshold = 0.3) {
    const out = [];
    for (let p = 0; p < this.nPairs; p++) {
      const f = this.pF[p]; if (f <= 0) continue;
      const s = f * this.pB[p]; if (s < threshold) continue;
      out.push({ i: this.pI[p], j: this.pJ[p], order: this.pN[p], strength: s });
    }
    return out;
  }
  /* Connected fragments using bonds with strength > 0.5. Returns {comp: Int32Array, list: [[indices]]}. */
  fragments(threshold = 0.5) {
    const N = this.N, parent = new Int32Array(N);
    for (let i = 0; i < N; i++) parent[i] = i;
    const find = i => { while (parent[i] !== i) { parent[i] = parent[parent[i]]; i = parent[i]; } return i; };
    for (let p = 0; p < this.nPairs; p++) {
      if (this.pF[p] * this.pB[p] > threshold) { const a = find(this.pI[p]), b = find(this.pJ[p]); if (a !== b) parent[a] = b; }
    }
    const comp = new Int32Array(N), groups = new Map();
    for (let i = 0; i < N; i++) { const r = find(i); if (!groups.has(r)) groups.set(r, []); groups.get(r).push(i); }
    const list = [...groups.values()];
    list.forEach((g, k) => { for (const i of g) comp[i] = k; });
    return { comp, list };
  }
  formulaOf(indices) {
    const cnt = {}; let charge = 0;
    for (const i of indices) { const s = ELEMENTS[this.type[i]].sym; cnt[s] = (cnt[s] || 0) + 1; charge += this.formal[i]; }
    const keys = Object.keys(cnt);
    let order;
    if (cnt.C) order = ['C', 'H', ...keys.filter(k => k !== 'C' && k !== 'H').sort()];
    else order = keys.sort();
    let s = order.filter(k => cnt[k]).map(k => k + (cnt[k] > 1 ? cnt[k] : '')).join('');
    if (charge) s += (Math.abs(charge) > 1 ? Math.abs(Math.round(charge)) : '') + (charge > 0 ? '+' : '−');
    return s;
  }
  /* A fragment is a radical if any atom has unused valence: valence − Σ(bond order) ≥ 0.5. */
  isRadical(indices) {
    if (!this._boSum || this._boSum.length < this.N) this._boSum = new Float64Array(this.cap);
    if (this._boStamp !== this.stepCount + ':' + this.N + ':' + this.Epot) {
      const s = this._boSum; s.fill(0, 0, this.N);
      for (let p = 0; p < this.nPairs; p++) { const f = this.pF[p]; if (f > 0) { const v = f * this.pN[p]; s[this.pI[p]] += v; s[this.pJ[p]] += v; } }
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
      pos: this.pos.slice(0, n3), vel: this.vel.slice(0, n3), frc: this.frc.slice(0, n3),
      type: this.type.slice(0, N), formal: this.formal.slice(0, N), val: this.val.slice(0, N), lp: this.lp.slice(0, N),
      pinned: this.pinned.slice(0, N), ids: this.ids.slice(0, N), cos0: this.cos0.slice(0, N), bo: Float64Array.from(bo)
    };
  }
  restore(s) {
    if (s.N > this.cap) this._alloc(Math.max(s.N, 2 * this.cap));
    this.N = s.N; this.time = s.time; this.stepCount = s.stepCount; this.rngState = s.rng; this.nextId = s.nextId;
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
    this._rebuild(); this.computeForces(); this.frc.set(f); this.pN.set(pn.subarray(0, Math.min(pn.length, this.nPairs))); this.cos0.set(c0);
    this.needRebuild = true;
  }
  _checkpoint() {
    const s = this.snapshot();
    const bytes = s.pos.byteLength * 3 + s.N * 40 + s.bo.byteLength + 200;
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
    return { format: 'chem-playground/scene@1', box: this.box, T: this.T, tau: this.tau, thermostat: this.thermostat, time: this.time, atoms };
  }
}

return { Engine, ELEMENTS, BY_SYM, PAIR, TUNE, KB, KEU, BAR, valenceFor, lonePairs };
});
