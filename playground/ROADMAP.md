# Roadmap: a playground where textbook chemistry happens by itself

The goal: put molecules in a chamber, set temperature and pressure, and see what real chemistry would do,
at a pace a person can follow. Fast steps are slowed down so the mechanism can be watched; slow ones are
skipped with an honest clock saying how long they really took. It should be good enough to learn textbook
organic chemistry from, and loose enough to experiment in.

## Where it stands

- **Watching.** A reactive molecular-dynamics engine (bond-order Morse force field, continuous bond
  orders, energy bookkeeping) runs in the browser at about a picosecond of chemistry per second. Radical
  reactions, combustion, chlorination and polymerisation happen on their own when conditions allow.
- **Light.** A UV flash (L) splits halogen molecules the way a lamp does, so a chlorination starts
  the way it does on a bench: nothing in the dark, a chain after the flash.
- **Forecasting.** The hourglass tool measures the barrier of any reaction you point at on the
  engine's own energy surface, in a worker, and turns it into a half-life with a range and a temperature
  table (transition-state theory).
- **Skipping.** The clock can jump by a waiting time drawn from that rate; the reaction then plays out
  live from just past its barrier.
- **What happens next.** The chamber can be searched for every reaction it could undergo, ranked by
  chance, and skipped to the next one (kinetic Monte Carlo), once or over and over on its own. Each
  reaction is tried on a copy before it is played, so what the chamber shows is what was forecast.
  Radical chlorination of methane runs as a chain this way until the chlorine is used up: 20 steps,
  24 µs of chemistry, three minutes.

## What still has to happen, most important first

### 1. Barriers good to about 5 kJ/mol

First target: a hydrogen shared between Cl (or O) and C is held by about 60 kJ/mol of Coulomb
attraction at 1.8 Å, because the charges follow the geometric bond switch rather than how much bond has
actually formed. That makes the Cl···H···C complexes, the high Cl barriers and the stray CH₃Cl₂ /
C₂H₄Cl₃ species in long chains. Fixing it means tying charge transfer to the formed bond and refitting
every polar bond energy against the benchmark.

Every time on the clock depends exponentially on a barrier: 6 kJ/mol is a factor of 10 at room
temperature. The forecast tool now shows where the force field is wrong:

| | model | real |
| --- | --- | --- |
| ethane C–H | 382 | 423 |
| formaldehyde C–H (formyl radical) | 493 | 369 |
| OH + H₂ barrier | 36 | 15 |
| H + HCl barrier | 33 | 15 |
| Cl + CH₄ barrier | 31 | 11 |
| Cl + H₂ barrier | 37 | 19 |
| H + C₂H₆ barrier | 12 | 38 |
| CH₃ + C₂H₄ barrier | 9 | 31 |
| Cl + C₂H₄ barrier | 11 | 0 |
| O₃, NO bond energies | 32%, 12% low | |

Fixed on the way: vinyl C–H (412 → 462, real 465), acetylene C–H (421 → 556, real 556) and phenyl C–H
(351 → 472, real 473), by giving σ radicals and bent π bonds the cost they have in reality. The barrier
scans themselves were also corrected: two-way scans (OH + CH₄ was 96, now 17, real 15), and the pass
back from the products may no longer swap to an equivalent product.

What to do:
- `tests/benchmark.cjs` is the benchmark: 39 reactions and bonds (barriers rms 12.8 kJ/mol), about a
  minute. Grow it to 60 or more (NIST kinetics database, standard compilations).
- Fit the force field against that table automatically, rather than tuning one reaction at a time.
  The fit must also score the wells along each path, not only the barrier: a chlorine-only crossing
  stabilisation fixed all four chlorine barriers and every bond energy, but deepened the Cl···H···Cl
  and Cl···CH₄ complexes by 25–50 kJ/mol (README, What was tried).
- Where the chlorine barrier and its complex come from, measured on Cl + CH₄ with the Cl–H distance
  held: the C–H holds on unchanged to 2.0 Å (+48 kJ/mol), then flips within 0.2 Å to a shared
  Cl···H···C state with bond strengths 0.43 / 0.48 and charges Cl −0.26, H +0.33, at −12 kJ/mol. With
  the charge transfer switched off the same state is +35. So the barrier is the abrupt switch between
  two bonding patterns, and the complex is the shared state over-stabilised by charge transfer. The lead
  is to let a shared hydrogen pass its bond over gradually and to limit the charge a hydrogen held by
  two partners can take, then refit, rather than adding stabilisation.
- The errors fall into three families: hydrogen passing to or from O or Cl is 15–23 kJ/mol too high
  (the half-made crossings cost too much); every C–H next to another carbon is too weak (ethyl radical
  too stable), which makes H + C₂H₆ and CH₃ + C₂H₄ too easy; and the formyl radical lacks the C–O
  strengthening that makes aldehyde C–H weak (oxygen's valence caps C–O at a double bond).

### 2. The live simulation and the forecast must agree

They agree for methane at 2000 K, and the live chamber's barrier for H + H₂ matches the forecast. Two
disagreements remain:
- The live H + H₂ exchange runs about 40 times faster than both the forecast and the real reaction,
  because the model's hydrogen captures H₂ from about twice the real distance. This is a too-loose
  transition state, fixable in the force field once item 1's benchmark exists.
- Above about 2200 K ethane gives off H₂ within picoseconds in the live chamber, while its barrier
  (408–435 kJ/mol) says microseconds or longer. Not yet explained. The leading suspect is that the bond
  orders are a state the dynamics carries, not a function of the geometry, so the live system can reach
  states the static energy surface does not contain. The structural fix is to make bond orders an
  explicit, differentiable function of the atom positions (as ReaxFF does), with forces that include
  their gradient. That is a large change to the engine and would need the whole test suite and a refit.

### 3. Solution chemistry: the biggest missing piece for organic textbooks

Most textbook mechanisms (SN1, SN2, E1, E2, aldol, esterification, acid and base catalysis) involve ions
in a solvent. The engine has neutral atoms in a gas. Needed:
- **Charges that move.** Charge equilibration (QEq / electronegativity equalisation) so a leaving
  group can take its electrons with it and a carbocation can exist.
- **Proton transfer.** A proton must be able to hop from acid to base; today hydrogen only moves as an
  atom with its electron.
- **Solvent.** An implicit continuum (Born / generalised Born) to stabilise ions cheaply, and a way to
  add a shell of explicit water where a mechanism needs it (proton relays).
- **Calibration** against pKa values and the classic SN2 rate series (methyl > primary > secondary).

Without this, the forecast and fast-forward will work well for radical and gas-phase chemistry
(halogenation, combustion, cracking, radical polymerisation) and not at all for ionic chemistry.

### 4. A complete kinetics engine

- Running on by itself takes about 9 s a step for a handful of molecules (barriers it has measured
  are remembered). The hydrogen-passing complexes (item 1) are what stop long runs today.
- Search reactions between two closed-shell molecules (Diels–Alder, ene, [2+2]) and rearrangements
  inside one molecule.
- Replace the fixed prefactors with ones calculated from vibrational frequencies at the barrier and in
  the reactants (harmonic transition-state theory).
- Pressure fall-off for unimolecular and association steps (needed before hydrogen combustion can be
  forecast: H· + O₂ → HO₂ only wins at high pressure), tunnelling and zero-point energy for hydrogen
  transfers.
- A penalty for concerted additions that orbital symmetry forbids (two closed-shell molecules adding
  in one step, such as ethene + methane through a four-centre ring; the search skips these below 500 K).
- Large chambers: forecast once per distinct situation and reuse it, instead of rescanning every copy.

### 5. The teaching layer

- Curved-arrow mechanisms drawn from the bond orders that changed during a reaction.
- Energy diagrams for a whole sequence of steps, not only one.
- Ready-made textbook scenarios: four radical ones exist (Experiments tab: chlorination of methane,
  chlorine and ethene, methyl recombination, hydroxyl and methane). Next: combustion of hydrogen,
  radical polymerisation of ethene; SN2 and E2 once item 3 exists.

## Suggested order

1 and 2 first: without trustworthy barriers, nothing else on the clock can be trusted. Then 4 (running
on by itself), because that turns the playground into the "hours in minutes" experience for everything
radical. Then 3, which is the largest piece of work and what opens up most of an organic chemistry
course. The teaching layer can grow alongside.
