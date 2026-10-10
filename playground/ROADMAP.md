# Roadmap: a playground where textbook chemistry happens by itself

The goal: put molecules in a chamber, set temperature and pressure, and see what real chemistry would do,
at a pace a person can follow. Fast steps are slowed down so the mechanism can be watched; slow ones are
skipped with an honest clock saying how long they really took. It should be good enough to learn textbook
organic chemistry from, and loose enough to experiment in.

## Where it stands

- **Watching.** A reactive molecular-dynamics engine (bond-order Morse force field, continuous bond
  orders, energy bookkeeping) runs in the browser at about a picosecond of chemistry per second. Radical
  reactions, combustion, chlorination and polymerisation happen on their own when conditions allow.
- **Light.** A lamp (L) is a condition that stays on: Cl₂, Br₂, I₂ and F₂ split at their measured
  rates, so a chlorination starts the way it does on a bench: nothing in the dark, a chain under the
  lamp. While it is on, radicals meet at the rate they would in a lit flask, so the chains run their
  real length. Methane, ethene and hydrogen chlorination all run at room temperature.
- **Forecasting.** The hourglass tool measures the barrier of any reaction you point at on the
  engine's own energy surface, in a worker, and turns it into a half-life with a range and a temperature
  table (transition-state theory).
- **Skipping.** The clock can jump by a waiting time drawn from that rate; the reaction then plays out
  live from just past its barrier.
- **What happens next.** The chamber can be searched for every reaction it could undergo, ranked by
  chance, and skipped to the next one (kinetic Monte Carlo), once or over and over on its own. Each
  reaction is tried on a copy before it is played, so what the chamber shows is what was forecast.
  Under the lamp at room temperature methane goes CH₄ → CH₃Cl → CH₂Cl₂, ethene to
  1,2-dichloroethane (or 1,2-dibromoethane with Br₂), and H₂ + Cl₂ to HCl: one photon starts a chain
  that runs until the chlorine or the ethene is used up. In the dark, ethene and Cl₂ or Br₂ still add on
  the glass, slowly, as in a real flask. The search for the next step starts in the background while
  the last one plays, so most skips take 0.2–8 s (an occasional one still takes 15–30 s). After a
  skip, molecules left stuck to each other drift apart, as hours or nanoseconds of gas motion would make
  them; with chlorine in excess, the leftover Cl· goes on to 1,1,2-trichloroethane, as it does under a
  real lamp. Methane and O₂ at 1200 K start with O₂ + CH₄ → CH₃· + HO₂·. Once radicals are made, they
  meet at a real flask's rate, so chains are not cut short by the chamber's crowding. Running on stops
  by itself at an equilibrium or a loop and names it.
- **Starting.** The chamber starts empty, with no text: a tap places a hydrogen, + opens the atoms and a
  built-in set of molecules. On a phone everything works by touch: hold an atom to inspect it, pinch to
  zoom, drag to pan or move.

## What still has to happen, most important first

### 1. Barriers good to about 5 kJ/mol

Why single knobs fail here, measured on Cl + CH₄: the model's path is reactants (0) → an abrupt switch
of bonding pattern (+31 kJ/mol, the barrier) → a half-shared Cl···H···C complex (−20) → products (+8),
where the real path has the half-shared state as its top, at about +11. So the fix has to lower the
switch and raise the shared state at the same time. Four levers tried on 4 October each moved both the
same way: chlorine crossing stabilisation, chlorine Pauli softening, removing a bond's own charge from
its Coulomb term, and the chlorine–hydrogen saturation range (`kb` for Cl–H alone: 1.0 gives barrier
46 / well −15, 1.7 gives 8 / −55). The likely structural fix is a saturation function for hydrogen
that lets it share its bond gradually (no abrupt switch) together with a cost on the shared state
itself, fitted against the benchmark's barrier and well columns at once.

What it costs today: the live simulation still has the wrong barriers. The skip-ahead corrects them
with measured values per reaction class (README, Measured barriers), which is a patch over
the force field, not a fix of it: a fixed force field would let the corrections go.

Started (October 9, evening): the bond charges of H–Cl, H–Br and H–I are now scaled by 0.7, 0.6 and
0.5. The one electronegativity rule had made HCl twice as polar as it is (2.2 D against 1.08), and the
too-large charge on its hydrogen was most of the CH₃···HCl and Cl···H···C complex wells. Those wells
are now about half as deep (Cl + CH₄: 20 → 11 kJ/mol, Cl + C₂H₆: 24 → 13) for a 2 kJ/mol higher
Cl + CH₄ barrier; going all the way to the measured dipole (0.5) leaves 7 kJ/mol but raises that
barrier to 41, because the barriers were tuned with the overcharged hydrogen. The rest of this item
is still to do.

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
| CH₃OO–H (peroxide O–H) | ≈ 418 | 357 |
| CH₃O₂· + CH₄ → CH₃OOH + CH₃· | ΔE +21, barrier 23 | ΔE +82 |

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

- Running on by itself takes a few seconds a step for a handful of molecules on an ordinary
  multi-core machine (barriers it has measured are remembered, and hopeless routes are dropped early),
  and up to a minute once a polymer has grown past C₉. It stops by itself at an equilibrium or a loop.
  The hydrogen-passing complexes (item 1) are what still clutter long runs: CH₃···HCl at room
  temperature, HCl···Cl, NH₂···HCl.
  A freshly made CH₃O₂· sometimes falls apart again in the live simulation within a picosecond of
  its skip (the CH₃–OO bond is too weak while the product is still hot), and is remade on the next one.
- A ground-state O atom's abstraction from methane is placed with the oxygen swung round onto the
  carbon (the scan's path bends that way on this spinless surface), so the CH₃· and ·OH pair up at once
  and some trials fail; holding the path straight made the barrier worse (README, What was tried).
- Search reactions between two closed-shell molecules (Diels–Alder, ene, [2+2]) and rearrangements
  inside one molecule.
- Replace the fixed prefactors with ones calculated from vibrational frequencies at the barrier and in
  the reactants (harmonic transition-state theory).
- Pressure fall-off for unimolecular and association steps (needed before hydrogen combustion can be
  forecast: H· + O₂ → HO₂ only wins at high pressure), tunnelling and zero-point energy for hydrogen
  transfers.
- A penalty in the force field for concerted additions that orbital symmetry forbids (two closed-shell
  molecules adding in one step through a four-centre ring), and for a ground-state (triplet) O or S
  atom inserting into a σ bond. The search already refuses to count both; the live simulation can
  still do them, because the model carries no spin.
- Skip each scan whose best possible rate cannot matter. In ethene + Cl₂ under the lamp, the state
  C₂H₄Cl· + C₂H₄Cl₂ + Cl· costs ten quick scans of 1.6–5 s each (about 7 s of waiting on four cores),
  though Cl· taking Cl from C₂H₄Cl₂ is about 100 kJ/mol uphill and could never compete with the two
  radicals pairing up. A bound from bond energies before scanning (rate ≤ collision rate × e^(−ΔH/RT))
  would drop those; the catch is that the engine's own Cl bond energies are not yet trustworthy enough
  to prune on.
- Free-form mixtures from the molecules list that still misbehave: benzene + Cl₂ under the lamp (the Cl·
  adducts of benzene come apart again in the live simulation, so the chain wanders instead of heading for
  hexachlorocyclohexane); H₂ at 1000 K (each H₂ → 2 H· skip is undone at once, because the two atoms meet
  again in a box this small and the live simulation joins them without a third body); H₂ + O₂ at 1200 K
  (the radicals shuttle H· + H₂O ⇌ H₂ + ·OH, since in a lit or burning flask they meet each other a
  million times less often than in the chamber, so water forms slowly); NH₃ + HCl (no ionic pairing, so
  no ammonium chloride); a spark in H₂ + O₂ at room temperature (the Kelvin bath touching every atom takes
  the spark's heat away within a few hundred femtoseconds, so the chain never branches; a real flame heats
  itself. A spark that heats a small region and holds the bath back there for a while would let a flame
  front run).
- Ethene + O₂ when hot has no first step: the scan of a vinyl hydrogen meeting O₂ ends in a concerted
  C₂H₃O· + ·OH (256 kJ/mol) instead of the clean transfer to C₂H₃· + HO₂·, and two closed-shell molecules
  may only meet as a clean abstraction, so the pair is dropped and the forecast falls back on splitting O₂
  or ethene, 10⁴ years at 900 K where the real first step takes hours. The transfer channel needs to be
  held to its products along the path itself (the C₂H₃ end already bonds to the far oxygen while the
  hydrogen moves, so relaxing the end state is not enough), and given the measured 237 kJ/mol.
- Live speed. A step of ethene + Cl₂ (16 atoms) costs about 230 µs: two force evaluations (every O–H,
  H–H and C–H stretch asks for two sub-steps) of about 110 µs each, so the chamber plays about 4 ps a
  second on a desktop and 0.2–2 ps on a phone, against the 20 ps a second that 1× stands for. The time
  goes evenly into the many passes of the force field (angle screening, insertion screening, bond orders,
  pair energies); no single loop dominates. Worth trying: one pass over the pairs that computes what
  several passes now recompute, and holding the angle and insertion screens for a few sub-steps.
  Moving the helpers' return values from closure variables onto one object took 10–18 % off a step
  (V8 boxed every double written to a closure slot: about 20 kB of garbage a step) with bit-identical
  trajectories. A step still allocates about 75 kB, most of it inside computeForces and the closures
  that _torsions, _angScreen and _updateBondOrders create on every call; those are the next place to look.
- Large chambers: forecast once per distinct situation and reuse it, instead of rescanning every copy.
- Large molecules: every scan point relaxes every atom of both partners, so a C₉ radical meeting ethene
  costs 25–30 s a candidate. Scanning only the atoms within two bonds of the reacting pair, capped with
  hydrogens, would make it independent of chain length; the catch is placing the result back, since
  the rest of the chain must follow the moved region without straining the cut bonds.

### 5. The teaching layer

- Curved-arrow mechanisms drawn from the bond orders that changed during a reaction.
- Energy diagrams for a whole sequence of steps: the last fourteen skips are drawn as one diagram in
  the chamber's corner. Next: label the species at each level, and include steps that happened live.
- Ready-made textbook scenarios: six radical ones are verified (built by hand from the molecules list: chlorination of
  methane, chlorine and ethene, radical polymerisation of ethene, methyl recombination, hydroxyl with
  methane in air, hydrogen and chlorine). Next: combustion of hydrogen (at the chamber's 100 bar it runs
  through HO₂· and H₂O₂, the high-pressure route, with some exotic H₂O₃ on the way); SN2 and E2 once
  item 3 exists.

## Suggested order

1 and 2 first: without trustworthy barriers, nothing else on the clock can be trusted. Then 4 (running
on by itself), because that turns the playground into the "hours in minutes" experience for everything
radical. Then 3, which is the largest piece of work and what opens up most of an organic chemistry
course. The teaching layer can grow alongside.
