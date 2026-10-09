# Chem Playground

A reactive molecular-dynamics sandbox at true scale. Atoms and molecules live in a 3D slab, 1.2 nm deep by default, which you view from above, so reactions cannot drift away along z. The manual view zooms from 0.05 nm to 15 nm. Chamber dimensions extend to 50 nm; Fit can zoom farther out to show the complete box. In the sandbox, bonds form and break from an experimental potential energy model.

Open `playground/` on the site. An empty chamber shows the ready-made experiments as tiles, each
with a picture of what goes in; one click loads it, and the button to press next glows (☀ for light,
then the double arrow for the next reaction). Single atoms come from the dock, molecules from **+**
or from [Nomenclature](../nomenclature.html) → **Playground**.

## Thermostats

Three ways to hold a temperature, chosen in **Environment** or from the temperature gauge. The
difference between them is who gets touched and when.

### Building a molecule

The Kelvin stat is the default, because it is the one that lets a molecule be built. A bond
forming in open space releases its binding energy on the spot; a boundary heater has no way to
take that away from the middle of the chamber, so the molecule you just made blows itself apart.
That is correct physics — recombination genuinely needs a third body to carry the energy off —
and useless as a default. The Kelvin stat is that third body, everywhere.

Two things follow, and both are worth knowing before building anything:

- **Add atoms one at a time, keeping like with like apart.** Four loose hydrogens near each other
  pair into H₂ long before they reach a carbon, and H₂ will not then add to it. Starting from CH₂
  and dragging in one lone hydrogen at a time gives CH₃ and then CH₄; starting from a carbon and
  four adjacent hydrogens gives CH₂ + H₂ and stops there.
- **A mixture that sits still is usually right.** Reactants are mostly metastable. Use the Spark.

### Kelvin stat

Every unpinned atom is coupled to a bath at the target temperature by an exact Ornstein–Uhlenbeck (Langevin) velocity update: `v ← c·v + σ·ξ`, with `c = e^(−dt/τ)` and fluctuation–dissipation noise `σ² = (1 − c²)·kT/m`. Each atom relaxes toward the bath on its own, so a molecule that has just formed sheds its new bond's energy without any other atom being touched, and the sample samples the canonical ensemble with its natural kinetic-energy fluctuations.

The coupling time is `tau` (default **250 fs**) for hydrogen and scales with mass, `τ_atom = tau · m / 1.008`, the way friction from a light buffer gas does — about 3 ps for carbon, 9 ps for chlorine. The temperature is therefore **held on average, not exactly**: a reaction that releases heat faster than this drains it warms the sample until the bath catches up, and a runaway chain can outpace it entirely. Lower `tau` to clamp harder, at the cost of more strongly damped dynamics.

Setting a new temperature rescales velocities once, immediately (`_kelvinSet`); startup from rest initializes velocities. All energy the bath adds or removes is recorded in `kelvinWork`. Pinned atoms stay fixed, and work accounting survives replay. See [Langevin thermostats](https://docs.lammps.org/fix_langevin.html).

Only the pointer's own contribution is excluded from all of this, and that is the dragged
cluster's centre-of-mass velocity — the coherent part the servo imposes. A dragged molecule still
vibrates, still rotates, and still has to give up the energy a new bond releases. Exempting those
atoms wholesale meant two hydrogens dragged together formed H₂, kept the 436 kJ/mol that released,
and tore straight back apart.

### Wall heater

A **spatial wall thermostat**, not a global velocity rescale:

- Heater targets are **15–350 °C**. Changing the target does not jump the sample or wall temperature.
- The wall follows `dT/dt = (T_target - T_wall) / tau`; its response time is adjustable in **simulated ps**, default 10 ps. This specifies a nanoscale thermal reservoir, not the heating time of a particular laboratory appliance.
- Only atoms within **0.2 nm of a wall** couple to the reservoir. The coupling vanishes quadratically at the inner edge. An exact Ornstein–Uhlenbeck velocity update combines drag and Gaussian noise with fluctuation–dissipation variance. Interparticle forces carry energy through the interior.
- Effective wall heat capacity is 1000 k_B. Heat given to the sample is subtracted from the wall; heater work and sample heat are recorded separately and included in deterministic snapshots. This finite reservoir approximation is not an explicit solid-wall atomistic model.
- Against a void wall the heater cannot reach its setpoint: it only touches the boundary layer, which is the same place the void takes from. At 400 K with a temperature void the sample settles near 310 K and keeps wandering. Use the Kelvin stat when the number matters.

### What counts as temperature

An atom held by the pointer is being driven from outside, so its motion is not thermal. It is
left out of `dof()`, `temperature()` and `thermalKinetic()`, out of every thermostat's rescale,
out of the void-temperature cap and out of the boundary bath, and it is never itself rescaled —
the pointer owns it while you hold it, and it counts again the moment you let go.

Dragging is a servo on the whole cluster, not a spring on one atom: the pointer drives the
dragged molecule's centre of mass with damping and a speed limit, distributed mass-weighted so
every atom takes the same acceleration and the molecule feels no internal stress. A stiff spring
on a single atom reached its force cap as soon as the pointer was a few ångström ahead, tore the
bond, and flung the pieces at kilometres per second — which then counted as the sample's heat.

Counting it was why dragging one molecule stopped every other one. A stat that holds the total
kinetic energy saw the drag as an enormous excess and scaled the whole chamber down to
compensate, while the dragged atom, re-accelerated by the tweezer every step, was the only thing
still moving. Measured over 3 ps of dragging with twenty argon atoms, the kinetic energy of the
rest of the chamber went 90.2 → 0.0 under the Kelvin stat and 78.7 → 0.0 under a temperature
void; it now goes 90.2 → 94.8 and 78.7 → 94.7, with the gauge reading the fluid rather than the
drag. With no thermostat the chamber still heats as you stir it, which is the work you are doing.

### Off

Editing temperature immediately rescales kinetic energy to the requested temperature, including
an initialization from rest. Pinned atoms remain fixed. Subsequent dynamics have no thermostat.

Isolated molecule conditioning retains a separate CSVR sampling bath.

[Spatial Langevin thermostat methodology](https://docs.lammps.org/fix_langevin.html). This engine uses its own exact OU update, not the LAMMPS integration algorithm.

## Chamber boundaries

The **Conditions** panel sets the chamber geometry and boundary.

- **Solid** reflects atomic centres at all six faces during the integration drift. Normal velocity reverses; tangential velocity and kinetic energy are preserved. Corner and multiple-face crossings are handled. Contour tails may extend beyond the face because the contours are schematic atomic envelopes.
- **Forcefield** uses a conservative harmonic potential, 5 kJ/mol/Å², starting 1.2 Å inside the face. Atoms can penetrate this soft boundary before turning back.
### Void walls

A reflecting chamber hands back everything it receives, so energy a sample releases as it
settles comes straight off the walls again. A void wall makes the chamber open in one chosen
respect. None of them can add energy or reverse a velocity; all are deliberate sinks, not
models of vacuum.

- **Temperature** radiates heat away the instant it appears, everywhere — not only where an atom
  happens to touch a wall. A sample that heats itself, through friction, a reaction or work done
  on it, never ends a step hotter than the temperature it was set to. One-sided: it only removes,
  so a cold chamber stays cold and nothing here can drive the sample. The whole sample is scaled
  at once, so no bond is ever pulled harder at one end than the other.
- **Velocity** stops an atom dead. Not damped, not reflected: every component is zeroed at the
  moment of contact, and the atom is then left alone — free to be moved by whatever else acts on
  it. The momentum it delivered is still wall stress and is still reported.
- **Pressure** voids the pressure, not the motion. Atoms arrive and rebound exactly as they do
  in a reflecting chamber — a voided chamber's trajectories are identical to a closed one's, atom
  for atom — but the wall records nothing, so the gauge reads zero with the chamber exactly the
  size it was. It removes no energy; it is a reading that is let go, not a sink.

Velocity acts **only on contact, and only on motion driving into a face the atom has actually
reached.** Measuring contact as a shell reaching inward from the face was wrong twice over: an
atom gliding past a wall it never touched was stopped, and an atom approaching one was halted
short of it, hanging in mid-air. It also clamped an atom inside the shell every step while its
bonded partners outside it kept moving, which stretches the bond and drags the molecule — a wall
that appeared to shove things for no reason. A molecule resting near a wall now behaves exactly
as it does in a reflecting chamber.

The stop happens at the end of a sub-step, after the integrator has finished with the atom, and the
atom is not moved: the face's own spring holds it, and gives back exactly what it stored. Stopping it
halfway through the step, or putting it back on the face, did work nobody recorded. A CH₄ thrown into
the wall at 3 km/s lost or gained 20–40 kJ/mol that `voidHeat` never saw; it now balances to 0.4. Only
an atom more than `WALL_FOLD` outside (left there by an edit or a resize) is still carried back.

**The wall is measured, not set.** `wallMeasured` is the kinetic temperature of the fluid lying
within `wallSkin` of a face — the wall is whatever the sample against it is, and a heater only
sets what it *aims* for. With a temperature void the heater stands down entirely, because it
cannot warm a wall that radiates everything away, and the wall simply reports the fluid; that is
what the gauge and the Environment panel show. To hold an open chamber at a fixed temperature,
use the **Kelvin stat**, which keeps replacing what the void removes.

Removed energy is recorded in `voidHeat`. Rejected integration attempts restore this ledger, so
an absorbed collision is counted once. Measured over 40,000 steps with hydrogen and oxygen let go
at 300 K: a closed chamber heats itself past 600 K, while the same chamber with a temperature
void never exceeds 300.00 K and radiates 26,000 kJ/mol away.

### Pressure

Wall pressure is normal momentum flux (solid collisions) or normal reaction force (soft fields),
divided by total wall area. That passive reading is always reported. **Hold** turns on a
Berendsen-style controller: every 100 steps the chamber breathes toward the target with
`μ = 1 + gain·(P − P_target)/span`, clamped to 0.4 % per adjustment, across width and height only
— the slab depth stays what you set. Unpinned fragments move with the walls without scaling their internal bond lengths. Pinned fragments remain fixed and limit compression. This heuristic controller is not a validated NPT ensemble or an explicit piston; its mechanical work is not included in the energy ledger.

`node playground/tests/thermal-walls.cjs` checks all three thermostats, including exactness,
replacement of void losses, pinned atoms and bit-exact replay.

`node playground/tests/bounds.cjs` checks six-face containment, energy conservation, collision
pressure, multiple crossings, all three void channels at both wall kinds, pressure control and
deterministic replay. A boundary you drag is a user edit, not a simulated piston.

## Interface and performance preferences

The compact settings panel includes persistent contour quality, grid, atom-label, interface-motion and compute-priority controls. Reduced motion follows the operating system by default. Display quality changes rendering resolution only; the physical step stays 1 fs. The actual simulation rate remains visible when CPU-bound.

The graphite interface is organized around one bottom control rail. It combines manipulation tools, the last-used element, an **Add** drawer for atoms and imported molecules, and playback. The clock opens sample measurements and composition; the scale bar fits the chamber. Target temperature remains directly editable, while the pressure reading opens chamber measurements and dimensions. No separate inventory chips, Data/Fit row, or molecule-tray button occupies the canvas.

Secondary actions stay near their controls: right-click Add for molecules, the current element to change it, temperature or pressure for Conditions, the scale for Display, and the clock for speed presets. Every secondary action also has a left-click or keyboard route. Text selection and dragging are disabled on interface chrome; value editors and the command input retain normal editing.

Settings tabs crossfade inside the same surface; dropdowns expand and collapse in place; panels retain their contents during exit. Motion uses a shared, non-bouncing curve with 340–600 ms transitions. Closing content becomes inert immediately. Reduced motion disables both CSS transitions and the new scripted transitions. Rendering and simulation timing are independent of interface animation.

## Time

- Every step is exactly **1 fs**: velocity Verlet in a render-decoupled accumulator (`while (acc >= 1) { step(); acc -= 1 }`), drawn with interpolation between the last two steps.
- **1.0× = 20 000 steps/s = 20 ps of simulated time per real second.** The scale is exponential below 1×: 0.1× is one step per second. Above 1× it is linear, e.g. 2× = 40 000 steps/s. Type any value into the speed field, including ones outside the slider range.
- When a frame cannot fit all the steps, the rate readout says **CPU-bound** and shows the actual speed.
- A step that runs into extreme curvature splits itself into up to 16 sub-steps. This happens when a hydrogen is squeezed between competing bonds, or in a hard collision at thousands of kelvin. Either the curvature from the previous step triggers it, or a per-step energy check does. The step still advances exactly 1 fs.
- **Step back** is exact. The engine keeps a checkpoint every 64 steps and replays deterministically, bit-identical, because the RNG state is part of the checkpoint. Hold the button to rewind.

## Igniting a reaction

A mixture of reactants usually sits there doing nothing, and that is correct rather than broken.
H₂ and O₂ are metastable: they do not react by direct collision at any temperature you can reach
in a browser, because the initiation step is expensive. Real combustion runs as a *chain*, and
the chain needs a radical to start it. On a bench that comes from a spark, a flame, a hot
surface or a photon.

The **Spark** tool is that starter. Click, and the atoms within 2.2 Å get a short outward burst
of directed motion — 25,000 K worth — which is enough to pull a bond apart and set two radicals
loose. The energy is real, is added to the ledger, and is not a scripted reaction: what happens
afterwards is ordinary dynamics.

Because the chamber holds only a handful of atoms, a thermostat here is instantaneous and global
and would erase a spark on the step it landed. So a spark opens a 2 ps **ignition window** during
which the stat and the temperature void stand back. That window is part of the saved state and
replays exactly. It is the one place where a temperature void can be briefly exceeded, and
deliberately so: a spark is hotter than its surroundings, which is the entire point of one.

Spark and Break bond start the chamber running if it was paused. Reactions appear in the
list at the bottom left only once the new molecules have lasted 100 fs: a freshly made bond vibrates
through the bonding threshold a few times, and listing every crossing used to show C₂H₄ + H· → C₂H₅·
followed at once by C₂H₅· → C₂H₄ + H·, which never happened. A cluster in which an atom holds more
bonds than its valence (a hydrogen half-shared between a carbon and a chlorine) is listed as the
encounter it is, split at its weakest bond, CH₃···HCl rather than "CH₄Cl".

Measured, three H₂ and two O₂ in a 1.6 nm chamber at 300 K, 150 ps:

| | left alone | after one spark |
|---|---|---|
| no thermostat | 3 H₂ + 2 O₂ | 2 H₂ + 2 HO + O₂, sample runs to 2300 K |
| Kelvin stat | 3 H₂ + 2 O₂ | **H₂ + 2 H₂O + O₂**, held back at 300 K |

That is the whole of hydrogen combustion: metastable until lit, chain-branching once started,
and the stat carrying the heat away afterwards. Adding a lone H atom from the dock does the same
job more gently — at 300 K it stalls at HO₂, the real chain-terminating step, and only branches
to water once the chamber is hot.

## Valence sharing

An atom cannot give a whole bond to two neighbours at once. The coordination each bond sees is
therefore scaled back where the raw claims exceed the valence:

    k_i = 1 − share·(1 − V_i / Z_i)     for Z_i > V_i,  k_i = 1 otherwise
    C_ij = (Z_i − f_ij)·k_i             the competing coordination bond ij is judged against

with `share = 0.5`. `k` is exactly 1 for every atom at or under its valence — every equilibrium
structure — so bond lengths, bond energies and thermochemistry are untouched by construction, and
only the half-made crossing of a reaction changes.

Without it, a hydrogen mid-handover contributed a full 1.000 to *both* its old and its new
partner, reaching coordination 2.0 on a valence of 1. Each bond then saw a full competitor, kept
0.120 of its attraction, and the total collapsed to **0.24 of one bond** — an artificial wall of
roughly 330 kJ/mol across every reaction where an atom changes partners. With sharing the total
through a CH₂ + H₂ crossing runs 0.98 → 0.84 → 0.99 → 0.81, which is bonding roughly conserved,
as it should be.

Measured against the fitted set, `share = 0.5` changes no barrier at all: H + H₂ 45, H + CH₄ 40,
F + H₂ 12, Cl + H₂ 50, H + Cl₂ −9, OH + H₂ 45, 2 H₂ + O₂ −486, deepest spurious complex −3 — all
identical to `share = 0`. H₂ + O₂ still sits unreacted through 100 ps at 300 K. Pushing it to
`share = 1` does start to bite: H + CH₄ falls to 32, F + H₂ to 3, the handover total overshoots
past one bond, and a 36 kJ/mol spurious H₃ complex appears — so half is where it stays.

Forces remain the exact gradient of the energy; `regression.cjs` checks analytic against finite
differences at `share` 0, 0.5 and 1 on a deliberately over-coordinated geometry.

**It does not make CH₂ + H₂ react.** That barrier's remainder is set by the saturation curve
itself — at a perfectly shared handover each bond still keeps only 0.758, about 105 kJ/mol short
of one bond — and that curve is the same number that gives H + H₂ its 40 kJ/mol and keeps H₂ and
O₂ apart. Flattening it for carbenes would flatten it for everything. Singlet CH₂ inserts because
of orbital structure this model has no way to express, the same absence that makes O₂ need a
hand-patch.

## Force field

A custom, experimental bond-order model inspired by reactive force-field ideas. It is not an implementation of Tersoff, Brenner, or ReaxFF, and its fitted examples do not establish general chemical accuracy.

| Term | Form | Data |
| --- | --- | --- |
| Covalent | Morse `De[e^(−2y) − 2·b·e^(−y)]`, with re, De and a interpolated in the continuous bond order n | Measured bond lengths, dissociation energies and stretch force constants (≈80 pairs); Pauling/Pyykkö estimates fill the rest |
| Valence used | Each bond uses `s = e^(−1.3·a·Δr)` of both partners' valence — it decays like the Morse attraction it still provides — times its multiplicity | — |
| Saturation | `b = p_i·p_j`, `p = 1/(1 + 0.6x + 10.8x⁶)`, x = excess valence if this bond were fully formed | Fitted to the barriers below while keeping closed-shell dimers non-sticky |
| Evans–Polanyi | Excess valence is weighted by `√(De_competing / De_new)`, so a stronger incoming bond displaces a weaker one more easily. An atom at its normal valence is unaffected | Fitted to F + H₂ and H + Cl₂ |
| Screening | Two atoms bonded to a common neighbour do not compete for each other's valence, unless they are bonded themselves | — |
| Bond order | Spare valence shared between neighbours → C=C, C≡C, O=O, N≡N, CO₂, aromatic 1.5. Each atom first offers its spare valence in proportion to what its partners can take, then re-divides it `pairIter` times toward partners that offer back (a damped Sinkhorn pairing), so a partner with nowhere else to go wins: butadiene alternates 1.84 / 1.16 (1.351 / 1.498 Å, exp. 1.338 / 1.467) and comes out 9.7 kJ/mol more stable s-trans, while symmetric systems — benzene, allyl, CO₂ — are a fixed point and stay exactly as they were. It follows the same long-ranged, screened bond order the saturation uses, so a partner weakens a π bond *while* its own bond forms; a partner claims valence only insofar as it can bond at all, so a radical frees the π bond and a saturated molecule drifting past does not | Radical addition to alkenes |
| Angular saturation | A saturated centre cannot take a partner side-on to one of its σ bonds: the Morse attraction of a newcomer n at centre c is screened by each arm c–k it approaches at less than ≈ 70°, weighted by how σ-like that arm is (geometric π measure), and the screen lifts for a radical newcomer, for ring partners, for a carbene (≥ 1.2 free valence counted with long-reach bonds) at any corner of the triangle, for metals and ionic pairs, and when neither c nor n has any free valence. Contacts from behind an arm (≥ 160°) are untouched, so collinear transfer barriers do not change | CH₅ and five-coordinate carbon wells removed; every barrier and atomization energy unchanged |
| π blocking | A π bond blocks an incoming partner at 0.85 of a σ bond, being weaker and more polarizable. O₂'s π counts 0.78 (triplet O₂ is a diradical) and an O–O bond with one unpaired oxygen gains 0.45 order (the three-electron bond of HO₂·) | Cl + ethene; O₂ and HO₂ chemistry |
| Angles | VSEPR, θ₀ a smooth function of the continuous steric number | 109.5°, 107°, 104.5°, 120°, 180° |
| σ over-coordination | Every attracting contact of a carbon, nitrogen or oxygen is counted as a σ partner by how much Morse attraction it still gives (weak contacts down to a fifth of a bond, 1-3 neighbours excluded, π order not counted). Holding more than valence + 0.3 costs 200·excess² kJ/mol. Stretched bonds count as little used valence, so a hot carbon could otherwise hold five bonds' worth of attraction | Equilibrium molecules and every fitted barrier unchanged; five-coordinate carbon probe +162 kJ/mol |
| σ radicals | A carbon with π bonds, at most two σ partners and valence left over after both (vinyl, ethynyl: the unpaired electron sits in an sp or sp² orbital, not in the π system) costs 50 kJ/mol (two σ partners) to 135 (one), switched smoothly on each of the three. Closed-shell molecules, CH₃, ethyl, allyl and carbenes are untouched | Vinyl C–H 412 → 462 kJ/mol (exp. 465), acetylene C–H 421 → 556 (exp. 556); everything fitted unchanged |
| Bent π | A bent carbon with two σ partners has one p orbital out of the plane, so more than one π bond's worth on it costs 200·(π − 1)² kJ/mol, fading out between 143° and 162°. A phenyl radical otherwise borrowed a second π bond from both neighbours (1.89 each) and its C–H came out 122 kJ/mol too weak. A cap on the bond order itself was tried and left free valence that hydrogens smeared over (What was tried); a cost leaves the bond orders alone | Phenyl C–H 351 → 472 kJ/mol (exp. 473); acetylene and allene, being linear, are untouched; hot benzene no longer finds states below its own minimum |
| Three-membered rings | An angle whose two outer atoms are themselves bonded bends at 6 % of the usual cost (`ring3` 0.94, switched by that bond): cyclopropane's bonds bow outward, so its 60° angles do not cost what bending tetrahedral carbon to 60° would. Without it cyclopropane was 611 kJ/mol underbound | Cyclopropane atomization 3403 kJ/mol (exp. 3401); C₄–C₆ rings, not fitted, within 3 % (cyclobutane +2.3 %) |
| π torsion | Across a bond of order n > 1: `(n − 1)²·(De₂ − De₁)·⟨sin²φ·sin²θ₁·sin²θ₂⟩ / 0.5625`, averaged over the substituent pairs and switched with every bond involved. It fades out as either end gains a third substituent (2 → 2.5), since that centre is turning sp³. Twisting breaks the π bond, so the barrier is the π-bond energy itself; the angle factors keep it smooth through linear geometries | Ethene 90° twist: 248 kJ/mol unrelaxed (exp. ≈ 272) |
| Non-bonded | Shielded Lennard-Jones (UFF); its Pauli wall fades where the Morse term already repels, and between atoms that can still bond. Shifted-force Coulomb between saturating bond-polarisation charges | UFF, Pauling electronegativity |
| Charge scale | 0.38 e per unit electronegativity difference per bond, calibrated on condensed water rather than the gas-phase dimer — the same deliberate over-polarisation the TIP3P family uses. Covalent bonds, NaCl and HCl are unaffected | Water cohesion, O–O distance |
| Walls | Specular hard reflection or optional soft harmonic field; pressure from normal momentum transfer | Idealized boundaries |
| Wall reservoir | Local OU coupling at the boundary; finite heat capacity, gradual heater response | Model parameters, not a calibrated apparatus |

Every force is the exact gradient of the energy, including the many-body saturation, screening, charge and
VSEPR terms; `physics-check.cjs` verifies this numerically (errors ≈1e-6 kJ/mol/Å). The only lagged
quantity is the bond multiplicity n, which relaxes over about 12 fs.

Because n relaxes on its own clock, it changes the energy without any force doing work, and finer
sub-steps cannot remove that. The engine pays for it instead. n is held fixed through every sub-step
of a step, so each sub-step integrates one energy surface, and relaxes once at the end of the step;
a pair moves only when that would shift it more than 0.001 from the value last paid for, so a calm
scene almost never pays. When n does move, the energy change at those final positions is taken from (or
given to) the kinetic energy of the atoms whose bonds changed, about their own centre of mass, so
momentum is untouched, and the step's last half-kick uses the forces of the n it was integrated with.
Relaxing n inside every sub-step instead, with sub-threshold moves left unpaid until they added up,
let n change at one geometry and be paid for at another: a lone ethane at about 1200 K, insulated,
gained 12–14 kJ/mol in 10 ps that way (now within 2.5), and the gain shrank with the sub-step rather
than vanishing, which is how it was found. If they and their
bonded neighbours cannot afford it, n holds where it was until they can, unless the pointer is
dragging something: then the hand pays, and the cost goes into `servoWork` like any other work it
does. Without that, Cl· dragged onto cold ethene stalled at 3 Å, because the π bond starts to give
way about 14 kJ/mol before the new C–Cl bond pays it back. An atom placed less than
50 fs ago settles its bonds unpaid, because a dropped, drawn or reloaded molecule starts from single bonds
and its π bonds forming is not a reaction; changes made by `minimize` are not paid either.

Before this, an isolated cyclohexene at 1200 K gained 1000–1400 kJ/mol in 2 ps and tore itself apart,
and benzene at 1500 K gained 300–1000 kJ/mol on the seeds where it broke up. Now an intact molecule holds
its energy to about 5 kJ/mol, and one that breaks up stays within about 40. No bookkeeping is needed while
n is still, so a scene at room temperature runs at the same speed.

## Validation

`node playground/tests/physics-check.cjs` — structure, energy conservation, determinism:

| Check | Model | Experiment |
| --- | --- | --- |
| H–H, H₂ bond energy | 0.741 Å, 436 kJ/mol | 0.741 Å, 436 |
| H₂O O–H, angle | 0.959 Å, 104.2° | 0.958 Å, 104.5° |
| CH₄ C–H, angle, atomisation | 1.087 Å, 109.47°, 1665 kJ/mol | 1.09 Å, 109.47°, 1652 |
| C–C / C=C / C≡C | 1.55 / 1.34 / 1.20 Å | 1.54 / 1.34 / 1.20 |
| Benzene C–C | 1.42 Å | 1.39 |
| O=O, N≡N energies | 498, 945 kJ/mol | 498, 945 |
| Water dimer | −26.2 kJ/mol, O···O 2.885 Å | −21 (De), 2.91 Å |
| Methane dimer | −3.8 kJ/mol | −2.2 |
| NVE energy drift, 2 ps | < 2 kJ/mol | — |
| Step back ×2 | bit-identical | — |

`node playground/tests/reactions.cjs [--dyn]` — reaction barriers from constrained collinear scans,
reaction energies, and whole mixtures run for tens of picoseconds:

| Elementary step | Barrier, model | Barrier, lit. | ΔE model | ΔH lit. |
| --- | --- | --- | --- | --- |
| H + H₂ → H₂ + H | 43 | 40 | 0 | 0 |
| H + CH₄ → H₂ + CH₃ | 41 | 50 | −15 | +3 |
| H + Cl₂ → HCl + Cl | barrierless | 8 | −189 | −189 |
| F + H₂ → HF + H | 15 | 4 | −131 | −130 |
| O + H₂ → OH + H | 44 | 37 | −27 | +8 |
| OH + H₂ → H₂O + H | 44 | 17 | −29 | −62 |
| Cl + H₂ → HCl + H | 58 | 23 | +5 | +4 |
| H + O₂ → HO₂· | ≈0 | 0 | −185 | −205 |
| H + O₂ → OH + O | 29 | 70 | +35 | +70 |

| Overall reaction | ΔE model | ΔH lit. |
| --- | --- | --- |
| 2 H₂ + O₂ → 2 H₂O | −486 | −484 |
| H₂ + Cl₂ → 2 HCl | −184 | −185 |
| H₂ + F₂ → 2 HF | −543 | −546 |
| N₂ + 3 H₂ → 2 NH₃ | −116 | −92 |
| CH₄ + 2 O₂ → CO₂ + 2 H₂O | −755 | −802 |
| 2 CO + O₂ → 2 CO₂ | −1062 | −566 |
| C₂H₄ + H₂ → C₂H₆ | −116 | −136 |
| Na + Cl → NaCl | −410 | −412 (bond energy) |

Mixtures that react on their own (40 ps each, `--dyn`):

- Na + Cl₂ → NaCl already at 300 K; complete conversion at 1000 K.
- H· + Cl₂ → HCl + Cl· at 300 K.
- H₂ + Cl₂ → HCl at 3000 K; H₂ + O₂ → H₂O at 3000 K (H· initiation visible from 1500 K).
- CH₄ + O₂ at 3500 K goes through CH₂O, OH· and H₂O — the real combustion intermediates.
- Scattered H and O atoms build H₂, O₂, OH· and H₂O.

**Why a mixture may appear inactive.** At the default target rate, one wall-clock second covers only 20 picoseconds (often less when CPU-bound). A reaction may need activation, a catalyst, solvent, or a mechanism this model cannot represent. Lack of a reaction here is not evidence that a real mixture is unreactive. Heating the sandbox or adding radicals explores this model only.

**Limits.** A classical reactive model, not quantum chemistry. Spin is absent: O₂ is patched to behave
like the triplet diradical it is, but the general case is not. Some barriers are off by 20–40 kJ/mol
(Cl + H₂ and OH + H₂ too high, H + O₂ → OH + O too low), so *relative rates are qualitative*. CO is
modelled with a double instead of a triple bond, so carbon-monoxide energetics are wrong. Water condenses
but freezes far too cold, and π stacking is over-bound (benzene dimer ≈ −18 kJ/mol against ≈ −10).
Transition metals, hypervalent geometry, tunnelling, excited states and solvent chemistry are not
represented.

## Nomenclature bridge

**Playground** in Nomenclature exports `chem-playground/molecule@1`. It contains atoms in Å with every implicit hydrogen made explicit, bonds and orders, and formal charges. Wedge and dash bonds lift the wide end to ±0.9 Å in z. The molecule is delivered through a same-origin `BroadcastChannel` and a shared `localStorage` inbox. If no Playground tab is open, it opens one with the molecule in the URL hash. You can also paste the JSON into the Playground.

In the Playground, a molecule is first **conditioned** in a spherical cell:

- **Room**: a 298 K sampling bath; nominal confinement size estimated from 1 atm.
- **Zero K**: classical geometry minimisation and zero initial velocity; no quantum zero-point motion.
- **Custom**: temperature and a nominal pressure used to estimate confinement size. This cell has no barostat and does not establish bulk pressure. A 2 ps preview is sampling, not proof of equilibration.

Then click to place, drag to throw, Q/E to rotate, or Shift-click to place several. **Place with zero initial velocity** disables both internal and translational velocities, including a drag launch; forces can still accelerate the molecule after placement. Otherwise conditioned internal motion and Maxwell molecular translation are retained.

## The console

The apps icon in the gauges opens a window that sits **on** the scene rather than over it. Opening it holds the clock where it stands and closing hands time back exactly where it was left; nothing in the simulation changes in between. The field keeps drawing behind it, so a change made in **Environment** shows in the chamber while it is still being made. The window can be dragged by its title bar and remembers where it was put.

Three apps: **Environment** (bounds, void walls, chamber size, pressure, thermostat, and a live account of what the boundary has absorbed), **About** (the model notes below), and **Keybinds**. Environment states itself in drawings rather than prose: the diagram at the top runs an atom into the boundary that is actually configured, so switching to a forcefield visibly turns it earlier and switching on a temperature void flattens the leg it leaves on.

## Seeing the third dimension

The view is orthographic and looks straight down z by default. **Right-drag empty space** turns the slab about its own centre for as long as the button is held; releasing it swings back to face you. While the view is turned, the chamber is drawn as a wireframe with its near edges bright, and every pointer gesture — dragging an atom, placing one, the heat brush, the tweezer — works in the plane you are looking at, because pointer positions are mapped back through the same rotation. Chamber-edge resizing is disabled while the view is turned, since the faces are no longer screen-aligned.

Value fields inside the console do not take the scroll wheel; the panel scrolls instead.

## The atom inspector

Right-click or double-tap an atom with the hand tool. **Orbitals** displays the actual converged molecular-orbital coefficients, with orbital index, energy, α/β spin, xy/xz/yz slice and phase/probability controls. Gold and blue indicate opposite wavefunction phases, not charges. A plane passing through a node can be empty; switch planes to see the lobes. A finite minimal basis can have no unoccupied orbital in a spin channel. A plotted empty orbital is not an occupied electron. Probability brightness is gamma-scaled for visibility and the displayed square is a spatial slice, not a projection or a particle path.

The calculation stays in a cancellable worker and reuses converged orbitals when changing views. Tests check wavefunction normalization, antibonding nodes and reconstruction of density from occupied α/β orbitals.


Right-click an atom. The card is pinned beside it with a leader line drawn on the field, and follows the atom as the camera moves.

**Electron cloud** shows that one atom's share of the fragment's calculated density, by Hirshfeld's stockholder rule:

    w_A(r) = rho_A_free(|r - R_A|) / sum_B rho_B_free(|r - R_B|),    rho_A(r) = w_A(r) rho_mol(r)

The free-atom references are ground-state atoms in the same STO-3G basis, spherically averaged over 26 directions onto a radial table. The molecular density being divided is a real self-consistent Hartree–Fock solution for the whole fragment at the geometry it has at that instant, so the shape shown is that atom as the neighbours within this isolated fragment and the present geometry have made it. Surrounding fragments are omitted. **Density** switches to the undivided density. Both are computed in one worker pass, so the toggle is instant.

Hirshfeld shares sum exactly to the molecular density and no share can exceed it; `node playground/tests/quantum.cjs` checks both.

## Light

**Lamp** (sun button in the top bar, or `L`): a condition, like temperature, that stays on until you
turn it off. It shines 330 nm light at 10¹⁷ photons cm⁻² s⁻¹. Cl₂, Br₂, I₂ and F₂ absorb it with
their measured cross-sections (Cl₂ 2.6·10⁻¹⁹ cm², so each Cl₂ splits about once in 40 s); nothing
else in the chamber absorbs at that wavelength. A split bond comes apart with the photon's energy as
motion along it (362 kJ/mol for Cl₂: 243 to break it, the rest as the speed of the two atoms). In a
real gas the two atoms fly about 100 Å apart before they slow down, more than the chamber is wide, so
one of them is set down in a clear spot elsewhere in the chamber. With the
lamp on, photolysis is one more reaction the skip can pick: ethene and Cl₂ sit unchanged in the dark
(the clean gas-phase reaction is far too slow; the fast dark addition in a flask runs through ions on
the glass or in a solvent, which the engine does not have), and with the lamp on `J` goes Cl₂ → 2 Cl·,
Cl· + C₂H₄ → C₂H₄Cl·, C₂H₄Cl· + Cl₂ → C₂H₄Cl₂ + Cl·, and round again.

**The chamber under the lamp is a window into a lit flask.** A chamber of six molecules with two
chlorine atoms in it has one radical for every three molecules; a real lit flask has about one in a
million, because radicals are made slowly and destroyed as soon as two meet. That ratio is what lets a
real chain run thousands of times before it stops, and in the small chamber the radicals would mostly
find each other instead. So while the lamp is on, two radicals meet at the rate they would in the flask:
the steady radical concentration √(R/k_t), from the lamp's own photolysis rate R and a typical
radical–radical rate constant k_t = 3·10⁻¹¹ cm³ s⁻¹, instead of the chamber's own. Everything else
(radical + molecule, molecule alone) keeps the chamber's concentrations. With the lamp off, or with
nothing left that absorbs, radicals meet at the chamber's rate, so radicals you place yourself still
pair up in picoseconds.

## Measured barriers

The force field gets most barriers within 10 kJ/mol (the scorecard below), but a few reaction classes
are far enough off to change which reaction comes next, and the barrier the search finds for them
also moves by ±10 kJ/mol with the angle the two molecules meet at. For those classes the skip-ahead
turns the measured barrier into the rate instead of the one it found, and the measured prefactor
where the model's own (a tenth of the collision rate) is far from it. The reverse of each reaction gets
the measured barrier minus the measured reaction enthalpy, with the same prefactor, so both directions
follow real thermochemistry. A hydrogen taken from a radical (Cl· + C₂H₄Cl· → C₂H₃Cl + HCl) is
not in any class: radicals meeting have no barrier. Neither is a hydrogen on a double-bonded carbon
(vinyl C–H, 465 kJ/mol), which keeps the model's own barrier. The live simulation and the reaction played after a skip are unchanged.

| class | measured, kJ/mol (search found) | prefactor |
| --- | --- | --- |
| Cl· adds to C=C | Cl· + C₂H₄: 0 (10); back: 75 (249, the scan cannot find it) | ×3 |
| H· adds to C=C | H· + C₂H₄: 9 (0) | ×0.2 |
| alkyl adds to C=C | CH₃· + C₂H₄: 31 (9) | ×0.01 |
| Cl· + H–H | Cl· + H₂: 18 (36); H· + HCl: 14 (31) | ×0.45 |
| Cl· + H–CH₃ | Cl· + CH₄: 11 (25); CH₃· + HCl: 4 (18) | ×0.2 |
| Cl· + H–C beside a halogen | Cl· + C₂H₄Cl₂, CH₃Cl, CH₂Cl₂: 10 (8) | |
| H· + H–O | H· + H₂O: 76 (97) | |
| H· + H–C | H· + C₂H₆: 38 (11); H· + CH₄: 50 (41) | |
| CH₃· + H–C | CH₃· + C₂H₆: 45 (6) | ×0.01 |
| OH· + H–C | OH· + C₂H₆: 9 (0); OH· + CH₄: 15 (24) | ×0.2, ×0.05 |
| ROO· + H–C | CH₃O₂· + CH₄: 90 (29) | |
| OH· + H–OOR | OH· + CH₃OOH: 0 (2); back, CH₃O₂· + H₂O: 139 (78) | ×0.1 |
| Br· adds to C=C | Br· + C₂H₄: 0 (25), back 30 | |
| Br· + H–H, H–C | Br· + H₂: 82 (72); Br· + CH₄: 75 (75); Br· + C₂H₆: 57 (26) | |

A radical adding to a double bond, for example, has a tight transition state and a real prefactor
near 10⁻¹³ cm³ s⁻¹, a hundred times below a typical collision.

With them the classic chains run at room temperature as they do on a bench: methane and Cl₂ under the
lamp give CH₃Cl and then CH₂Cl₂, H₂ and Cl₂ give HCl, ethene and Cl₂ give 1,2-dichloroethane, and
ethene and Br₂ give 1,2-dibromoethane. The
underlying errors (ethyl C–H 41 kJ/mol too weak, peroxide O–H too strong, Cl···H complexes too deep)
are still in the force field and in the live simulation; ROADMAP has them.

## Forecasting and skipping ahead

Most reactions take far longer than a browser can simulate: at about a picosecond of chemistry per
second of screen time, one real second would take tens of thousands of years. A reaction is also not a
slow creep: a molecule rattles in its well for a long time and then crosses in about 100 fs. So the slow
part is never simulated. It is calculated, and only the crossing itself is shown.

**Forecast** (toolbar hourglass, `W`): click two atoms. The chamber pauses, and the answer appears
between them in the scene: the half-life, what they become, and the energy along the way. Click either
atom again (or Enter) to skip ahead to it; the reaction is first tried on three copies, as with `J`, and
a red ✕ instead of a skip means it did not end as forecast. Bonded atoms ask when that bond breaks;
atoms of different molecules ask when the two react. A copy of just those molecules is pulled
through the reaction in a worker while every other atom relaxes (a relaxed scan; the live chamber is not
touched). Several routes are tried and the lowest barrier kept:

- joining the two atoms, with the attacker brought in along its current line, from the target's open
  side, and from both faces of a flat (sp²) centre, turned so its own free valence leads;
- transfers, where one of the target's existing bonds lets go, pulled along the difference between the
  bond forming and the bond breaking and approached collinearly from behind the bond that breaks;
- breaking, by stretching the bond 3.2 Å.

Each route is scanned both ways, reactants to products and back, and the lower energy at each point
is kept: pulling one way only lets a bond hang on past the point where it should give way and snap later,
which put OH + CH₄ at 96 kJ/mol where the live chamber (and reality) says about 15. The pass back from the products only
keeps states whose bonds all exist in the reactants or the products: pulled back, the hydrogen of
H + Cl₂ grabbed the other chlorine and gave HCl + Cl at the reactants' end of the path. A route that ends with
an atom holding more than half a bond over its valence (a hydrogen stuck between two partners) is not
counted as a product.

A barrier is never taken lower than the reaction's own energy: a scan pulls the products only a little
way apart, so an uphill route used to report less than it costs to finish it (CH₃· + ethene → propene
+ H· read 16 kJ/mol for a reaction 24 kJ/mol uphill, and came up as a stray H· in polymerisation runs).

The barrier becomes a half-life by transition-state theory, rate = A·e^(−barrier/RT): A ≈ 10^15.5 s⁻¹
for a bond simply coming apart, kT/h for other unimolecular steps, and for two molecules the
hard-sphere collision rate × 0.1 times the partner's number density in the chamber. The scene shows the
half-life and the energy along the path; the uncertainty behind it is ±10 kJ/mol on the barrier and
×/÷10 on A.

`tests/benchmark.cjs` scans 39 textbook reactions and bonds and prints the scorecard: barriers rms
12.8 kJ/mol over 23 reactions, bonds rms 33 kJ/mol over 16 (14 without the formyl outlier). Selected rows:

| Reaction | Forecast | Real |
| --- | --- | --- |
| H + H₂ → H₂ + H | 43 | 40 |
| H + CH₄ → H₂ + CH₃ | 46 | 50 |
| OH + CH₄ → H₂O + CH₃ | 17 | 15 |
| OH + H₂ → H₂O + H | 36 | 15 |
| H + HCl → H₂ + Cl | 33 | 15 |
| Cl + CH₄ → HCl + CH₃ | 31 | 11 |
| Cl + H₂ → HCl + H | 37 | 19 |
| F + H₂ → HF + H | 7 | 4 |
| H + C₂H₆ → H₂ + C₂H₅ | 12 | 38 |
| Cl + C₂H₄ → C₂H₄Cl | 11 | ≈ 0 |
| CH₃ + C₂H₄ → C₃H₇ | 9 | 31 |
| CH₄ → CH₃ + H | 439 | 439 |
| C₂H₆ → C₂H₅ + H | 382 | 423 |
| CH₃–Cl | 357 | 350 |
| HCO–H (formaldehyde) | 493 | 369 |

The scorecard also gives the deepest complex along each path. Hydrogen taken by a heavy atom passes
through one bound by 20–27 kJ/mol (Cl + CH₄, Cl + C₂H₆, OH + NH₃) where the real ones are
bound by a few at most; these are the CH₄Cl-type complexes the chlorination chain runs into.

The errors fall into three families, which is what a fit has to target:
- a hydrogen passing to or from oxygen or chlorine (OH + H₂, H + H₂O, H + HCl, Cl + H₂, Cl + CH₄)
  is 15–23 kJ/mol too high: the half-made O···H···H and Cl···H···C crossings cost too much;
- anything that makes an ethyl radical is too easy because ethyl's C–H is 41 kJ/mol too weak (H +
  C₂H₆ 26 too low, CH₃ + C₂H₄ 21 too low);
- the formyl radical HCO· gets none of the stabilisation that makes formaldehyde's C–H so weak in
  reality: oxygen's valence of 2 caps the C–O bond at a double bond, where in HCO· it strengthens
  toward a triple. Aldehyde chemistry needs this before it can be forecast.

Where the forecast and reality disagree, it is the force field's barrier that is wrong, and the
forecast shows it. Its H + H₂ rate at 1500 K, 1.7×10⁻¹² cm³/s, is close to the measured 1.4×10⁻¹².

Reactant and product energies are relaxed from a copy nudged by up to 0.03 Å, because a perfectly
symmetric start stays symmetric through every relaxation: water formed by OH + H₂ on one line stayed
linear, 136 kJ/mol above bent water. The scans themselves are not nudged: even 0.005 Å changed where
some reactions are committed (Cl· + CH₄ played out in 4 of 6 trials instead of 6), and a live chamber
is never perfectly symmetric anyway. A hand-built collinear start can still read a barrier too high
(OH + H₂ on one line: 69, against 36 slightly off the line).

**Skipping ahead**: the clock jumps by a waiting time drawn at random from the
forecast rate, the molecules are placed at the first point clearly past the top of the barrier with a
thermal push along the path, and the live simulation carries out the reaction itself. The skipped time
is shown beside the simulated time, is part of saved scenes and undo, and Ctrl+Z puts everything back.
The skipped time floats up from the reacting pair with the reaction's own energy curve drawn under it
(the climb to the barrier, marked with a dot, and what the products give back), and the first 0.8 ps after a skip run at 0.5× (about
80 fs a second), so the bond-making itself can be watched, and the speed then returns to what it was (kept-going chains stay at
full speed).
The reacting pair is placed where its own scan put it, in a copy that held only those two molecules;
any other molecule found within 2.4 Å of them in the real chamber is moved straight out, whole, first
(over the time that was skipped it would have drifted anyway). Without that, a pair could land on a
bystander, and an HCl struck that way once fell apart into H· and Cl· at room temperature. A bystander
that cannot be pushed straight out, because a wall is in the way, is set down whole in the emptiest
place in the chamber instead, and the reacting pair itself is kept inside the walls. Before that, one
caught between a pair and a wall was left 0.9 Å away and blown apart: HCl into H· and Cl·, a water
into H· and ·OH, both at room temperature, which is how this was found.

Where a path climbs its barrier, drops into a well and then climbs a second, lower bump before the
products (the Cl···H···C wells, above), the pair is placed past that last bump instead, provided the
well in between is deeper than 8 kJ/mol. Placed past the first top, Cl· + CH₃Cl stayed stuck as a
CH₃Cl₂ complex in 11 of 30 seeded trials; past the last bump, in none. The barrier and the rate are
still those of the highest top; only the starting point moves.

A reaction with no barrier that joins two molecules into one (CH₃· + CH₃·, OH· + OH· → H₂O₂) starts
further along its path, at 60 % of the way down instead of 30 %: from the shallower point the two
hydroxyls bounced apart in 5 of 6 trials, and from the deeper one they join in 6 of 6.
A skipped reaction's heat is taken away the way a gas would, by later collisions: 80, 160, 300, 500 and
800 fs after the skip, each piece that is already one of the expected products is set back to the
chamber temperature, and at 800 fs whatever is left is too. Pieces that are still together wait. When
everything was cooled at 80 fs whatever its state, H· + Cl₂ → HCl + Cl· (188 kJ/mol downhill, the
products flying apart) was stopped mid-separation and dropped into the model's Cl···H···Cl well in
half of all trials at 500 K. A trial also counts as done when the products are there but still
touching, which is how the feed already reads such a cluster (HCl···Cl·): the over-coordinated atom's
weakest link is set aside and what remains is compared. In a crowded 500 K chamber H· + Cl₂ now passes
12 of 12 trials instead of 5–6, and a hydrogen–chlorine chain no longer stalls on repeated ✕.

Played this way, Cl· + ethene, CH₃· + ethene and H· + CH₄ complete in 9 of 9 browser trials and in
at least 3 of 4 seeded runs each in the tests.

**Next reaction** (`J`, or the double-arrow button beside step-forward): the chamber is paused
and searched for everything that could happen: a radical meeting any atom of another molecule (the radical is always
the one that attacks), a π bond meeting the hydrogens of another closed-shell molecule (above 500 K only: below it
that four-centre addition, about 170 kJ/mol here and symmetry-forbidden in reality, never wins), and every
distinct bond coming apart (at most 24
candidates, nearest first). While it searches, each candidate pair is a dashed line in the scene, and
it turns solid as its barrier comes in: brighter and thicker the easier the reaction. Each is forecast
in a pool of workers, quickly at first and then properly for the three lowest barriers (in the quick
pass a route is dropped as soon as it climbs more than max(40 kJ/mol, 14·RT) above the best route
already found for that pair, a million times slower, which takes about a fifth off a growing polymer's
search; asking about a pair by hand still reports every route), and counted as
often as it occurs (methane's four hydrogens make four times the rate). Then it draws which reaction comes first and when from those rates, as a
real flask would (kinetic Monte Carlo), and plays it in the live simulation; press J again for the one
after. A stable chamber says so in the scene instead ("stable · nothing for …"), and a reaction that
does not play out as forecast on a trial copy is marked with a red ✕ and not skipped, and the chamber
then runs on, so the next J meets molecules that have moved instead of failing the same way again.
A lone oxygen or sulfur atom is a triplet in its ground state, and a triplet atom cannot slip into a
σ bond of a closed-shell molecule in one step (only the excited O(¹D) does that): the model has no
spin, and without this rule it made methanol from O + CH₄ at 42 kJ/mol, as easily as the real
abstraction. Such single-product routes are not counted, so O + CH₄ → CH₃· + ·OH (38 kJ/mol, real
about 40), O + H₂ → ·OH + H· (40, real 37) and O + H₂O → 2 ·OH; additions to π bonds and to radicals
are untouched. A radical knocking a hydrogen off a carbon in one step (substitution at carbon,
CH₄ + O· → CH₃O· + H·) is not tried either: the model put it at its bare reaction energy, 53 kJ/mol,
where the real backside attack costs far more and never competes with taking the hydrogen; at a double
bond the real route is addition and then loss of the hydrogen, which the search finds as two steps.
Leaving these out took about a third off a growing polymer's search.
When a scan that hands an atom over ends with its two products already apart, they are relaxed apart
too; relaxed together, CH₃· and ·OH from O + CH₄ paired up into methanol and the scan reported that
second step as the reaction. A trial now also counts if the products were there at 40 fs or at any
cooling point, since a radical pair born in contact may pair up straight afterwards (the feed then shows
both steps), but not if it ended back where it started: a pair that crossed the barrier and came
straight back over it did not react, and a skip that changes nothing is worse than an honest ✕.
Two lone atoms joining into a diatomic (Cl· + Cl· → Cl₂) count at a thousandth of their meeting rate:
on their own they fly apart again within one vibration, and in a real gas at 1 atm another molecule
arrives to carry the energy off only about once in a thousand such meetings.
Every distinct outcome a scan finds counts with its own rate, not only the lowest: Cl· + ethene lists
addition, abstraction and substitution, and CH₃· + C₂H₃· both recombination and disproportionation.
Reactions that only swap identical partners are left out of the list. If the molecules a forecast was
made for have changed by the time you skip, it refuses and asks for a new forecast.

After two or more skips the chamber's lower-left corner draws the last fourteen as one energy diagram,
each step starting where the one before ended: a chlorination chain reads as small climbs and large
drops, energy going down step by step. Loading an experiment or emptying the chamber starts it afresh.

**Keep going** (`Shift+J`, or double-click the double arrow, which then pulses; `J` or Esc stops it)
repeats this: skip to the next reaction, watch it for 3 ps, search again, until nothing more can happen
on a human timescale. A reaction that does not play out on its trial copies does not end the run: the
chamber is left to move for 3 ps and searched again, and only three such misses in a row stop it.
It also stops when the chamber has only gone round two or three states over its last six searches,
each seen at least twice, and says so in the scene ("equilibrium · CH₃Cl ⇌ CH₃· + Cl·", "going round in
circles · …"): at 1000 K a methyl chloride otherwise came apart and re-formed for as long as anyone
watched, each round skipping minutes, and Cl· + 2 ethene at 500 K settled into adding and losing an
ethene, a ceiling-temperature balance. A real chain never trips this, because it uses its reactants up
and no state comes back. Three misses in a row stop it too, and it says that as well. `J` still steps
on by hand. One Cl· in a chamber of three CH₄ and three Cl₂ at 298 K ran the textbook chain
by itself: Cl· + CH₄ → CH₃· + HCl, CH₃· + Cl₂ → CH₃Cl + Cl·, then CH₃Cl → CH₂Cl· → CH₂Cl₂, with the
reverse CH₃· + HCl → CH₄ + Cl· in between; 30 reactions and 323 µs of chemistry in about six minutes.
It also produced CH₄Cl and CH₃Cl₂ complexes, which are the force field's over-coordination wells
(above), not chemistry.

In Cl· + ethene + methane at 298 K it gives the addition to the double bond 99.6% and methane's
hydrogen 0.4% (in reality about 99.9% and 0.1%); taking ethene's own hydrogen needs 46 kJ/mol. (Before the σ-radical
term, a vinyl C–H 50 kJ/mol too weak made that abstraction look downhill and it came out first.) in a 600 K mixture of H·, CH₃·, ethene and methane,
H· + ethene 72%, CH₃· + H· 26%, CH₃· + ethene 1%. Two methanes at room temperature: nothing on any
human timescale.

**Experiments** (tiles in the middle of an empty chamber, and a third tab in the add-atoms library):
ready-made textbook chambers, each loading its own conditions (temperature, and the lamp for the ones
marked ☀): radical chlorination of methane (room temperature, lamp on: CH₃Cl, then CH₂Cl₂), chlorine
adding across ethene (lamp on: 1,2-dichloroethane by the chain), radical polymerisation of ethene
(CH₃· → C₃H₇· → C₅H₁₁·, about one step every 0.1 ms at 350 K with the measured 31 kJ/mol barrier, against 0.3 ms from the measured rate), methyl
radicals recombining, hydroxyl with methane in air (OH· + CH₄ → H₂O + CH₃·, then CH₃· + O₂ → CH₃O₂·,
the first two steps by which the atmosphere removes methane; after them the model's peroxide O–H, about
60 kJ/mol too strong, lets CH₃O₂· take a hydrogen from water within hours, which in reality is far
uphill and never happens), and hydrogen with chlorine, the classic photochemical chain (Cl· + H₂ → HCl + H·,
H· + Cl₂ → HCl + Cl·) at room temperature until only HCl is left. Loading one is undoable; `J` then
skips to whatever happens next.

The search also knows two-bond swaps: two atoms of one molecule joining while each lets go of a
partner, such as H₂ leaving ethane (1,2-elimination, 400 kJ/mol; 1,1 to a carbene, 403). A route whose
products, relaxed apart, would sit more than 60 kJ/mol above where the scan left them is a bound
complex on its way somewhere else (a CH₂ half inserted into methane is ethane in the making), not those
products, and is not counted.

### How far the forecast and the live simulation agree

The forecast measures the model's energy surface; the live chamber moves on it. Checked against each other:

- **CH₄ at 2000 K**: the forecast says days; in the live chamber it stays whole for 30 ps in 8 of 8 runs.
  The splits seen earlier at that temperature were the fragment count briefly treating a stretched C–H
  as broken.
- **H + H₂**: an Arrhenius fit of exchanges counted in the live chamber (1000, 1500, 2500 K) gives a
  barrier of 44–46 kJ/mol, the same as the forecast, but a prefactor about 40 times the collision
  estimate, so the live chamber exchanges about 40 times faster than both the forecast and the real
  reaction. The model's H captures H₂ from roughly twice the real distance.
- **Ethane above about 2200 K**: the live chamber eliminates H₂ within picoseconds; the forecast's
  barrier for that (1,2-elimination 400 kJ/mol, measured with two-way four-centre scans) says
  microseconds or longer. The cause is not yet found. A C–C bond order above 1 in ethane is not a cheap
  shortcut either: at ethane's own geometry it costs 64 kJ/mol at 1.2 and 425 at 1.5. It is not the thermostat (the statistics match equipartition at 800–2000 K, and under CSVR it
  still happens within picoseconds: 4 of 6 runs in 20 ps) and not energy accounting (insulated runs hold to 2 kJ/mol). The bond orders are implicated —
  with them frozen ethane survives 20 ps at 2500 K — but freezing them also stops ethylene from forming
  its π bond, so that test alone does not prove it. A nudged-elastic-band search for a lower path
  fell apart (images flew to separated atoms through the force field's cutoffs) and was not kept.

Repeating the search remembers every barrier it has measured, by reaction type and temperature,
and only scans what is new; the reaction it skips to is played from the geometry its full scan found
(scanned afresh only if the search had it from memory), and is tried on a copy first. If it does not
end as forecast on the copy, a red ✕ marks the pair and nothing is skipped (ending in another of the
same molecules' known outcomes counts: OH· + OH· may become H₂O + O· or H₂O₂, and which one is the
dynamics' choice); it does not fall back to a less
likely reaction, which would bend the odds. Two artefacts the search showed up earlier are gone:
C₂H₃· + H₂ used to come out as H₂ adding whole to the vinyl radical (C₂H₅·, barrierless); it is now the
real abstraction, C₂H₃· + H₂ → C₂H₄ + H·, over 27 kJ/mol (real about 40, the reaction energy right),
since the scan back from the products may no longer swap to an equivalent product. And a hot C₂H₅
formed by H· + ethene lost its hydrogen again within picoseconds in the live chamber.
That second one is now handled the way a real gas handles it: a freshly formed molecule carries its
whole reaction energy, and in a real flask collisions with the surrounding gas take it away within
nanoseconds, while the chamber holds too few molecules to do that. The product is therefore brought
back to the chamber's temperature at 80, 160 and 300 fs after the skip, standing in for those
collisions. Before a skip is played, the chosen event is also tried on three copies at once (400 fs
each); if fewer than two of those end in one of its known outcomes, nothing is skipped.

The search used to scan ethene's carbon reaching for the H atom rather than the H atom reaching for
ethene, which gave a poorer starting geometry: H· + ethene then completed in 6 of 10 runs. Radicals
are now always the ones that attack, and it completes in 6 of 6 seeded runs and 5 of 5 browser trials.

### What it does not do yet

- A random sweep (eight chambers from twelve small molecules and radicals at 298, 600 and 1200 K, three
  skips each) ran without errors or impossible fragments, and the trial copies turned away two products
  the force field should not make (HClO₂ from O₂ + HCl, and Cl· + NH₃ ending in a complex). One got
  through: at 1200 K O₂ + NH₃ formed H₃N–O₂, because nitrogen may take a valence of five (for nitro
  groups and N-oxides) and the cost of doing so is too low for an amine and O₂.

- Kept going on its own, one Cl· with three CH₄ and three Cl₂ runs the chain until the Cl₂ is used
  up (20 steps in three minutes, 24 µs of chemistry, ending in CH₃Cl and HCl, then Cl· and CH₃Cl
  swapping a hydrogen back and forth) and then stops: hydrogen passing between chlorine and carbon
  can end in a Cl···H···C complex 20–27 kJ/mol deep in the model (the scorecard's well column), and
  when the trial run ends there a red ✕ marks it and nothing is skipped. The fix belongs in the force
  field (ROADMAP, item 1). Whether the molecules are still the ones the search saw is judged from the
  same rebuilt copy the search used, so a borderline complex no longer reads as a change.

- The search takes 10–60 s for a handful of molecules and is capped at 24 candidates, so a crowded
  chamber is only partly searched.
- Two closed-shell molecules reacting through their π systems (Diels–Alder, ene), and rearrangements
  inside one molecule, are not searched; the hourglass can still be pointed at them.
- Solution chemistry: ions, proton transfer and solvent are not in the engine, so most textbook
  organic mechanisms (SN1/SN2, E1/E2, acid catalysis) cannot be forecast or watched yet.
- Pressure-dependent rates (fall-off), tunnelling and zero-point energy are left out. Fall-off is what
  keeps hydrogen combustion from being forecast yet: H· + O₂ → HO₂ wins on barrier, but at high
  temperature HO₂ needs a third molecule to carry away its energy, and without that H· + O₂ → OH· + O·
  (the branching step, which the search does list, at 7% at 1500 K) should win.
- Concerted four-centre additions that orbital symmetry forbids (H₂ or a C–H adding whole across a
  double bond, O₂ inserting whole into a C–H) are not penalised by the force field: a bond-order model
  has no orbital phases to forbid them. The search therefore never counts them: when two closed-shell
  molecules meet through a π bond and a hydrogen, only outcomes that make as many molecules as went in
  (an abstraction) are kept, and below 500 K such meetings are not searched at all. Methane and O₂ at
  2000 K used to start with O₂ inserting into a C–H (CH₃OOH in one step); now the first step is a C–H
  bond breaking, and the radicals take it from there. The live simulation itself can still do it.

## Temperature readings

The gauge shows a running mean of the kinetic temperature with the measured spread beside it. The instantaneous value genuinely wanders — a sample of N atoms has relative fluctuations of about `sqrt(2/3N)`, so twenty atoms swing by tens of kelvin from one instant to the next. That is the sample being small, not friction or a numerical fault. Hovering the gauge gives the mean, the spread, the instantaneous value and the expected fluctuation for the current atom count.

## Controls

Press `?` for every shortcut. Click a key to rebind it; conflicts move automatically, and bindings persist. Letter bindings follow the physical key, so they work on any keyboard layout. `Ctrl+K` opens a command palette that also understands `500 K`, `80 °C` and `2x`.

## Additional checks and model limitations

Run `node playground/tests/regression.cjs` for force-query invariance, charge conservation, variable-radius force gradients, hot rewind across checkpoints, pinned atoms, failure handling, and malformed imports.

Bond multiplicity relaxes as an internal heuristic variable; the energy its relaxation would create or destroy is settled against the kinetic energy of the atoms involved (see Force field). The heat bath exchanges energy; the engine also counts safety velocity clamps in extreme collisions. Neither behavior should be mistaken for isolated, rigorously conservative dynamics. The UI reports those clamps when they occur.

The camera is orthographic: depth changes shading, not apparent atomic radii, and turning the view changes nothing physical. Chamber contours are summed Gaussians, not electron-density calculations; the inspector's cloud is a real Hartree–Fock calculation and the only quantum result in the app. Playback speed is a nonlinear UI setting, not a multiplier of physical real time. With any void wall enabled, energy is deliberately removed and tallied; the chamber is then an open system by design.

Imported molecules now receive Maxwell–Boltzmann center-of-mass translation at their conditioning temperature. Previously, removing this motion for preview alignment and never restoring it left molecules vibrating in place. A drag-to-throw overrides that translation; 0 K placements remain motionless. Existing scenes can use **Resample thermal motion** (in the temperature menu, live observations, or **Shift+T**) to draw fresh thermal velocities. This is an explicit state edit and can be undone.


## Current physics work

The scripted reaction notebook has been removed. No reaction recipes select products or animate an imposed trajectory.

`node playground/tests/thermal-walls.cjs` checks spatial isolation, heater response, heat accounting, velocity statistics, direct temperature edits, setpoint limits and wall-state rewind.

### Radical addition to π bonds — fixed

`node playground/tests/radical-addition.cjs` reports the constrained Cl· approach to ethene. The entrance
barrier was about 128 kJ/mol because the π bond only released once the newcomer was bonded, while the
newcomer could not bond until the π released. Bond order now responds at the range where bonding actually
begins, and a π bond blocks a newcomer less than a σ bond does. The constrained path is now downhill
(≈ −1 kJ/mol at the top), and chlorine radicals chlorinate ethene at 300 K in dynamics, including radical
oligomerisation. Every single-bond barrier in the table above is numerically unchanged by this.
Reference: [Cl + ethylene potential-energy surface](https://doi.org/10.1021/jp001221u).

### Condensed water — partly fixed

With the calibrated charge scale, a relaxed water cluster has a cohesive energy of **−42 kJ/mol per
molecule** (experiment ≈ −44, the enthalpy of vaporisation) at a mean nearest O–O distance of **2.86 Å**
(experiment 2.8). Water therefore condenses: at 273–330 K a cluster stays together as a hydrogen-bonded
liquid droplet with about 2 hydrogen bonds per molecule, and it boils between 350 and 500 K.

It does **not** freeze at 273 K. Cooling turns the droplet glassy only below roughly 150 K, and the network
carries ~2–2.5 hydrogen bonds per molecule where real water has 3.5–4. The missing ingredient is acceptor
directionality: the model places one negative charge on the oxygen, with no lone-pair geometry, so the
tetrahedral network that makes ice is never strongly preferred. Water models that reproduce ice
(TIP4P/Ice, TIP5P) add off-atom charge sites for exactly this reason. A small cluster also melts far below
bulk ice in reality, so part of the gap is finite size. See
[water models including TIP4P/Ice](https://docs.lammps.org/stable/Howto_tip4p.html).

### Over-coordination — mostly fixed

`node playground/tests/overcoordination.cjs` pushes an H onto the back of CH₄. Real CH₅ does not exist.
It used to be bound by about 118 kJ/mol, and five-coordinate carbon from CH₃ + CH₄ by about 84; both
scans are now repulsive everywhere (deepest point about +36 and +23 kJ/mol). Every barrier in the table
above and every atomization energy is unchanged, and CH₂ still inserts into H₂.

The cause was a mismatch between two curves. A C–H stretched by 0.3 Å still collects about 81 % of its
Morse attraction, but `satF` counts it as only 0.68 of a used valence. The stretched hydrogens then
looked open-shell and bonded to each other or to a newcomer side-on. The angular saturation term
(force-field table) removes that side-on attraction instead of changing the saturation curve, so the
fitted barriers are not touched. Together with the bond-order pairing it is the main cost of a step: measured
against the engine before both, a 14-atom scene runs at 0.68 of the old speed, 96 atoms of cyclohexene
at 0.58 and 72 atoms of benzene at 0.54.

A leaving H could also draw partial attraction from two or three carbons at once, bridging them 2.0–2.3 Å
out while using almost none of its own valence. Contacts down to a fifth of a bond (`ANG_G0` 0.12) now
act as arms of the screen too, with two rules that keep ordinary chemistry intact: an arm screens only a
contact no stronger than itself (`ANG_DM`, `ANG_DU`), so the nearer of two weak contacts wins instead of
both cancelling, and a 1-3 contact held by a shared neighbour (an H on one carbon 2.1 Å from the next)
is no arm at all. Cl· + ethene keeps an entrance barrier of about 10 kJ/mol (8.5 before).

Hot hydrocarbons, 20 ps in a bath (`major` mixtures), before → after. The two ethylene rows are from a
26 Å box, three seeds each; a 20 Å box packed the molecules into each other and every result from it
was a crash, not chemistry.

| Mixture | Before | After | Real, on picoseconds |
| --- | --- | --- | --- |
| cyclohexene, 1200 K, alone | loses an H in 0.2–1.6 ps | intact 3 ps, 4 of 4 seeds | intact |
| benzene + H₂, 700 K | 19 H·, fragments | no reaction | no reaction |
| butadiene + ethylene, 900 K | atomised | no reaction | no reaction (Diels–Alder too slow) |
| ethylene + 2 CH₃·, 450 K | not rerun | the two CH₃· recombine to C₂H₆, ethylene untouched (3 of 3) | mostly recombination; CH₃· + C₂H₄ has a 33 kJ/mol barrier |
| ethylene alone, 450 K | not rerun | no reaction (3 of 3) | no reaction |
| benzene, 1500 K, 4 molecules | all broken up | 2 of 4 intact on one seed; broken up on two others | intact |
| cyclohexene, 1500–3000 K | broken up | still loses an H within 0.2–0.6 ps | intact (ns–µs) |

Still open: above about 1500 K hydrocarbons still shed H far too fast. A C–H stretched to 1.4 Å costs
35–45 kJ/mol in cyclohexene against 76 in methane and about 85 for a real C–H, and the relaxed scan lets
the allylic H walk to the other end of the allyl unit (a suprafacial [1,3] shift, forbidden in reality)
for under 160 kJ/mol.

What was tried, so it is not repeated:

- `kb = 1.0` (valence use decays like the attraction) closes both wells but roughly doubles every
  transfer barrier, because `c1`, `c4` and `kb` were fitted together.
- A joint refit of `kb`, `c1`, `c4`, `cap`, `mu`, `kbo`, `k` and `oo3e` needs `c1` ≈ 0.14, which makes a
  small excess on hydrogen almost free: CH₂ + H₂ then stalls as a CH₂·H₂ complex and H₂ dissociates
  among noble gases at 3000 K. With `c1` held at 0.6, closing the wells needs `cap` ≥ 0.19, and the
  barriers rise past 75 kJ/mol.
- Screening every 1-3 Morse attraction removes the H–H part, but the minimiser finds the next well
  (an H capping three stretched hydrogens).
- Gating the screen on the newcomer's free valence alone blocks CH₂ insertion; it has to be gated at all
  three corners, counting free valence with long-reach bonds so a stretched C–H cannot fake a carbene.
- `kb` 1.1–1.2 with the half-made-valence stabilisation `tsStab` 8–22 to win the barriers back: the
  barriers fit about as well (H + H₂ 40), but `tsStab` also rewards the 0.15 spare valence every alkene
  and aromatic carbon carries (benzene 24 kJ/mol too stable, allylic and phenyl C–H weaker), and hot
  combustion then grew CH₅ and CH₆ fragments. A threshold (`tsStab` acting only above 0.2 spare) removes
  the π side effect but fits the barriers no better than the current set.
- The same half-made-valence stabilisation for chlorine alone (25 kJ/mol at half a bond) fits every
  chlorine barrier (Cl + CH₄ 28 → 16, real 11; Cl + H₂ 34 → 21, real 19; H + HCl 31 → 22, real 15;
  Cl + C₂H₄ 11 → 0, real 0) and leaves every bond energy unchanged, but it is also the energy of the
  complexes the chlorine paths pass through, which it deepens: Cl···H···Cl from −22 to −71 kJ/mol,
  Cl₃ from −4 to −36, Cl···CH₄ from −19 to −44. Those would trap atoms in the live chamber. For oxygen
  (10) it also made O₂ 14 kJ/mol too stable, because O₂'s bond orders do not sit at whole valence.
  The barrier and the well have to be separated before a crossing term can be fitted.
- Leaving a bond's own polarisation charge out of the Coulomb term between its two atoms. The term was
  switched off as (1 − f) while the charges grow with f, so a half-made polar bond attracted itself
  (f²(1 − f), largest at about two thirds formed). Removing that keeps every equilibrium molecule and
  bond energy exactly as before and makes the chlorine complexes shallower (Cl + CH₄ −20 → −15, Cl +
  C₂H₆ −24 → −17 kJ/mol), but the same attraction was carrying the crossings: Cl + CH₄ rises 31 → 40,
  OH + CH₄ 17 → 36, OH + NH₃ 1 → 21, and the benchmark goes from 12.8 to 14.4 kJ/mol rms. Adding the
  chlorine/oxygen crossing stabilisation on top lowers the barriers again but deepens the complexes
  with them (Cl + CH₄ −37, Cl···H···Cl −43): in this bond-order model the top of the barrier and the
  complex are the same half-shared state, so no per-atom term separates them. The fix has to change
  how a hydrogen hands its bond from one partner to the other.
- Taking the Cl···H···C complex apart term by term (Cl–H held at 1.65 Å, everything else relaxed). Two
  things hold it together. First, the two ends bond to each other across the hydrogen: C and Cl, 3.3 Å
  apart, reach a bond order of 0.52 and gain about 22 kJ/mol, where real chemistry has the opposite
  (in X···H···Y the electrons on X and Y cannot pair: the triplet repulsion of bond-energy/bond-order
  theory). Second, the C–H at 1.6 Å has already left the 1.33–1.64 Å window that switches the Coulomb
  term off inside a bond, so the hydrogen's +0.3 and the carbon's −0.4 attract at full strength: about
  66 kJ/mol once the Cl–H range is widened. Both were tried as terms with exact forces. Suppressing
  the end-to-end bond through a monovalent bridge, and fading the Coulomb term out with the bond's
  saturation instead of the short switch (with the 1-3 exclusion made consistent, so bond energies
  stay put), removes the complex: the shared state goes from −23 to +18 kJ/mol. Every
  hydrogen transfer to or from oxygen or nitrogen rises with it, because that same attraction was
  carrying those crossings: OH + CH₄ 17 → 50, OH + NH₃ 1 → 46, H + H₂O 99 → 106. With the Cl–H range
  widened (`kb` 1.7) to win the chlorine barriers back, the benchmark ends at 14–17 kJ/mol rms
  against 12.8. Same verdict as above, now with the cause measured: the fix is a hydrogen that hands
  its bond over continuously (n_XH + n_HY ≈ 1 along the path) and carries charge in proportion to the
  bond it has, refitted as a whole.
- Keeping a hydrogen hand-over near a straight line during the scan (a one-sided restraint that lets
  X···H···Y bend to about 120°, after O + CH₄ was seen to swing the oxygen round onto the carbon) left
  the benchmark where it was (12.8 kJ/mol rms) but dropped O + CH₄ from 50 to 23 kJ/mol (real about
  40) and passed fewer of its trials, so it was not kept.
- Scanning a radical meeting a closed-shell molecule from the radical's side only (posing the radical
  onto the molecule, not also the molecule onto the radical) takes the quick search for a growing
  polymer from 13.4 to 8.4 s and the benchmark from 12.8 to 12.5 kJ/mol rms (CH₃ + C₂H₄ 9 → 14, CH₃ + H₂
  35 → 39), but OH + CH₄ rises from 17 to 24 (real 15): the reverse approach is what finds that
  crossing, and the hydroxyl experiment's first step, now at the measured rate, would come 17 times too
  slowly. Both approaches are kept.
- Capping π bonding by the p orbitals an atom has left (bent two-coordinate C, N, O get one) raised the
  phenyl C–H from 351 to 389 kJ/mol, but the valence kept out of π went into extra σ partners in hot,
  bending fragments: CH₄ + O₂ at 3500 K grew C₃H₁₂O₄ clumps. A narrower version (carbon only, and only
  while its two σ partners are bent below about 145°) leaves every fitted energy unchanged and lifts the
  phenyl C–H to 421 and butadiene's inner C–H from 275 to 390, but a stripped ring then keeps one spare
  valence on each carbon, and hydrogens spread over three to six of those at 2.1 Å bind almost a whole
  C–H bond each (Morse at 2.1 Å still gives 0.2–0.3 of the attraction, `satF` counts 0.03–0.1 of it):
  C₆ with six hydrogens hovering over it came out 890 kJ/mol below benzene, and one benzene at 1500 K
  ran straight into it. Without the cap the same geometry relaxes to C₆ + 3 H₂. What works instead is
  to leave the bond orders alone and charge for the second π bond on a bent carbon (Bent π, force-field
  table): phenyl C–H 472, and hot benzene no longer finds those states.
- Softening the Pauli wall between open-shell partners (`pauliOpen` 0.3 / 0.6) lowers Cl + CH₄ from 31
  to 27 / 23 kJ/mol, but lowers every other radical barrier with it (H + CH₄ 40 / 36, CH₃ + ethene 1 / 0)
  and worsens the bond energies; the remaining Cl excess is not a general radical effect. A
  chlorine-only version (0.45 wherever a chlorine meets an open-shell partner) lowered every chlorine
  barrier without deepening the complexes (Cl + CH₄ 31 → 26, Cl + H₂ 37 → 31, H + HCl 33 → 28,
  Cl· + ethene 11 → 2; benchmark 12.8 → 11.6 rms) and made Cl· add to ethene five times as often in
  live runs. It was not kept: the Cl· + ethene path became flat and bumpy, so a skipped addition placed
  at its commit point bounced off (1 of 6 trials instead of 6 of 6), and Cl· + CH₄ fell into the CH₄Cl
  complex (2 of 6). At 0.3 the gain is too small to matter (Cl· + ethene 9).
- Screening weak (below 0.12 of a bond) attractions harder (`ANG_B0` 0.04) stiffens a stretched C–H to
  the methane value, but makes Cl· + ethene's entrance barrier 13 kJ/mol and hot mixtures no better.

Impossible-fragment scans of hot mixtures (5 ps, 22 Å box) found a carbon in CH₄ + O₂ at 3500 K holding
about five bonds' worth of attraction (a double bond, a C–C stretched to 1.9 Å and three C–H at 1.45 Å,
which `satF` counts as only 3.3 used valence) and fragments such as C₂H₇. The σ over-coordination term
(force-field table) now keeps carbon near four: the worst carbon in CH₄ + O₂, C₂H₄ + O₂ and C₆H₆ + O₂ is
4.1–4.2. Still open: in CH₄ + O₂ at 3500 K the excess moves to hydrogen (an H holding 2.3 bonds' worth,
fragments such as H₄O), which the term does not cover, because a hydrogen half-way through a transfer
legitimately holds two partners; and ethane at 2500 K still sheds most of its hydrogen within 5 ps. H₂ + O₂, N₂ + H₂, NH₃, H₂O, CO₂ + H₂ and CH₂O at 2500–3000 K stay
clean.
