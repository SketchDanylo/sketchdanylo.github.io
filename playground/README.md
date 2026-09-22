# Chem Playground

A reactive molecular-dynamics sandbox at true scale. Atoms and molecules live in a 3D slab, 1.2 nm deep by default, which you view from above, so reactions cannot drift away along z. The manual view zooms from 0.05 nm to 15 nm. Chamber dimensions extend to 50 nm; Fit can zoom farther out to show the complete box. In the sandbox, bonds form and break from an experimental potential energy model.

Open `playground/` on the site. Molecules come only from [Nomenclature](../nomenclature.html) → **Playground**; single atoms come from the dock.

## Thermostats

Three ways to hold a temperature, chosen in **Environment** or from the temperature gauge. The
difference between them is who gets touched and when.

### Kelvin stat

All unpinned velocities receive the same factor `sqrt(K_target / K)` each step, after void-wall absorption. This fixes the **total kinetic temperature**, not each atom's energy. Startup from rest initializes velocities and records the injected energy in `kelvinWork`.

This is a numerical isokinetic control, not a physical wall bath. It preserves instantaneous speed ratios but changes relative velocities, suppresses kinetic-energy fluctuations, and can alter reaction dynamics. An exact temperature readout does not establish accurate reaction rates or canonical sampling. Pinned atoms stay fixed, and work accounting survives replay. See [velocity-rescaling thermostats](https://docs.lammps.org/fix_temp_rescale.html).

It shares heat out as well as fixing the total. A single global factor preserves whatever ratio
the last event happened to leave, which sounds harmless and is not: a molecule that has just
formed carries the whole of its new bond's energy, gets scaled hard for it, and then stays frozen
while the rest of the chamber runs warm — nothing in a uniform rescale can ever warm it again.
Each atom is first pulled gently toward its own share, `kelvinMix` per fs, and the exact total is
imposed after. Measured on an H₂ formed in a chamber of argon, per-atom energy relative to the
spectators: **0.48 before, 0.97 after.** Set `kelvinMix: 0` for the old pure rescale.

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

**The wall is measured, not set.** `wallMeasured` is the kinetic temperature of the fluid lying
within `wallSkin` of a face — the wall is whatever the sample against it is, and a heater only
sets what it *aims* for. With a temperature void the heater stands down entirely, because it
cannot warm a wall that radiates everything away, and the wall simply reports the fluid; that is
what the gauge and the Environment panel show. To hold an open chamber at a fixed temperature,
use the **Kelvin stat**, which replaces what the void removes on the same tick.

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
| Bond order | Spare valence shared between neighbours → C=C, C≡C, O=O, N≡N, CO₂, aromatic 1.5. It follows the same long-ranged, screened bond order the saturation uses, so a partner weakens a π bond *while* its own bond forms; a partner claims valence only insofar as it can bond at all, so a radical frees the π bond and a saturated molecule drifting past does not | Radical addition to alkenes |
| π blocking | A π bond blocks an incoming partner at 0.85 of a σ bond, being weaker and more polarizable. O₂'s π counts 0.78 (triplet O₂ is a diradical) and an O–O bond with one unpaired oxygen gains 0.45 order (the three-electron bond of HO₂·) | Cl + ethene; O₂ and HO₂ chemistry |
| Angles | VSEPR, θ₀ a smooth function of the continuous steric number | 109.5°, 107°, 104.5°, 120°, 180° |
| Non-bonded | Shielded Lennard-Jones (UFF); its Pauli wall fades where the Morse term already repels, and between atoms that can still bond. Shifted-force Coulomb between saturating bond-polarisation charges | UFF, Pauling electronegativity |
| Charge scale | 0.38 e per unit electronegativity difference per bond, calibrated on condensed water rather than the gas-phase dimer — the same deliberate over-polarisation the TIP3P family uses. Covalent bonds, NaCl and HCl are unaffected | Water cohesion, O–O distance |
| Walls | Specular hard reflection or optional soft harmonic field; pressure from normal momentum transfer | Idealized boundaries |
| Wall reservoir | Local OU coupling at the boundary; finite heat capacity, gradual heater response | Model parameters, not a calibrated apparatus |

Every force is the exact gradient of the energy, including the many-body saturation, screening, charge and
VSEPR terms; `physics-check.cjs` verifies this numerically (errors ≈1e-6 kJ/mol/Å). The only lagged
quantity is the bond multiplicity n, which relaxes over about 12 fs.

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

## Temperature readings

The gauge shows a running mean of the kinetic temperature with the measured spread beside it. The instantaneous value genuinely wanders — a sample of N atoms has relative fluctuations of about `sqrt(2/3N)`, so twenty atoms swing by tens of kelvin from one instant to the next. That is the sample being small, not friction or a numerical fault. Hovering the gauge gives the mean, the spread, the instantaneous value and the expected fluctuation for the current atom count.

## Controls

Press `?` for every shortcut. Click a key to rebind it; conflicts move automatically, and bindings persist. Letter bindings follow the physical key, so they work on any keyboard layout. `Ctrl+K` opens a command palette that also understands `500 K`, `80 °C` and `2x`.

## Additional checks and model limitations

Run `node playground/tests/regression.cjs` for force-query invariance, charge conservation, variable-radius force gradients, hot rewind across checkpoints, pinned atoms, failure handling, and malformed imports.

Bond multiplicity relaxes as an internal heuristic variable. Fixed-state force gradient checks do not prove total energy conservation during reactions. The heat bath exchanges energy; the engine also counts safety velocity clamps in extreme collisions. Neither behavior should be mistaken for isolated, rigorously conservative dynamics. The UI reports those clamps when they occur.

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
