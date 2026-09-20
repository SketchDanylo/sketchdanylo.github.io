# Chem Playground

A reactive molecular-dynamics sandbox at true scale. Atoms and molecules live in a 3D slab, 1.2 nm deep by default, which you view from above, so reactions cannot drift away along z. The manual view zooms from 0.05 nm to 15 nm. Chamber dimensions extend to 50 nm; Fit can zoom farther out to show the complete box. In the sandbox, bonds form and break from an experimental potential energy model.

Open `playground/` on the site. Molecules come only from [Nomenclature](../nomenclature.html) → **Playground**; single atoms come from the dock.

## Boundary thermostat

The laboratory uses a **spatial wall thermostat**, not a global velocity rescale:

- Heater targets are **15–350 °C**. Changing the target does not jump the sample or wall temperature.
- The wall follows `dT/dt = (T_target - T_wall) / tau`; its response time is adjustable in **simulated ps**, default 10 ps. This specifies a nanoscale thermal reservoir, not the heating time of a particular laboratory appliance.
- Only atoms within **0.2 nm of a wall** couple to the reservoir. The coupling vanishes quadratically at the inner edge. An exact Ornstein–Uhlenbeck velocity update combines drag and Gaussian noise with fluctuation–dissipation variance. Interparticle forces carry energy through the interior.
- Effective wall heat capacity is 1000 k_B. Heat given to the sample is subtracted from the wall; heater work and sample heat are recorded separately and included in deterministic snapshots. This finite reservoir approximation is not an explicit solid-wall atomistic model.
- With the thermostat **off**, editing temperature immediately rescales kinetic energy to the requested temperature, including an initialization from rest. Pinned atoms remain fixed. Subsequent dynamics have no thermostat.
- Isolated molecule conditioning retains a separate CSVR sampling bath.

[Spatial Langevin thermostat methodology](https://docs.lammps.org/fix_langevin.html). This engine uses its own exact OU update, not the LAMMPS integration algorithm.

## Bounds and the void wall

The **Environment** app in the console sets what the chamber face does.

- **Bounds — Solid**: the default stiff harmonic wall, 60 kJ/mol/Å², beginning exactly at the face. An atom bounces on contact.
- **Bounds — Forcefield**: a softer harmonic, 6 kJ/mol/Å², shifted 3 Å inside the face and continuing outward without limit. An atom is turned around gradually and can lean past the face, but a finite energy can never escape a harmonic that keeps growing.
- **Void wall — Temperature**: an absorbing boundary. Beyond the face, kinetic energy drains with a 40 fs time constant that reaches full strength 2 Å out. The energy is tallied in `voidHeat` and never returned, so the sample cools at its edges. This is an open system: total energy is deliberately not conserved.
- **Void wall — Pressure**: the face still applies its restoring force — atoms are pushed back exactly as before — but the impulse is booked to `voidForce` instead of `wallForce`, so the pressure gauge reads what an open chamber would read.

Both void options can be on at once, and both are part of snapshots, of the deterministic replay used by step-back, and of the saved scene. `node playground/tests/bounds.cjs` checks containment, energy conservation with no void selected, drainage only outside the face, absorbed pressure, and replay.

## Time

- Every step is exactly **1 fs**: velocity Verlet in a render-decoupled accumulator (`while (acc >= 1) { step(); acc -= 1 }`), drawn with interpolation between the last two steps.
- **1.0× = 20 000 steps/s = 20 ps of simulated time per real second.** The scale is exponential below 1×: 0.1× is one step per second. Above 1× it is linear, e.g. 2× = 40 000 steps/s. Type any value into the speed field, including ones outside the slider range.
- When a frame cannot fit all the steps, the rate readout says **CPU-bound** and shows the actual speed.
- A step that runs into extreme curvature splits itself into up to 16 sub-steps. This happens when a hydrogen is squeezed between competing bonds, or in a hard collision at thousands of kelvin. Either the curvature from the previous step triggers it, or a per-step energy check does. The step still advances exactly 1 fs.
- **Step back** is exact. The engine keeps a checkpoint every 64 steps and replays deterministically, bit-identical, because the RNG state is part of the checkpoint. Hold the button to rewind.

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
| Walls | Soft harmonic container; the normal force on the walls is the displayed pressure | — |
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
| Water dimer | −18 kJ/mol, O···O 2.93 Å | −21 (De), 2.91 Å |
| Methane dimer | −3.8 kJ/mol | −2.2 |
| NVE energy drift, 2 ps | < 2 kJ/mol | — |
| Step back ×2 | bit-identical | — |

`node playground/tests/reactions.cjs [--dyn]` — reaction barriers from relaxed minimum-energy paths,
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
| CH₄ + 2 O₂ → CO₂ + 2 H₂O | −756 | −802 |
| C₂H₄ + H₂ → C₂H₆ | −127 | −136 |
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

Three apps: **Environment** (bounds, void wall, chamber size, thermostat, and a live account of what the boundary has absorbed), **About** (the model notes below), and **Keybinds**.

## Seeing the third dimension

The view is orthographic and looks straight down z by default. **Right-drag empty space** turns the slab about its own centre; **double right-click** faces it again. While the view is turned, the chamber is drawn as a wireframe with its near edges bright, and every pointer gesture — dragging an atom, placing one, the heat brush, the tweezer — works in the plane you are looking at, because pointer positions are mapped back through the same rotation. Chamber-edge resizing is disabled while the view is turned, since the faces are no longer screen-aligned.

## The atom inspector

Right-click an atom. The card is pinned beside it with a leader line drawn on the field, and follows the atom as the camera moves.

**Electron cloud** shows that one atom's share of the fragment's calculated density, by Hirshfeld's stockholder rule:

    w_A(r) = rho_A_free(|r - R_A|) / sum_B rho_B_free(|r - R_B|),    rho_A(r) = w_A(r) rho_mol(r)

The free-atom references are ground-state atoms in the same STO-3G basis, spherically averaged over 26 directions onto a radial table. The molecular density being divided is a real self-consistent Hartree–Fock solution for the whole fragment at the geometry it has at that instant, so the shape shown is that atom as its neighbours — bonded and not — and the present thermal geometry have made it. **Whole fragment** switches back to the undivided density. Both are computed in one worker pass, so the toggle is instant.

Hirshfeld shares sum exactly to the molecular density and no share can exceed it; `node playground/tests/quantum.cjs` checks both.

## Temperature readings

The gauge shows a running mean of the kinetic temperature with the measured spread beside it. The instantaneous value genuinely wanders — a sample of N atoms has relative fluctuations of about `sqrt(2/3N)`, so twenty atoms swing by tens of kelvin from one instant to the next. That is the sample being small, not friction or a numerical fault. Hovering the gauge gives the mean, the spread, the instantaneous value and the expected fluctuation for the current atom count.

## Controls

Press `?` for every shortcut. Click a key to rebind it; conflicts move automatically, and bindings persist. Letter bindings follow the physical key, so they work on any keyboard layout. `Ctrl+K` opens a command palette that also understands `500 K`, `80 °C` and `2x`.

## Additional checks and model limitations

Run `node playground/tests/regression.cjs` for force-query invariance, charge conservation, variable-radius force gradients, hot rewind across checkpoints, pinned atoms, failure handling, and malformed imports.

Bond multiplicity relaxes as an internal heuristic variable. Fixed-state force gradient checks do not prove total energy conservation during reactions. The heat bath exchanges energy; the engine also counts safety velocity clamps in extreme collisions. Neither behavior should be mistaken for isolated, rigorously conservative dynamics. The UI reports those clamps when they occur.

The camera is orthographic: depth changes shading, not apparent atomic radii, and turning the view changes nothing physical. Chamber contours are summed Gaussians, not electron-density calculations; the inspector's cloud is a real Hartree–Fock calculation and the only quantum result in the app. Playback speed is a nonlinear UI setting, not a multiplier of physical real time. With a void wall selected the chamber is an open system by design, and total energy is not conserved.

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
