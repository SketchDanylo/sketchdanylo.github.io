# Chem Playground

A reactive molecular-dynamics sandbox at true scale. Atoms and molecules live in a 3D slab, 1.2 nm deep by default, which you view from above, so reactions cannot drift away along z. The view zooms from 0.05 nm to 5 nm. Bonds are never scripted: they form and break from the potential energy surface.

Open `playground/` on the site. Molecules come only from [Nomenclature](../nomenclature.html) → **Playground**; single atoms come from the dock.

## Time

- Every step is exactly **1 fs**: velocity Verlet in a render-decoupled accumulator (`while (acc >= 1) { step(); acc -= 1 }`), drawn with interpolation between the last two steps.
- **1.0× = 20 000 steps/s = 20 ps of simulated time per real second.** The scale is exponential below 1×: 0.1× is one step per second. Above 1× it is linear, e.g. 2× = 40 000 steps/s. Type any value into the speed field, including ones outside the slider range.
- When a frame cannot fit all the steps, the rate readout says **CPU-bound** and shows the actual speed.
- A step that runs into extreme curvature splits itself into up to 16 sub-steps. This happens when a hydrogen is squeezed between competing bonds, or in a hard collision at thousands of kelvin. Either the curvature from the previous step triggers it, or a per-step energy check does. The step still advances exactly 1 fs.
- **Step back** is exact. The engine keeps a checkpoint every 64 steps and replays deterministically, bit-identical, because the RNG state is part of the checkpoint. Hold the button to rewind.

## Force field

A compact bond-order reactive force field in the Tersoff / Brenner / ReaxFF family, written for this app
and fitted to measured data.

| Term | Form | Data |
| --- | --- | --- |
| Covalent | Morse `De[e^(−2y) − 2·b·e^(−y)]`, with re, De and a interpolated in the continuous bond order n | Measured bond lengths, dissociation energies and stretch force constants (≈80 pairs); Pauling/Pyykkö estimates fill the rest |
| Valence used | Each bond uses `s = e^(−1.3·a·Δr)` of both partners' valence — it decays like the Morse attraction it still provides — times its multiplicity | — |
| Saturation | `b = p_i·p_j`, `p = 1/(1 + 0.6x + 10.8x⁶)`, x = excess valence if this bond were fully formed | Fitted to the barriers below while keeping closed-shell dimers non-sticky |
| Evans–Polanyi | Excess valence is weighted by `√(De_competing / De_new)`, so a stronger incoming bond displaces a weaker one more easily. An atom at its normal valence is unaffected | Fitted to F + H₂ and H + Cl₂ |
| Screening | Two atoms bonded to a common neighbour do not compete for each other's valence, unless they are bonded themselves | — |
| Bond order | Spare valence shared between neighbours → C=C, C≡C, O=O, N≡N, CO₂, aromatic 1.5. O₂'s π bond counts as only 0.78 of oxygen's valence (triplet O₂ is a diradical), and an O–O bond with one unpaired oxygen gains 0.45 order (the three-electron bond of HO₂·) | O₂ and HO₂ chemistry |
| Angles | VSEPR, θ₀ a smooth function of the continuous steric number | 109.5°, 107°, 104.5°, 120°, 180° |
| Non-bonded | Shielded Lennard-Jones (UFF); its Pauli wall fades where the Morse term already repels, and between atoms that can still bond. Shifted-force Coulomb between saturating bond-polarisation charges | UFF, Pauling electronegativity |
| Walls | Soft harmonic container; the normal force on the walls is the displayed pressure | — |
| Heat bath | Bussi CSVR thermostat (canonical); switch it off for an isolated box | — |

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

**Why nothing happens at room temperature.** This is real kinetics, not a bug. A second of simulation
covers tens of nanoseconds, while a hydrogen–oxygen mixture at 298 K takes years to ignite by itself.
To see a textbook reaction, either raise the temperature (1500–4000 K), or start the chain yourself by
placing single atoms — radicals react immediately.

**Limits.** A classical reactive model, not quantum chemistry. Spin is absent: O₂ is patched to behave
like the triplet diradical it is, but the general case is not. Some barriers are off by 20–40 kJ/mol
(Cl + H₂ and OH + H₂ too high, H + O₂ → OH + O too low), so *relative rates are qualitative*. CO is
modelled with a double instead of a triple bond, so carbon-monoxide energetics are wrong. Transition
metals, hypervalent geometry, tunnelling, excited states and solvent chemistry are not represented.

## Nomenclature bridge

**Playground** in Nomenclature exports `chem-playground/molecule@1`. It contains atoms in Å with every implicit hydrogen made explicit, bonds and orders, and formal charges. Wedge and dash bonds lift the wide end to ±0.9 Å in z. The molecule is delivered through a same-origin `BroadcastChannel` and a shared `localStorage` inbox. If no Playground tab is open, it opens one with the molecule in the URL hash. You can also paste the JSON into the Playground.

In the Playground, a molecule is first **conditioned** in a spherical cell:

- **Room**: 298 K and 1 atm.
- **Absolute zero**: energy-minimised and motionless.
- **Custom**: any T and P. The cell radius follows the ideal-gas volume per molecule, so only high pressure squeezes it.

Then you place it: click to drop it, drag to throw it (the arrow shows speed and kinetic energy), Q/E to rotate, Shift-click to place several.

## Controls

Press `?` for every shortcut. Click a key to rebind it; conflicts move automatically, and bindings persist. Letter bindings follow the physical key, so they work on any keyboard layout. `Ctrl+K` opens a command palette that also understands `500 K`, `80 °C` and `2x`.
