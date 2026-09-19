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

It is a compact bond-order reactive force field in the Tersoff / Brenner / ReaxFF family. It was written for this app and fitted to measured data.

| Term | Form | Data |
| --- | --- | --- |
| Covalent | Morse `De[e^(−2y) − 2·b·e^(−y)]`, with re, De and a interpolated in the continuous bond order n | Measured bond lengths, dissociation energies and stretch force constants (≈80 pairs). Pauling/Pyykkö estimates fill the rest |
| Saturation | `b = p_i·p_j`, `p = 1/(1 + 1.3x + 20x⁴)`, x = valence already used + 1 − valence | Tuned to the H + H₂ and H + CH₄ barriers and the CH₄ dimer |
| Bond order | Spare valence shared between neighbours → C=C, C≡C, O=O, N≡N, CO₂, aromatic 1.5 | — |
| Angles | VSEPR, with θ₀ a smooth function of the continuous steric number | 109.5°, 107°, 104.5°, 120°, 180° |
| Non-bonded | Shielded Lennard-Jones (UFF) with its Pauli part switched off between atoms that can bond. Shifted-force Coulomb between saturating bond-polarisation charges | UFF, Pauling electronegativity |
| Walls | Soft harmonic container. The normal force on the walls is the displayed pressure | — |
| Heat bath | Bussi CSVR thermostat (canonical); switch it off for an isolated box | — |

Every force is the exact gradient of the energy, including the many-body saturation, screening, charge and VSEPR terms. The test suite checks this numerically. The only lagged quantity is the bond multiplicity n, which relaxes over about 12 fs.

## Validation

`node playground/tests/physics-check.cjs` (add `--strict` to fail on any miss):

| Check | Model | Experiment |
| --- | --- | --- |
| H–H, H₂ bond energy | 0.741 Å, 436 kJ/mol | 0.741 Å, 436 |
| H₂O O–H, angle | 0.960 Å, 104.4° | 0.958 Å, 104.5° |
| CH₄ C–H, angle, atomisation | 1.089 Å, 109.47°, 1654 kJ/mol | 1.09 Å, 109.47°, 1652 |
| C–C / C=C / C≡C | 1.55 / 1.34 / 1.20 Å | 1.54 / 1.34 / 1.20 |
| Benzene C–C | 1.42 Å | 1.39 |
| O=O, N≡N energies | 498, 945 kJ/mol | 498, 945 |
| Water dimer | −16 kJ/mol, O···O 2.98 Å | −21 (De), 2.91 Å |
| Methane dimer | −2.7 kJ/mol | −2.2 |
| H + H₂ → H₂ + H barrier | 38 kJ/mol | ≈ 40 |
| H + CH₄ → H₂ + CH₃ barrier | 52 kJ/mol | ≈ 50 |
| NVE energy drift, 2 ps | < 2 kJ/mol | — |

Emergent chemistry you can reproduce: H and O atoms form H₂, O₂, OH· and H₂O. Hot Na + Cl nucleates NaCl clusters. At 3000–6000 K bonds break, and the reaction feed logs every event.

**Limits.** This is a classical reactive model, not quantum chemistry. Spin states are not modelled; O₂ is closed-shell here, not the real triplet. CO is double-bonded rather than triple-bonded. Transition metals and hypervalent geometries are approximate. Radical additions to π bonds have barriers that are too high. Elements without bonding parameters cannot be imported.

## Nomenclature bridge

**Playground** in Nomenclature exports `chem-playground/molecule@1`. It contains atoms in Å with every implicit hydrogen made explicit, bonds and orders, and formal charges. Wedge and dash bonds lift the wide end to ±0.9 Å in z. The molecule is delivered through a same-origin `BroadcastChannel` and a shared `localStorage` inbox. If no Playground tab is open, it opens one with the molecule in the URL hash. You can also paste the JSON into the Playground.

In the Playground, a molecule is first **conditioned** in a spherical cell:

- **Room**: 298 K and 1 atm.
- **Absolute zero**: energy-minimised and motionless.
- **Custom**: any T and P. The cell radius follows the ideal-gas volume per molecule, so only high pressure squeezes it.

Then you place it: click to drop it, drag to throw it (the arrow shows speed and kinetic energy), Q/E to rotate, Shift-click to place several.

## Controls

Press `?` for every shortcut. Click a key to rebind it; conflicts move automatically, and bindings persist. Letter bindings follow the physical key, so they work on any keyboard layout. `Ctrl+K` opens a command palette that also understands `500 K`, `80 °C` and `2x`.
