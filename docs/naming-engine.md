# Offline nomenclature engine

`nomenclature.html` uses `src/js/native-namer.js`, a pure JavaScript structure-to-name engine, through `src/js/naming-worker.js`. It ships with the page, needs no API key, and performs no network requests. PubChem is an optional, explicitly selected database lookup. A database miss is not treated as a chemical error.

The output is a **systematic connectivity name**, not a certified preferred IUPAC name (PIN). The editor does not specify tetrahedral stereochemistry, and the engine does not assign R/S, E/Z, D/L or alpha/beta descriptors. The previous coordinate-based CIP approximation has been removed from naming and export because it could falsely imply stereochemical certainty. English names are generated directly; German names remain automatic translations.

## What changed

- Nitrogen and oxygen linkages are traversed recursively. Substituted amino groups retain their attached carbon chains and rings, even when the amine is not the principal group.
- Every successful name accounts for every heavy atom in the input. Unnamed branches, multiply connected substituents, unrecognized charges and unsupported environments return an explicit unsupported result, never the name of just the recognized fragment.
- Functional groups inside branches, ester alcohol/acyl fragments and amide N-substituents are named recursively.
- Parent selection considers principal groups, chain length, unsaturation and locants. Tied numberings are resolved independently of atom insertion order. Chains are never created by cutting open rings.
- Long carbon stems and multiplying prefixes are generated through 100 carbons instead of relying on the old 30-entry array.
- Bounded graph search, memoization, a worker and an 8-second UI deadline keep large requests from hanging the editor. Graph edits cancel pending work and invalidate previous results. Moving atoms does not change connectivity names.
- Formula-only special-name matching has been removed. Ethanol and dimethyl ether remain distinct despite sharing a formula.

## Coverage

Supported combinations include acyclic saturated/unsaturated hydrocarbons; halogens; nitro groups; alcohols; ethers; thiols and thioethers; aldehydes; ketones; carboxylic acids and carboxylates; nitriles; amines; amides; esters with matching ester substituents; acyl halides; acyclic anhydrides; sulfonic acids and sulfonates; and neutral orthophosphate acids/esters.

The ring rules cover carbon monocyclic parents and substituents, benzene, multiple independently linked rings, saturated all-carbon bicycles (including simple fused/bridged bicycles) and two-ring spiro systems. A fixed-numbering library covers 30 heterocyclic parents, including pyridine, pyrimidine, imidazole, furan, thiophene, piperidine, piperazine, morpholine, oxolane and oxane. Carbonyl derivatives of supported saturated heterocycles can also be named.

Biochemical structures are handled through these same structural rules, not by guessing a biological name from a formula. The regression corpus includes 18 amino-acid skeletons, peptides through a 12-residue example, open-chain and cyclic sugars, a glycosidic disaccharide, fatty esters, and glycerol phosphate. These are connectivity names: the engine cannot distinguish a specific stereoisomer of a sugar or amino acid without stereochemical information.

Current explicit limits:

- 1024 heavy atoms, 100 atoms per carbon parent, 64 nested substituent levels, and a bounded search budget. These are ceilings, not guarantees that every structure below them is supported.
- Arbitrary fused aromatic or heterocyclic systems, larger polycycles, steroids, porphyrins, general nucleotides, polymers, coordination/organometallic compounds, most charged nitrogen environments, polyphosphates, mixed ester substituents, isotope labels and stereochemical naming remain outside coverage.
- Multiple disconnected fragments and salts need a separate naming method.
- Unattached explicit hydrogens are rejected before creating the heavy-atom view. Attached explicit hydrogens are folded into hydrogen counts.
- A successful full-atom coverage check prevents omitted fragments; it is not by itself a proof that every IUPAC rule has been implemented.

## Screenshot regressions

1. `CC(NCOCO)CCC(Cl)C` → `({[(5-chlorohexan-2-yl)amino]methyl}oxy)methanol`
2. `CC(NC1CCCCC1)C(O)C(=O)C` → `4-(cyclohexylamino)-3-hydroxypentan-2-one`

The second example previously lost all six cyclohexyl carbons from the name. Both examples now pass independent name-to-structure reconstruction.

## Validation

Run the dependency-free engine suite:

```sh
node --test tests/native-namer.test.cjs
```

It tests expected names, invalid/unsupported structures, atom accounting, unchanged inputs, search limits, long chains and five shuffled atom/bond orders per positive fixture. Tests deliberately include same-formula isomers and nonalternating rings.

Run browser checks with Playwright available on the Node module search path:

```sh
node tests/naming-browser.cjs
```

The default browser channel is `msedge`; override `BROWSER_CHANNEL` for another installed Chromium channel. Optional `NAMING_SCREENSHOT` saves a screenshot. The script hosts the page temporarily on loopback, tests workers, both screenshot regressions, no automatic PubChem request, copy, failed optional lookup, narrow-screen name wrapping, edit cancellation, unsupported output and worker fallback.

For independent chemical reconstruction, install [OPSIN](https://github.com/dan2097/opsin/releases/tag/2.9.0) and RDKit **as development tools**, then set `OPSIN_JAR` to the CLI JAR, `CHEM_PYTHON` to a Python executable with RDKit, and `PYTHONPATH` if necessary:

```sh
node tests/opsin-roundtrip.cjs
```

The script generates names locally, sends them to a local OPSIN process, and compares the resulting structures with the input SMILES using RDKit canonical SMILES. OPSIN is a name-to-structure parser; it is not used to generate names or as an online runtime fallback. At implementation time all 115 positive fixtures passed this independent check, and all 124 engine tests passed. This verifies structural equivalence on the corpus, not universal IUPAC/PIN compliance.

## Rule references

- [IUPAC Blue Book, P-14/P-15: locants, numerical prefixes and name construction](https://iupac.qmul.ac.uk/BlueBook/P1.html).
- [P-2: parent hydrides, ring systems and heterocyclic parents](https://iupac.qmul.ac.uk/BlueBook/P2.html).
- [P-5, especially P-57: substituent naming and complex prefixes](https://iupac.qmul.ac.uk/BlueBook/P5.html).
- [P-6: amines, ethers, hydroxy compounds, carbonyl compounds and acid derivatives](https://iupac.qmul.ac.uk/BlueBook/P6.html).
- [OPSIN documentation: independent conversion from systematic names to structures](https://www.ebi.ac.uk/opsin/).
