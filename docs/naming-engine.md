# Offline nomenclature engine

`nomenclature.html` is the self-contained deployable engine: styles, JavaScript, the worker source and the RDKit WebAssembly runtime are embedded in this single file. It needs no API key or external runtime assets. PubChem is an optional, explicitly selected database lookup. A database miss is not treated as a chemical error. The readable development sources remain in `src/`; they are not runtime dependencies of the generated HTML.

## Single-file build

Run `node scripts/build-nomenclature.cjs` after editing `src/nomenclature.template.html` or the modules under `src/js/`. The builder embeds the runtime and its BSD license, uses a Blob worker, and instantiates WebAssembly directly from embedded bytes. System fonts avoid external stylesheet/font requests. Only the resulting `nomenclature.html` is required for publication or local use; it is approximately 10 MB. It also works opened directly from disk with the network disabled. Do not edit the generated HTML by hand.

The output is a **systematic name with specified stereochemistry**, not a certified preferred IUPAC name (PIN). The engine assigns R/S and E/Z from explicit configuration labels or wedge/dash bonds. It also derives D/L and alpha/beta annotations for supported amino-acid and sugar units. Unmarked drawings remain stereochemically unspecified. English names are generated directly; German names remain automatic translations.

## Stereochemistry

- Right-click an atom (tap on touch devices) to set R or S, or choose Unspecified. Right-click a single bond for solid/dashed wedge; its **narrow endpoint** is the stereocenter. The reverse button swaps the narrow endpoint. Use one marked bond per center. Right-click a double bond to set E or Z explicitly.
- R/S and E/Z labels specify configuration independently of the layout. Wedges use the drawn geometry and the narrow endpoint. Moving a wedge-connected drawing invalidates its previous result; unmarked coordinates never imply configuration. Tidy converts validated wedges into equivalent R/S labels before rearranging coordinates. Undo/redo and saved structures preserve the annotations.
- `stereochemistry.js` uses RDKit's accurate CIP labeler, validates valence and hydrogen counts, and checks that every requested assignment survives serialization. Contradictory or nonstereogenic labels produce an error. Unspecified stereogenic units produce a notice. `naming-pipeline.js` maps raw atom indices into the hydrogen-folded graph and integrates descriptors into the corresponding parent/substituent names. A conservative descriptor-count check rejects names that omit assigned units.
- `biochemical-stereo.js` measures relative configuration against chemical reference ligands. It does **not** equate L with S: L-cysteine is correctly recognized with an R alpha carbon. D/L labels describe individual alpha-amino-acid units, including peptide residues, rather than assigning an unverified retained name to an entire molecule.
- For ordinary unbranched aldoses through six carbons and 2-ketoses through seven, sugar D/L uses the highest-numbered configurational atom. Alpha/beta uses the anomeric/reference relationship in five- or six-membered O-rings. Complete configurations identify standard aldose families through hexoses and ketoses through hexoses; incomplete configurations yield only the justified series annotation. Substituted glycosides are labeled as configured sugar units. No optical rotation is inferred.
- Optional PubChem lookup receives validated **isomeric SMILES**. A stereo-validation error disables lookup rather than exporting a stripped, potentially different stereoisomer.

RDKit `@rdkit/rdkit` **2026.03.6** is vendored under `src/vendor/rdkit/` with its BSD-3-Clause license for reproducible builds. Its WebAssembly binary is about 7.3 MB before base64 encoding and initializes locally on first naming. The generated HTML supports HTTP(S) and direct `file:` loading. No CDN or online RDKit service is required.

## What changed

- Nitrogen and oxygen linkages are traversed recursively. Substituted amino groups retain their attached carbon chains and rings, even when the amine is not the principal group.
- Every successful name accounts for every heavy atom in the input. Unnamed branches, multiply connected substituents, unrecognized charges and unsupported environments return an explicit unsupported result, never the name of just the recognized fragment.
- Functional groups inside branches, ester alcohol/acyl fragments and amide N-substituents are named recursively.
- Parent selection considers principal groups, chain length, unsaturation and locants. Tied numberings are resolved independently of atom insertion order. Chains are never created by cutting open rings.
- Long carbon stems and multiplying prefixes are generated through 100 carbons instead of relying on the old 30-entry array.
- Bounded graph search, memoization, a worker and a 20-second UI deadline keep large requests from hanging the editor. Graph edits cancel pending work and invalidate previous results. A browser blocking workers uses the main-thread fallback for at most 256 heavy atoms. Moving atoms does not change connectivity or explicitly labeled R/S/E/Z configurations; wedge geometry is significant.
- Formula-only special-name matching has been removed. Ethanol and dimethyl ether remain distinct despite sharing a formula.

## Coverage

Supported combinations include acyclic saturated/unsaturated hydrocarbons; halogens; nitro groups; alcohols; ethers; thiols and thioethers; aldehydes; ketones; carboxylic acids and carboxylates; nitriles; amines; amides; esters with matching ester substituents; acyl halides; acyclic anhydrides; sulfonic acids and sulfonates; and neutral orthophosphate acids/esters.

The ring rules cover carbon monocyclic parents and substituents, benzene, multiple independently linked rings, saturated all-carbon bicycles (including simple fused/bridged bicycles) and two-ring spiro systems. A fixed-numbering library covers 30 heterocyclic parents, including pyridine, pyrimidine, imidazole, furan, thiophene, piperidine, piperazine, morpholine, oxolane and oxane. Carbonyl derivatives of supported saturated heterocycles can also be named.

Biochemical structures are handled through these same structural rules, not by guessing a biological name from a formula. The regression corpus includes 18 amino-acid skeletons, peptides through a 12-residue connectivity example, open-chain and cyclic sugars, a glycosidic disaccharide, fatty esters, and glycerol phosphate. Additional configured fixtures cover amino acids, dipeptides, sugar anomers, a glycoside, esters and nested chiral substituents. A specific stereoisomer requires explicit stereochemical information.

Current explicit limits:

- 1024 atoms including explicit hydrogens for stereo validation, 100 atoms per carbon parent, 64 nested substituent levels, and a bounded search budget. These are ceilings, not guarantees that every structure below them is supported.
- Arbitrary fused aromatic or heterocyclic systems, larger polycycles, steroids, porphyrins, general nucleotides, polymers, coordination/organometallic compounds, most charged nitrogen environments, polyphosphates, mixed ester substituents and isotope naming remain outside coverage.
- Supported stereochemical naming focuses on tetrahedral carbon centers and ordinary E/Z double bonds. Axial/planar/helical chirality, general heteroatom stereochemistry, enhanced stereo groups, racemate notation, and comprehensive pseudoasymmetric/interdependent configurations are not guaranteed. Sugar annotations do not cover general deoxy/amino/branched/multiply bridged sugars or complex multiple configurational prefixes. Unsupported assignments are rejected; these limits must not be described as universal stereochemical coverage.
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
node --test tests/stereochemistry.test.cjs
```

It tests expected names, invalid/unsupported structures, atom accounting, unchanged inputs, search limits, long chains and five shuffled atom/bond orders per positive fixture. Tests deliberately include same-formula isomers and nonalternating rings.

Run browser checks with Playwright available on the Node module search path:

```sh
node tests/naming-browser.cjs
node tests/standalone-naming.cjs
```

The default browser channel is `msedge`; override `BROWSER_CHANNEL` for another installed Chromium channel. Optional `NAMING_SCREENSHOT` and `STEREO_SCREENSHOT` save screenshots. The script hosts the page temporarily on loopback, tests workers, both screenshot regressions, no automatic PubChem request, copy, failed optional lookup, narrow-screen name wrapping, edit cancellation, unsupported output and worker fallback. Stereo checks exercise R/S and E/Z controls, D/L and alpha/beta results, undo/redo, isomeric export, saved structures and tidy preservation.

The standalone check opens the generated file directly with all HTTP(S) requests blocked. It verifies embedded worker/WASM operation, configured molecules, biochemical annotations, and the fallback when workers are blocked. The HTTP browser suite also asserts that no JavaScript or WASM dependencies are fetched.

For independent chemical reconstruction, install [OPSIN](https://github.com/dan2097/opsin/releases/tag/2.9.0) and RDKit **as development tools**, then set `OPSIN_JAR` to the CLI JAR, `CHEM_PYTHON` to a Python executable with RDKit, and `PYTHONPATH` if necessary:

```sh
node tests/opsin-roundtrip.cjs
node tests/stereo-roundtrip.cjs
```

The scripts generate names locally, send them to a local OPSIN process, and compare the resulting structures with the input SMILES using RDKit canonical isomeric SMILES. OPSIN is a name-to-structure parser; it is not used to generate runtime names or as an online fallback. Validation covers 115 connectivity fixtures plus 46 configured reference structures, with 124 core tests and 52 stereo tests. Stereo references are independently built from names in `tests/stereo-reference-names.json`; regenerate the checked-in structures with `node tests/build-stereo-fixtures.cjs` using the same development dependencies. Tests include L-cysteine, meso structures, mixed R/S and E/Z, D/L sugar anomers, peptides, chiral ester fragments, explicit H, wedge inversion, conflicting specifications and reordered atoms. This verifies structural equivalence on the corpus, not universal IUPAC/PIN compliance.

## Rule references

- [IUPAC Blue Book, P-14/P-15: locants, numerical prefixes and name construction](https://iupac.qmul.ac.uk/BlueBook/P1.html).
- [P-2: parent hydrides, ring systems and heterocyclic parents](https://iupac.qmul.ac.uk/BlueBook/P2.html).
- [P-5, especially P-57: substituent naming and complex prefixes](https://iupac.qmul.ac.uk/BlueBook/P5.html).
- [P-6: amines, ethers, hydroxy compounds, carbonyl compounds and acid derivatives](https://iupac.qmul.ac.uk/BlueBook/P6.html).
- [OPSIN documentation: independent conversion from systematic names to structures](https://www.ebi.ac.uk/opsin/).
- [IUPAC/IUBMB amino-acid configuration, 3AA-3](https://iupac.qmul.ac.uk/AminoAcid/AA3t5.html).
- [Carbohydrate D/L configurations, 2-Carb-4](https://iupac.qmul.ac.uk/2carb/03n04.html).
- [Carbohydrate alpha/beta relationships, 2-Carb-6](https://iupac.qmul.ac.uk/2carb/06n07.html).
- [RDKit MinimalLib implementation of stereo tags and CIP assignment](https://github.com/rdkit/rdkit/blob/Release_2026_03/Code/MinimalLib/minilib.cpp).
