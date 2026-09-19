const assert = require('node:assert/strict');
// Deliberately small test-only Kekule SMILES reader. Fixtures are independent
// structures, not calls to the naming algorithm's own group detection.
function graph(smiles) {
  const atoms = [], stack = [], rings = new Map();
  let at = -1, order = 1;
  const bond = (a, b, o) => { atoms[a].nb.push({ n: b, o }); atoms[b].nb.push({ n: a, o }); };
  const tokens = smiles.match(/\[[^\]]+\]|Cl|Br|[A-Z][a-z]?|[()=#.]|%\d\d|\d/g);
  assert.equal(tokens.join(''), smiles, 'Unsupported test SMILES: ' + smiles);
  for (const token of tokens) {
    if (token === '(') { stack.push(at); continue; }
    if (token === ')') { at = stack.pop(); continue; }
    if (token === '=') { order = 2; continue; }
    if (token === '#') { order = 3; continue; }
    if (token === '.') { at = -1; continue; }
    if (/^%?\d/.test(token)) {
      if (rings.has(token)) { const [other, previousOrder] = rings.get(token); bond(at, other, Math.max(order, previousOrder)); rings.delete(token); }
      else rings.set(token, [at, order]);
      order = 1; continue;
    }
    const el = token.match(/[A-Z][a-z]?/)[0];
    const charge = token.includes('+') ? 1 : token.includes('-') ? -1 : 0;
    const h = token.startsWith('[') ? Number(token.match(/H(\d*)/)?.[1] || (token.includes('H') ? 1 : 0)) : null;
    const next = atoms.length;
    atoms.push({ el, charge, h, nb: [] });
    if (at >= 0) bond(at, next, order);
    at = next; order = 1;
  }
  assert.equal(rings.size, 0);
  atoms.forEach(a => {
    if (a.h !== null) return;
    const sum = a.nb.reduce((s, e) => s + e.o, 0);
    const vals = a.el === 'C' ? [4] : a.el === 'N' ? [3] : a.el === 'O' ? [2] : a.el === 'S' ? [2, 4, 6] : [1];
    a.h = Math.max(0, (vals.find(v => v >= sum) || 0) - sum);
  });
  return atoms;
}
function permute(g, seed) {
  let state = seed;
  const random = () => (state = (Math.imul(state, 1664525) + 1013904223) >>> 0) / 4294967296;
  const order = g.map((_, i) => i);
  for (let i = order.length - 1; i > 0; i--) { const j = Math.floor(random() * (i + 1)); [order[i], order[j]] = [order[j], order[i]]; }
  const idx = new Map(order.map((old, next) => [old, next]));
  return order.map(i => ({ ...g[i], nb: [...g[i].nb].reverse().map(e => ({ n: idx.get(e.n), o: e.o })) }));
}
const fixtures = [
  ['methane', 'C', 'methane'], ['ethane', 'CC', 'ethane'], ['ethene', 'C=C', 'ethene'],
  ['ethyne', 'C#C', 'ethyne'], ['branch', 'CC(C)CC', '2-methylbutane'],
  ['alcohol', 'CCO', 'ethanol'], ['ketone', 'CC(=O)C', 'propan-2-one'],
  ['acid', 'CCC(=O)O', 'propanoic acid'], ['dicarboxylic acid', 'OC(=O)CC(=O)O', 'propanedioic acid'],
  ['aldehyde', 'CCC=O', 'propanal'], ['dialdehyde', 'O=CCC=O', 'propanedial'],
  ['nitrile', 'CCC#N', 'propanenitrile'], ['amide', 'CC(N)=O', 'ethanamide'],
  ['N-methylamide', 'CC(=O)NC', 'N-methylethanamide'],
  ['N,N-dimethylamide', 'CC(=O)N(C)C', 'N,N-dimethylethanamide'],
  ['secondary amine', 'CCNCC', 'N-ethylethanamine'],
  ['tertiary amine', 'CCN(CC)CC', 'N,N-diethylethanamine'],
  ['ether', 'COCC', 'methoxyethane'], ['diol', 'OCCO', 'ethane-1,2-diol'],
  ['thiol', 'CCS', 'ethanethiol'], ['thioether', 'CSCC', '(methylsulfanyl)ethane'],
  ['carboxylate', 'CC(=O)[O-]', 'ethanoate'],
  ['nitro', 'C[N+](=O)[O-]', 'nitromethane'], ['acyl chloride', 'CC(=O)Cl', 'ethanoyl chloride'],
  ['sulfonic acid', 'CS(=O)(=O)O', 'methanesulfonic acid'],
  ['ester', 'CC(=O)OCC', 'ethyl ethanoate'],
  ['branched acyl ester', 'CC(C)C(=O)OC', 'methyl 2-methylpropanoate'],
  ['halogenated acyl ester', 'ClCC(=O)OCC', 'ethyl 2-chloroethanoate'],
  ['hydroxy ester', 'CC(O)C(=O)OC', 'methyl 2-hydroxypropanoate'],
  ['dimethyl diester', 'COC(=O)CC(=O)OC', 'dimethyl propanedioate'],
  ['cyclohexane', 'C1CCCCC1', 'cyclohexane'], ['benzene', 'C1=CC=CC=C1', 'benzene'],
  ['phenol', 'OC1=CC=CC=C1', 'phenol'], ['benzoic acid', 'OC(=O)C1=CC=CC=C1', 'benzoic acid'],
  ['cyclohexanol', 'OC1CCCCC1', 'cyclohexan-1-ol'],
  ['phenyl ethanol', 'OCCC1=CC=CC=C1', '2-phenylethanol'],
  ['two isolated rings', 'C1CCCCC1C2CCCCC2', 'cyclohexylcyclohexane'],
  ['three isolated rings', 'C1CCCCC1N(C2CCCCC2)C3CCCCC3', 'N,N-dicyclohexylcyclohexan-1-amine'],
  ['pyridine', 'N1=CC=CC=C1', 'pyridine'], ['pyrimidine', 'N1=CN=CC=C1', 'pyrimidine'],
  ['imidazole', 'N1C=NC=C1', '1H-imidazole'], ['furan', 'O1C=CC=C1', 'furan'],
  ['piperidine', 'N1CCCCC1', 'piperidine'], ['morpholine', 'O1CCNCC1', 'morpholine'],
  ['N-methylpiperidine', 'CN1CCCCC1', '1-methylpiperidine'],
  ['nicotinic acid skeleton', 'OC(=O)C1=CC=CN=C1', 'pyridine-3-carboxylic acid'],
  ['glycine', 'NCC(=O)O', '2-aminoethanoic acid'],
  ['alanine', 'CC(N)C(=O)O', '2-aminopropanoic acid'],
  ['valine', 'CC(C)C(N)C(=O)O', '2-amino-3-methylbutanoic acid'],
  ['leucine', 'CC(C)CC(N)C(=O)O', '2-amino-4-methylpentanoic acid'],
  ['serine', 'NC(CO)C(=O)O', '2-amino-3-hydroxypropanoic acid'],
  ['cysteine', 'NC(CS)C(=O)O', '2-amino-3-sulfanylpropanoic acid'],
  ['methionine', 'CSCCC(N)C(=O)O', '2-amino-4-(methylsulfanyl)butanoic acid'],
  ['phenylalanine', 'NC(CC1=CC=CC=C1)C(=O)O', '2-amino-3-phenylpropanoic acid'],
  ['lysine', 'NCCCCC(N)C(=O)O', '2,6-diaminohexanoic acid'],
  ['aspartic acid', 'NC(CC(=O)O)C(=O)O', '2-aminobutanedioic acid'],
  ['proline', 'OC(=O)C1CCCN1', 'pyrrolidine-2-carboxylic acid'],
  ['glycerol', 'OCC(O)CO', 'propane-1,2,3-triol'],
  ['open-chain hexose', 'O=CC(O)C(O)C(O)C(O)CO', '2,3,4,5,6-pentahydroxyhexanal'],
  ['sugar ring', 'OCC1OC(O)C(O)C(O)C1O', '6-(hydroxymethyl)oxane-2,3,4,5-tetrol'],
  ['citric acid skeleton', 'OC(=O)CC(O)(CC(=O)O)C(=O)O', '3-carboxy-3-hydroxypentanedioic acid'],
  ['glycylglycine', 'NCC(=O)NCC(=O)O', '2-[(2-aminoethanoyl)amino]ethanoic acid'],
  ['image 1', 'CC(NCOCO)CCC(Cl)C', '({[(5-chlorohexan-2-yl)amino]methyl}oxy)methanol'],
  ['image 2', 'CC(NC1CCCCC1)C(O)C(=O)C', '4-(cyclohexylamino)-3-hydroxypentan-2-one'],
  ['water', 'O', 'water'], ['ammonia', 'N', 'ammonia'],
  ['anhydride', 'CC(=O)OC(=O)C', 'ethanoic anhydride'],
  ['mixed anhydride', 'CCC(=O)OC(=O)C', 'ethanoic propanoic anhydride'],
  ['branched anhydride', 'CC(C)C(=O)OC(=O)C(C)C', '2-methylpropanoic anhydride'],
  ['fused bicyclic', 'C1CCC2CCCCC2C1', 'bicyclo[4.4.0]decane'],
  ['bridged bicyclic', 'C1CC2CCC1C2', 'bicyclo[2.2.1]heptane'],
  ['spiro bicycle', 'C1CCC2(CC1)CCCC2', 'spiro[4.5]decane'],
  ['phosphoric acid', 'OP(=O)(O)O', 'phosphoric acid'],
  ['phosphate ester', 'COP(=O)(O)O', 'methyl dihydrogen phosphate'],
  ['trialkyl phosphate', 'COP(=O)(OC)OC', 'trimethyl phosphate'],
  ['glycerol phosphate', 'OCC(O)COP(=O)(O)O', '(2,3-dihydroxypropyl) dihydrogen phosphate'],
  ['piperidone', 'O=C1CCCCN1', 'piperidin-2-one'],
  ['histidine skeleton', 'NC(CC1=CN=CN1)C(=O)O', '2-amino-3-(1H-imidazol-5-yl)propanoic acid'],
  ['isoleucine', 'CCC(C)C(N)C(=O)O', '2-amino-3-methylpentanoic acid'],
  ['threonine', 'CC(O)C(N)C(=O)O', '2-amino-3-hydroxybutanoic acid'],
  ['glutamic acid', 'OC(=O)CCC(N)C(=O)O', '2-aminopentanedioic acid'],
  ['tyrosine skeleton', 'NC(CC1=CC=C(O)C=C1)C(=O)O', '2-amino-3-(4-hydroxyphenyl)propanoic acid'],
  ['asparagine', 'NC(CC(N)=O)C(=O)O'],
  ['glutamine', 'NC(CCC(N)=O)C(=O)O'],
  ['12-residue peptide', 'NCC(=O)'.repeat(12) + 'O'],
  ['glycosidic disaccharide', 'OCC1OC(OC2OC(CO)C(O)C(O)C2O)C(O)C(O)C1O'],
  ['fatty ester', 'CCCCCCCCCCCCCCCC(=O)OCC(O)CO'],
  ['31-carbon chain', 'C'.repeat(31), 'hentriacontane'],
  ['100-carbon chain', 'C'.repeat(100), 'hectane'],
  ['pyridazine', 'N1N=CC=CC=1', 'pyridazine'],
  ['pyrazine', 'N1=CC=NC=C1', 'pyrazine'],
  ['pyrrole', 'N1C=CC=C1', '1H-pyrrole'],
  ['thiophene', 'S1C=CC=C1', 'thiophene'],
  ['oxazole', 'O1C=NC=C1', '1,3-oxazole'],
  ['thiazole', 'S1C=NC=C1', '1,3-thiazole'],
  ['piperazine', 'N1CCNCC1', 'piperazine'],
  ['oxolane', 'O1CCCC1', 'oxolane'],
  ['oxirane', 'O1CC1', 'oxirane'],
  ['unsaturated ketone', 'CC(=O)C=CC', 'pent-3-en-2-one'],
  ['en-yne locants', 'C=CCC#C', 'pent-1-en-4-yne'],
  ['principal chain before unsaturation', 'CCC(=C)CCC', '3-methylidenehexane'],
  ['alphabetical numbering tie', 'CCC(Cl)C(Br)CC', '3-bromo-4-chlorohexane'],
  ['same-formula isomer ethanol', 'CCO', 'ethanol'],
  ['same-formula isomer dimethyl ether', 'COC', 'methoxymethane'],
  ['carbon dioxide', 'O=C=O', 'carbon dioxide'],
  ['hydrogen peroxide', 'OO', 'hydrogen peroxide'],
  ['dioxygen', 'O=O', 'dioxygen'],
  ['dinitrogen', 'N#N', 'dinitrogen'],
  ['hydrogen chloride', 'Cl', 'hydrogen chloride'],
  ['sulfuric acid', 'OS(=O)(=O)O', 'sulfuric acid'],
  ['hydrogen sulfate', 'OS(=O)(=O)[O-]', 'hydrogen sulfate'],
  ['sulfate', '[O-]S(=O)(=O)[O-]', 'sulfate'],
  ['nitric acid', 'O[N+](=O)[O-]', 'nitric acid'],
  ['nitrate', '[O-][N+](=O)[O-]', 'nitrate'],
  ['numbered N-substituents', 'CNC(=O)CCC(=O)NC', 'N^1,N^4-dimethylbutanediamide'],
];

module.exports = { graph, permute, fixtures };
