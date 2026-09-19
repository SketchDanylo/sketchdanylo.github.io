/* Optional independent chemistry check. OPSIN is a name-to-structure parser,
 * NOT a runtime naming dependency. Run with OPSIN_JAR, CHEM_PYTHON and
 * PYTHONPATH (RDKit) configured. All processing is local.
 */
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const { spawnSync } = require('node:child_process');
const { fixtures, graph } = require('./naming-fixtures.cjs');
const { name } = require('../src/js/native-namer.js');
if (!process.env.OPSIN_JAR) throw new Error('Set OPSIN_JAR to a local OPSIN CLI JAR.');
const dir = fs.mkdtempSync(path.join(os.tmpdir(), 'naming-roundtrip-'));
const examples = fixtures.map(([label, smiles]) => {
  const result = name(graph(smiles));
  if (result.status !== 'ok') throw new Error(label + ': ' + result.reason);
  return { label, smiles, name: result.name };
});
fs.writeFileSync(path.join(dir, 'names.txt'), examples.map(x => x.name).join('\n') + '\n');
fs.writeFileSync(path.join(dir, 'fixtures.json'), JSON.stringify(examples));
const opsin = spawnSync('java', ['-jar', process.env.OPSIN_JAR, '-o', 'smi', path.join(dir, 'names.txt'), path.join(dir, 'structures.smi')], { encoding: 'utf8' });
if (opsin.status) throw new Error(opsin.stderr);
if (opsin.stderr) process.stderr.write(opsin.stderr);
const compare = spawnSync(process.env.CHEM_PYTHON || 'python', ['-c', `
import json,sys
from pathlib import Path
from rdkit import Chem
p=Path(sys.argv[1])
fixtures=json.loads((p/'fixtures.json').read_text())
lines=(p/'structures.smi').read_text().splitlines()
assert len(lines)==len(fixtures), (len(lines),len(fixtures))
failures=[]
for x,actual in zip(fixtures,lines):
    expected_mol=Chem.MolFromSmiles(x['smiles'])
    actual_mol=Chem.MolFromSmiles(actual) if actual.strip() else None
    if actual_mol is None or expected_mol is None or Chem.MolToSmiles(expected_mol)!=Chem.MolToSmiles(actual_mol):
        failures.append(dict(x,actual=actual))
print(json.dumps({'checked':len(fixtures),'failures':failures},indent=2))
sys.exit(bool(failures))
`, dir], { encoding: 'utf8' });
process.stdout.write(compare.stdout || '');
process.stderr.write(compare.stderr || '');
process.exitCode = compare.status ?? 1;
