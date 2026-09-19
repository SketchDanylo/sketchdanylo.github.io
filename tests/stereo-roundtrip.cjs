// Independent name-to-structure verification, including isomeric identity.
const fs=require('node:fs'),os=require('node:os'),path=require('node:path'),{spawnSync}=require('node:child_process');
const P=require('../src/js/naming-pipeline.js'),fixtures=require('./stereo-fixtures.json');
function heavy(m){const g=m.atoms.map((a,i)=>({el:a.el,h:a.h,charge:a.charge,orig:i,nb:[]}));m.bonds.forEach(b=>{g[b.a].nb.push({n:b.b,o:b.order});g[b.b].nb.push({n:b.a,o:b.order});});return g;}
(async()=>{
 const rows=[];for(const f of fixtures){const r=await P.name(heavy(f.model),f.model);if(r.status!=='ok')throw Error(f.reference+': '+r.reason);rows.push({...f,name:r.name});}
 const dir=fs.mkdtempSync(path.join(os.tmpdir(),'stereo-roundtrip-'));
 fs.writeFileSync(path.join(dir,'names.txt'),rows.map(x=>x.name).join('\n')+'\n');fs.writeFileSync(path.join(dir,'rows.json'),JSON.stringify(rows));
 const p=spawnSync('java',['-jar',process.env.OPSIN_JAR,'-o','smi',path.join(dir,'names.txt'),path.join(dir,'out.smi')],{encoding:'utf8'});if(p.status)throw Error(p.stderr);
 const check=spawnSync(process.env.CHEM_PYTHON||'python',['-c',`
import json,sys
from pathlib import Path
from rdkit import Chem
p=Path(sys.argv[1]);rows=json.loads((p/'rows.json').read_text());lines=(p/'out.smi').read_text().splitlines();failures=[]
assert len(rows)==len(lines)
for row,smi in zip(rows,lines):
 m=Chem.MolFromSmiles(smi) if smi else None
 if m is None or Chem.MolToSmiles(m)!=Chem.MolToSmiles(Chem.MolFromSmiles(row['smiles'])):failures.append(dict(reference=row['reference'],name=row['name'],actual=smi))
print(json.dumps(dict(checked=len(rows),failures=failures),indent=2));sys.exit(bool(failures))
`,dir],{encoding:'utf8'});process.stdout.write(check.stdout);process.stderr.write(check.stderr);process.exitCode=check.status;
})();
