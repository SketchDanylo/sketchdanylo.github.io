// Optional regeneration from independently parsed reference names (OPSIN 2.9).
const fs=require('node:fs'),path=require('node:path'),os=require('node:os'),{spawnSync}=require('node:child_process');
const names=require('./stereo-reference-names.json'),dir=fs.mkdtempSync(path.join(os.tmpdir(),'stereo-references-'));
fs.writeFileSync(path.join(dir,'names.txt'),names.map(x=>x[0]).join('\n')+'\n');
const parsed=spawnSync('java',['-jar',process.env.OPSIN_JAR,'-o','smi',path.join(dir,'names.txt'),path.join(dir,'out.smi')],{encoding:'utf8'});
if(parsed.status)throw Error(parsed.stderr);
const smiles=fs.readFileSync(path.join(dir,'out.smi'),'utf8').trim().split(/\r?\n/);
if(smiles.length!==names.length||smiles.some(s=>!s))throw Error(parsed.stderr);
fs.writeFileSync(path.join(dir,'inputs.json'),JSON.stringify(names.map((n,i)=>({reference:n[0],series:n[1],kind:n[2],anomeric:n[3],family:n[4],smiles:smiles[i]}))));
const py=spawnSync(process.env.CHEM_PYTHON||'python',['-c',`
import json,sys
from rdkit import Chem
from rdkit.Chem import rdDepictor,rdCIPLabeler
rows=json.load(open(sys.argv[1],encoding='utf-8'))
for row in rows:
 m=Chem.MolFromSmiles(row['smiles']); Chem.Kekulize(m,clearAromaticFlags=True); rdDepictor.Compute2DCoords(m); rdCIPLabeler.AssignCIPLabels(m)
 c=m.GetConformer(); atoms=[]; bonds=[]
 for a in m.GetAtoms():
  p=c.GetAtomPosition(a.GetIdx()); d=dict(el=a.GetSymbol(),charge=a.GetFormalCharge(),h=a.GetTotalNumHs(),x=p.x*50,y=-p.y*50)
  if a.HasProp('_CIPCode'): d['cip']=a.GetProp('_CIPCode')
  atoms.append(d)
 for b in m.GetBonds():
  d=dict(a=b.GetBeginAtomIdx(),b=b.GetEndAtomIdx(),order=int(b.GetBondTypeAsDouble()))
  if str(b.GetStereo()) in ('STEREOE','STEREOTRANS'): d['stereo']='E'
  if str(b.GetStereo()) in ('STEREOZ','STEREOCIS'): d['stereo']='Z'
  bonds.append(d)
 row['model']=dict(atoms=atoms,bonds=bonds)
json.dump(rows,open(sys.argv[2],'w'),indent=2)
`,path.join(dir,'inputs.json'),path.join(__dirname,'stereo-fixtures.json')],{encoding:'utf8'});
if(py.status)throw Error(py.stderr);
console.log('Built '+names.length+' independent reference structures.');
