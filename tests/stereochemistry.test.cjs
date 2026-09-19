const test=require('node:test'),assert=require('node:assert/strict');
const S=require('../src/js/stereochemistry.js'),B=require('../src/js/biochemical-stereo.js'),P=require('../src/js/naming-pipeline.js');
const fixtures=require('./stereo-fixtures.json');
function heavy(model){
  const graph=model.atoms.map((a,i)=>({el:a.el,charge:a.charge,h:a.h,orig:i,nb:[]}));
  model.bonds.forEach(b=>{graph[b.a].nb.push({n:b.b,o:b.order});graph[b.b].nb.push({n:b.a,o:b.order});});return graph;
}
for(const f of fixtures)test(f.reference,async()=>{
  const input=JSON.stringify(f.model),s=await S.assign(f.model);
  assert.equal(s.ok,true,s.reason);
  for(const [i,a]of f.model.atoms.entries())if(a.cip)assert.equal(s.atoms.find(x=>x.atom===i)?.label,a.cip);
  for(const b of f.model.bonds)if(b.stereo)assert.equal(s.bonds.find(x=>x.a===b.a&&x.b===b.b||x.b===b.a&&x.a===b.b)?.label,b.stereo);
  const bio=B.describe(f.model,s);
  if(f.kind){const label=bio.labels.find(x=>x.kind===f.kind);assert(label,JSON.stringify(bio));assert.equal(label.series,f.series);if(f.kind==='sugar'){assert.equal(label.anomeric,f.anomeric);assert.equal(label.family,f.family);}}
  const named=await P.name(heavy(f.model),f.model);assert.equal(named.status,'ok',named.reason);
  assert.equal(JSON.stringify(f.model),input,'mutated input');
});
test('unspecified centers and alkenes remain unspecified',async()=>{
  for(const index of [0,9,19]){const m=structuredClone(fixtures[index].model);m.atoms.forEach(a=>delete a.cip);m.bonds.forEach(b=>delete b.stereo);
    const s=await S.assign(m);assert(s.ok,s.reason);assert.equal(s.atoms.length+s.bonds.length,0);assert(s.warnings.length);assert.equal(B.describe(m,s).labels.length,0);
  }
});
test('rejects labels on nonstereogenic atoms and bonds',async()=>{
  const m=structuredClone(fixtures[9].model);m.atoms[0].cip='R';assert.equal((await S.assign(m)).ok,false);
  delete m.atoms[0].cip;m.bonds[0].stereo='E';assert.equal((await S.assign(m)).ok,false);
});
test('wedge inversion, reversed endpoint and conflicting annotations',async()=>{
  const m=structuredClone(fixtures[0].model),center=m.atoms.findIndex(a=>a.cip);m.atoms.forEach(a=>delete a.cip);
  const bond=m.bonds.find(b=>b.a===center||b.b===center);if(bond.b===center)[bond.a,bond.b]=[bond.b,bond.a];
  bond.stereo='wedge';const a=await S.assign(m);assert(a.ok,a.reason);bond.stereo='hash';const b=await S.assign(m);assert(b.ok,b.reason);assert.notEqual(a.atoms[0].label,b.atoms[0].label);
  m.atoms[center].cip=a.atoms[0].label;assert.equal((await S.assign(m)).ok,false,'conflicting wedge accepted');delete m.atoms[center].cip;
  [bond.a,bond.b]=[bond.b,bond.a];assert.equal((await S.assign(m)).ok,false,'nonstereogenic wedge endpoint accepted');
});
test('explicit hydrogen retains R/S and biochemical reference parity',async()=>{
  const m=structuredClone(fixtures[0].model),center=m.atoms.findIndex(a=>a.cip),i=m.atoms.length;
  m.atoms[center].h=0;m.atoms.push({el:'H',charge:0,h:0,x:m.atoms[center].x+20,y:m.atoms[center].y+25});m.bonds.push({a:center,b:i,order:1});
  const s=await S.assign(m);assert(s.ok,s.reason);assert.equal(s.atoms.find(x=>x.atom===center).label,'S');assert.equal(B.describe(m,s).labels[0].series,'L');
});
test('atom and bond reordering preserve names and biochemical assignments',async()=>{
  for(const index of [0,2,8,11,19,35,39,41,44]){
    const m=structuredClone(fixtures[index].model),n=m.atoms.length;
    const reordered={atoms:m.atoms.toReversed(),bonds:m.bonds.toReversed().map(b=>({...b,a:n-1-b.a,b:n-1-b.b}))};
    const a=await P.name(heavy(m),m),b=await P.name(heavy(reordered),reordered);
    assert.equal(b.status,'ok',b.reason);assert.equal(a.name,b.name);assert.deepEqual(a.biochemical.labels.map(x=>[x.kind,x.series,x.anomeric,x.family]),b.biochemical.labels.map(x=>[x.kind,x.series,x.anomeric,x.family]).reverse());
  }
});
test('all connectivity regression structures survive stereo validation',async()=>{
  const{fixtures:rows,graph}=require('./naming-fixtures.cjs');
  for(const[label,smiles]of rows){const g=graph(smiles),m={atoms:g.map(a=>({el:a.el,h:a.h,charge:a.charge})),bonds:g.flatMap((a,i)=>a.nb.filter(e=>e.n>i).map(e=>({a:i,b:e.n,order:e.o})))};
    const r=await P.name(g,m);assert.equal(r.status,'ok',label+': '+r.reason);assert.equal(r.stereo.atoms.length+r.stereo.bonds.length,0,label+' guessed stereochemistry');
  }
});
module.exports={heavy};
