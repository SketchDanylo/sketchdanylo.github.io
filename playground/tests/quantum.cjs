const assert=require('node:assert/strict');
for(const file of ['basis','linalg','integrals','eri','scf','uhf']) require('../vendor/qchem/'+file+'.js');
const Q=require('../quantum-core.js');
Q.installBasis(require('../vendor/qchem/sto-3g-h-ar.json'));
const near=(a,b,tol)=>assert.ok(Math.abs(a-b)<tol,`${a} vs ${b}`);
const atom=(Z,x=0,y=0,z=0)=>({Z,x,y,z});
let count=0;
function test(label,run){run();console.log('PASS',label);count++;}
let h2;
test('H2 RHF/STO-3G energy matches textbook reference',()=>{
  h2=Q.compute({atoms:[atom(1,0,0,-.3704),atom(1,0,0,.3704)]});
  near(h2.scf.E,-1.1167,5e-4); near(h2.electronTrace,2,1e-9);
});
test('Water RHF/STO-3G energy matches the reference geometry',()=>{
  const w=Q.compute({atoms:[atom(8),atom(1,.758130,0,.635742),atom(1,-.758130,0,.635742)]});
  near(w.scf.E,-74.9659,.002); near(w.electronTrace,10,1e-8);
  near(w.scf.mulliken.reduce((a,b)=>a+b,0),0,1e-8);
});
test('Real-space H2 density integrates to the molecular electron count',()=>{
  let total=0;const h=.2,values=new Float64Array(2);
  for(let x=-6;x<=6;x+=h)for(let y=-6;y<=6;y+=h)for(let z=-6;z<=6;z+=h) total+=Q.densityAt(h2,x,y,z,values)*h**3;
  near(total,2,.006);
});
test('Density responds to bond geometry and preserves spatial symmetry',()=>{
  const stretched=Q.compute({atoms:[atom(1,0,0,-1),atom(1,0,0,1)]});
  assert.ok(Q.densityAt(h2,0,0,0)>Q.densityAt(stretched,0,0,0)*2);
  near(Q.densityAt(h2,.2,.3,.4),Q.densityAt(h2,-.2,-.3,-.4),1e-12);
});
test('Chlorine radical uses an open-shell molecular calculation',()=>{
  const cl=Q.compute({atoms:[atom(17)],mult:2});
  near(cl.electronTrace,17,1e-7);near(cl.scf.S2,.75,1e-7);
  assert.equal(cl.scf.uhf,true);assert.ok(Number.isFinite(Q.densityAt(cl,.4,0,0)));
});
test('Invalid charges, spins, unsupported elements and overlapping nuclei fail explicitly',()=>{
  for(const input of [{atoms:[atom(1)],charge:.3},{atoms:[atom(1)],mult:1},{atoms:[atom(35)]},{atoms:[atom(1),atom(1)]}]) assert.throws(()=>Q.compute(input));
});
console.log(`${count} quantum checks passed. References: Szabo/Ostlund H2 and upstream canonical STO-3G water geometry. These are implementation checks, not general chemical-accuracy claims.`);
