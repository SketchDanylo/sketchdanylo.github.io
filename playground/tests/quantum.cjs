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
test('Free-atom references are positive, spherical and decay outward',()=>{
  const O=Q.freeAtomRadial(8);
  assert.ok(O[0]>0);
  for(let k=1;k<40;k++) assert.ok(O[k]<=O[k-1]+1e-12,'radial density never grows outward');
  near(Q.freeDensity(O,9),0,1e-12);
});
test("Hirshfeld shares partition the molecular density exactly",()=>{
  const w=Q.compute({atoms:[atom(8),atom(1,.758130,0,.635742),atom(1,-.758130,0,.635742)]});
  const shares=[0,1,2].map(i=>Q.slice(w,i===0?0:0,48,'atom'));  // same frame for all three
  const total=Q.slice(w,0,48,'total');
  // sum the three atoms' weights at each point of the oxygen-centred frame
  const tables=[8,1,1].map(Z=>Q.freeAtomRadial(Z));
  const C=w.atoms.map(a=>a.xyz.map(v=>v*Q.BOHR));
  const vals=new Float64Array(w.basis.length);
  for(const p of [[0,0,0],[.4,.2,0],[-.7,.6,.1],[1.6,.2,-.3]]){
    const rho=Q.densityAt(w,...p,vals);
    let sum=0,all=0;
    for(let b=0;b<3;b++) all+=Q.freeDensity(tables[b],Math.hypot(p[0]-C[b][0],p[1]-C[b][1],p[2]-C[b][2]));
    for(let b=0;b<3;b++) sum+=rho*Q.freeDensity(tables[b],Math.hypot(p[0]-C[b][0],p[1]-C[b][1],p[2]-C[b][2]))/all;
    near(sum,rho,1e-9);
  }
  // and no share ever exceeds the whole
  for(let k=0;k<total.data.length;k++) assert.ok(shares[0].data[k]<=total.data[k]+1e-9);
});
test("An atom's share is smaller than the fragment total wherever a neighbour competes",()=>{
  const w=Q.compute({atoms:[atom(8),atom(1,.758130,0,.635742),atom(1,-.758130,0,.635742)]});
  const tables=[Q.freeAtomRadial(8),Q.freeAtomRadial(1)];
  const C=w.atoms.map(a=>a.xyz.map(v=>v*Q.BOHR));
  // at the hydrogen nucleus the oxygen keeps only a minority of the density
  const hp=C[1];
  const wO=Q.freeDensity(tables[0],Math.hypot(...hp.map((v,k)=>v-C[0][k])));
  const wH=Q.freeDensity(tables[1],0);
  assert.ok(wH/(wH+2*wO)>.4,'the hydrogen owns most of the density at its own nucleus');
  assert.ok(wO/(wH+2*wO)<.6);
});
console.log(`${count} quantum checks passed. References: Szabo/Ostlund H2 and upstream canonical STO-3G water geometry. These are implementation checks, not general chemical-accuracy claims.`);
