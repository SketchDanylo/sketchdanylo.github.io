/* Molecular density from a self-consistent Hartree–Fock calculation.
 * Coordinates into this adapter: Å. Integral engine: atomic units.
 * This is an inspector calculation, independent of the classical MD forces. */
(function (root) {
'use strict';
const BOHR = .529177210903;
const SYMBOLS = ['', 'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P','S','Cl','Ar'];
function installBasis(data) {
  const A = root.App;
  A.tr = key => ({ 'err.rhf.conv': 'Restricted Hartree–Fock did not converge.', 'err.uhf.conv': 'Unrestricted Hartree–Fock did not converge.', 'err.uhf.mult': 'Spin multiplicity is incompatible with the electron count.' }[key] || key);
  for (const [z, atom] of Object.entries(data.elements)) {
    A.BASIS_TABLES['STO-3G'][z] = atom.electron_shells.flatMap(shell => shell.angular_momentum.map((l, i) => ({ l, e: shell.exponents.map(Number), c: shell.coefficients[i].map(Number) })));
    A.SYMBOLS[z] = SYMBOLS[z];
  }
}
function compute(input, progress = () => {}) {
  if (!input || !Array.isArray(input.atoms) || !input.atoms.length || input.atoms.length > 30) throw Error('Choose a molecule with 1–30 atoms.');
  const A = root.App, charge = input.charge ?? 0;
  if (!Number.isInteger(charge)) throw Error('The molecular charge must be an integer.');
  const atoms = input.atoms.map(a => {
    if (!Number.isInteger(a.Z) || a.Z < 1 || a.Z > 18) throw Error('This quantum basis currently supports H–Ar. No substitute cloud is drawn for other elements.');
    if (![a.x,a.y,a.z].every(Number.isFinite)) throw Error('Non-finite atomic position.');
    return { Z: a.Z, xyz: [a.x/BOHR, a.y/BOHR, a.z/BOHR] };
  });
  const nelec = atoms.reduce((s, a) => s + a.Z, 0) - charge;
  const mult = input.mult || (nelec % 2 ? 2 : 1);
  if (!Number.isInteger(mult) || mult < 1 || nelec < 1 || mult - 1 > nelec || (nelec - mult + 1) % 2) throw Error('Spin multiplicity is incompatible with the electron count.');
  for (let i = 0; i < atoms.length; i++) for (let j = i+1; j < atoms.length; j++) {
    if (Math.hypot(...atoms[i].xyz.map((v,k) => v-atoms[j].xyz[k])) < .15) throw Error('Nuclei are too close for a stable calculation.');
  }
  const basis = A.buildBasis(atoms, 'STO-3G');
  if (basis.length > 48) throw Error('This inspector currently allows 48 basis functions per calculation.');
  if ((nelec + mult - 1) / 2 > basis.length) throw Error('Electron count exceeds the basis capacity.');
  progress('One-electron integrals');
  const ints = A.integrals.oneElectron(basis, atoms);
  const eri = A.eri.computeERI(basis, f => progress('Electron repulsion · ' + Math.round(f*100) + '%'));
  progress('Solving molecular orbitals');
  const scf = mult > 1 ? A.uhf.runUHF(atoms,basis,ints,eri,charge,mult) : A.scf.runRHF(atoms,basis,ints,eri,charge);
  const electronTrace = A.linalg.trace2(scf.D, ints.S, basis.length);
  if (!scf.converged || !Number.isFinite(scf.E) || Math.abs(electronTrace-nelec) > 1e-5) throw Error('The calculation failed its electron-count check.');
  return { atoms, basis, scf, mult, electronTrace, method: mult > 1 ? 'UHF / STO-3G' : 'RHF / STO-3G' };
}
// rho(r) = sum_mu,nu D_mu,nu chi_mu(r) chi_nu(r), in electrons / Å³.
function densityAt(result, x, y, z, values = new Float64Array(result.basis.length)) {
  x /= BOHR; y /= BOHR; z /= BOHR;
  const basis = result.basis, n = basis.length, D = result.scf.D;
  for (let i = 0; i < n; i++) {
    const f = basis[i], dx=x-f.center[0], dy=y-f.center[1], dz=z-f.center[2], r2=dx*dx+dy*dy+dz*dz;
    let v=0; for (let k=0;k<f.exps.length;k++) v+=f.coefs[k]*Math.exp(-f.exps[k]*r2);
    values[i] = v * dx**f.l * dy**f.m * dz**f.n;
  }
  let rho=0;
  for (let i=0;i<n;i++) {
    rho+=D[i*n+i]*values[i]*values[i];
    for(let j=0;j<i;j++) rho+=2*D[i*n+j]*values[i]*values[j];
  }
  return Math.max(0,rho/BOHR**3);
}
// Individual molecular orbital amplitude in Å^(-3/2); its square integrates to one.
function orbitalAt(result, orbital, x, y, z, spin='alpha') {
  const n=result.basis.length, C=spin==='beta' && result.scf.uhf ? result.scf.CB : result.scf.C;
  if(!Number.isInteger(orbital)||orbital<0||orbital>=n) throw Error('Invalid orbital index.');
  x/=BOHR; y/=BOHR; z/=BOHR;
  let psi=0;
  for(let i=0;i<n;i++) {
    const f=result.basis[i],dx=x-f.center[0],dy=y-f.center[1],dz=z-f.center[2],r2=dx*dx+dy*dy+dz*dz;
    let v=0; for(let k=0;k<f.exps.length;k++) v+=f.coefs[k]*Math.exp(-f.exps[k]*r2);
    psi+=C[i*n+orbital]*v*dx**f.l*dy**f.m*dz**f.n;
  }
  return psi/BOHR**1.5;
}
function orbitalSlice(result, selected, orbital, plane='xy', spin='alpha', resolution=160) {
  const axes={xy:[0,1],xz:[0,2],yz:[1,2]}[plane];
  if(!axes) throw Error('Invalid slice plane.');
  const origin=result.atoms[selected].xyz.map(v=>v*BOHR);
  const extent=Math.max(2.5,...result.atoms.map(a=>Math.hypot(...a.xyz.map((v,k)=>v*BOHR-origin[k]))+2));
  const data=new Float32Array(resolution**2); let max=0;
  for(let y=0;y<resolution;y++) for(let x=0;x<resolution;x++) {
    const p=origin.slice(); p[axes[0]]+=((x+.5)/resolution*2-1)*extent; p[axes[1]]+=((y+.5)/resolution*2-1)*extent;
    const v=orbitalAt(result,orbital,...p,spin); data[y*resolution+x]=v; max=Math.max(max,Math.abs(v));
  }
  return {data,resolution,extent,origin,max,mode:'orbital',plane,axes,orbital,spin};
}
/* One atom's share of the molecular cloud, by Hirshfeld's stockholder rule:
 *   w_A(r) = rho_A_free(|r - R_A|) / sum_B rho_B_free(|r - R_B|),   rho_A(r) = w_A(r) rho_mol(r)
 * The molecular density already carries every neighbour, near and far, at the geometry it was
 * given, so the share drawn for one atom is that atom as its surroundings have made it.
 * The free-atom references are spherically averaged ground-state atoms in the same basis. */
const FREE_MULT=[0,2,1,2,1,2,3,4,3,2,1,2,1,2,3,4,3,2,1];
const FREE_MAX=8, FREE_N=320;                       // Å, radial samples
const freeCache=new Map();
function freeAtomRadial(Z) {
  if (freeCache.has(Z)) return freeCache.get(Z);
  const atom=compute({atoms:[{Z,x:0,y:0,z:0}],mult:FREE_MULT[Z]||(Z%2?2:1)});
  // spherical average over the 26 directions of a 3x3x3 stencil
  const dirs=[];
  for(let a=-1;a<=1;a++) for(let b=-1;b<=1;b++) for(let c=-1;c<=1;c++) {
    if(!a&&!b&&!c) continue;
    const n=Math.hypot(a,b,c); dirs.push([a/n,b/n,c/n]);
  }
  const table=new Float64Array(FREE_N+1), values=new Float64Array(atom.basis.length);
  for(let k=0;k<=FREE_N;k++) {
    const r=k/FREE_N*FREE_MAX;
    let sum=0;
    for(const d of dirs) sum+=densityAt(atom,d[0]*r,d[1]*r,d[2]*r,values);
    table[k]=sum/dirs.length;
  }
  freeCache.set(Z,table);
  return table;
}
function freeDensity(table, r) {
  if (r >= FREE_MAX) return 0;
  const t=r/FREE_MAX*FREE_N, i=t|0, f=t-i;
  return table[i]+(table[i+1]-table[i])*f;
}
function hirshfeldWeights(result) { return result.atoms.map(a=>freeAtomRadial(a.Z)); }
/* mode 'atom' draws the selected atom's share; 'total' draws the whole fragment's density. */
function slice(result, selected, resolution = 160, mode = 'atom') {
  const atoms=result.atoms, origin=atoms[selected].xyz.map(v=>v*BOHR);
  const extent=Math.max(2, ...atoms.map(a=>Math.hypot(a.xyz[0]*BOHR-origin[0],a.xyz[1]*BOHR-origin[1])+1.5));
  const data=new Float32Array(resolution*resolution), values=new Float64Array(result.basis.length);
  const tables=mode==='atom'?hirshfeldWeights(result):null;
  const centres=atoms.map(a=>a.xyz.map(v=>v*BOHR));
  let max=0;
  for(let y=0;y<resolution;y++) for(let x=0;x<resolution;x++) {
    const px=origin[0]+((x+.5)/resolution*2-1)*extent, py=origin[1]+((y+.5)/resolution*2-1)*extent, pz=origin[2];
    let rho=densityAt(result,px,py,pz,values);
    if (tables) {
      let mine=0, all=0;
      for(let b=0;b<centres.length;b++) {
        const c=centres[b], d=freeDensity(tables[b],Math.hypot(px-c[0],py-c[1],pz-c[2]));
        all+=d; if(b===selected) mine=d;
      }
      rho = all>1e-12 ? rho*(mine/all) : 0;
    }
    data[y*resolution+x]=rho; max=Math.max(max,rho);
  }
  return {data,resolution,extent,origin,max,mode};
}
root.ChemQuantum={installBasis,compute,densityAt,orbitalAt,orbitalSlice,slice,freeAtomRadial,freeDensity,BOHR};
if(typeof module==='object') module.exports=root.ChemQuantum;
})(globalThis);
