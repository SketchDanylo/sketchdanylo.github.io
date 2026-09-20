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
function slice(result, selected, resolution = 160) {
  const atoms=result.atoms, origin=atoms[selected].xyz.map(v=>v*BOHR);
  const extent=Math.max(2, ...atoms.map(a=>Math.hypot(a.xyz[0]*BOHR-origin[0],a.xyz[1]*BOHR-origin[1])+1.5));
  const data=new Float32Array(resolution*resolution), values=new Float64Array(result.basis.length);
  let max=0;
  for(let y=0;y<resolution;y++) for(let x=0;x<resolution;x++) {
    const rho=densityAt(result,origin[0]+((x+.5)/resolution*2-1)*extent,origin[1]+((y+.5)/resolution*2-1)*extent,origin[2],values);
    data[y*resolution+x]=rho; max=Math.max(max,rho);
  }
  return {data,resolution,extent,origin,max};
}
root.ChemQuantum={installBasis,compute,densityAt,slice,BOHR};
if(typeof module==='object') module.exports=root.ChemQuantum;
})(globalThis);
