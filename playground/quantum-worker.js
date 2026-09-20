/* Runs only local, vendored numerical code. No molecular data leaves the device. */
importScripts('vendor/qchem/basis.js','vendor/qchem/linalg.js','vendor/qchem/integrals.js','vendor/qchem/eri.js','vendor/qchem/scf.js','vendor/qchem/uhf.js','quantum-core.js');
const ready = fetch('vendor/qchem/sto-3g-h-ar.json').then(r => {
  if (!r.ok) throw Error('Could not load the quantum basis.');
  return r.json();
}).then(ChemQuantum.installBasis);
self.onmessage = async ({data}) => {
  try {
    await ready;
    const result = ChemQuantum.compute(data.input, status => self.postMessage({status}));
    self.postMessage({status:'Sampling electron density'});
    const slice=ChemQuantum.slice(result,data.selected);
    self.postMessage({result:{method:result.method,scf:result.scf,mult:result.mult,electronTrace:result.electronTrace,slice}},[slice.data.buffer]);
  } catch(error) { self.postMessage({error:error.message}); }
};
