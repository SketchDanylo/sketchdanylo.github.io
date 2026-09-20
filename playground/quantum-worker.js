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
    self.postMessage({status:'Free-atom references'});
    const atomSlice=ChemQuantum.slice(result,data.selected,160,'atom');
    self.postMessage({status:'Sampling electron density'});
    const totalSlice=ChemQuantum.slice(result,data.selected,160,'total');
    self.postMessage({result:{method:result.method,scf:result.scf,mult:result.mult,electronTrace:result.electronTrace,atomSlice,totalSlice}},
      [atomSlice.data.buffer,totalSlice.data.buffer]);
  } catch(error) { self.postMessage({error:error.message}); }
};
