/* Runs only local, vendored numerical code. No molecular data leaves the device. */
importScripts('vendor/qchem/basis.js?v=20260922-dropfix5','vendor/qchem/linalg.js?v=20260922-dropfix5','vendor/qchem/integrals.js?v=20260922-dropfix5','vendor/qchem/eri.js?v=20260922-dropfix5','vendor/qchem/scf.js?v=20260922-dropfix5','vendor/qchem/uhf.js?v=20260922-dropfix5','quantum-core.js?v=20260922-dropfix5');
const ready = fetch('vendor/qchem/sto-3g-h-ar.json').then(r => {
  if (!r.ok) throw Error('Could not load the quantum basis.');
  return r.json();
}).then(ChemQuantum.installBasis);
let calculation=null;
self.onmessage = async ({data}) => {
  try {
    await ready;
    if(data.view) {
      if(!calculation) throw Error('Calculate a fragment first.');
      const v=data.view, grid=ChemQuantum.orbitalSlice(calculation,v.selected,v.orbital,v.plane,v.spin);
      self.postMessage({grid,requestId:data.requestId},[grid.data.buffer]); return;
    }
    const result = calculation = ChemQuantum.compute(data.input, status => self.postMessage({status}));
    self.postMessage({status:'Free-atom references'});
    const atomSlice=ChemQuantum.slice(result,data.selected,160,'atom');
    self.postMessage({status:'Sampling electron density'});
    const totalSlice=ChemQuantum.slice(result,data.selected,160,'total');
    self.postMessage({result:{method:result.method,scf:result.scf,mult:result.mult,electronTrace:result.electronTrace,atomSlice,totalSlice}},
      [atomSlice.data.buffer,totalSlice.data.buffer]);
  } catch(error) { self.postMessage({error:error.message}); }
};
