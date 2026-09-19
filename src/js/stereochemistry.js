/* RDKit 2026.03.6 (BSD-3-Clause), bundled locally. No network service.
 * Explicit wedges/hash bonds and requested R/S or E/Z configurations only.
 * Plain 2D layout is never interpreted as a stereochemical specification.
 */
(function(root,factory){
  const api=factory();
  if(typeof module==='object'&&module.exports) module.exports=api;
  else root.Stereochemistry=api;
})(typeof globalThis!=='undefined'?globalThis:this,function(){
  'use strict';
  let runtime;
  async function init(){
    if(!runtime){
      if(typeof module==='object'&&module.exports) runtime=require('../vendor/rdkit/dist/RDKit_minimal.js')();
      else {
        const embedded=globalThis.NAMING_WASM_BASE64||(typeof document!=='undefined'?document.getElementById('naming-wasm')?.textContent.trim():null);
        if(embedded){
          const bytes=Uint8Array.from(atob(embedded),c=>c.charCodeAt(0));
          // This MinimalLib build exposes instantiateWasm, but does not consume
          // the generic Emscripten wasmBinary option. Instantiate bytes directly.
          runtime=new Promise((resolve,reject)=>{
            initRDKitModule({instantiateWasm(imports,receive){
              WebAssembly.instantiate(bytes,imports).then(result=>receive(result.instance)).catch(reject);
              return {};
            }}).then(resolve,reject);
          });
        }else runtime=initRDKitModule({locateFile:file=>new URL('../vendor/rdkit/dist/'+file,
          typeof document!=='undefined'?new URL('src/js/stereochemistry.js',document.baseURI).href:self.location.href).href});
      }
    }
    return runtime;
  }
  const key=(a,b)=>Math.min(a,b)+':'+Math.max(a,b);
  const strip=s=>s.replace(/[()]/g,'');
  function molblock(input){
    const {atoms,bonds}=input;
    const lines=['Nomenclature stereo','  NativeNamer       2D','','  0  0  0     0  0            999 V3000',
      'M  V30 BEGIN CTAB',`M  V30 COUNTS ${atoms.length} ${bonds.length} 0 0 0`,'M  V30 BEGIN ATOM'];
    atoms.forEach((a,i)=>{
      if(!/^[A-Z][a-z]?$/.test(a.el)) throw Error('Invalid atom element.');
      const x=Number.isFinite(a.x)?a.x/50:i*1.5,y=Number.isFinite(a.y)?-a.y/50:0;
      lines.push(`M  V30 ${i+1} ${a.el} ${x.toFixed(5)} ${y.toFixed(5)} 0 0${a.charge?' CHG='+a.charge:''}${a.isotope?' MASS='+a.isotope:''}`);
    });
    lines.push('M  V30 END ATOM','M  V30 BEGIN BOND');
    bonds.forEach((b,i)=>{
      if(!atoms[b.a]||!atoms[b.b]||b.a===b.b||![1,2,3].includes(b.order)) throw Error('Invalid stereo graph.');
      if(b.stereo&&!['wedge','hash','E','Z'].includes(b.stereo)) throw Error('Unknown stereochemical bond annotation.');
      if(b.stereo&&((['wedge','hash'].includes(b.stereo)&&b.order!==1)||(['E','Z'].includes(b.stereo)&&b.order!==2))) throw Error('Stereo annotation does not match the bond order.');
      const cfg=b.stereo==='wedge'?1:b.stereo==='hash'?3:b.order===2?2:0;
      lines.push(`M  V30 ${i+1} ${b.order} ${b.a+1} ${b.b+1}${cfg?' CFG='+cfg:''}`);
    });
    lines.push('M  V30 END BOND','M  V30 END CTAB','M  END','');
    return lines.join('\n');
  }
  function materialize(R,data){
    // JSON preserves ordered tetrahedral parity. Generate fresh coordinates so
    // explicit alkene assignments survive RDKit's stereo-cleaning pass without
    // borrowing any stereochemistry from the user's unmarked 2D coordinates.
    const copy=JSON.parse(JSON.stringify(data));
    const molecule=copy.molecules[0]; delete molecule.conformers; molecule.extensions=[];
    const jsonMol=R.get_mol(JSON.stringify(copy),JSON.stringify({assignStereo:false,removeHs:false}));
    if(!jsonMol) throw Error('The stereochemical graph could not be interpreted.');
    try {
      const m=R.get_mol(jsonMol.get_v3Kmolblock(),JSON.stringify({removeHs:false}));
      if(!m) throw Error('The stereochemical graph failed validation.');
      return m;
    } finally {jsonMol.delete();}
  }
  function extract(m){
    const t=JSON.parse(m.get_stereo_tags()), json=JSON.parse(m.get_json());
    return {atoms:t.CIP_atoms.filter(x=>/^[RSrs]$/.test(strip(x[1]))).map(x=>({atom:x[0],label:strip(x[1])})),
      bonds:t.CIP_bonds.filter(x=>/^[EZ]$/.test(strip(x[2]))).map(x=>({a:x[0],b:x[1],label:strip(x[2])})),
      unspecified:t.CIP_atoms.filter(x=>strip(x[1])==='?').map(x=>x[0]),json};
  }
  async function assign(input){
    let m;
    try {
      if(input.atoms.length>1024) throw Error('Stereochemistry is limited to 1024 atoms per request.');
      const R=await init();
      const initial=molblock(input);
      m=R.get_mol(initial,JSON.stringify({removeHs:false}));
      if(!m) throw Error('Invalid valence or stereochemical drawing.');
      let data=JSON.parse(m.get_json()); m.delete();m=null;
      const mol=data.molecules[0];
      if(mol.atoms.length!==input.atoms.length) throw Error('Atom mapping changed during stereo parsing.');
      // Hydrogen and charge agreement avoids silently reinterpreting the graph.
      input.atoms.forEach((a,i)=>{
        const rd=mol.atoms[i];
        if(a.h!==undefined && a.el!=='H' && (rd.impHs??data.defaults.atom.impHs)!==a.h)
          throw Error('Hydrogen count does not match the structure at atom '+(i+1)+'.');
      });
      mol.bonds.forEach(b=>{delete b.stereo;delete b.stereoAtoms;});
      const wedges=input.bonds.filter(b=>b.stereo==='wedge'||b.stereo==='hash');
      const centers=new Set(wedges.map(b=>b.a));
      for(const center of centers){
        const marks=wedges.filter(b=>b.a===center);
        if(marks.length>1) throw Error('Use one wedge or dashed bond per stereocenter to avoid conflicting specifications.');
        if(!['cw','ccw'].includes(mol.atoms[center].stereo)) throw Error('The wedge at atom '+(center+1)+' is ambiguous or not stereogenic.');
      }
      const specifiedAtoms=input.atoms.map((a,i)=>({atom:i,label:a.cip})).filter(a=>a.label);
      for(const a of specifiedAtoms){
        if(!['R','S','r','s'].includes(a.label)) throw Error('Invalid atom configuration.');
        if(!centers.has(a.atom)) mol.atoms[a.atom].stereo='cw';
      }
      const specifiedBonds=input.bonds.filter(b=>b.stereo==='E'||b.stereo==='Z');
      for(const b of specifiedBonds){
        const refA=input.bonds.find(e=>(e.a===b.a||e.b===b.a)&&e!==b);
        const refB=input.bonds.find(e=>(e.a===b.b||e.b===b.b)&&e!==b);
        if(!refA||!refB) throw Error('This double bond has no E/Z configuration.');
        const target=mol.bonds.find(e=>key(...e.atoms)===key(b.a,b.b));
        target.stereo='trans';target.stereoAtoms=[refA.a===b.a?refA.b:refA.a,refB.a===b.b?refB.b:refB.a];
      }
      m=materialize(R,data);let result=extract(m);m.delete();m=null;
      for(const a of specifiedAtoms){
        const found=result.atoms.find(x=>x.atom===a.atom);
        if(!found) throw Error('Atom '+(a.atom+1)+' is not an assignable tetrahedral stereocenter.');
        if(found.label!==a.label){
          if(centers.has(a.atom)) throw Error('The wedge and requested configuration conflict at atom '+(a.atom+1)+'.');
          mol.atoms[a.atom].stereo=mol.atoms[a.atom].stereo==='cw'?'ccw':'cw';
        }
      }
      for(const b of specifiedBonds){
        const found=result.bonds.find(x=>key(x.a,x.b)===key(b.a,b.b));
        if(!found) throw Error('This double bond is not an assignable E/Z stereogenic unit.');
        if(found.label!==b.stereo)mol.bonds.find(e=>key(...e.atoms)===key(b.a,b.b)).stereo='cis';
      }
      m=materialize(R,data);result=extract(m);
      for(const a of specifiedAtoms) if(result.atoms.find(x=>x.atom===a.atom)?.label!==a.label)
        throw Error('Interdependent stereocenters could not be resolved consistently.');
      for(const b of specifiedBonds) if(result.bonds.find(x=>key(x.a,x.b)===key(b.a,b.b))?.label!==b.stereo)
        throw Error('The requested double-bond configurations could not be resolved consistently.');
      for(const center of centers) if(!result.atoms.some(x=>x.atom===center)) throw Error('The marked atom is not a stereogenic center.');
      // Check possible but unassigned alkenes by trying a temporary stereo tag.
      const unspecifiedBonds=[];
      for(const b of input.bonds.filter(b=>b.order===2&&!b.stereo)){
        const neighbors=i=>input.bonds.filter(e=>(e.a===i||e.b===i)&&e!==b).map(e=>e.a===i?e.b:e.a);
        const left=neighbors(b.a),right=neighbors(b.b);
        if(!left.length||!right.length||left.length>2||right.length>2) continue;
        const probe=JSON.parse(JSON.stringify(data));
        const bond=probe.molecules[0].bonds.find(e=>key(...e.atoms)===key(b.a,b.b));
        bond.stereo='trans';bond.stereoAtoms=[left[0],right[0]];
        let candidate;
        try{candidate=materialize(R,probe);if(extract(candidate).bonds.some(e=>key(e.a,e.b)===key(b.a,b.b))) unspecifiedBonds.push({a:b.a,b:b.b});}
        finally{candidate?.delete();}
      }
      const warnings=[];
      if(result.unspecified.length||unspecifiedBonds.length) warnings.push('Some stereogenic atoms or double bonds are unspecified; no configuration was guessed.');
      return {ok:true,...result,unspecifiedBonds,warnings,smiles:m.get_smiles(),molblock:m.get_v3Kmolblock(),version:R.version()};
    } catch(error){return {ok:false,reason:error.message||'Stereochemical assignment failed.'};}
    finally{m?.delete();}
  }
  // Orientation against a semantic ligand order, independent of CIP priority.
  // +1 / -1 invert with any exchange; an implicit H is represented by -1.
  function orderedParity(result,center,ligands){
    const mol=result.json?.molecules?.[0], a=mol?.atoms?.[center];
    if(!a||!['cw','ccw'].includes(a.stereo)) return null;
    const neighbors=mol.bonds.filter(b=>b.atoms.includes(center)).map(b=>b.atoms[0]===center?b.atoms[1]:b.atoms[0]);
    if(neighbors.length===3&&(a.impHs??result.json.defaults.atom.impHs)===1) neighbors.push(-1);
    if(neighbors.length!==4||ligands.length!==4||new Set(ligands).size!==4||ligands.some(x=>!neighbors.includes(x))) return null;
    const perm=ligands.map(x=>neighbors.indexOf(x));let sign=a.stereo==='cw'?1:-1;
    for(let i=0;i<4;i++)for(let j=i+1;j<4;j++)if(perm[i]>perm[j])sign=-sign;
    return sign;
  }
  return {init,assign,orderedParity,molblock};
});
