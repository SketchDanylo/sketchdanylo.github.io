(function(root,factory){
  if(typeof module==='object'&&module.exports) module.exports=factory(require('./native-namer.js'),require('./stereochemistry.js'),require('./biochemical-stereo.js'));
  else root.NamingPipeline=factory(root.NativeNamer,root.Stereochemistry,root.BiochemicalStereo);
})(typeof globalThis!=='undefined'?globalThis:this,function(N,S,B){
  'use strict';
  async function name(graph,model){
    const stereo=await S.assign(model);
    if(!stereo.ok) return {status:'unsupported',reason:stereo.reason,stereoError:true};
    const indices=new Map(graph.map((a,i)=>[a.orig??i,i]));
    const mapped={atoms:stereo.atoms.map(a=>({...a,atom:indices.get(a.atom)})),
      bonds:stereo.bonds.map(b=>({...b,a:indices.get(b.a),b:indices.get(b.b)}))};
    const result=N.name(graph,{stereo:mapped});
    return {...result,isomericSMILES:stereo.smiles,stereo:{atoms:stereo.atoms,bonds:stereo.bonds,
      unspecified:stereo.unspecified,unspecifiedBonds:stereo.unspecifiedBonds,warnings:stereo.warnings},
      biochemical:B.describe(model,stereo)};
  }
  return {name};
});
