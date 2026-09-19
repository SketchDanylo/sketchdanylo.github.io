/* IUPAC/IUBMB 3AA-3, 2-Carb-4 and 2-Carb-6. Relative configurations
 * are measured against semantic ligands, never inferred from R/S alone.
 */
(function(root,factory){
  if(typeof module==='object'&&module.exports) module.exports=factory(require('./stereochemistry.js'));
  else root.BiochemicalStereo=factory(root.Stereochemistry);
})(typeof globalThis!=='undefined'?globalThis:this,function(S){
  'use strict';
  const ALDO={3:{'+':'glyceraldehyde'},4:{'++':'erythrose','-+':'threose'},
    5:{'+++':'ribose','-++':'arabinose','+-+':'xylose','--+':'lyxose'},
    6:{'++++':'allose','-+++':'altrose','+-++':'glucose','--++':'mannose','++-+':'gulose','-+-+':'idose','+--+':'galactose','---+':'talose'}};
  const KETO={4:{'+':'erythrulose'},5:{'++':'ribulose','-+':'xylulose'},
    6:{'+++':'allulose','-++':'fructose','+-+':'sorbose','--+':'tagatose'}};
  function describe(input,stereo){
    if(!stereo.ok) return {labels:[],notes:[]};
    const atoms=input.atoms, nb=atoms.map(()=>[]),labels=[],notes=[];
    input.bonds.forEach(b=>{nb[b.a].push({n:b.b,o:b.order});nb[b.b].push({n:b.a,o:b.order});});
    const el=i=>atoms[i]?.el, heavy=i=>nb[i].filter(e=>el(e.n)!=='H');
    const h=i=>nb[i].find(e=>el(e.n)==='H')?.n??((atoms[i].h===1)?-1:null);
    const parity=(i,ligands)=>S.orderedParity(stereo,i,ligands);
    const carbox=i=>el(i)==='C'&&nb[i].some(e=>el(e.n)==='O'&&e.o===2)
      &&nb[i].some(e=>['O','N'].includes(el(e.n))&&e.o===1);
    for(let i=0;i<atoms.length;i++){
      if(el(i)!=='C'||h(i)===null) continue;
      const ns=heavy(i).filter(e=>el(e.n)==='N'&&e.o===1);
      const cs=heavy(i).filter(e=>el(e.n)==='C'&&e.o===1);
      if(ns.length!==1||cs.length!==2||cs.filter(e=>carbox(e.n)).length!==1) continue;
      const co=cs.find(e=>carbox(e.n)).n,side=cs.find(e=>!carbox(e.n)).n;
      const sign=parity(i,[ns[0].n,co,side,h(i)]);
      if(sign) labels.push({kind:'amino-acid',atom:i,series:sign===-1?'L':'D',text:(sign===-1?'L':'D')+' amino-acid configuration at atom '+(i+1)+' (α carbon)'});
    }
    // Carbon components of ordinary unbranched aldoses/ketoses. Ring oxygens,
    // including the anomeric bridge, do not break the underlying carbon chain.
    const seen=new Set();
    for(let start=0;start<atoms.length;start++){
      if(el(start)!=='C'||seen.has(start))continue;
      const component=[],stack=[start];
      while(stack.length){const i=stack.pop();if(seen.has(i))continue;seen.add(i);component.push(i);nb[i].filter(e=>el(e.n)==='C').forEach(e=>stack.push(e.n));}
      if(component.length<3||component.length>7)continue;
      const carbonNB=i=>nb[i].filter(e=>el(e.n)==='C');
      if(component.some(i=>carbonNB(i).length>2||carbonNB(i).some(e=>e.o!==1)))continue;
      const ends=component.filter(i=>carbonNB(i).length===1);
      if(ends.length!==2)continue;
      for(const first of ends){
        const chain=[];let prev=-1,current=first;
        while(current!==undefined){chain.push(current);const next=carbonNB(current).find(e=>e.n!==prev)?.n;prev=current;current=next;}
        const oxy=i=>heavy(i).filter(e=>el(e.n)==='O');
        const carbonylIndex=chain.findIndex(i=>oxy(i).some(e=>e.o===2));
        const acetalIndex=chain.findIndex(i=>oxy(i).filter(e=>e.o===1).length===2);
        const cyclic=carbonylIndex<0&&acetalIndex>=0;
        const ai=cyclic?acetalIndex:carbonylIndex;
        if(![0,1].includes(ai)||(ai===0&&component.length>6))continue;
        const anomer=chain[ai], anomerO=oxy(anomer);
        if((ai===0&&h(anomer)===null)||(ai===1&&(atoms[anomer].h||nb[anomer].some(e=>el(e.n)==='H'))))continue;
        let ringO,exoO,closure;
        if(cyclic){
          ringO=anomerO.find(e=>heavy(e.n).length===2&&heavy(e.n).some(x=>chain.includes(x.n)&&x.n!==anomer))?.n;
          if(ringO===undefined)continue;
          closure=heavy(ringO).find(e=>e.n!==anomer).n;
          const ringSize=chain.indexOf(closure)-ai+2;
          if(![5,6].includes(ringSize))continue;
          exoO=anomerO.find(e=>e.n!==ringO)?.n;
        }
        // Restrict to oxygenated monosaccharide skeletons; deoxy/amino/branched
        // and multiply bridged sugars need additional nomenclature rules.
        if(chain.some((i,k)=>{
          if(heavy(i).some(e=>!['C','O'].includes(el(e.n))))return true;
          if(k===ai)return false;
          const os=oxy(i);return os.length!==1||os[0].o!==1
            ||(k>0&&k<chain.length-1&&h(i)===null);
        }))continue;
        const positions=chain.slice(ai+1,-1);
        if(!positions.length)continue;
        const signs=positions.map(i=>{const k=chain.indexOf(i);return parity(i,[oxy(i)[0].n,chain[k-1],chain[k+1],h(i)]);});
        const ref=positions.at(-1),refSign=signs.at(-1);
        if(!refSign){notes.push('Sugar configuration at atom '+(ref+1)+' is unspecified; D/L and its anomeric relationship were not guessed.');break;}
        const series=refSign===1?'D':'L';
        let anomeric=null;
        if(cyclic){
          const sign=parity(anomer,[exoO,ringO,chain[ai+1],ai===0?h(anomer):chain[0]]);
          if(sign)anomeric=sign===refSign?'α':'β';
          else notes.push('Anomeric atom '+(anomer+1)+' is unspecified; no α/β assignment was made.');
        }
        const pattern=signs.map(x=>x===null?'?':x*refSign>0?'+':'-').join('');
        const family=(ai===0?ALDO:KETO)[chain.length]?.[pattern];
        const ringSize=cyclic?chain.indexOf(closure)-ai+2:0;
        const form=cyclic?(ringSize===5?'furanose':'pyranose'):'';
        let sugar=family?(cyclic?family.replace(/se$/,'')+form:family):(ai===0?'aldose':'ketose')+(cyclic?' '+form:'');
        const substituted=chain.some(i=>oxy(i).some(e=>e.n!==ringO&&heavy(e.n).length>1));
        const text=(anomeric?anomeric+'-':'')+series+'-'+sugar+(substituted?' configuration (substituted sugar unit)':family?'':' series')+' · reference atom '+(ref+1);
        labels.push({kind:'sugar',series,anomeric,family:family||null,reference:ref,anomer:cyclic?anomer:null,atoms:chain,text});break;
      }
    }
    return {labels,notes};
  }
  return {describe};
});
