/* Atom inspector. A frozen molecular geometry is calculated in a cancellable worker. */
(function(root){
'use strict';
const ISOTOPE={H:1,He:4,Li:7,Be:9,B:11,C:12,N:14,O:16,F:19,Ne:20,Na:23,Mg:24,Al:27,Si:28,P:31,S:32,Cl:35,Ar:40,K:39,Ca:40,Fe:56,Cu:63,Zn:64,Br:79,Kr:84,I:127,Xe:132};
const ATOMIC_MULT=[0,2,1,2,1,2,3,4,3,2,1,2,1,2,3,4,3,2,1];
class AtomInspector {
  constructor(actions){
    this.actions=actions;this.worker=null;this.input=null;this.i=-1;
    this.panel=document.createElement('aside');this.panel.className='atom-inspector';this.panel.inert=true;
    this.panel.setAttribute('role','dialog');this.panel.setAttribute('aria-label','Atom inspector');
    this.panel.innerHTML=`<div class="inspector-head"><div><p class="eyebrow" id="inspectEyebrow">UNDER THE SURFACE</p><h2 id="inspectName"></h2></div>
      <div class="inspector-acts">
        <button class="iconbtn danger" id="inspectRemove" aria-label="Remove this atom" title="Remove this atom"><svg viewBox="0 0 16 16"><path d="M3.2 4.6h9.6M6.4 4.6V3.3h3.2v1.3M4.5 4.6l.6 8.1h5.8l.6-8.1M6.6 7v3.4M9.4 7v3.4"/></svg></button>
        <button class="iconbtn" id="inspectClose" aria-label="Close atom inspector" title="Close"><svg viewBox="0 0 16 16"><path d="M4.6 4.6l6.8 6.8M11.4 4.6l-6.8 6.8"/></svg></button>
      </div></div>
      <div class="atom-identity"><span id="inspectSymbol"></span><div id="inspectIdentity"></div><canvas id="nucleusCanvas" width="96" height="96" aria-label="Schematic reference nucleus, enlarged"></canvas></div><div class="nucleus-legend" id="nucleusLegend"></div>
      <p class="nucleus-note" id="nucleusNote"></p><dl class="atom-properties" id="atomProperties"></dl>
      <section class="quantum-section"><div class="quantum-heading"><h3>Electron cloud</h3><span id="quantumMethod">Molecular calculation</span></div>
      <div class="seg-mini cloud-seg" id="cloudSeg"><button data-v="atom" aria-pressed="true">This atom</button><button data-v="total" aria-pressed="false">Whole fragment</button></div>
      <div class="density-view"><canvas id="electronCanvas" width="320" height="280" aria-label="Calculated electron density in a plane through the selected atom"></canvas><div id="quantumStatus" role="status"></div><span class="density-plane" id="densityPlane">xy slice through this nucleus</span></div>
      <div class="quantum-controls"><label for="quantumSpin">Spin multiplicity</label><select id="quantumSpin"><option value="1">1 · singlet</option><option value="2">2 · doublet</option><option value="3">3 · triplet</option><option value="4">4 · quartet</option><option value="5">5 · quintet</option><option value="6">6 · sextet</option><option value="7">7 · septet</option></select><button id="quantumRecalc" class="ghost">Recalculate ↻</button></div>
      <div class="quantum-results" id="quantumResults"></div><details class="quantum-notes"><summary>What is calculated</summary><p>Self-consistent Hartree–Fock orbitals in a minimal STO-3G basis. The color shows electron density in electrons/Å³ on a logarithmic scale. It is a spatial slice, not electron paths. Divide density by the electron count for the probability density of a randomly chosen electron.</p><p><b>This atom</b> is the atom's own share of that cloud, by Hirshfeld's stockholder rule: every point is split between the atoms in proportion to what free atoms would put there. The molecular density it divides already carries every neighbour, bonded or not, so the shape you see is this atom as its surroundings and the present thermal geometry have made it.</p><p>The calculation uses this entire isolated fragment at the inspected geometry and selected spin. It omits surrounding molecules and electron correlation. A converged solution need not be the electronic ground state. These orbitals do not supply the MD forces.</p><p>Available for H–Ar, up to 30 atoms and 48 basis functions. <a href="vendor/qchem/LICENSE" target="_blank" rel="noopener">Numerical core: W1neSkin · MIT ↗</a> · <a href="https://www.basissetexchange.org/" target="_blank" rel="noopener">Basis Set Exchange ↗</a></p></details></section>
      <div class="inspector-footer"><span id="inspectTime"></span></div>`;
    document.body.append(this.panel);
    this.el=id=>this.panel.querySelector('#'+id);
    this.el('inspectClose').onclick=()=>this.close();
    this.el('inspectRemove').onclick=()=>{const id=this.atomId;this.close();actions.remove(id);};
    this.el('quantumSpin').onchange=()=>this.calculate();
    this.cloudMode='atom';
    this.el('cloudSeg').querySelectorAll('button').forEach(b=>b.onclick=()=>{
      this.cloudMode=b.dataset.v;
      this.el('cloudSeg').querySelectorAll('button').forEach(x=>x.setAttribute('aria-pressed',x.dataset.v===this.cloudMode));
      this.el('densityPlane').textContent=this.cloudMode==='atom'?"this atom's share":'whole fragment';
      if(this.slices) this.drawDensity(this.slices[this.cloudMode]);
    });
    this.el('quantumRecalc').onclick=()=>this.calculate();
    this.panel.addEventListener('keydown',e=>{if(e.key==='Escape'){e.stopPropagation();this.close();}});
  }
  /* The card belongs to the atom it describes: it is placed beside it on screen and kept there
     while the camera moves, and the field draws the leader line that ties the two together. */
  place(sx,sy){
    const p=this.panel, w=p.offsetWidth||300, h=p.offsetHeight||360, gap=26;
    let x=sx+gap, y=sy-h/2;
    if(x+w>innerWidth-12) x=sx-gap-w;
    x=Math.max(12,Math.min(innerWidth-w-12,x));
    y=Math.max(12,Math.min(innerHeight-h-12,y));
    p.style.left=x+'px'; p.style.top=y+'px';
    this.anchor=[x<sx?x+w:x, Math.max(y+18,Math.min(y+h-18,sy))];
  }
  close(){this.worker?.terminate();this.worker=null;cancelAnimationFrame(this.nucleusRaf);this.panel.classList.remove('open');this.panel.inert=true;this.i=-1;this.slices=null;this.anchor=null;this.actions.focus?.();}
  open(i){
    const e=this.actions.engine,E=root.ChemEngine.ELEMENTS,el=E[e.type[i]];
    this.actions.pause();this.worker?.terminate();this.atomId=e.ids[i];this.i=i;
    const fragments=e.fragments(),list=fragments.list[fragments.comp[i]],a3=3*i;
    this.selected=list.indexOf(i);
    this.input={atoms:list.map(j=>({Z:E[e.type[j]].Z,x:e.pos[3*j]-e.pos[a3],y:e.pos[3*j+1]-e.pos[a3+1],z:e.pos[3*j+2]-e.pos[a3+2]})),charge:list.reduce((s,j)=>s+e.formal[j],0)};
    const electrons=this.input.atoms.reduce((s,a)=>s+a.Z,0)-this.input.charge;
    let mult=electrons%2?2:1;
    if(list.length===1&&this.input.charge===0) mult=ATOMIC_MULT[el.Z]||mult;
    if(list.length===2&&this.input.charge===0&&this.input.atoms.every(a=>a.Z===8)) mult=3;
    this.el('quantumSpin').value=String(mult);
    this.el('inspectName').textContent=el.name;
    this.el('inspectSymbol').textContent=el.sym;this.el('inspectSymbol').style.color=el.color;
    this.el('inspectIdentity').textContent=`Z ${el.Z} · ${e.formulaOf(list)}\nAtom ${e.ids[i]} · ${list.length} atom${list.length===1?'':'s'} in fragment`;
    const A=ISOTOPE[el.sym];
    this.el('nucleusNote').textContent=`${el.sym}-${A} reference nucleus · enlarged schematic. The MD uses average atomic masses; its isotopes are unspecified.`;
    const v=Math.hypot(...e.vel.subarray(a3,a3+3));
    const rows=[['Protons',el.Z],['Neutrons in reference isotope',A-el.Z],['Molecular electrons',electrons],['MD partial charge',`${e.q[i]>=0?'+':''}${e.q[i].toFixed(3)} e`],['Formal charge',`${e.formal[i]} e`],['Atomic mass',`${e.mass[i].toFixed(3)} u`],['Speed',`${(v*100).toFixed(3)} km/s`],['Kinetic energy',`${(.5*e.mass[i]*v*v*1e4).toFixed(3)} kJ/mol`]];
    this.el('atomProperties').replaceChildren(...rows.map(([key,value])=>{const d=document.createElement('div'),dt=document.createElement('dt'),dd=document.createElement('dd');dt.textContent=key;dd.textContent=value;d.append(dt,dd);return d;}));
    this.electronCount=electrons;this.startNucleus(el.Z,A);
    this.el('inspectTime').textContent=`Held at ${(e.time/1000).toFixed(3)} ps · sample ${e.temperature().toFixed(0)} K`;
    this.el('inspectEyebrow').textContent=`ATOM ${e.ids[i]} · THIS INSTANT'S GEOMETRY`;
    this.cloudMode='atom';
    this.el('cloudSeg').querySelectorAll('button').forEach(x=>x.setAttribute('aria-pressed',x.dataset.v==='atom'));
    this.el('densityPlane').textContent="this atom's share";
    this.panel.inert=false;this.panel.classList.add('open');
    this.place(...(this.actions.screenOf?.(i)||[innerWidth/2,innerHeight/2]));
    this.el('inspectClose').focus();this.calculate();
  }
  /* Nucleus artwork: A nucleons packed in a ball, lit and depth-sorted, turning slowly.
     Placement is schematic — a nucleus has no fixed particle arrangement. */
  startNucleus(Z,A){
    cancelAnimationFrame(this.nucleusRaf);
    const canvas=this.el('nucleusCanvas'),dpr=Math.min(2,devicePixelRatio||1),W=canvas.clientWidth||96,H=canvas.clientHeight||96;
    canvas.width=W*dpr;canvas.height=H*dpr;
    const ctx=canvas.getContext('2d');ctx.setTransform(dpr,0,0,dpr,0,0);
    // deterministic ball packing: Fibonacci shells so nucleons never overlap oddly
    const nucleons=[];
    for(let k=0;k<A;k++){
      const t=(k+.5)/A,z=1-2*t,rho=Math.sqrt(Math.max(0,1-z*z)),theta=k*2.399963229;
      const shell=Math.cbrt((k+.5)/A);
      nucleons.push({x:rho*Math.cos(theta)*shell,y:rho*Math.sin(theta)*shell,z:z*shell,proton:false});
    }
    // assign protons evenly through the packing so the render reads as a mixed nucleus
    nucleons.map((n,i)=>({n,i})).sort((a,b)=>((a.i*2654435761)%A)-((b.i*2654435761)%A)).forEach((e,rank)=>e.n.proton=rank<Z);
    const R=Math.min(W,H)*.5-6,scale=R/Math.max(1,Math.cbrt(A)*.62),rad=Math.max(2.4,scale*.58);
    const reduce=matchMedia('(prefers-reduced-motion: reduce)').matches;
    const frame=now=>{
      const a=reduce?.6:now/2600,ca=Math.cos(a),sa=Math.sin(a),tilt=.42,ct=Math.cos(tilt),st=Math.sin(tilt);
      ctx.clearRect(0,0,W,H);
      const cx=W/2,cy=H/2;
      const halo=ctx.createRadialGradient(cx,cy,R*.2,cx,cy,R*1.5);
      halo.addColorStop(0,'rgba(255,190,130,.14)');halo.addColorStop(.55,'rgba(120,170,220,.07)');halo.addColorStop(1,'rgba(0,0,0,0)');
      ctx.fillStyle=halo;ctx.beginPath();ctx.arc(cx,cy,R*1.5,0,7);ctx.fill();
      const view=nucleons.map(n=>{
        const x1=n.x*ca-n.z*sa,z1=n.x*sa+n.z*ca,y1=n.y*ct-z1*st,z2=n.y*st+z1*ct;
        return{x:cx+x1*scale,y:cy+y1*scale,z:z2,proton:n.proton};
      }).sort((p,q)=>p.z-q.z);
      for(const p of view){
        const depth=(p.z+1)/2,r=rad*(.82+.22*depth),shade=.55+.45*depth;
        const g=ctx.createRadialGradient(p.x-r*.4,p.y-r*.45,r*.08,p.x,p.y,r);
        if(p.proton){g.addColorStop(0,`rgba(255,${Math.round(214*shade)},${Math.round(170*shade)},1)`);g.addColorStop(.55,`rgba(${Math.round(235*shade)},${Math.round(124*shade)},${Math.round(64*shade)},1)`);g.addColorStop(1,`rgba(${Math.round(120*shade)},${Math.round(48*shade)},${Math.round(22*shade)},1)`);}
        else{g.addColorStop(0,`rgba(${Math.round(226*shade)},${Math.round(240*shade)},255,1)`);g.addColorStop(.55,`rgba(${Math.round(140*shade)},${Math.round(176*shade)},${Math.round(206*shade)},1)`);g.addColorStop(1,`rgba(${Math.round(56*shade)},${Math.round(80*shade)},${Math.round(104*shade)},1)`);}
        ctx.beginPath();ctx.arc(p.x,p.y,r,0,7);ctx.fillStyle=g;ctx.fill();
        ctx.beginPath();ctx.arc(p.x-r*.34,p.y-r*.38,r*.26,0,7);ctx.fillStyle=`rgba(255,255,255,${.28*depth+.08})`;ctx.fill();
      }
      if(!reduce) this.nucleusRaf=requestAnimationFrame(frame);
    };
    this.nucleusRaf=requestAnimationFrame(frame);
    this.el('nucleusLegend').innerHTML=`<span class="pt p"><i></i>${Z} proton${Z===1?'':'s'}</span><span class="pt n"><i></i>${A-Z} neutron${A-Z===1?'':'s'}</span><span class="pt e"><i></i>${this.electronCount} electron${this.electronCount===1?'':'s'}</span>`;
  }
  calculate(){
    if(!this.input)return;this.worker?.terminate();
    const status=this.el('quantumStatus');status.hidden=false;status.className='calculating';status.textContent='Preparing calculation';
    this.el('quantumResults').textContent='';this.el('quantumMethod').textContent='HF / STO-3G';this.slices=null;
    this.el('electronCanvas').getContext('2d').clearRect(0,0,320,280);
    const worker=this.worker=new Worker('quantum-worker.js');
    worker.onmessage=({data})=>{
      if(this.worker!==worker)return;
      if(data.status){status.textContent=data.status;return;}
      status.className='';
      if(data.error){status.textContent=data.error;worker.terminate();this.worker=null;return;}
      status.hidden=true;
      this.slices={atom:data.result.atomSlice,total:data.result.totalSlice};
      this.drawDensity(this.slices[this.cloudMode]);
      const {scf,method,electronTrace}=data.result;
      this.el('quantumMethod').textContent=method;
      const q=scf.mulliken[this.selected];
      this.el('quantumResults').textContent=`${scf.nelec} electrons · ${scf.iterations} iterations · E = ${scf.E.toFixed(6)} Eh\nMulliken charge here: ${q>=0?'+':''}${q.toFixed(3)} e${scf.uhf?' · ⟨S²⟩ '+scf.S2.toFixed(3):''}\nElectron-count check: ${electronTrace.toFixed(6)}`;
      worker.terminate();this.worker=null;
    };
    worker.onerror=event=>{console.error('Quantum worker:',event.message,event.filename,event.lineno);status.className='';status.textContent='Quantum worker failed. Recalculate to retry.';worker.terminate();this.worker=null;};
    worker.postMessage({input:{...this.input,mult:Number(this.el('quantumSpin').value)},selected:this.selected});
  }
  drawDensity(grid){
    if(!grid)return;
    const canvas=this.el('electronCanvas'),ctx=canvas.getContext('2d'),n=grid.resolution;
    ctx.clearRect(0,0,canvas.width,canvas.height);
    const off=document.createElement('canvas');off.width=off.height=n;const c=off.getContext('2d'),im=c.createImageData(n,n);
    const norm=Math.log1p(grid.max/.002);
    for(let k=0;k<grid.data.length;k++){
      const t=Math.log1p(grid.data[k]/.002)/norm;
      im.data[4*k]=Math.round(8+178*t*t);im.data[4*k+1]=Math.round(19+218*t);im.data[4*k+2]=Math.round(30+205*t);im.data[4*k+3]=255;
    }
    c.putImageData(im,0,0);ctx.imageSmoothingEnabled=true;ctx.drawImage(off,20,0,280,280);
    // Atomic nuclei are marks laid over the density, not orbiting electrons.
    const factor=140/grid.extent;
    this.input.atoms.forEach((a,i)=>{const x=160+(a.x-grid.origin[0])*factor,y=140+(a.y-grid.origin[1])*factor;ctx.strokeStyle=i===this.selected?'#ffba66':'#cdeaff88';ctx.lineWidth=1;ctx.beginPath();ctx.arc(x,y,i===this.selected?5:2,0,2*Math.PI);ctx.stroke();});
    ctx.fillStyle='#a4c3d5';ctx.font='9px monospace';ctx.textAlign='right';
    ctx.fillText(`${(2*grid.extent).toFixed(2)} Å across · log ρ · peak ${grid.max.toFixed(1)} e/Å³`,300,263);
  }
}
root.AtomInspector=AtomInspector;
})(window);
