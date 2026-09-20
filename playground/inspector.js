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
    this.panel.innerHTML=`<div class="inspector-head"><div><p class="eyebrow">UNDER THE SURFACE</p><h2 id="inspectName"></h2></div><button class="x" id="inspectClose" aria-label="Close atom inspector">×</button></div>
      <div class="atom-identity"><span id="inspectSymbol"></span><div id="inspectIdentity"></div><canvas id="nucleusCanvas" width="100" height="100" aria-label="Schematic reference nucleus, enlarged"></canvas></div>
      <p class="nucleus-note" id="nucleusNote"></p><dl class="atom-properties" id="atomProperties"></dl>
      <section class="quantum-section"><div class="quantum-heading"><h3>Electron density</h3><span id="quantumMethod">Molecular calculation</span></div>
      <div class="density-view"><canvas id="electronCanvas" width="320" height="280" aria-label="Calculated molecular electron density in a plane through the selected atom"></canvas><div id="quantumStatus" role="status"></div><span class="density-plane">xy slice through this nucleus</span></div>
      <div class="quantum-controls"><label for="quantumSpin">Spin multiplicity</label><select id="quantumSpin"><option value="1">1 · singlet</option><option value="2">2 · doublet</option><option value="3">3 · triplet</option><option value="4">4 · quartet</option><option value="5">5 · quintet</option><option value="6">6 · sextet</option><option value="7">7 · septet</option></select><button id="quantumRecalc" class="ghost">Recalculate ↻</button></div>
      <div class="quantum-results" id="quantumResults"></div><details class="quantum-notes"><summary>What is calculated</summary><p>Self-consistent Hartree–Fock orbitals in a minimal STO-3G basis. The color shows electron density in electrons/Å³ on a logarithmic scale. It is a spatial slice, not electron paths. Divide density by the electron count for the probability density of a randomly chosen electron.</p><p>The calculation uses this entire isolated fragment at the inspected geometry and selected spin. It omits surrounding molecules and electron correlation. A converged solution need not be the electronic ground state. These orbitals do not supply the MD forces.</p><p>Available for H–Ar, up to 30 atoms and 48 basis functions. <a href="vendor/qchem/LICENSE" target="_blank" rel="noopener">Numerical core: W1neSkin · MIT ↗</a> · <a href="https://www.basissetexchange.org/" target="_blank" rel="noopener">Basis Set Exchange ↗</a></p></details></section>
      <div class="inspector-footer"><span id="inspectTime"></span><button class="ghost" id="inspectRemove">Remove atom</button></div>`;
    document.body.append(this.panel);
    this.el=id=>this.panel.querySelector('#'+id);
    this.el('inspectClose').onclick=()=>this.close();
    this.el('inspectRemove').onclick=()=>{const id=this.atomId;this.close();actions.remove(id);};
    this.el('quantumSpin').onchange=()=>this.calculate();
    this.el('quantumRecalc').onclick=()=>this.calculate();
    this.panel.addEventListener('keydown',e=>{if(e.key==='Escape'){e.stopPropagation();this.close();}});
  }
  close(){this.worker?.terminate();this.worker=null;this.panel.classList.remove('open');this.panel.inert=true;this.i=-1;this.actions.focus?.();}
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
    this.drawNucleus(el.Z,A);
    this.el('inspectTime').textContent=`Paused snapshot · ${(e.time/1000).toFixed(3)} ps`;
    this.panel.inert=false;this.panel.classList.add('open');this.el('inspectClose').focus();this.calculate();
  }
  drawNucleus(Z,A){
    const ctx=this.el('nucleusCanvas').getContext('2d');ctx.clearRect(0,0,100,100);
    const dots=[];for(let k=0;k<A;k++){const z=1-2*(k+.5)/A,r=Math.sqrt(1-z*z),theta=k*2.399963;dots.push({x:50+24*r*Math.cos(theta),y:50+24*r*Math.sin(theta),z,p:(k*Z)%A<Z});}
    // Exactly Z proton symbols; placement is schematic, not a nuclear wavefunction.
    const order=[...dots].sort((a,b)=>a.x-b.x);order.forEach((d,k)=>d.p=k<Z);
    dots.sort((a,b)=>a.z-b.z);
    for(const d of dots){const r=6.5*(.85+.15*d.z),g=ctx.createRadialGradient(d.x-2,d.y-2,0,d.x,d.y,r);g.addColorStop(0,d.p?'#ffd0a0':'#c8e4fa');g.addColorStop(1,d.p?'#a45635':'#466d8b');ctx.beginPath();ctx.arc(d.x,d.y,r,0,Math.PI*2);ctx.fillStyle=g;ctx.fill();if(d.p&&A<45){ctx.fillStyle='#542a19';ctx.font='8px sans-serif';ctx.textAlign='center';ctx.fillText('+',d.x,d.y+3);}}
  }
  calculate(){
    if(!this.input)return;this.worker?.terminate();
    const status=this.el('quantumStatus');status.hidden=false;status.className='calculating';status.textContent='Preparing calculation';
    this.el('quantumResults').textContent='';this.el('quantumMethod').textContent='HF / STO-3G';
    this.el('electronCanvas').getContext('2d').clearRect(0,0,320,280);
    const worker=this.worker=new Worker('quantum-worker.js');
    worker.onmessage=({data})=>{
      if(this.worker!==worker)return;
      if(data.status){status.textContent=data.status;return;}
      status.className='';
      if(data.error){status.textContent=data.error;worker.terminate();this.worker=null;return;}
      status.hidden=true;this.drawDensity(data.result.slice);
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
    const canvas=this.el('electronCanvas'),ctx=canvas.getContext('2d'),n=grid.resolution;
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
    ctx.fillStyle='#a4c3d5';ctx.font='9px monospace';ctx.textAlign='right';ctx.fillText(`${(2*grid.extent).toFixed(2)} Å across · log ρ`,300,263);
  }
}
root.AtomInspector=AtomInspector;
})(window);
