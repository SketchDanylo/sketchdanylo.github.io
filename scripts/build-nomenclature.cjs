/* Build the deployable, self-contained HTML from readable development sources. */
const fs=require('node:fs'),path=require('node:path'),assert=require('node:assert/strict');
const root=path.resolve(__dirname,'..');
const read=file=>fs.readFileSync(path.join(root,file),'utf8').replace(/\r\n/g,'\n');
const files=['src/js/native-namer.js','src/vendor/rdkit/dist/RDKit_minimal.js','src/js/stereochemistry.js','src/js/biochemical-stereo.js','src/js/naming-pipeline.js'];
const runtime=files.map(file=>'\n/* '+file+' */\n'+read(file)).join('\n');
const worker=read('src/js/naming-worker.js').replace(/^importScripts\([^\n]*\);\n/m,'');
assert(!/<\/script/i.test(runtime+worker),'Embedded source contains a closing script tag');
const wasm=fs.readFileSync(path.join(root,'src/vendor/rdkit/dist/RDKit_minimal.wasm')).toString('base64');
const license=read('src/vendor/rdkit/LICENSE');
let html=read('src/nomenclature.template.html');
const scripts=files.map(file=>'<script src="'+file+'"></script>').join('\n');
assert(html.includes(scripts),'Runtime script block changed; update the bundler');
html=html.replace(scripts,()=>
  '<!-- Self-contained build. Regenerate with node scripts/build-nomenclature.cjs.\nRDKit 2026.03.6 third-party license:\n'+license+'-->\n'+
  '<script id="naming-wasm" type="application/octet-stream">'+wasm+'</script>\n'+
  '<script id="naming-runtime">'+runtime+'\n</script>\n'+
  '<script id="naming-worker-source" type="text/plain">'+worker+'</script>');
// System fonts keep both canvas and controls usable with no network access.
html=html.replace(/^<link[^\n]*fonts\.googleapis\.com[^\n]*\n/gm,'').replace(/\bInter\b/g,'Arial');
assert(!/<script[^>]+\bsrc=|<link[^>]+\bhref=/i.test(html),'External script or stylesheet remains');
fs.writeFileSync(path.join(root,'nomenclature.html'),html);
console.log('Built nomenclature.html: '+fs.statSync(path.join(root,'nomenclature.html')).size.toLocaleString()+' bytes; JS, WASM, worker and styles embedded.');
