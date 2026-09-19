const{chromium}=require('playwright'),assert=require('node:assert/strict'),path=require('node:path'),{pathToFileURL}=require('node:url');
const fixtures=require('./stereo-fixtures.json');
(async()=>{
  const browser=await chromium.launch({channel:process.env.BROWSER_CHANNEL||'msedge',headless:true});
  try{
    const context=await browser.newContext(),page=await context.newPage(),errors=[],requests=[];
    page.on('pageerror',e=>errors.push(e.message));page.on('request',r=>requests.push(r.url()));
    await context.route(/^https?:/,route=>route.abort());
    await page.goto(process.env.STANDALONE_URL||pathToFileURL(path.resolve(__dirname,'../nomenclature.html')).href);
    for(const index of [0,2,9,10,19,22,39,41]){
      await page.evaluate(async m=>{atoms=m.atoms;bonds=m.bonds;afterChange();await runNaming();},structuredClone(fixtures[index].model));
      assert(await page.locator('.the-name').count()>0,await page.locator('#resultsBody').innerText());
      if(fixtures[index].series)assert((await page.locator('#resultsBody').innerText()).includes(fixtures[index].series));
    }
    await page.evaluate(()=>{window.Worker=class{constructor(){throw Error('blocked')}};});
    await page.evaluate(async m=>{atoms=m.atoms;bonds=m.bonds;afterChange();await runNaming();},structuredClone(fixtures[19].model));
    assert.match(await page.locator('#resultsBody').innerText(),/α-D-glucopyranose/);
    assert.deepEqual(errors,[]);assert(!requests.some(url=>/^https?:/.test(url)),'Offline file attempted a network request');
    console.log('Standalone HTML passes with all HTTP(S) blocked: file loading, embedded worker/WASM, R/S, E/Z, biochemical configurations and no-worker fallback.');
  }finally{await browser.close();}
})().catch(e=>{console.error(e);process.exitCode=1;});
