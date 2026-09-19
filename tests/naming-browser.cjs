const { chromium } = require('playwright');
const assert = require('node:assert/strict');
const fs = require('node:fs');
const path = require('node:path');
const http = require('node:http');
const { graph } = require('./naming-fixtures.cjs');
const root = path.resolve(__dirname, '..');
const server = http.createServer((req, res) => {
  const file = path.resolve(root, '.' + new URL(req.url, 'http://localhost').pathname);
  if (!file.startsWith(root + path.sep)) { res.writeHead(403).end(); return; }
  fs.readFile(file, (err, data) => {
    if (err) { res.writeHead(404).end(); return; }
    res.setHeader('Content-Type', file.endsWith('.js') ? 'text/javascript' : 'text/html'); res.end(data);
  });
});
(async () => {
  await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
  const origin = 'http://127.0.0.1:' + server.address().port;
  const browser = await chromium.launch({ channel: process.env.BROWSER_CHANNEL || 'msedge', headless: true });
  try {
    const context = await browser.newContext({ viewport: { width: 1440, height: 900 }, permissions: ['clipboard-read', 'clipboard-write'] });
    const page = await context.newPage(), errors = [], requests = [];
    page.on('pageerror', error => errors.push(error.message));
    page.on('request', request => requests.push(request.url()));
    await page.route('https://pubchem.ncbi.nlm.nih.gov/**', route => route.fulfill({ status: 404, body: '{}' }));
    await page.goto(origin + '/nomenclature.html');
    const load = async (smiles, positions) => {
      await page.evaluate(({ g, positions }) => {
        atoms = g.map((a, i) => ({ el: a.el, charge: a.charge, manualCharge: !!a.charge,
          x: positions?.[i]?.[0] ?? 160 + (i % 7) * 75, y: positions?.[i]?.[1] ?? 180 + Math.floor(i / 7) * 85 }));
        bonds = g.flatMap((a, i) => a.nb.filter(e => e.n > i).map(e => ({ a: i, b: e.n, order: e.o })));
        skeletal = true; afterChange(); fitView();
      }, { g: graph(smiles), positions });
    };
    await load('CC(NC1CCCCC1)C(O)C(=O)C', [
      [250, 450], [315, 410], [315, 335], [250, 295], [250, 220], [185, 182],
      [120, 220], [120, 295], [185, 332], [380, 450], [380, 525], [445, 410], [445, 335], [510, 450]
    ]);
    await page.locator('#nameBtn').click();
    await page.getByText('4-(cyclohexylamino)-3-hydroxypentan-2-one', { exact: true }).waitFor();
    assert.match(await page.locator('#resultsBody').innerText(), /14 \/ 14 heavy atoms/);
    assert(requests.some(url => url.includes('naming-worker.js')), 'Worker did not load');
    assert(!requests.some(url => url.includes('pubchem')), 'Offline naming sent a network lookup');
    await page.getByRole('button', { name: 'Copy EN', exact: true }).click();
    assert.equal(await page.evaluate(() => navigator.clipboard.readText()), '4-(cyclohexylamino)-3-hydroxypentan-2-one');
    if (process.env.NAMING_SCREENSHOT) await page.screenshot({ path: process.env.NAMING_SCREENSHOT, fullPage: true });
    await page.getByRole('button', { name: 'Look up on PubChem' }).click();
    await page.getByText(/Database lookup unavailable/).waitFor();
    assert.equal(await page.locator('.the-name').first().innerText(), '4-(cyclohexylamino)-3-hydroxypentan-2-one');
    await page.locator('#drawerClose').click();
    await load('CC(NCOCO)CCC(Cl)C');
    await page.locator('#nameBtn').click();
    await page.getByText('({[(5-chlorohexan-2-yl)amino]methyl}oxy)methanol', { exact: true }).waitFor();
    assert.match(await page.locator('#resultsBody').innerText(), /12 \/ 12 heavy atoms/);
    await page.setViewportSize({ width: 390, height: 844 });
    assert(await page.locator('.the-name').first().evaluate(el => el.scrollWidth <= el.clientWidth + 1), 'Long name overflowed on mobile');
    await page.evaluate(async () => {
      const pending = runNaming(); atoms.push({ el: 'C', x: 0, y: 0, charge: 0 }); afterChange(); await pending;
    });
    assert.match(await page.locator('#resultsBody').innerText(), /Structure changed/);
    await load('C1=CC=C2C=CC=CC2=C1');
    await page.evaluate(() => runNaming());
    assert.match(await page.locator('#resultsBody').innerText(), /Not yet supported offline/);
    assert.equal(await page.locator('.the-name').count(), 0);
    assert.deepEqual(errors, []);
    // The no-worker path also works for browsers restricting workers/file URLs.
    await page.evaluate(() => { window.Worker = class { constructor() { throw new Error('blocked'); } }; });
    await load('CCO'); await page.evaluate(() => runNaming());
    assert.equal(await page.locator('.the-name').first().innerText(), 'ethanol');
    await page.evaluate(async () => { const job=runNaming(); atoms[0].x+=10; await job; });
    assert.equal(await page.locator('.the-name').first().innerText(), 'ethanol');
    await page.evaluate(async () => { atoms.push({el:'H',x:0,y:0,charge:0}); afterChange(); await runNaming(); });
    assert.match(await page.locator('#resultsBody').innerText(), /Unattached hydrogen/);
    assert.equal(await page.locator('.the-name').count(), 0);
    await load('CCO');
    await page.evaluate(async () => {
      atoms.push({el:'H',x:0,y:0,charge:0}); bonds.push({a:2,b:3,order:1}); afterChange(); await runNaming();
    });
    assert.equal(await page.locator('.the-name').first().innerText(), 'ethanol');
    console.log('Browser checks passed: worker, offline naming, both screenshots, copy, failed optional lookup, mobile wrapping, edit cancellation, unsupported structures, worker fallback.');
  } finally { await browser.close(); server.close(); }
})().catch(error => { console.error(error); server.close(); process.exitCode = 1; });
