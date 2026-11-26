/* Chem Puzzles offline cache. VERSION is stamped by tools/generate.cjs. */
const VERSION = 'cp-49a74baef0';
const FILES = ['./', './index.html', './puzzles.js?v=49a74baef0', './namer-de.js?v=49a74baef0'];
self.addEventListener('install', e => { e.waitUntil(caches.open(VERSION).then(c => c.addAll(FILES)).then(() => self.skipWaiting())); });
self.addEventListener('activate', e => { e.waitUntil(caches.keys().then(ks => Promise.all(ks.filter(k => k.startsWith('cp-') && k !== VERSION).map(k => caches.delete(k)))).then(() => self.clients.claim())); });
self.addEventListener('fetch', e => {
  const req = e.request; if (req.method !== 'GET' || new URL(req.url).origin !== location.origin) return;
  // pages: network first so updates arrive; data: cache first (the URL carries the version)
  if (req.mode === 'navigate' || req.url.endsWith('.html') || req.url.endsWith('/OCpuzzles/')) {
    e.respondWith(fetch(req).then(r => { const copy = r.clone(); caches.open(VERSION).then(c => c.put(req, copy)); return r; }).catch(() => caches.match(req, { ignoreSearch: true }).then(r => r || caches.match('./index.html'))));
  } else {
    e.respondWith(caches.match(req).then(r => r || fetch(req).then(res => { const copy = res.clone(); caches.open(VERSION).then(c => c.put(req, copy)); return res; })));
  }
});
