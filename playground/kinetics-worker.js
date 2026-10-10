importScripts('engine.js?v=20261009-build122', 'kinetics.js?v=20261009-build122');
let last = null;
onmessage = ev => {
  const m = ev.data, K = self.ChemKinetics;
  try {
    if (m.type === 'verify') postMessage({ id: m.id, result: K.verifyEvent(m.scene, m.event, m.products, m.T, m.seeds || 3, 400, m.first || 1) });
    else if (m.type === 'scan') {
      if (!last || last.key !== m.key) last = { key: m.key, e: K.fromScene(m.scene), cache: new Map(), store: new Map() };
      postMessage({ id: m.id, result: K.scanCandidate(last.e, m.c, m.T, m.quick, last.cache, last.store) });
    } else postMessage({ id: m.id, result: K.study(m.scene, m.i, m.j) });
  } catch (err) { postMessage({ id: m.id, result: { ok: false, reason: 'the calculation failed: ' + err.message } }); }
};
