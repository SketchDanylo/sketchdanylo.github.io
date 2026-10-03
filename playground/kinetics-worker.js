importScripts('engine.js?v=20261003-build52', 'kinetics.js?v=20261003-build52');
onmessage = ev => {
  const { id, scene, i, j } = ev.data;
  try { postMessage({ id, result: self.ChemKinetics.study(scene, i, j) }); }
  catch (err) { postMessage({ id, result: { ok: false, reason: 'the calculation failed: ' + err.message } }); }
};
