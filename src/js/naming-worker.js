'use strict';
importScripts('native-namer.js');
self.onmessage = ({ data }) => {
  try { self.postMessage({ id: data.id, result: NativeNamer.name(data.graph) }); }
  catch (_) { self.postMessage({ id: data.id, result: { status: 'unsupported', reason: 'The local engine could not complete this structure. No partial name was returned.' } }); }
};
