'use strict';
importScripts('../vendor/rdkit/dist/RDKit_minimal.js','native-namer.js','stereochemistry.js','biochemical-stereo.js','naming-pipeline.js');
self.onmessage = async ({ data }) => {
  try { self.postMessage({ id: data.id, result: await NamingPipeline.name(data.graph,data.model) }); }
  catch (_) { self.postMessage({ id: data.id, result: { status: 'unsupported', reason: 'The local engine could not complete this structure. No partial name was returned.' } }); }
};
