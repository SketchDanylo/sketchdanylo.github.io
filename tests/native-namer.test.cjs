const { test } = require('node:test');
const assert = require('node:assert/strict');
const { name, stem } = require('../src/js/native-namer.js');

const { graph, permute, fixtures } = require('./naming-fixtures.cjs');
for (const [label, smiles, expected] of fixtures) test(label, () => {
  const g = graph(smiles), result = name(g);
  assert.equal(result.status, 'ok', result.reason);
  if (expected) assert.equal(result.name, expected);
  assert.equal(result.coveredAtoms, g.length);
  for (let seed = 1; seed <= 5; seed++) assert.equal(name(permute(g, seed)).name, result.name, 'Atom order changed the name');
});
for (const [label, smiles] of [
  ['unsaturated fused rings', 'C1=CC=C2C=CC=CC2=C1'],
  ['charged amine', 'C[NH3+]'], ['separate fragments', 'C.C'], ['unsupported phosphorus', 'CP(=O)(O)O'],
  ['invalid carbon valence', 'C(C)(C)(C)(C)C'], ['nonalternating six-membered ring', 'C1=C=CC=CC1'],
]) test('honest unsupported: ' + label, () => {
  const result = name(graph(smiles));
  if (label === 'nonalternating six-membered ring') {
    assert.notEqual(result.name, 'benzene');
  } else { assert.equal(result.status, 'unsupported'); assert.equal(result.name, undefined); }
});
test('long chains and multipliers', () => {
  for (const n of [31, 40, 50, 99, 100]) {
    const g = graph('C'.repeat(n));
    assert.equal(name(g).name, stem(n) + 'ane');
    assert.equal(name(permute(g, n)).name, stem(n) + 'ane');
  }
  assert.equal(stem(31), 'hentriacont');
  assert.equal(name(graph('C'.repeat(101))).status, 'unsupported');
});
test('bounded search returns no guessed name', () => {
  const result = name(graph('CC(C)CC(C)CC(C)CC'), { maxWork: 10 });
  assert.equal(result.status, 'unsupported'); assert.equal(result.name, undefined);
});
test('input is immutable', () => {
  const g = graph('CC(NC1CCCCC1)C(O)C(=O)C'), before = JSON.stringify(g);
  name(g); assert.equal(JSON.stringify(g), before);
});

