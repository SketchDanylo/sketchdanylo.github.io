/* Shared, testable boundaries for Nomenclature imports and playback units. */
(function (root, factory) {
  if (typeof module === 'object' && module.exports) module.exports = factory();
  else root.ChemProtocol = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';
  function stepsPerSecond(speed) {
    if (!Number.isFinite(speed) || speed <= 0) return 0;
    return speed <= 1 ? Math.pow(20000, (speed - 0.1) / 0.9) : 20000 * speed;
  }
  function validMolecule(m, elements) {
    if (!m || m.format !== 'chem-playground/molecule@1' || (m.units != null && m.units !== 'angstrom')) return false;
    if (typeof m.id !== 'string' || !/^[\w-]{1,100}$/.test(m.id) || typeof m.formula !== 'string' || m.formula.length > 500) return false;
    if (m.name != null && (typeof m.name !== 'string' || m.name.length > 2000)) return false;
    if (m.smiles != null && (typeof m.smiles !== 'string' || m.smiles.length > 10000)) return false;
    if (!Array.isArray(m.atoms) || !m.atoms.length || m.atoms.length > 500 || !Array.isArray(m.bonds) || m.bonds.length > 2000) return false;
    if (!m.atoms.every(a => a && Object.hasOwn(elements, a.el) && [a.x, a.y, a.z ?? 0].every(v => Number.isFinite(v) && Math.abs(v) <= 10000) && Number.isInteger(a.charge ?? 0) && Math.abs(a.charge ?? 0) <= 8)) return false;
    const seen = new Set();
    return m.bonds.every(b => {
      if (!b || !Number.isInteger(b.a) || !Number.isInteger(b.b) || b.a < 0 || b.b < 0 || b.a >= m.atoms.length || b.b >= m.atoms.length || b.a === b.b || ![1, 1.5, 2, 3].includes(b.order)) return false;
      const key = Math.min(b.a, b.b) + ':' + Math.max(b.a, b.b);
      if (seen.has(key)) return false;
      seen.add(key); return true;
    });
  }
  return { stepsPerSecond, validMolecule };
});
