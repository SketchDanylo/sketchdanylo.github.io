/* Building a molecule by hand is the thing this sandbox is for, and for a long time it did not
   work: the moment a bond formed under the pointer, the step fabricated more than a thousand
   kJ/mol and the molecule that had just formed blew itself apart. Two things were wrong. The
   energy check that would have caught it was switched off whenever the pointer held an atom,
   because the servo does real work and the check could not tell that work from an error. And a
   newly formed molecule in a nearly empty chamber had nothing to hand its bond energy to.
   These lock both down. */
const assert = require('node:assert/strict');
const { Engine } = require('../engine.js');
let count = 0;
function test(name, run) { run(); console.log('PASS', name); count++; }
const total = e => e.Epot + e.kinetic();
const products = e => {
  const c = {};
  for (const f of e.fragments().list) { const s = e.formulaOf(f); c[s] = (c[s] || 0) + 1; }
  return Object.keys(c).sort().map(k => (c[k] > 1 ? c[k] + ' ' : '') + k).join(' + ');
};
function methylene(e, cx, cy) {
  e.addAtom('C', cx, cy, 0, { thermal: false });
  e.addAtom('H', cx - .55, cy - .9, 0, { thermal: false });
  e.addAtom('H', cx + .55, cy - .9, 0, { thermal: false });
  e.touch(); e.minimize(600, .2); e.refresh();
}
/* The pointer closing on a target, exactly as the field does it: the point it aims for is always
   a little ahead of where the held atom is now. `to` is read every step, because someone watching
   the screen follows the molecule they are aiming at rather than a fixed spot on the glass. */
function dragOnto(e, h, to, fs) {
  e.tweezer = { i: h, x: e.pos[3 * h], y: e.pos[3 * h + 1] };
  const n = Math.round(fs / e.dt);
  for (let s = 0; s < n; s++) {
    const [tx, ty] = typeof to === 'function' ? to() : to;
    e.tweezer.x = e.pos[3 * h] + (tx - e.pos[3 * h]) * .02;
    e.tweezer.y = e.pos[3 * h + 1] + (ty - e.pos[3 * h + 1]) * .02;
    e.step();
  }
  e.tweezer = null;
}
const justAbove = (e, c, gap) => () => [e.pos[3 * c], e.pos[3 * c + 1] + gap];
function reactiveDrag(opts = {}) {
  const e = new Engine({ width: 36, height: 20, depth: 12, T: 300, thermostat: false, ...opts });
  methylene(e, 18, 10);
  const h = e.addAtom('H', 18, 14.5, 0, { thermal: false });
  e.touch(); e.refresh();
  return { e, h, E0: total(e) };
}

test('A bond formed under the pointer does not fabricate energy', () => {
  const { e, h, E0 } = reactiveDrag();
  dragOnto(e, h, [18, 11.1], 1500);
  const drift = total(e) - E0 - e.servoWorkTotal;
  // before the fix this was +1100 kJ/mol and the molecule came apart into four loose atoms
  assert.ok(Math.abs(drift) < 60, `drift ${drift.toFixed(1)} kJ/mol`);
  // a four-atom molecule swings its energy between motion and stretch every few fs, so one frame's
  // kinetic temperature says nothing: average over the vibrations, and check it held together
  let T = 0; for (let i = 0; i < 500; i++) { e.step(); T += e.temperature() / 500; }
  assert.equal(products(e), 'CH3');
  assert.ok(T < 3500, `left the sample at ${T.toFixed(0)} K on average`);
});

test('The outcome of a reactive drag no longer depends on the step size', () => {
  const at = dt => { const { e, h, E0 } = reactiveDrag(); e.dt = dt; dragOnto(e, h, [18, 11.1], 1200); return total(e) - E0; };
  const coarse = at(1), fine = at(0.1);
  assert.ok(Math.abs(coarse - fine) < 60, `1 fs gives ${coarse.toFixed(0)}, 0.1 fs gives ${fine.toFixed(0)}`);
});

test('The servo\'s work is measured, and it is all the energy a non-reactive drag adds', () => {
  const e = new Engine({ width: 36, height: 20, depth: 12, T: 300, thermostat: false });
  const c = e.addAtom('C', 10, 10, 0, { thermal: false });
  e.addAtom('H', 9.45, 9.1, 0, { thermal: false }); e.addAtom('H', 10.55, 9.1, 0, { thermal: false });
  e.touch(); e.minimize(600, .2); e.refresh();
  const E0 = total(e);
  dragOnto(e, c, [26, 10], 3000);                 // right across the chamber, nothing to react with
  assert.equal(products(e), 'CH2');
  const unexplained = total(e) - E0 - e.servoWorkTotal;
  assert.ok(Math.abs(unexplained) < 1, `${unexplained.toFixed(3)} kJ/mol unaccounted`);
  assert.ok(e.servoWorkTotal > 0, 'the drag did no work at all');
});

test('An insulated chamber keeps all the heat a reaction releases', () => {
  /* There is no switch for a third body. What carries a new bond's energy away is whatever is
     really there: other molecules, the walls, the surrounding bath. With the thermostat off the
     sample is insulated, and two H atoms that bond leave all 436 kJ/mol in the chamber. */
  const e = new Engine({ width: 24, height: 20, depth: 12, T: 300, seed: 48, thermostat: false });
  e.addAtom('H', 10, 10, 0, { thermal: true }); e.addAtom('H', 12.4, 10.1, 0, { thermal: true });
  e.touch(); e.refresh();
  const E0 = total(e);
  for (let i = 0; i < 12000; i++) e.step();
  assert.ok(Math.abs(total(e) - E0) < 15, `energy went from ${E0.toFixed(1)} to ${total(e).toFixed(1)} kJ/mol`);
});

test('Two H atoms alone cannot keep their bond; in a bath they can', () => {
  // H + H needs a third body to take the energy: in an empty insulated box there is none
  const run = opts => {
    let made = 0, met = 0;
    for (let s = 0; s < 12; s++) {
      const e = new Engine({ width: 24, height: 20, depth: 12, T: 300, seed: 31 + s * 17, ...opts });
      e.addAtom('H', 10, 10, 0, { thermal: true }); e.addAtom('H', 14, 10.3, 0, { thermal: true });
      e.touch(); e.refresh();
      let close = false;
      for (let i = 0; i < 12000; i++) { e.step(); if (Math.hypot(e.pos[0] - e.pos[3], e.pos[1] - e.pos[4], e.pos[2] - e.pos[5]) < 1.2) close = true; }
      if (close) met++;
      if (products(e) === 'H2') made++;
    }
    return { made, met };
  };
  assert.ok(run({ thermostat: false }).made <= 4, 'an empty insulated chamber should not hold a fresh bond together');
  // in a bath the atoms diffuse rather than fly, so not every pair meets in 12 ps; every pair that does, stays
  const bath = run({ thermostatMode: 'kelvin' });
  assert.ok(bath.met >= 7 && bath.made === bath.met, `bath: ${bath.made} bonded of ${bath.met} that met`);
});

test('A molecule can be built by hand whichever way the temperature is held', () => {
  /* The pointer brings a hydrogen to a methyl radical and lets go. The new C–H bond releases
     ~440 kJ/mol. The Kelvin bath takes it within a few hundred fs. Behind a wall heater or insulated, the fresh
     methane keeps it for a while and now and then throws the H back off, as a real one would. */
  const built = opts => {
    let made = 0;
    for (let s = 0; s < 8; s++) {
      const e = new Engine({ width: 30, height: 20, depth: 12, T: 300, seed: 12 + s * 23, ...opts });
      e.addAtom('C', 12, 10, 0, { thermal: false });
      e.addAtom('H', 13.08, 10, 0, { thermal: false });
      e.addAtom('H', 11.46, 10.93, 0, { thermal: false });
      e.addAtom('H', 11.46, 9.07, 0, { thermal: false });
      e.touch(); e.minimize(800, .1); e.refresh(); e.thermalize(300);
      const h = e.addAtom('H', 12, 14.5, 0, { thermal: false });
      e.touch(); e.refresh();
      dragOnto(e, h, justAbove(e, 0, 1.1), 1500);
      for (let i = 0; i < 4000; i++) e.step();
      if (products(e) === 'CH4') made++;
    }
    return made;
  };
  for (const [opts, need] of [[{ thermostatMode: 'kelvin' }, 6], [{ thermostatMode: 'wall', wallT: 300 }, 4], [{ thermostat: false }, 4]]) {
    const label = opts.thermostat === false ? 'off' : opts.thermostatMode, n = built(opts);
    assert.ok(n >= need, `${label} built methane only ${n} times in 8`);
  }
});

console.log(count + ' reactive-drag checks passed.');
