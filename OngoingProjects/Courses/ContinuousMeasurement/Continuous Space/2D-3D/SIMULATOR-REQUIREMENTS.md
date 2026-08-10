# What a continuous-measurement simulator should support, and whether we need full 3D

A plain-language report, grounded in a 2024-26 literature audit (see `related-papers/`). It answers
two questions: do we need to model motion in full three dimensions, and what is missing from a
simulator that wants to be useful to the people actually doing these experiments and the theorists
modeling them.

## The short answer on 3D

There are two separate worlds, and they give opposite answers.

When the particle's state keeps a simple bell-curve shape (a particle in an ordinary trap, gently
watched), the thing you track is just a small table of numbers: where it is, how fast it is moving,
and how uncertain you are about each. That table grows only gently as you add directions: one
direction, three directions, or all six directions of a spinning particle cost almost nothing extra.
In this world, "full 3D" is essentially free, and our engine's tracker already handles any number of
directions. So here the answer is yes, we need it, and it is already handled.

When the state takes a complicated shape (a particle that can sit in either of two wells, a strong
nonlinearity, a genuine superposition of two places), the small table is no longer enough and you
have to carry the full wave on a grid. A grid in three dimensions is genuinely expensive: the cost
grows like the cube of the resolution along each axis, so a finely-resolved 3D box is tens of
megabytes and about a tenth of a second per step, and a finely-resolved four-or-more-dimensional box
is simply out of reach. So in this world, full 3D is the costly thing.

The saving grace is that the complicated shape almost always lives in **one** direction. A double
well is a double well along a single axis; a nonlinearity sits in one coordinate; the other
directions stay simple bell curves that ride along. So you very rarely need a full 3D grid. You need
a small grid for the one or two hard directions, glued to the cheap table-of-numbers for all the
easy ones.

**Verdict.** Do not build a full 3D wave-on-a-grid as the main tool. Keep the 3D grid core available
for the rare case where all three directions are hard at once (we have it, and it is checked).
Invest instead in two things: the cheap any-number-of-directions tracker (mostly done), and the
hybrid that couples a small grid for the hard direction to the simple table for the rest. That
hybrid is what the physics actually asks for, and it is not built yet.

There is a little genuinely three-dimensional physics that a two-dimensional model cannot show, and
it is worth one demonstration rather than a rebuild: a relay, where watching one direction tells you
about a third only through a middle one, with a delay you can watch build up; and a symmetry where a
round trap watched along one line leaves a whole flat plane of motion invisible, not just a single
line.

## What the field actually runs (so the simulator has to match it)

The audit makes the working reality concrete, and it is not the textbook picture of one particle,
one detector, white noise.

- **The object is several coupled directions, not one.** Levitated-particle experiments now cool and
  watch three position directions at once, and often three spinning directions on top, all in a
  single setup (2205.10193 cools all six).
- **The directions are watched together and their readouts leak into each other.** The signal in a
  directional-force sensor lives entirely in the correlation between the x and y readouts; the
  single-direction spectra are blind to it (2307.06765). Real detectors also cross-talk, where the
  reading of one axis contaminates another, and recent schemes turn that leakage into a control knob
  (2506.17172). A model that treats each direction on its own is structurally unable to reproduce any
  of this.
- **The detector always misses some of the light, and that missing fraction is the deciding knob.**
  Whether the watched motion ends up quieter than the quantum floor or noisier than it is set almost
  entirely by how much of the scattered light the detector actually catches (efficiency around a
  third to nine-tenths in current work; 2505.10157, and the parent survey's Magrini anchor sitting
  just above the floor).
- **The loop is closed.** The estimate is fed back to push on the particle, by electronics, by a
  delayed light field, or all-optically; the delay and the extra noise the loop injects are part of
  the physics, not details (2506.21341, 2505.10157). Some schemes even skip the estimator and act on
  the raw record (2603.06370), which is worth knowing so we scope the tracker to where it earns its
  place.
- **The noise is realistic and structured.** Gas collisions and photon recoil heat the particle,
  the trap damps it, the feedback light carries phase noise with its own frequency shape, and the
  whole setup drifts from shot to shot. None of this is white noise.
- **Verification looks at the future of the record, not only the present.** The sharpest estimate of
  where the particle was uses the part of the record that came afterward, a step beyond real-time
  tracking that a recent experiment demonstrated on a vibrating resonator (2510.16754).

The software the community reaches for sets the performance floor: the standard toolkit runs many
noise histories at once, on a graphics card, in a form you can differentiate so control recipes can
be tuned automatically (QuTiP 5, 2412.04705; dynamiqs). And the method literature adds one useful
warning: the total cost of this kind of simulation is not one number but three, the memory per run,
the time per run, and the number of runs you must average, and clever choices trade one against the
other rather than reducing all three (2606.13779).

## What should be added, in priority order

1. **Close the loop (feedback).** Let the estimate push back on the trap, with a realistic delay.
   Every cooling and stabilization experiment is closed-loop. We have the estimator; we do not yet
   act on it.

2. **Model the light the detector misses (imperfect detection).** This turns the tracked state fuzzy
   in a specific, quantifiable way, and it is the single knob that separates the headline "below the
   quantum floor" result from "above" it. Today this lives only in our small tracker, not in the
   grid, because a grid of a single pure wave cannot represent a fuzzy state; fixing it means either
   a heavier grid that carries the fuzziness directly, at small sizes, or an extra pretend detector
   that no one reads.

3. **Put in realistic noise.** Gas heating and damping, photon recoil, structured (frequency-shaped)
   feedback and phase noise, and shot-to-shot drift in the setup. White noise is the exception in the
   lab, not the rule.

4. **Couple the directions, with their own strengths and efficiencies, and surface the cross-
   readouts.** The correlation between two directions' readouts is the actual measured signal in
   directional sensing, and it is also what lets watching one direction narrow another. Most of this
   is already in the tracker; the cross-readouts need to be exposed as outputs, and per-direction
   detection efficiency needs to be a knob.

5. **Add smoothing.** Use the later part of the record to sharpen the earlier estimate. This is now a
   standard verification step and gives the best-possible reconstruction.

6. **Package time-dependent recipes.** A trap that changes in time, a measurement strength that ramps,
   sudden frequency jumps (the expansion protocols that dominate the current experimental literature).
   The single step already allows this; it needs to be packaged as named protocols rather than
   hand-assembled each time.

7. **Keep the complicated-shape capability sharp where it matters.** Double wells, cubic traps, and
   the preparation of cats and number states are exactly where the grid earns its place over the
   small tracker, because the tracker cannot describe two lumps at once. This is the grid's reason to
   exist alongside the tracker.

8. **Let the truth and the estimator disagree (model mismatch).** Real experiments never know their
   parameters exactly, and the estimate degrades when the model is wrong. The recent machine-learning
   filters exist entirely to be robust to this; at minimum the simulator should let you feed one
   model's record into a different model's estimator and watch what happens.

9. **Make it fast.** Averaging over noise means running many independent histories, which is trivially
   parallel; this should run batched, on a graphics card, and ideally in a differentiable form so
   control recipes can be optimized rather than guessed. This is the baseline the standard tools
   already meet.

## Where our engine stands against this list

Already in place: the any-number-of-directions tracker (cheap, and checked at up to three directions
with everything left as symbols); the full grid for hard shapes in one, two, or three directions;
exact-shape checks against the cases with known answers; the physics where coupling makes an
unwatched direction visible; and the dark-direction physics where a symmetry hides a whole mode.

Missing, in the order above: feedback actuation, imperfect detection on the grid, realistic and
structured noise, exposed cross-readouts and per-direction efficiency, smoothing, packaged time-
dependent protocols, model mismatch, and batched/graphics-card/differentiable execution. The single
highest-value addition is the hybrid architecture (item 7 crossed with the 3D verdict above): a
small grid for the one hard direction coupled to the cheap table for the rest, which is what turns
this from a 1D-and-2D demonstration into the tool the levitated-particle and membrane experiments
would actually run.

Files saved for future use are in `related-papers/` (eleven new sources plus four from the search
sweep, TeX extracted, titles verified, manifest in `related-papers/README.md`); the parent corpus in
`../related-papers/` holds the experimental frontier this report builds on.
