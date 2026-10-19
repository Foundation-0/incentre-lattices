# program.md -- autoresearch loop for the folding octahedra packing question

QUESTION. Cube cells of side 2 tile space; each holds 6 dots at c +- k e_i (k = 2 - sqrt2), ball radius
R = sqrt2 - 1; bars of length 2R join the 12 octahedron edges per cell and each +e_i dot to the -e_i dot of
the next cell. The framework folds. What is the MAXIMUM packing density reachable by a continuous motion that
keeps every bar length and never overlaps balls (lattice may shear), and what is the end structure?
Known: cubic-lattice [111] checkerboard fold reaches 27*pi*sqrt2/256 = 0.4686 at 60 deg (6 contacts per ball).

## Files
- `verify.py` LOCKED trusted checker (sha256 in `verify.sha256`). NEVER edit it. `run.sh` refuses to run if it changes,
  and checks again after search.py and before it logs or keeps a result (against `verify.sha256` and, in a git checkout,
  against HEAD). The hash lock is a tripwire against accidental edits, not a sandbox: search.py runs with full write
  access to this folder, so a result in `best/` counts only after it is replayed with the pinned checker from a clean
  checkout.
- `search.py` the ONLY file you edit. It must write `out/path.npz` (keys m, Ls (T,3,3), Xs (T,6m^3,3)) and
  keep a one-line `DESC = "..."` describing the current idea.
- `run.sh` one experiment (3300 s cap): verify hash, run search.py, verify hash and path, append a line to
  `results.tsv`, copy to `best/` (path.npz, search.py, verify.json) if valid AND max density beat the best so far.
- `results.tsv` the log: timestamp, search hash, valid, final density, max density, end contacts, flex count, m, exit, desc.
- `search.py.best` the last accepted version of search.py. `notes.md` your ideas, findings, dead ends.
- Metric: `max_density` of a VALID path (scalar, higher is better). Invalid = worth nothing.

## The loop (budget: 60 iterations by default; one iteration = one ./run.sh)
1. Read `results.tsv`, `notes.md`, `best/verify.json` (and out/search.log for the last run).
2. Form ONE hypothesis (a change to search.py with an expected effect). Write it in `notes.md` first.
3. Edit `search.py` (update DESC). Run `./run.sh`.
4. KEEP if the line is `valid=true` and max_density improved on the best: `cp search.py search.py.best`.
   Otherwise DISCARD: `cp search.py.best search.py` (the failed line stays in results.tsv, that is the record).
   A run that is valid, equal density but gives a simpler or faster search.py may also be kept; say so in notes.md.
5. Repeat. Never run two experiments at once (they share out/).

## Rules
- Never edit verify.py, verify.sha256, test_verify.py (run `.venv/bin/python test_verify.py` if in doubt).
- Never claim a density that verify.py rejected. Quote numbers only from results.tsv / best/verify.json.
- `first_order_jammed` / `flex_bars_plus_contacts = 0` is only a FIRST-ORDER statement (tests are linear). Say
  "first-order rigid", never "proven maximal". Second-order flexes and other branches may exist.
- Different supercells constrain the answer differently: results for one m are lower bounds on what a periodic
  structure of that size can reach. Try m = 2, 3, 4 (checkerboard needs even m; the engine does not).
  Note m=3,4 are slow (dense least squares); lower HX cost.
- search.py stops after MAX_STEPS stage-B steps or STALL_KICKS kicks without progress; WALL_LIMIT (3000 s) is a safety
  limit that makes the run INCOMPLETE (exit 2, no path). A full m=2 run takes about 4 minutes; larger m may need a
  lower MAX_STEPS to finish inside WALL_LIMIT (notes.md iter3: about 200 steps in 240 s at m=4). Keep WALL_LIMIT
  below run.sh's 3300 s cap.
- No commits, no posting. Everything stays in this folder.

## Idea menu (non-exhaustive)
- Different starts: skip the cubic fold (FOLD_END=0), stop it early (e.g. 20 or 40 deg), other seeds (SEED).
- Different supercells (M = 1, 2, 3, 4), different symmetry-breaking: add random tangent kicks at the start,
  perturb the start by a small random null-space motion.
- Different tilt systems: replace stage A by another analytic family (all octahedra tilt about [111] in phase,
  columns alternating, [100]-axis tilts, mixed axes); any start must lie on the manifold (the engine retracts).
- Sheared lattices: let the LP run free (default). Try objective variants (e.g. maximise a different functional
  than det L, penalise anisotropy) to steer to other basins, then finish with the pure density descent.
- Contact sliding: change CONTACT_MARGIN, step caps HX/HL, tie-break weight eps, kick strategy (random cone
  directions vs biased), retraction details.
- Random restarts of the whole pipeline with many seeds; keep the best valid.
- When two different runs end at the same density, compare lattices and contact graphs: same structure?

## Stopping rule
Stop when ANY holds: 60 iterations are used; or the last 12 iterations did not improve max_density by > 1e-4
(plateau); or the best end state is first-order rigid and 10+ varied restarts (seeds, m, starts) all land at or
below it. On stopping, characterise the best end state and write it up in notes.md:
 - lattice parameters of the END state (lengths, angles, det, density) and, from the dot positions, the
   primitive cell (is it trigonal / cubic / ...?); the number of distinct ball environments;
 - contact graph: degree per ball, coordination sequence (BFS shells over the contact graph in the
   infinite periodic structure), contact-graph vs bar graph, cycle statistics;
 - compare with known packings: fcc/hcp (0.7405, 12 contacts), bcc (0.6802), sc (0.5236), diamond (0.3401),
   and with RCSR nets (name, vertex symbol, coordination sequence, e.g. "sms" 6/3/h2 for the 6-contact end);
 - the first-order verdict from best/verify.json and what remains open.
