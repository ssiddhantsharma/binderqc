# Contributing to binderqc

Thanks for taking an interest. This is a small, deliberately dependency-light
package, so the bar for changes is mostly about keeping it that way.

## Reporting a problem

Open an issue: https://github.com/ssiddhantsharma/binderqc/issues

For a wrong or surprising number, the useful report has three things — the
structure (or a PDB id), the exact command, and what you expected instead. The
numbers are all deterministic given coordinates, so a reproducer almost always
pins it immediately.

If a metric looks wrong, it is worth checking [METHODS.md](METHODS.md) first:
several columns are deliberately proxies (`n_hbonds`, `buried_unsat_polar` and
`interface_packing` especially), and a few depart from the source they are named
after in ways that are documented there.

## Asking a question

Open an issue with the `question` label, or start a discussion. Questions about
what a number means are welcome and usually mean METHODS.md needs to be clearer —
those get treated as documentation bugs.

## Contributing code

1. Fork, branch, and open a pull request against `main`.
2. `pip install -e ".[test]"`, then `pytest` and `ruff check .`. CI runs both.
3. Keep the dependency set as it is — biotite, numpy, pandas. A new third-party
   dependency needs a reason in the PR; a new *model* dependency almost certainly
   belongs in the caller instead (`grippability_consensus` takes a learned score
   as an argument rather than invoking one, for exactly this reason).

### If you add or change a metric

- Document it in `METHODS.md`: the formula, the source, and **where it departs
  from that source**. The departures are the part people need.
- Add a test that would fail if the implementation were wrong, not only one that
  checks the sign. `test_constant_tables.py` is the pattern — it is
  mutation-checked, so a transposed digit in a published scale fails the suite.
- If it is a port or reimplementation of a published method, say so and cite it.
  Where an external reference exists, a correlation check against it is worth far
  more than an assertion about the output's range: `tests/pisa_correctness.py`
  does this for interface area against PDBePISA.

### Exposure conventions

Any new structural metric has to declare whether it reads exposure on the
**complex** or on the **free** binder. Mixing the two silently is the failure mode
this package is most careful about — see
[METHODS.md](METHODS.md#the-two-exposure-conventions).

## Scope

binderqc computes physical properties of a finished complex. It does not fold,
score confidence, or design. Proposals that would add a model dependency, require
a GPU, or need network access at scoring time are out of scope — not because they
are bad ideas, but because being CPU-only and model-agnostic is what makes the
panel comparable across predictors.

## Licence

MIT. By contributing you agree your contribution is licensed under it.
