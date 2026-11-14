# The golden regression battery

The battery runs a fixed set of fits and computations with fixed data and seeds, saves everything they return, and compares two installs of RiskMap task by task. A refactor that is meant to change nothing (class C in the release plan) must leave every task identical bit for bit. A change that moves fits must move only the tasks it should, and the tolerance mode of `compare.R` then says by how much. The package's tests check behaviour to tolerances and cannot see a change in the last bits, which is why the battery exists beside them.

It moved here from the MCML paper's repository (`theory/family-object-golden.R` and `-compare.R`), where it checked every PR of #253. It is listed in `.Rbuildignore`, so it is not in the tarball and CRAN never runs it. It runs before a merge, not in CI, since fits reproduce bit for bit only on one machine type.

## Files

| File | What it is |
|---|---|
| `run.R` | The tasks and the runner. One output file per task. |
| `compare.R` | The comparison of two arms. |
| `groups.txt` | The tasks in seven lines, for parallel runs. |
| `data/surveys.rds` | Two surveys from the MCML paper's studies, saved exactly (see below). |
| `hpc/push.sh` | Copies a commit's archive to Wayland and checks its md5 sums (run on the laptop). |
| `hpc/submit.sh` | Submits the install and the run of one arm (run on the Wayland login node). |
| `hpc/install.sh`, `hpc/run.sh` | The two Slurm jobs that `submit.sh` submits. |

## The tasks

`Rscript --vanilla tools/golden/run.R --list` lists them.

- **F01 to F29, fits.** Binomial, Poisson and Gaussian fits with and without a nugget, an `re()` term, households, held nuggets, integrated and sampled noise, named links and links given as functions, an imperfect test, and the edge cases (no positives, separation, sparse and nearly all positive surveys). A fit task keeps the fit, `check_fit()` and the printed fit, summary and check.
- **D tasks, downstream of five fits** (logit, Poisson, probit, integrated nugget, Gaussian). Each refits its base fit and keeps one of `predict()` and `simulate_glgpm()`, `parametric_bootstrap()`, `profile()`, `check_fit(extra = )`, or `assess_prediction()` with `assess_simulation()`.
- **L01 to L03, low-level.** Links, observation terms, likelihoods and engines evaluated at fixed inputs through internals, and the adaptive MCML driver at a bound. They reach internals through `asNamespace("RiskMap")` and follow whichever entry points the commit has, so older commits can still be run.
- **`--smoke`.** Thirteen of them (fits of each family, a nugget, an `re()` term, a named link, both saved surveys, two predictions and the low-level tasks), about a minute on one Wayland core.

Each output file is a record with the task's status, its value (what is compared), the time and garbage-collection time, the environment (commit, RiskMap version and path, R, BLAS, host, thread settings, the versions of every loaded package) and `fits`, the fits the task made at its top level, which the tolerance mode reads. A failed task is saved with its message and compared like any other.

## What makes two arms comparable

Bits reproduce only when everything but the RiskMap commit is the same.

- **One BLAS and OpenMP thread.** `run.R` stops unless `OPENBLAS_NUM_THREADS` and `OMP_NUM_THREADS` are 1. Use `Rscript --vanilla`, because a user's `~/.Renviron` that sets these variables overrides the shell's values, and `--vanilla` is what skips it.
- **One machine type.** Laptop arms compare only with laptop arms, and Wayland arms only with Wayland arms of one node type (64core256G, with comp06 excluded, which `submit.sh` sets).
- **The same R, container and dependencies.** On Wayland every arm runs in the same container and loads its dependencies from the same libraries. `compare.R` prints the environments and warns when package versions differ.
- **An installed package.** Install each commit from a clean `git archive` with `R CMD INSTALL` into a library of its own. Never use `devtools::load_all()`, which compiles without optimisation and reuses object files.
- **The same battery.** Run both arms with the battery of the newer commit, so that they have the same tasks. On Wayland that is the third argument of `submit.sh`.

## An arm on Wayland

The battery's home is `~/storage/riskmap-golden` (`GOLDEN_HOME`).

| Folder | What it holds |
|---|---|
| `src/<sha>` | The clean archive of a commit, unpacked by `push.sh`. |
| `library/<sha>` | RiskMap installed from it, with `<sha>.install.log`. |
| `library/deps` | Imports of RiskMap that the shared library `~/storage/riskmap-sim-env/library` lacks. `install.sh` adds any new import here, so that every arm loads the same versions. |
| `arms/<arm>` | The outputs of an arm, one `.rds` file per task. |
| `logs` | The jobs' logs. |

From a RiskMap checkout on the laptop, copy the commits, then submit from the Wayland login node.

```bash
tools/golden/hpc/push.sh <sha> <base-sha>
```

```bash
timeout 120 ssh fronter2@wayland-2022.hec.lancaster.ac.uk '~/storage/riskmap-golden/src/<sha>/tools/golden/hpc/submit.sh <sha> <arm>'
```

`submit.sh <commit> <arm> [<battery> [<tasks>]]` installs the commit unless its library exists and runs the battery after the install. Installs run one at a time (all are named `golden-install`, with `--dependency=singleton`), so two sessions that submit arms of one commit never write its library at once. `<battery>` is the commit whose `tools/golden` runs, by default the commit itself. For the base of a PR, give the PR's commit, so that both arms have the PR's tasks. `<tasks>` is `all` (the default, one job of about 20 minutes), `smoke`, a comma-separated list, or `groups` for an array over the lines of `groups.txt`, whose elements take only a few minutes each, so keep it for when the wait matters. The install log should report six compiled source files.

Check the jobs with `squeue -u fronter2` and `sacct -j <job>`, and confirm the arm by counting its files rather than by the job's state. An element that timed out or failed can be resubmitted as it is, since tasks with an output are skipped. Bring the arm back with `rsync -a` into a folder outside the checkout and compare it there.

```bash
rsync -a fronter2@wayland-2022.hec.lancaster.ac.uk:storage/riskmap-golden/arms/<arm>/ <arms>/<arm>/
```

## An arm on the laptop

Install the commit from a clean archive into a library of its own, then run the battery.

```bash
git archive --format=tar --prefix=RiskMap/ <sha> | tar -x -C <scratch>
```

```bash
R CMD INSTALL --preclean --library=<libraries>/<sha> <scratch>/RiskMap
```

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript --vanilla tools/golden/run.R --lib=<libraries>/<sha> --out=<arms>/<arm> --smoke
```

A full arm took 30 to 50 minutes in one process on the laptop. Run it detached, with its memory capped, for example `tmux new -d -s golden-<arm> 'systemd-run --user --scope -p MemoryMax=8G bash -c "ulimit -v 6000000; OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 Rscript --vanilla tools/golden/run.R --lib=... --out=... > <log> 2>&1"'`, or run the lines of `groups.txt` in two or three such processes.

## Comparing two arms

```bash
Rscript --vanilla tools/golden/compare.R --a=<arms>/<base> --b=<arms>/<pr> --report=<arms>/compare-<base>-<pr>.txt
```

Arm a is the reference. `compare.R` takes about five seconds and needs only sf. It exits with status 1 if any task fails the comparison.

- **Identity**, the default. A task passes when its value and its printed text are identical in both arms, after timings, environments and the environments of functions are removed (functions are compared by their deparsed text). This is the test for a class C change.
- **Declared paths**, `--declared=<file or list>`. For a PR that changes part of a fit on purpose, for example by adding a field. The paths it declares are removed from both arms before the comparison and reported as same, differs or only in one arm, and everything else must still be identical. Each entry is a regular expression matched against paths as the report prints them, optionally preceded by a regular expression for the task names it applies to. Entries are lines of a file, where `#` starts a comment, or a comma-separated list.

  ```
  # A version stamp on every fit, and its line in the printed summaries
  ^/fit/version$
  ^F ^/text/2$
  ```

- **Tolerance**, `--tolerance=<k>`. For a PR that changes the optimiser's path or the random streams, so that every number may move. A task passes when
  - the statuses agree, and the values have the same structure (the same paths, classes and lengths, the same missing values),
  - every fit the task recorded has the same outcome, certification, number of stages and parameters in both arms, its estimates move by at most k standard errors of arm a (or of arm b for a parameter that has none in a), on the scale of `fit$estimate` (log variances and scale), its standard errors change by at most a fraction `--se-tolerance` and its log-likelihood by at most `--loglik-tolerance` (both k by default), and an estimate or standard error that is finite in one arm is finite in the other,
  - strings and the printed text agree once their numbers are masked,
  - and its other numbers differ only if one of its fits moved.

  A parameter without a standard error in either arm (at a boundary, for example) and a change in the proposals of the stages are reported as notes, but a fit whose estimates move with no standard error at all fails. Declared paths are removed from the recorded fits too (at `/fit` for a fit task and `/base` for a downstream task), so a declared field alone does not make a fit count as moved. Downstream numbers that follow from a moved fit are counted, with the largest relative difference, but not judged. The fits of D tasks are recorded only by this `run.R`, so take both arms with it. State the tolerances in the PR and why they fit the change.

- **`--allow-new`** lets arm b have tasks that arm a lacks, for an arm with new tasks against an older reference. They are listed and not compared. A task missing from arm b always fails, as does a comparison in which no task is in both arms.

Each task gets one line, `identical`, `DIFFERENT value`, `DIFFERENT text`, `MISSING in a`, `STATUS done against failed`, `within tolerance` or `OUTSIDE tolerance`. A difference is followed by up to twelve paths, such as `/fit/estimate/phi` or `/fit@optimizer/mcml_outcome` for an attribute, each with its largest absolute and relative difference, and by the first differing line of each printed text. The environments of the two arms come next (commit, RiskMap, R, BLAS, hosts, and the packages whose versions differ), then the names removed as timings and a last line such as `Tasks: 57 compared, 0 differ, all identical: TRUE`.

## The reference

The reference for `main` at 6010772 is `arms/kernel-a2b307e-wayland` on Wayland (a2b307e has the tree of 6010772), taken by the MCML paper's battery on a 64core256G node. `arms/golden-e1fd248`, taken by this battery as #314 merged it, is identical to it in all 57 tasks and also records the fits of the downstream tasks and the package versions, so it is the better reference for the tolerance mode. A PR that changes no number keeps the reference. A PR that moves numbers or adds tasks takes a new one at its merge commit, noted here.

## Adding a task

- Add it to `fits`, `downstream` or the low-level list of `run.R` with fixed seeds, and to a line of `groups.txt`. Keep it to a minute or two on one core where possible.
- Never change the call, data or seed of an existing task, since arms are compared by task name. Give a changed task a new name instead.
- Keep data in the script or, when it cannot be simulated, small in `data/`.
- Take a new reference that includes the task, by running the base commit with the new battery.

## The saved surveys

`data/surveys.rds` holds two surveys exactly as the MCML paper's repository generated them, since regenerating them through BLAS on another machine could change their last bits and with them F26 to F28. Each entry has the data frame, its EPSG code and its source.

- `n05_2`, for F27 and F28, is survey n05-2 of `nugget/datasets.R` (seed 64010000), 200 sites with 20 tested each, a spatial variance of 3, a nugget of 0.15 and a prevalence of 3%.
- `households`, for F26, is survey r01-1 of `integrated/re-setup.R` (seed 73001001), 795 households at 200 locations with a household effect of variance 0.15 and a prevalence of 3%.
