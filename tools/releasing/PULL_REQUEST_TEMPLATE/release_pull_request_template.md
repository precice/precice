## How to work with this template

* [ ] assign a release manager, who takes care of the process
* [ ] assign each point below to a responsible person, before you continue. Use `@member`.

Only the release manager should update this post (even tickboxes, due to race conditions in editing). Everybody else should comment on the PR with the progress.


## Pre-PR steps

* [ ] Look over [PRs](https://github.com/precice/precice/pulls?q=is%3Apr+no%3Amilestone+is%3Aclosed) and [Issues](https://github.com/precice/precice/issues?q=is%3Aissue+no%3Amilestone+is%3Aclosed) without an assigned version. (all)
* [ ] Look over entries in [`docs/changelog`](https://github.com/precice/precice/blob/develop/docs/changelog) (all)
   * Add missing entries, if necessary
   * Fix wording and tense
* [ ] Make sure you have the latest `develop` and `main` branches locally.
* [ ] Merge `main` to `develop` ( This should result in no commits )
* [ ] Check code base w.r.t code formatting (run `pre-commit run -va`)
* [ ] Update the list of operating systems for the package generation in `.github/workflows/release.yml`
* [ ] Create branch `release-vX.Y.Z` from develop: `git switch -c release-vX.Y.Z develop`
* [ ] Bump the version: `tools/releasing/bumpversion.sh MAJOR.MINOR.PATCH`
* [ ] Look over the updated `CHANGELOG.md` of the release branch (all)
   * Check for merged lines
   * Add things, if necessary
   * Fix wording and tense
   * Sort the entries lexicographically
   * Merge or remove entries that are fully contained in the release
   * Extract summary
* [ ] Verify the version changes in:
   * [ ] `CHANGELOG`
   * [ ] `CMakeLists.txt`
   * [ ] `tools/releasing/packaging/debian/changelog`
* [ ] Commit the version bump: `git commit -m "Bump version to X.Y.Z"`
* [ ] Push the release branch to the precice repository: `git push -u upstream release-vX.Y.Z`
* Prepare independent releases
   * [ ] [Python bindings](https://github.com/precice/python-bindings/blob/develop/docs/ReleaseGuide.md)
   * [ ] (if necessary!) [MATLAB bindings](https://github.com/precice/matlab-bindings/blob/develop/docs/ReleaseGuide.md)
   * [ ] (if necessary!) [JULIA bindings](https://github.com/precice/PreCICE.jl)

## Step by step guide

* [ ] Open PR from `release-vX.Y.Z` to `main` (use [this template](https://github.com/precice/precice/blob/develop/tools/releasing/PULL_REQUEST_TEMPLATE/release_pull_request_template.md))
* [ ] Trigger the system tests using the `trigger-system-tests` label ([`release_test` suite](https://github.com/precice/tutorials/blob/develop/tools/tests/tests.yaml)). After any force-push, remove and add the label again.
* [ ] Do any additional regression tests using the release branch (specific revision) _list below :arrow_down:_ (all)
* [ ] Fix potential problems in develop (all)
* [ ] Rebase the release branch on develop to pull in fixes
* [ ] Draft release notes
* [ ] Write a draft "blog post" on [Discourse](https://precice.discourse.group/)
* [ ] Update documentation (all)
   * [ ] Update [XML configuration reference](https://github.com/precice/precice.github.io/blob/master/_includes/xmlreference.md)
   * [ ] Update version in [precice/precice.github.io](https://github.com/precice/precice.github.io) `_config.yml`
* [ ] Approve the PR with at least two reviews (all)
* [ ] Merge PR to `main`: `git merge --no-ff -m "Release vX.Y.Z" release-vX.Y.Z`
* [ ] Create an annotated tag on `main`: `git tag -a -m "preCICE version X.Y.Z" vX.Y.Z main`
* [ ] Verify the tag: `git describe --tags main`. It should be exactly `vX.Y.Z`
* [ ] Switch to `develop` and merge `main` back into it: `git merge --no-ff -m "Merge release back"`
* [ ] Verify the tag on develop: `git describe --tags develop`. It should start with `vX.Y.Z-1-` (i.e. tag plus the merge commit).
* [ ] Triple check that you haven't messed anything up. You can always discard local changes using `git reset --hard upstream BRANCH` or by cloning the precice repository again and start from scratch.
* [ ] Push `main` and the `vX.Y.Z` tag: `git push upstream main`, `git push upstream v3.3.0`
* [ ] Push `develop`: `git push upstream develop`
* [ ] Wait for the release pipeline
  * [ ] [To create a new draft release on GitHub](https://github.com/precice/precice/releases)
  * [ ] To automatically generate packages for the latest Debian and the two latest Ubuntu LTS versions.
* [ ] Write release text
* [ ] Publish the GitHub release


## Regression Tests

Use the following branches:
* precice `release-vX.Y.Z`
* pyprecice `python-bindings-vX.Y.Z.1`
* matlab-bindings `matlab-bindings-vX.Y.Z.1`
* rest `main`(`master`)

Run all these tests manually on your system. If you succeed, please write a comment with the revisions of the components that you used below. Example: https://github.com/precice/precice/pull/507#issuecomment-530432289 and update the table.

Tests covered by the system tests: see `release_test` in [`tests.yaml`](https://github.com/precice/tutorials/blob/develop/tools/tests/tests.yaml) (and the respective job summary).

| State | Success | Failure | Skipped |
| --- | --- | --- | --- |
| Write | `:o:` | `:x:` | `:fast_forward:` |
| Read | :o: | :x: | :fast_forward: |

| State | Tester | Test |
| --- | --- | --- |
| | | [quickstart](https://github.com/precice/tutorials/tree/master/quickstart) fluid-openfoam - solid-cpp |
| | | [perpendicular-flap](https://github.com/precice/tutorials/tree/master/perpendicular-flap) fluid-openfoam - solid-dune |
| | | [perpendicular-flap](https://github.com/precice/tutorials/tree/master/perpendicular-flap) fluid-nutils - solid-calculix |
| | | [flow-over-heated-plate](https://github.com/precice/tutorials/tree/master/flow-over-heated-plate) fluid-openfoam - solid-openfoam parallel |
| | | [flow-over-heated-plate](https://github.com/precice/tutorials/tree/master/flow-over-heated-plate) fluid-openfoam - solid-fenics parallel |
| | | [flow-over-heated-plate](https://github.com/precice/tutorials/tree/master/flow-over-heated-plate) fluid-openfoam - solid-nutils |
| | | [flow-over-heated-plate-nearest-projection](https://github.com/precice/tutorials/tree/master/flow-over-heated-plate-nearest-projection) fluid-openfoam - solid-openfoam |
| | | [flow-over-heated-plate-steady-state](https://github.com/precice/tutorials/tree/master/flow-over-heated-plate-steady-state) fluid-openfoam - solid-codeaster |
| | | [heat-exchanger](https://github.com/precice/tutorials/tree/master/heat-exchanger) fluid-(inner+outer)-openfoam - solid-calculix |
| | | [partitioned-elastic-beam](https://github.com/precice/tutorials/tree/master/partitioned-elastic-beam) dirichlet-calculix - neumann-calculix |
| | | [partitioned-heat-conduction](https://github.com/precice/tutorials/tree/master/partitioned-heat-conduction) fenics - nutils |
| | | [partitioned-heat-conduction-complex](https://github.com/precice/tutorials/tree/master/partitioned-heat-conduction-complex) fenics- fenics |
| | | [partitioned-pipe](https://github.com/precice/tutorials/tree/master/partitioned-pipe) fluid1-openfoam-pimplefoam - fluid2-openfoam-sonicliquidfoam |
| | | [elastic-tube-3d](https://github.com/precice/tutorials/tree/master/elastic-tube-3d) fluid-openfoam - solid-calculix |
| | | MATLAB / MATLAB [ODEs](https://github.com/precice/matlab-bindings/tree/develop/tutorial) |
| | | Solverdummy [Fortran module](https://github.com/precice/fortran-module/tree/develop/examples/solverdummy) |
| | | Solverdummy [Python](https://github.com/precice/python-bindings/tree/develop/solverdummy) |
| | | Solverdummy [MATLAB](https://github.com/precice/matlab-bindings/tree/develop/solverdummy) |
| | | Solverdummy [Julia](https://github.com/precice/PreCICE.jl/tree/develop/solverdummy) |
| | | Alya |
| | | SuperMUC |

## Post-release

* [ ] Update version specific documentation
* [ ] Flag [Arch Linux AUR package](https://aur.archlinux.org/packages/precice) and dependants as out-of-date.
* [ ] Update Spack recipe
* [ ] Update `pyprecice` Spack
* [ ] Update Website:
    * [ ] Bump version in [`_config.yml`](https://github.com/precice/precice.github.io/blob/master/_config.yml)
    * [ ] Update the [XML reference](https://github.com/precice/precice.github.io/blob/master/_includes/xmlreference.md) using `binprecice md`
    * [ ] Look over the [Roadmap](https://www.precice.org/fundamentals-roadmap.html) and update entries.

### Release new version for bindings (to ensure compatibility with newest preCICE version)

- [ ] [Fortran module](https://github.com/precice/fortran-module/compare/master...develop)
- [ ] [MATLAB bindings](https://github.com/precice/matlab-bindings/blob/develop/docs/ReleaseGuide.md)
- [ ] [python bindings](https://github.com/precice/python-bindings/blob/develop/docs/ReleaseGuide.md)
- [ ] [Julia bindings](https://github.com/precice/PreCICE.jl)
- [ ] [Rust bindings](https://github.com/precice/rust-bindings)

### (only if breaking changes) Open PRs or issues `develop -> main` for all adapters

- [ ] [calculix-adapter](https://github.com/precice/calculix-adapter/compare/master...develop)
- [ ] [code_aster-adapter](https://github.com/precice/code_aster-adapter/compare/master...develop)
- [ ] [comsol-adapter](https://github.com/precice/comsol-adapter/compare/master...develop)
- [ ] [dealii-adapter](https://github.com/precice/dealii-adapter/compare/master...develop)
- [ ] [fenics-adapter](https://github.com/precice/fenics-adapter/compare/master...develop)
- [ ] [fluent-adapter](https://github.com/precice/fluent-adapter/compare/master...develop)
- [ ] lsdyna-adapter?
- [ ] mbdyn-adapter?
- [ ] [openfoam-adapter](https://github.com/precice/openfoam-adapter/compare/master...develop)
- [ ] [su2-adapter](https://github.com/precice/su2-adapter/compare/master...develop)

### (only if breaking changes) Open PRs or issues `develop -> main` for all other tools

- [ ] [aste](https://github.com/precice/aste/compare/master...develop)
- [ ] [elastictube1d](https://github.com/precice/elastictube1d/compare/master...develop)
- [ ] [tutorials](https://github.com/precice/tutorials/compare/master...develop)

### Marketing

* [ ] Finalize post on [Discourse](https://precice.discourse.group/)
* [ ] Write on [Matrix](https://matrix.to/#/#precice_lobby:gitter.im?web-instance[element.io]=app.gitter.im)
* [ ] CFD Online:
     * [ ] [News](https://www.cfd-online.com/Forum/news.cgi/form/0)
     * [ ] [Forum](https://www.cfd-online.com/Forums/main/261393-precice-releases.html)
* [ ] NADigest
* [ ] [NAFEMS](https://www.nafems.org/mynafems/submitnews/) (needs account, appears in [Upcoming Industry Events](https://www.nafems.org/events/industry-events/))
* [ ] Post on [Bluesky](https://bsky.app/profile/precice.org)
* [ ] Post on [Mastodon](https://fosstodon.org/@precice)
* [ ] Post in LinkedIn. Relevant places:
    * [ ] [preCICE company](https://www.linkedin.com/company/precice)
    * [ ] [OpenFOAM group](https://www.linkedin.com/groups/1920608/)
    * [ ] [HPC group](https://www.linkedin.com/groups/87791/)
    * [ ] Personal

See the [advertising a workshop](https://precice.org/precice-workshop-organizing.html#advertising) section for additional ideas.

### Misc

* [ ] Update the [PR template](https://github.com/precice/precice/blob/develop/tools/releasing/PULL_REQUEST_TEMPLATE/release_pull_request_template.md)

To open a new PR with this template, use this [PR template query](https://github.com/precice/precice/compare/new?template=release_pull_request_template.md)
