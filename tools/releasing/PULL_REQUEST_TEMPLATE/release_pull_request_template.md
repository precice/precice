## How to work with this template

* [ ] assign a release manager, who takes care of the process
* [ ] assign each point below to a responsible person, before you continue. Use `@member`.

Only the release manager should update this post (even tickboxes, due to race conditions in editing). Everybody else should comment on the PR with the progress.


## Pre-PR steps

* [ ] Look over [PRs](https://github.com/precice/precice/pulls?q=is%3Apr+no%3Amilestone+is%3Aclosed) and [Issues](https://github.com/precice/precice/issues?q=is%3Aissue+no%3Amilestone+is%3Aclosed) without an assigned version. (all)
* [ ] Look over entries in [`docs/changelog`](https://github.com/precice/precice/blob/develop/docs/changelog)) (all)
   * Add missing entries, if necessary
   * Fix wording and tense
* [ ] Make sure you have the latest `develop` and `main` branches locally.
* [ ] Merge `main` to `develop` ( This should result in no commits )
* [ ] Check code base w.r.t code formatting (run [`precice/tools/formatting/check-format`](https://github.com/precice/precice/blob/develop/tools/formatting/check-format)) and reformat if required (run [`precice/tools/formatting/format-all`](https://github.com/precice/precice/blob/develop/tools/formatting/format-all))
* [ ] Update the list of operating systems for the package generation in `.github/workflows/release.yml`
* [ ] Create branch `release-vX.Y.Z` from develop. If needed, `git rebase develop`.
* [ ] Run `tools/releasing/bumpversion.sh MAJOR.MINOR.PATCH` to bump the version
* [ ] Look over the updated `CHANGELOG.md` of the release branch (all)
   * Check for merged lines
   * Add things, if necessary
   * Fix wording and tense
   * Sort the entries lexicographically
   * Extract summary
* [ ] Verify the version changes in:
   * [ ] `CHANGELOG`
   * [ ] `CMakeLists.txt`
   * [ ] `tools/releasing/packaging/debian/changelog`
* [ ] Commit the version bump
* [ ] Push the release branch to the precice repository
* Prepare independent releases
   * [ ] [Python bindings](https://github.com/precice/python-bindings/blob/develop/docs/ReleaseGuide.md)
   * [ ] (if necessary!) [MATLAB bindings](https://github.com/precice/matlab-bindings/blob/develop/docs/ReleaseGuide.md)
   * [ ] (if necessary!) [JULIA bindings](https://github.com/precice/PreCICE.jl)

## Step by step guide

* [ ] Open PR from `release-vX.Y.Z` to `main` (use [this template](https://github.com/precice/precice/blob/develop/tools/releasing/PULL_REQUEST_TEMPLATE/release_pull_request_template.md))
* [ ] Do regression tests using the release branch (specific revision) _list below :arrow_down:_ (all)
* [ ] Fix potential problems in develop (all)
* [ ] Rebase the release branch on develop (all)
* [ ] Commit the version bump
* [ ] Draft message to mailing list
* [ ] Write a draft "blog post" on [Discourse](https://precice.discourse.group/)
* [ ] Update documentation (all)
   * [ ] Update [XML configuration reference](https://github.com/precice/precice.github.io/blob/master/_includes/xmlreference.md)
   * [ ] Update version in [precice/precice.github.io](https://github.com/precice/precice.github.io):
      * `_config.yml`
      * `_data/sidebars/docs_sidebar.yml`
* [ ] Approve the PR with at least two reviews (all)
* [ ] Merge PR to `main` ( use `git merge --no-ff release-vX.Y.Z` )
* [ ] Tag release on `main` `vX.Y.Z` and verify by running `git describe --tags`
* [ ] Merge back to `develop` and verify by running `git describe --tags`
* [ ] Triple check that you haven't messed anything up. (You can always discard local changes)
* [ ] Push `main` and push the `vX.Y.Z` tag
* [ ] Push `develop`
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

| State | Success | Failure | Skipped |
| --- | --- | --- | --- |
| Write | `:o:` | `:x:` | `:fast_forward:` |
| Read | :o: | :x: | :fast_forward: |

| State | Tester | Test |
| --- | --- | --- |
| | | [quickstart](https://github.com/precice/tutorials/tree/master/quickstart) fluid-openfoam - solid-cpp |
| | | [perpendicular-flap](https://github.com/precice/tutorials/tree/master/perpendicular-flap) fluid-openfoam - solid-dune |
| | | [perpendicular-flap](https://github.com/precice/tutorials/tree/master/perpendicular-flap) fluid-nutils - solid-calculix |
| | | [perpendicular-flap](https://github.com/precice/tutorials/tree/master/perpendicular-flap) fluid-openfoam - solid-dealii |
| | | [perpendicular-flap](https://github.com/precice/tutorials/tree/master/perpendicular-flap) fluid-su2 - solid-fenics |
| | | [multiple-perpendicular-flaps](https://github.com/precice/tutorials/tree/master/multiple-perpendicular-flaps) fluid-openfoam - solid-(left+right)-dealii |
| | | [flow-over-heated-plate](https://github.com/precice/tutorials/tree/master/flow-over-heated-plate) fluid-openfoam - solid-openfoam serial + parallel |
| | | [flow-over-heated-plate](https://github.com/precice/tutorials/tree/master/flow-over-heated-plate) fluid-openfoam - solid-fenics serial + parallel |
| | | [flow-over-heated-plate](https://github.com/precice/tutorials/tree/master/flow-over-heated-plate) fluid-openfoam - solid-nutils |
| | | [flow-over-heated-plate-nearest-projection](https://github.com/precice/tutorials/tree/master/flow-over-heated-plate-nearest-projection) fluid-openfoam - solid-openfoam |
| | | [flow-over-heated-plate-steady-state](https://github.com/precice/tutorials/tree/master/flow-over-heated-plate-steady-state) fluid-openfoam - solid-codeaster |
| | | [heat-exchanger](https://github.com/precice/tutorials/tree/master/heat-exchanger) fluid-(inner+outer)-openfoam - solid-calculix |
| | | [partitioned-elastic-beam](https://github.com/precice/tutorials/tree/master/partitioned-elastic-beam) dirichlet-calculix - neumann-calculix |
| | | [partitioned-heat-conduction](https://github.com/precice/tutorials/tree/master/partitioned-heat-conduction) fenics - nutils |
| | | [partitioned-heat-conduction-complex](https://github.com/precice/tutorials/tree/master/partitioned-heat-conduction-complex) fenics- fenics |
| | | [partitioned-pipe](https://github.com/precice/tutorials/tree/master/partitioned-pipe) fluid1-openfoam-pimplefoam - fluid2-openfoam-sonicliquidfoam |
| | | [elastic-tube-1d](https://github.com/precice/tutorials/tree/master/elastic-tube-1d) fluid-cpp - solid-python |
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

### (only if breaking changes) Open PRs or issues `develop -> main` for all adapters

- [ ] [athlet-adapter](https://github.com/precice/athlet-adapter/compare/master...develop)
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
- [ ] [run.precice](https://github.com/precice/run.precice-frontend/compare/master...develop)
- [ ] [tutorials](https://github.com/precice/tutorials/compare/master...develop)

### Marketing

* [ ] Finalize post on [Discourse](https://precice.discourse.group/)
* [ ] Write on [Matrix](https://matrix.to/#/#precice_lobby:gitter.im?web-instance[element.io]=app.gitter.im)
* [ ] Send announcement to the [mailing list](https://mailman.informatik.uni-stuttgart.de/mailman/listinfo/precice)
* [ ] [CFD-Online](https://www.cfd-online.com/Forum/news.cgi/form/0)
* [ ] NADigest
* [ ] Post on [Twitter](https://twitter.com/preCICE_org) (additionally to the automatic)
* [ ] Post on [ResearchGate](https://www.researchgate.net/project/preCICE)
* [ ] Post in [LinkedIn Group](https://www.linkedin.com/groups/9073912/)
* [ ] Submit a short article to the [Quartl](https://www.in.tum.de/en/i05/further-activities/quartl/)


### Misc

* [ ] Update the [PR template](https://github.com/precice/precice/blob/develop/tools/releasing/PULL_REQUEST_TEMPLATE/release_pull_request_template.md)

To open a new PR with this template, use this [PR template query](https://github.com/precice/precice/compare/new?template=release_pull_request_template.md)
