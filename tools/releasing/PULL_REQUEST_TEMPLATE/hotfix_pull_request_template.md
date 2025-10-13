## How to work with this template

* [ ] assign a release manager, who takes care of the process
* [ ] assign each point below to a responsible person, before you continue. Use `@member`.

Only the release manager should update this post (even tickboxes, due to race conditions in editing). Everybody else should comment on the PR with the progress.


## Pre-PR steps

* [ ] Make sure you have the latest `main` branch locally.
* [ ] Create branch `hotfix-vX.Y.Z` from `main`. If needed, `git rebase main`.
* [ ] Commit fixes to the hotfix branch
* [ ] Check code base w.r.t code formatting (run `pre-commit run -va`)
* [ ] Update the list of operating systems for the package generation in `.github/workflows/release.yml`
* [ ] Run `tools/releasing/bumpversion.sh MAJOR.MINOR.PATCH` to bump the version
* [ ] Look over the updated `CHANGELOG.md` of the hotfix branch (all)
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
* [ ] Push the hotfix branch to the precice repository
* Prepare independent hotfixes
   * [ ] [Python bindings](https://github.com/precice/python-bindings/blob/develop/docs/hotfixGuide.md)
   * [ ] (if necessary!) [MATLAB bindings](https://github.com/precice/matlab-bindings/blob/develop/docs/hotfixGuide.md)
   * [ ] (if necessary!) [JULIA bindings](https://github.com/precice/PreCICE.jl)


## Step by step guide

* [ ] Open PR from `hotfix-vX.Y.Z` to `main` (use [this template](https://github.com/precice/precice/blob/develop/tools/releasing/PULL_REQUEST_TEMPLATE/hotfix_pull_request_template.md))
* [ ] Do regression tests using the hotfix branch (specific revision) _list below :arrow_down:_ (all)
* [ ] Fix potential problems on the hotfix branch (all)
* [ ] Reorder the commits for the version bump to be the latest. `git rebase -i main`
* [ ] Draft message to mailing list
* [ ] Write a draft "blog post" on [Discourse](https://precice.discourse.group/)
* [ ] Approve the PR with at least two reviews (all)
* [ ] Merge PR to `main` ( use `git merge --no-ff hotfix-vX.Y.Z` )
* [ ] Tag hotfix on `main` `vX.Y.Z` and verify by running `git describe --tags`
* [ ] Merge `main` back to `develop` and verify by running `git describe --tags`
* [ ] Triple check that you haven't messed anything up. (You can always discard local changes)
* [ ] Push `main` and push the `vX.Y.Z` tag
* [ ] Push `develop`
* [ ] Wait for the release pipeline
  * [ ] [To create a new draft hotfix on GitHub](https://github.com/precice/precice/releases)
  * [ ] To automatically generate packages for the latest Debian and the two latest Ubuntu LTS versions.
* [ ] Write hotfix text
* [ ] Publish the GitHub hotfix


## Regression Tests

Use the following branches:
* precice `hotfix-vX.Y.Z`
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
- [ ] [MATLAB bindings](https://github.com/precice/matlab-bindings/blob/develop/docs/hotfixGuide.md)
- [ ] [python bindings](https://github.com/precice/python-bindings/blob/develop/docs/hotfixGuide.md)
- [ ] [Julia bindings](https://github.com/precice/PreCICE.jl)
- [ ] [Rust bindings](https://github.com/precice/rust-bindings)

### Marketing

* [ ] Finalize post on [Discourse](https://precice.discourse.group/)
* [ ] Write on [Matrix](https://matrix.to/#/#precice_lobby:gitter.im?web-instance[element.io]=app.gitter.im)
* [ ] CFD Online:
     * [ ] [News](https://www.cfd-online.com/Forum/news.cgi/form/0)
     * [ ] Forum (choose a central topic)
* [ ] NADigest
* [ ] [NAFEMS](https://www.nafems.org/mynafems/submitnews/) (needs account, appears in [Upcoming Industry Events](https://www.nafems.org/events/industry-events/))
* [ ] Post on [Twitter](https://twitter.com/preCICE_org) (additionally to the automatic)
* [ ] Post on [Mastodon](https://fosstodon.org/@precice)
* [ ] Post on [ResearchGate](https://www.researchgate.net/project/preCICE)
* [ ] Post in LinkedIn. Relevant places:
    * [ ] [preCICE group](https://www.linkedin.com/groups/9073912/)
    * [ ] [OpenFOAM group](https://www.linkedin.com/groups/1920608/)
    * [ ] [HPC group](https://www.linkedin.com/groups/87791/)
    * [ ] Personal
* [ ] reddit [r/cfd](https://www.reddit.com/r/CFD/)
* [ ] Submit a short article to the [Quartl](https://www.in.tum.de/en/i05/further-activities/quartl/)

### Misc

* [ ] Update the [PR template](https://github.com/precice/precice/blob/develop/tools/releasing/PULL_REQUEST_TEMPLATE/hotfix_pull_request_template.md)

To open a new PR with this template, use this [PR template query](https://github.com/precice/precice/compare/new?template=hotfix_pull_request_template.md)
