## How to work with this template

* [ ] assign a release manager, who takes care of the process
* [ ] assign each point below to a responsible person, before you continue. Use `@member`.

Only the release manager should update this post (even tickboxes, due to race conditions in editing). Everybody else should comment on the PR with the progress.


## Pre-PR steps

* [ ] Make sure you have the latest `master` branch locally.
* [ ] Create branch `hotfix-vX.Y.Z` from master. If needed, `git rebase master`.
* [ ] Commit fixes to the hotfix branch
* [ ] Check code base w.r.t code formatting (run [`precice/tools/formatting/check-format`](https://github.com/precice/precice/blob/develop/tools/formatting/check-format)) and reformat if required (run [`precice/tools/formatting/format-all`](https://github.com/precice/precice/blob/develop/tools/formatting/format-all))
* [ ] Run `tools/releasing/bumpversion.sh MAJOR.MINOR.PATCH` to bump the version
* [ ] Look over [`CHANGELOG.md`](https://github.com/precice/precice/blob/develop/CHANGELOG.md) (all)
   * Check for merged lines
   * Add things, if necessary
   * Fix wording and tense
   * Sort the entries lexicographically
   * Extract summary
* [ ] Verify the version changes in:
   * [ ] [CHANGELOG](https://github.com/precice/precice/blob/develop/CHANGELOG.md)
   * [ ] [CMakeLists.txt](https://github.com/precice/precice/blob/develop/CMakeLists.txt)
   * [ ] [debian changelog](https://github.com/precice/precice/blob/develop/tools/releasing/packaging/debian/changelog)
* [ ] Commit the version bump
* [ ] Push the hotfix branch to the precice repository
* Prepare independent hotfixs
   * [ ] [Python bindings](https://github.com/precice/python-bindings/blob/develop/docs/hotfixGuide.md)
   * [ ] (if necessary!) [MATLAB bindings](https://github.com/precice/matlab-bindings/blob/develop/docs/hotfixGuide.md)
   * [ ] (if necessary!) [JULIA bindings](https://github.com/precice/PreCICE.jl)


## Step by step guide

* [ ] Open PR from `hotfix-vX.Y.Z` to `master` (use [this template](https://github.com/precice/precice/blob/add_PR_template/.github/PULL_REQUEST_TEMPLATE/release_pull_request_template.md))
* [ ] Do regression tests using the hotfix branch (specific revision) _list below :arrow_down:_ (all)
* [ ] Fix potential problems on the hotfix branch (all)
* [ ] Reorder the commits for the version bump to be the latest. `git rebase -i master`
* [ ] Draft message to mailing list
* [ ] Write a draft "blog post" on [Discourse](https://precice.discourse.group/)
* [ ] Approve the PR with at least two reviews (all)
* [ ] Merge PR to master ( use `git merge --no-ff hotfix-vX.Y.Z` )
* [ ] Tag hotfix on master `vX.Y.Z` and verify by running `git describe --tags`
* [ ] Merge master back to develop and verify by running `git describe --tags`
* [ ] Triple check that you haven't messed anything up. (You can always discard local changes)
* [ ] Push master and push the `vX.Y.Z` tag
* [ ] Push develop
* [ ] Wait for the release pipeline
  * [ ] [To create a new draft hotfix on GitHub](https://github.com/precice/precice/releases)
  * [ ] To automatically generate packages for latest Debian and two latest Ubuntu LTS versions.
* [ ] Write hotfix text
* [ ] Publish the GitHub hotfix


## Regression Tests

Use the following branches:
* precice `hotfix-vX.Y.Z`
* pyprecice `python-bindings-vX.Y.Z.1`
* matlab-bindings `matlab-bindings-vX.Y.Z.1`
* rest `master`

Run all these tests manually on your system. If you succeed, please write a comment with the revisions of the components that you used below. Example: https://github.com/precice/precice/pull/507#issuecomment-530432289 and update the table.

| State | Success | Failure | Skipped |
| --- | --- | --- | --- |
| Write | `:o:` | `:x:` | `:fast_forward:` |
| Read | :o: | :x: | :fast_forward: |

| State | Tester | Test |
| --- | --- | --- |
| | | SU2 / CalculiX [flap_perp](https://github.com/precice/tutorials/tree/develop/FSI/flap_perp/SU2-CalculiX) |
| | | SU2 / FEniCS [flap_perp_2D](https://github.com/precice/tutorials/tree/develop/FSI/flap_perp_2D/SU2-FEniCS) |
| | | OpenFOAM / OpenFOAM [flow_over_plate](https://github.com/precice/openfoam-adapter/tree/master/tutorials/CHT/flow-over-plate) serial + parallel |
| | | OpenFOAM / OpenFOAM - NP mapping in OpenFOAM [flow_over_plate](https://github.com/precice/openfoam-adapter/tree/master/tutorials/CHT/flow-over-plate) |
| | | OpenFOAM / CalculiX FSI [flap perp](https://github.com/precice/tutorials/tree/develop/FSI/flap_perp/OpenFOAM-CalculiX) |
| | | OpenFOAM / CalculiX FSI - NP mapping in CalculiX [3D_Tube](https://github.com/precice/tutorials/tree/develop/FSI/3D_Tube/OpenFOAM-CalculiX) |
| | | OpenFOAM / CalculiX / OpenFOAM CHT [heat_exchanger](https://github.com/precice/tutorials/tree/develop/CHT/heat_exchanger/buoyantSimpleFoam-CalculiX) |
| | | OpenFOAM / deal.II [flap_perp_2D](https://github.com/precice/tutorials/tree/develop/FSI/flap_perp_2D/OpenFOAM-deal.II) (linear + non-linear, serial + parallel) |
| | | OpenFOAM / deal.II [flap_perp](https://github.com/precice/tutorials/tree/master/FSI/flap_perp/OpenFOAM-deal.II) |
| | | OpenFOAM / deal.II [FSI cylinderFlap_2D](https://github.com/precice/tutorials/tree/develop/FSI/cylinderFlap_2D/OpenFOAM-deal.II) |
| | | OpenFOAM / deal.II [FSI cylinderFlap](https://github.com/precice/tutorials/tree/develop/FSI/cylinderFlap_2D/OpenFOAM-deal.II) |
| | | OpenFOAM / FEniCS [flap_perp](https://github.com/precice/tutorials/tree/master/FSI/flap_perp/OpenFOAM-FEniCS) |
| | | OpenFOAM / FEniCS [flow-over-plate](https://github.com/precice/tutorials/tree/master/CHT/flow-over-plate/buoyantPimpleFoam-fenics) |
| | | OpenFOAM / FEniCS [cylinderFlap, only run first few minutes](https://github.com/precice/tutorials/tree/develop/FSI/cylinderFlap/OpenFOAM-FEniCS) |
| | | OpenFOAM / Nutils [flow-over-plate](https://github.com/precice/tutorials/tree/master/CHT/flow-over-plate/buoyantPimpleFoam-nutils) |
| | | FEniCS / FEniCS [partitioned-heat](https://github.com/precice/tutorials/tree/master/HT/partitioned-heat/fenics-fenics) |
| | | MATLAB / MATLAB [ODEs](https://github.com/precice/matlab-bindings/tree/develop/tutorial) |
| | | 1D-ElasticTube [C++](https://github.com/precice/elastictube1d/tree/develop/cxx) |
| | | 1D-ElasticTube [Python](https://github.com/precice/elastictube1d/tree/develop/python) |
| | | Solverdummy [Fortran module](https://github.com/precice/fortran-module/tree/develop/examples/solverdummy) |
| | | Solverdummy [Python](https://github.com/precice/python-bindings/tree/develop/solverdummy) |
| | | Solverdummy [MATLAB](https://github.com/precice/matlab-bindings/tree/develop/solverdummy) |
| | | Alya |
| | | ExaFSA: Ateles / FASTEST |
| | | SuperMUC |


## Post-release

* [ ] Flag [Arch Linux AUR package](https://aur.archlinux.org/packages/precice) and dependants as out-of-date.
* [ ] Update Spack recipe
* [ ] Update pyprecice Spack
* [ ] Update Website:
    * [ ] Bump version in [`_config.yml`](https://github.com/precice/precice.github.io/blob/master/_config.yml)
    * [ ] Update the [XML reference](https://github.com/precice/precice.github.io/blob/master/_includes/xmlreference.md) using `binprecice md`
    * [ ] Look over the [Roadmap](https://www.precice.org/fundamentals-roadmap.html) and update entries.

### release new version for bindings (to ensure compatibility with newest preCICE version)

- [ ] [Fortran module](https://github.com/precice/fortran-module/compare/master...develop)
- [ ] [MATLAB bindings](https://github.com/precice/matlab-bindings/blob/develop/docs/hotfixGuide.md)
- [ ] [python bindings](https://github.com/precice/python-bindings/blob/develop/docs/hotfixGuide.md)
- [ ] [Julia bindings](https://github.com/precice/PreCICE.jl)

### Marketing

* [ ] Finalize post on [Discourse](https://precice.discourse.group/)
* [ ] Write on [Gitter](https://gitter.im/precice/Lobby)
* [ ] Send announcement to the [mailing list](https://mailman.informatik.uni-stuttgart.de/mailman/listinfo/precice)
* [ ] [CFD-Online](https://www.cfd-online.com/Forum/news.cgi/form/0)
* [ ] NADigest
* [ ] Post on [Twitter](https://twitter.com/preCICE_org) (additionally to the automatic)
* [ ] Post on [ResearchGate](https://www.researchgate.net/project/preCICE)
* [ ] Post in [LinkedIn Group](https://www.linkedin.com/groups/9073912/)
* [ ] Submit a short article to the [Quartl](https://www.in.tum.de/en/i05/further-activities/quartl/)


### Misc

* [ ] Update the [PR template](https://github.com/precice/precice/blob/add_PR_template/.github/PULL_REQUEST_TEMPLATE/hotfix_pull_request_template.md)

To open a new PR with this template, use this [PR template query](https://github.com/precice/precice/compare/new?template=hotfix_pull_request_template.md)
