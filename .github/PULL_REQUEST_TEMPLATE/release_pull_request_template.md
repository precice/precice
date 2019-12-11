## How to work with this template

* [ ] assign a release manager, who takes care of the process
* [ ] assign each point below to a responsible person, before you continue. Use `@member`.

## Step by step guide
* [ ] Merge master to develop (No commits after the release on master)
* [ ] Create branch `release-N` from develop. If needed, `git rebase develop`.
* [ ] Open PR from `release-N` to `master` (use [this template](https://github.com/precice/precice/blob/add_PR_template/.github/PULL_REQUEST_TEMPLATE/release_pull_request_template.md))
* [ ] Look over [`CHANGELOG.md`](https://github.com/precice/precice/blob/develop/CHANGELOG.md) and add things if necessary and extract summary (all)
* [ ] Look over the Roadmap and update entries.
* [ ] Do regression tests using the release branch (specific revision) _list below :arrow_down:_ (all)
* [ ] Bump version in:
   * [ ] [CHANGELOG](https://github.com/precice/precice/blob/develop/CHANGELOG.md)
   * [ ] [CMakeLists.txt](https://github.com/precice/precice/blob/develop/CMakeLists.txt)
   * [ ] **???** [Python bindings](https://github.com/precice/python-bindings) 
* [ ] Draft message to mailing list
* [ ] Update documentation (all)
* [ ] Fix potential problems in develop (all)
* [ ] Approve the PR with at least two reviews (all)
* [ ] Merge PR to master 
* [ ] Tag release (on GitHub, make sure to select the release branch as target) and merge back to develop

## Regression Tests

Run all these tests manually on your system. If you succeed, please write a comment with the revisions of the components that you used below. Example: https://github.com/precice/precice/pull/507#issuecomment-530432289

* [ ] SU2/CalculiX
* [ ] OpenFOAM / OpenFOAM
* [ ] OpenFOAM / OpenFOAM - NP mapping in OpenFOAM
* [ ] OpenFOAM / OpenFOAM fluid coupling module 
* [ ] OpenFOAM / CalculiX FSI (flap perp)
* [ ] OpenFOAM / CalculiX FSI - NP mapping in CalculiX (3DTube)
* [ ] OpenFOAM / CalculiX / OpenFOAM CHT (heat exchanger)
* [ ] OpenFOAM / deal.II
* [ ] OpenFOAM / FEniCS [flap_perp](https://github.com/precice/tutorials/tree/master/FSI/flap_perp/OpenFOAM-FEniCS)
* [ ] OpenFOAM / FEniCS [flow-over-plate](https://github.com/precice/tutorials/tree/master/CHT/flow-over-plate/buoyantPimpleFoam-fenics)
* [ ] OpenFOAM / Nutils [flow-over-plate](https://github.com/precice/tutorials/tree/master/CHT/flow-over-plate/buoyantPimpleFoam-nutils)
* [ ] FEniCS / FEniCS [partitioned-heat](https://github.com/precice/tutorials/tree/master/HT/partitioned-heat/fenics-fenics)
* [ ] ExaFSA: Ateles / FASTEST
* [ ] Alya
* [ ] 1D-ElasticTube (C++)
* [ ] 1D-ElasticTube (Python)
* [ ] SuperMUC
* [ ] C++ Dummy
* [ ] C Dummy
* [ ] Fortran Dummy
* [ ] Fortran 2003 Dummy
* [ ] Python Dummy


## Post-release
* [ ] Generate packages
   * [ ] Ubuntu 18.04
   * [ ] Ubuntu 19.10
   * [ ] Arch Linux AUR Package
* [ ] Update Spack recipe
* [ ] Send email and do marketing
* [ ] Tweet
* [ ] Format the code base
* [ ] Update the [PR template](https://github.com/precice/precice/blob/add_PR_template/.github/PULL_REQUEST_TEMPLATE/release_pull_request_template.md)

To open a new PR with this template, use this [PR template query](https://github.com/precice/precice/compare/new?template=release_pull_request_template.md)
