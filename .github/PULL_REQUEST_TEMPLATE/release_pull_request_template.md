## How to work with this template

* [ ] assign a release manager, who takes care of the process
* [ ] assign each point below to a responsible person, before you continue. Use `@member`.

Only the release manager should update this post (even tickboxes, due to race conditions in editing). Everybody else should comment on the PR with the progress.

## Step by step guide
* [ ] Merge master to develop (No commits after the release on master)
* [ ] Create branch `release-N` from develop. If needed, `git rebase develop`.
* [ ] Open PR from `release-N` to `master` (use [this template](https://github.com/precice/precice/blob/add_PR_template/.github/PULL_REQUEST_TEMPLATE/release_pull_request_template.md))
* [ ] Look over [`CHANGELOG.md`](https://github.com/precice/precice/blob/develop/CHANGELOG.md) (all)
   * Add things, if necessary
   * Extract summary
* [ ] Look over the Roadmap and update entries.
* [ ] Do regression tests using the release branch (specific revision) _list below :arrow_down:_ (all)
* [ ] Bump version in:
   * [ ] [CHANGELOG](https://github.com/precice/precice/blob/develop/CHANGELOG.md)
   * [ ] [CMakeLists.txt](https://github.com/precice/precice/blob/develop/CMakeLists.txt)
   * [ ] (do independent release, if necessary!) [Python bindings](https://github.com/precice/python-bindings)
   * [ ] (do independent release, if necessary!) [MATLAB bindings](https://github.com/precice/matlab-bindings)
* [ ] Draft message to mailing list
* [ ] Update documentation (all)
  * [ ] Update markdown configuration reference in wiki
* [ ] Fix potential problems in develop (all)
* [ ] Approve the PR with at least two reviews (all)
* [ ] Merge PR to master 
* [ ] Tag release (on GitHub, make sure to select the release branch as target) and merge back to develop

## Regression Tests

Run all these tests manually on your system. If you succeed, please write a comment with the revisions of the components that you used below. Example: https://github.com/precice/precice/pull/507#issuecomment-530432289

* [ ] SU2 / CalculiX [flap_perp](https://github.com/precice/tutorials/tree/develop/FSI/flap_perp/SU2-CalculiX)
* [ ] OpenFOAM / OpenFOAM [flow_over_plate](https://github.com/precice/openfoam-adapter/tree/master/tutorials/CHT/flow-over-plate)
* [ ] OpenFOAM / OpenFOAM - NP mapping in OpenFOAM [flow_over_plate](https://github.com/precice/openfoam-adapter/tree/master/tutorials/CHT/flow-over-plate)
* [ ] OpenFOAM / OpenFOAM fluid coupling module **???**
* [ ] OpenFOAM / CalculiX FSI [flap perp](https://github.com/precice/tutorials/tree/develop/FSI/flap_perp/OpenFOAM-CalculiX)
* [ ] OpenFOAM / CalculiX FSI - NP mapping in CalculiX [3D_Tube](https://github.com/precice/tutorials/tree/develop/FSI/3D_Tube/OpenFOAM-CalculiX)
* [ ] OpenFOAM / CalculiX / OpenFOAM CHT [heat_exchanger](https://github.com/precice/tutorials/tree/develop/CHT/heat_exchanger/buoyantSimpleFoam-CalculiX)
* [ ] OpenFOAM / deal.II [flap_perp_2D](https://github.com/precice/tutorials/tree/develop/FSI/flap_perp_2D/OpenFOAM-deal.II)
* [ ] OpenFOAM / FEniCS [flap_perp](https://github.com/precice/tutorials/tree/master/FSI/flap_perp/OpenFOAM-FEniCS)
* [ ] OpenFOAM / FEniCS [flow-over-plate](https://github.com/precice/tutorials/tree/master/CHT/flow-over-plate/buoyantPimpleFoam-fenics)
* [ ] OpenFOAM / Nutils [flow-over-plate](https://github.com/precice/tutorials/tree/master/CHT/flow-over-plate/buoyantPimpleFoam-nutils)
* [ ] FEniCS / FEniCS [partitioned-heat](https://github.com/precice/tutorials/tree/master/HT/partitioned-heat/fenics-fenics)
* [ ] ExaFSA: Ateles / FASTEST
* [ ] Alya
* [ ] 1D-ElasticTube [C++](https://github.com/precice/elastictube1d/tree/develop/cxx)
* [ ] 1D-ElasticTube [Python](https://github.com/precice/elastictube1d/tree/develop/python)
* [ ] SuperMUC
* [ ] Solverdummy [C++](https://github.com/precice/precice/tree/develop/tools/solverdummies/cpp)
* [ ] Solverdummy [C](https://github.com/precice/precice/tree/develop/tools/solverdummies/c)
* [ ] Solverdummy [Fortran](https://github.com/precice/precice/tree/develop/tools/solverdummies/fortran)
* [ ] Solverdummy [Fortran 2003](https://github.com/precice/precice/tree/develop/tools/solverdummies/f2003)
* [ ] Solverdummy [Python](https://github.com/precice/python-bindings/tree/develop/solverdummy)
* [ ] Solverdummy [MATLAB](https://github.com/precice/matlab-bindings/tree/develop/solverdummy)


## Post-release
* [ ] Generate packages
   * [ ] Latest Ubuntu LTS
   * [ ] Latest Ubuntu
   * [ ] Arch Linux AUR Package
* [ ] Update Spack recipe
* [ ] Send email and do marketing
* [ ] Tweet
* [ ] Format the code base
* [ ] Update the [PR template](https://github.com/precice/precice/blob/add_PR_template/.github/PULL_REQUEST_TEMPLATE/release_pull_request_template.md)

To open a new PR with this template, use this [PR template query](https://github.com/precice/precice/compare/new?template=release_pull_request_template.md)
