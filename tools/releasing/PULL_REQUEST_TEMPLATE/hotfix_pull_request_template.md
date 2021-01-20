## How to work with this template

* [ ] assign a release manager, who takes care of the process
* [ ] assign each point below to a responsible person, before you continue. Use `@member`.

Only the release manager should update this post (even tickboxes, due to race conditions in editing). Everybody else should comment on the PR with the progress.

## Step by step guide
* [ ] Create branch `hotfix-N` from master.
* [ ] Apply your fixes on this branch
* [ ] Look over [`CHANGELOG.md`](https://github.com/precice/precice/blob/develop/CHANGELOG.md) (all)
   * Add things, if necessary
   * Extract summary
   * Fix wording and tense
   * Sort the entries lexicographically
* [ ] Check code base w.r.t code formatting (run [`precice/tools/formatting/check-format`](https://github.com/precice/precice/blob/develop/tools/formatting/check-format)) and reformat if required (run [`precice/tools/formatting/format-all`](https://github.com/precice/precice/blob/develop/tools/formatting/format-all))
* [ ] Run `tools/releasing/bumpversion.sh MAJOR.MINOR.PATCH` to bump the version
* [ ] Verify the version changes in:
   * [ ] [CHANGELOG](https://github.com/precice/precice/blob/develop/CHANGELOG.md)
   * [ ] [CMakeLists.txt](https://github.com/precice/precice/blob/develop/CMakeLists.txt)
   * [ ] [debian changelog](https://github.com/precice/precice/blob/develop/tools/releasing/packaging/debian/changelog)
   * [ ] (prepare independent release) [Python bindings](https://github.com/precice/python-bindings)
   * [ ] (prepare independent release, if necessary!) [MATLAB bindings](https://github.com/precice/matlab-bindings)
* [ ] Commit the version bump
* [ ] Open PR from `hotfix-N` to `master` (use [this template](https://github.com/precice/precice/blob/add_PR_template/.github/PULL_REQUEST_TEMPLATE/release_pull_request_template.md))
* [ ] Do regression tests using the hotfix branch (specific revision) _list below :arrow_down:_ (all)
* [ ] Fix potential problems in the hotfix branch (all)
* [ ] Draft message to mailing list
* [ ] Update documentation (all)
   * [ ] Update [XML configuration reference](https://github.com/precice/precice.github.io/blob/master/_includes/xmlreference.md)
   * [ ] Update version in [precice/precice.github.io](https://github.com/precice/precice.github.io):
      * `_config.yml`
      * `_data/sidebars/docs_sidebar.yml`
* [ ] Approve the PR with at least two reviews (all)
* [ ] Merge PR to master ( use `git merge --no-ff hotfix-N` )
* [ ] Tag hotfix on master `vN` and verify by running `git describe --tags`
* [ ] Merge back to develop and verify by running `git describe --tags`
* [ ] Push master and push the `vN` tag
* [ ] [Draft a new release on GitHub](https://github.com/precice/precice/releases/new)
* [ ] Generate packages and upload to the draft release
   * [ ] Latest Ubuntu LTS
   * [ ] Latest Ubuntu
* [ ] Publish the GitHub release

## Regression Tests

If necessary, run some these tests manually on your system. Please post the result of the test as a comment below and update the table.

| State | Success | Failure |
| --- | --- | --- |
| Write | `:o:` | `:x:` |
| Read | :o: | :x: |

| State | Tester | Test |
| --- | --- | --- |
| | | SU2 / CalculiX [flap_perp](https://github.com/precice/tutorials/tree/develop/FSI/flap_perp/SU2-CalculiX) |
| | | OpenFOAM / OpenFOAM [flow_over_plate](https://github.com/precice/openfoam-adapter/tree/master/tutorials/CHT/flow-over-plate) |
| | | OpenFOAM / OpenFOAM - NP mapping in OpenFOAM [flow_over_plate](https://github.com/precice/openfoam-adapter/tree/master/tutorials/CHT/flow-over-plate) |
| | | OpenFOAM / CalculiX FSI [flap perp](https://github.com/precice/tutorials/tree/develop/FSI/flap_perp/OpenFOAM-CalculiX) |
| | | OpenFOAM / CalculiX FSI - NP mapping in CalculiX [3D_Tube](https://github.com/precice/tutorials/tree/develop/FSI/3D_Tube/OpenFOAM-CalculiX) |
| | | OpenFOAM / CalculiX / OpenFOAM CHT [heat_exchanger](https://github.com/precice/tutorials/tree/develop/CHT/heat_exchanger/buoyantSimpleFoam-CalculiX) |
| | | OpenFOAM / deal.II [flap_perp_2D](https://github.com/precice/tutorials/tree/develop/FSI/flap_perp_2D/OpenFOAM-deal.II) |
| | | OpenFOAM / FEniCS [flap_perp](https://github.com/precice/tutorials/tree/master/FSI/flap_perp/OpenFOAM-FEniCS) |
| | | OpenFOAM / FEniCS [flow-over-plate](https://github.com/precice/tutorials/tree/master/CHT/flow-over-plate/buoyantPimpleFoam-fenics) |
| | | OpenFOAM / Nutils [flow-over-plate](https://github.com/precice/tutorials/tree/master/CHT/flow-over-plate/buoyantPimpleFoam-nutils) |
| | | FEniCS / FEniCS [partitioned-heat](https://github.com/precice/tutorials/tree/master/HT/partitioned-heat/fenics-fenics) |
| | | MATLAB / MATLAB [ODEs](https://github.com/precice/matlab-bindings/tree/develop/tutorial) |
| | | ExaFSA: Ateles / FASTEST |
| | | Alya |
| | | 1D-ElasticTube [C++](https://github.com/precice/elastictube1d/tree/develop/cxx) | 
| | | 1D-ElasticTube [Python](https://github.com/precice/elastictube1d/tree/develop/python) |
| | | SuperMUC |
| | | Solverdummy [Fortran 2003](https://github.com/precice/precice/tree/develop/tools/solverdummies/f2003) | 
| | | Solverdummy [Python](https://github.com/precice/python-bindings/tree/develop/solverdummy) |
| | | Solverdummy [MATLAB](https://github.com/precice/matlab-bindings/tree/develop/solverdummy) |


## Post-release

* [ ] (independent release) [Python bindings](https://github.com/precice/python-bindings)
* [ ] (independent release) [MATLAB bindings](https://github.com/precice/matlab-bindings)
* [ ] Update Arch Linux AUR Package
* [ ] Update Spack recipe
* [ ] Send email and do marketing
* [ ] Tweet
