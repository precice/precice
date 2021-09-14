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
| | | Alya |
| | | ExaFSA: Ateles / FASTEST |
| | | SuperMUC |

## Post-release

* [ ] (independent release) [Python bindings](https://github.com/precice/python-bindings)
* [ ] (independent release) [MATLAB bindings](https://github.com/precice/matlab-bindings)
* [ ] Update Arch Linux AUR Package
* [ ] Update Spack recipe
* [ ] Send email and do marketing
* [ ] Tweet
