## Step by step guide
* [ ] Merge master to develop (No commits after the release on master) @MakisH 
* [ ] Create branch `release-N` from develop. If needed, `git rebase develop`. @MakisH 
* [ ] Open PR from `release-N` to `master` (use [this template](https://github.com/precice/precice/blob/add_PR_template/.github/PULL_REQUEST_TEMPLATE/release_pull_request_template.md)) @MakisH 
* [ ] Look over [`CHANGELOG.md`](https://github.com/precice/precice/blob/develop/CHANGELOG.md) and add things if necessary and extract summary (all) @fsimonis 
* [ ] Do regression tests using the release branch (specific revision) _list below :arrow_down:_ (all)
* [ ] Bump version in CMakeLists.txt, SConstruct, and CHANGELOG (and Python bindings). @MakisH 
* [ ] Draft message to mailing list @uekerman @MakisH 
* [ ] Update documentation (all) @MakisH 
* [ ] Fix potential problems in develop (all)
* [ ] Merge develop to master @MakisH 
* [ ] Tag release (on GitHub) and merge back to develop @MakisH 

## Regression Tests

Run all these tests manually on your system. If you succeed, please write a comment with the revisions of the components that you used below. Example: https://github.com/precice/precice/pull/507#issuecomment-530432289

* [ ] SU2/CalculiX @MakisH 
* [ ] OpenFOAM / OpenFOAM @MakisH 
* [ ] OpenFOAM / OpenFOAM - NP mapping in OpenFOAM @MakisH
* [ ] OpenFOAM / OpenFOAM fluid coupling module @MakisH 
* [ ] OpenFOAM / CalculiX FSI (flap perp) @MakisH
* [ ] OpenFOAM / CalculiX FSI - NP mapping in CalculiX (3DTube) @KyleDavisSA 
* [ ] OpenFOAM / CalculiX / OpenFOAM CHT (heat exchanger) @MakisH
* [ ] OpenFOAM / deal.II @MakisH 
* [ ] OpenFOAM / FEniCS [flap_perp](https://github.com/precice/tutorials/tree/master/FSI/flap_perp/OpenFOAM-FEniCS) @BenjaminRueth 
* [ ] OpenFOAM / FEniCS [flow-over-plate](https://github.com/precice/tutorials/tree/master/CHT/flow-over-plate/buoyantPimpleFoam-fenics) @BenjaminRueth 
* [ ] OpenFOAM / Nutils [flow-over-plate](https://github.com/precice/tutorials/tree/master/CHT/flow-over-plate/buoyantPimpleFoam-nutils) @BenjaminRueth 
* [ ] FEniCS / FEniCS [partitioned-heat](https://github.com/precice/tutorials/tree/master/HT/partitioned-heat/fenics-fenics) @BenjaminRueth 
* [ ] ExaFSA: Ateles / FASTEST @atotoun 
* [ ] Alya @uekerman 
* [ ] 1D-ElasticTube (C++) @BenjaminRueth 
* [ ] 1D-ElasticTube (Python) @BenjaminRueth @KyleDavisSA 
* [ ] SuperMUC @fsimonis 
* [ ] C++ Dummy @fsimonis 
* [ ] C Dummy @fsimonis 
* [ ] Fortran Dummy @fsimonis 
* [ ] Fortran 2003 Dummy @fsimonis 
* [ ] Python Dummy @fsimonis 


## Post-release
* [ ] Generate packages @shkodm @fsimonis 
   * [ ] Ubuntu 18.04
   * [ ] Ubuntu 19.04
   * [ ] Arch Linux AUR Package
* [ ] Update Spack recipe @fsimonis 
* [ ] Send email and do marketing @uekerman 
   * [ ] Mention [xSDK](https://github.com/xsdk-project/xsdk-policy-compatibility/blob/master/precice-policy-compatibility.md)
* [ ] Tweet @MakisH 
* [ ] Format the code base @shkodm @fsimonis 
* [ ] Update the [PR template](https://github.com/precice/precice/blob/add_PR_template/.github/PULL_REQUEST_TEMPLATE/release_pull_request_template.md)?
