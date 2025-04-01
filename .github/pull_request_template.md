## Main changes of this PR


## Motivation and additional information

<!--
Short rational why preCICE needs this change. If this is already described in an issue a link to that issue (closes #123) is sufficient.
-->

## Author's checklist

* [ ] I used the [`pre-commit` hook](https://precice.org/dev-docs-dev-tooling.html#setting-up-pre-commit) to prevent dirty commits and used `pre-commit run --all` to format old commits.
* [ ] I added a changelog file with `make changelog` if there are user-observable changes since the last release.
* [ ] I added a test to cover the proposed changes in our test suite.
* [ ] For breaking changes: I documented the changes in the appropriate [porting guide](https://precice.org/couple-your-code-porting-overview.html).
* [ ] I stuck to C++17 features.
* [ ] I stuck to CMake version 3.22.1.
* [ ] I squashed / am about to squash all commits that should be seen as one.
* [ ] I ran the systemtests by adding the label `trigger-system-tests` (may be skipped if minor change)

## Reviewers' checklist

<!-- Tag people next to each point and add points for specific questions -->

* [ ] Does the changelog entry make sense? Is it formatted correctly?
* [ ] Do you understand the code changes?

<!-- add more questions/tasks if necessary -->
