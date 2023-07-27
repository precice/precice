## Main changes of this PR


## Motivation and additional information

<!--
Short rational why preCICE needs this change. If this is already described in an issue a link to that issue (closes #123) is sufficient.
-->

## Author's checklist

* [ ] I used the [`pre-commit` hook](https://precice.org/dev-docs-dev-tooling.html#setting-up-pre-commit) to prevent dirty commits and used `pre-commit run --all` to format old commits.
* [ ] I ran coverity on this PR and took care of defecty I may have introduced. (Use workflow dispatch [here](https://github.com/precice/precice/actions/workflows/coverity-scan.yml))
* [ ] I added a changelog file with `make changelog` if there are user-observable changes since the last release.
* [ ] I added a test to cover the proposed changes in our test suite.
* [ ] For breaking changes: I documented the changes in the appropriate [porting guide](https://precice.org/couple-your-code-porting-overview.html).
* [ ] I sticked to C++17 features.
* [ ] I sticked to CMake version 3.16.3.
* [ ] I squashed / am about to squash all commits that should be seen as one.

## Reviewers' checklist

<!-- Tag people next to each point and add points for specific questions -->

* [ ] Does the changelog entry make sense? Is it formatted correctly?
* [ ] Do you understand the code changes?

<!-- add more questions/tasks if necessary -->
