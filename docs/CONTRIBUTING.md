# How to collaborate within the preCICE project?

**Welcome to preCICE as a collaborator, we value your contribution!**

In order to get your contributions into the code base as smoothly as possible, please follow these contribution guidelines.

* Make sure you have a GitHub account.
* [Open an issue][newissue], assuming one does not already exist.
  * Clearly describe the issue including steps to reproduce when it is a bug.
  * Make sure you fill in the earliest version that you know has the issue.
* Fork the repository on GitHub.
* Create a feature branch, based on `develop`, from where you want to base your work. For simplicity, prefix the branch name either with `add-` or `fix-`.
* Make commits of logical units.
  Write [good commit messages][commit].
  Check for unnecessary whitespace with `git diff --check` before committing.
* Write tests to assure your feature works as expected and prevent it from getting broken in the future.
  See our example tests in `src/testing/tests/ExampleTests.cpp` and the [documentation of boost.test][boosttest] for more information.
* Run _all_ the tests to assure nothing else accidentally broke.
* Use the available [tools to format and check your contribution](https://precice.org/dev-docs-dev-tooling.html), in particular the pre-commit hook that will format the code following our [style guide][style].
* Submit a pull request to the repository in the preCICE organization and fill the provided template.
  See [Collaboration workflow with pull requests and issues][workflow] and [Creating a pull request][pullrequest] to get started.
* If applicable, add an entry for our `CHANGELOG.md` as a file `docs/changelog/123.md`, where `123` your Pull Request number.
  If you have [GitHub CLI](https://cli.github.com/) installed, you can use the script `createChangelog` in `tools/building`, or run `make changelog` to create the file.
  The content of the file is a markdown list of the entries using `*` as marker. Write 1 entry per line to simplify the merging process. Each entry should start with a verb in the past tense such as `Added`, `Fixed`, `Updated`, `Simplified`, `Improved`. This simplifies sorting the changelog and finding interesting content in it.

## Taking code from other projects
We believe in the power of Open Source or Free Software to share and reuse code from other projects. However, Free Software is not public domain, and not every code could be reused in every other project.

Please contact the maintainers before integrating non-trivial amount of code from other projects, so we can ensure the compatibility of licences. Same holds true for additional dependencies, libraries etc.

[newissue]: https://github.com/precice/precice/issues/new/choose
[boosttest]: https://www.boost.org/doc/libs/1_65_1/libs/test/doc/html/index.html
[commit]: http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html
[pullrequest]: https://help.github.com/articles/creating-a-pull-request
[style]: https://precice.org/dev-docs-dev-conventions.html
[workflow]: https://docs.github.com/en/pull-requests/collaborating-with-pull-requests
