# How to collaborate within the preCICE project?

**Welcome to preCICE as a collaborator, we value your contribution!**

In order to get your contributions into the code base as smoothly as possible, please follow these contribution guidelines.

* Make sure you have a GitHub account.
* Submit a ticket for your issue, assuming one does not already exist.
  * Clearly describe the issue including steps to reproduce when it is a bug.
  * Make sure you fill in the earliest version that you know has the issue.
* Fork the repository on GitHub.
* Create a topic branch, based on ```develop```, from where you want to base your work.
* Make commits of logical units.
* Check for unnecessary whitespace with `git diff --check` before committing.
* Write tests.
* Run _all_ the tests to assure nothing else was accidentally broken.
* Follow our [style guide][style].
* Write a [good commit message][commit].
* Submit a pull request to the repository in the preCICE organization.

## Taking code from other projects
We believe in the power of Open Source or Free Software Software to share and reuse code from other projects. However, Free Software is not public domain, and not every code could be reused in every other project.

Before taking non-trivial amount of code from other projects, check back with Benjamin (@uekerman) or Florian (@floli), so we can ensure the compatibility of licences. Same holds true for additional dependencies, libraries etc.

## Ressources
* [Collaboration workflow with pull requests and issues][workflow]
* [Creating a pull request][pullrequest]

[commit]: http://tbaggery.com/2008/04/19/a-note-about-git-commit-messages.html
[pullrequest]: https://help.github.com/articles/creating-a-pull-request
[style]: https://ipvs.informatik.uni-stuttgart.de/sgs/precice/docs/develop/conventions.html
[workflow]: https://help.github.com/categories/collaborating-with-issues-and-pull-requests
