## How to work with this template

* [ ] assign a release manager, who takes care of the process
* [ ] assign each point below to a responsible person, before you continue. Use `@member`.

Only the release manager should update this post (even tickboxes, due to race conditions in editing). Everybody else should comment on the PR with the progress.

## Pre-PR steps

* [ ] Make sure you have the latest `main` branch locally.
* [ ] Create branch `hotfix-vX.Y.Z` from `main`. If needed, `git rebase main`.
* [ ] Commit fixes to the hotfix branch
* [ ] Check code base w.r.t code formatting (run `pre-commit run -va`)
* [ ] Update the list of operating systems for the package generation in `.github/workflows/release.yml`
* [ ] Bump the version: `tools/releasing/bumpversion.sh MAJOR.MINOR.PATCH`
* [ ] Look over the updated `CHANGELOG.md` of the hotfix branch (all)
   * Check for merged lines
   * Add things, if necessary
   * Fix wording and tense
   * Sort the entries lexicographically
   * Extract summary
* [ ] Verify the version changes in:
   * [ ] `CHANGELOG`
   * [ ] `CMakeLists.txt`
   * [ ] `tools/releasing/packaging/debian/changelog`
* [ ] Commit the version bump: `git commit -m "Bump version to X.Y.Z"`
* [ ] Push the hotfix branch to the precice repository: `git push -u upstream hotfix-vX.Y.Z`

## Step by step guide

* [ ] Open PR from `hotfix-vX.Y.Z` to `main` (use [this template](https://github.com/precice/precice/blob/develop/tools/releasing/PULL_REQUEST_TEMPLATE/hotfix_pull_request_template.md))
* [ ] Trigger the system tests using the `trigger-system-tests` label ([`release_test` suite](https://github.com/precice/tutorials/blob/develop/tools/tests/tests.yaml)). After any force-push, remove and add the label again.
* [ ] Fix potential problems on the hotfix branch (all)
* [ ] Reorder the commits for the version bump to be the latest. `git rebase -i main`
* [ ] Draft release notes
* [ ] Write a draft "blog post" on [Discourse](https://precice.discourse.group/)
* [ ] Approve the PR with at least two reviews (all)
* [ ] Merge PR to `main` ( use `git merge --no-ff hotfix-vX.Y.Z` )
* [ ] Tag hotfix on `main` `vX.Y.Z` and verify by running `git describe --tags`
* [ ] Merge `main` back to `develop` and verify by running `git describe --tags`
* [ ] Triple check that you haven't messed anything up. You can always discard local changes using `git reset --hard upstream BRANCH` or by cloning the precice repository again and start from scratch.
* [ ] Push `main` and the `vX.Y.Z` tag: `git push upstream main`, `git push upstream v3.3.1`
* [ ] Push `develop`: `git push upstream develop`
* [ ] Wait for the release pipeline
  * [ ] [To create a new draft hotfix on GitHub](https://github.com/precice/precice/releases)
  * [ ] To automatically generate packages for the latest Debian and the two latest Ubuntu LTS versions.
* [ ] Write hotfix text
* [ ] Publish the GitHub hotfix

## Post-release

* [ ] Update version specific documentation
* [ ] Flag [Arch Linux AUR package](https://aur.archlinux.org/packages/precice) and dependants as out-of-date.
* [ ] Update Spack recipe
* [ ] Update Website:
    * [ ] Bump version in [`_config.yml`](https://github.com/precice/precice.github.io/blob/master/_config.yml)
    * [ ] Update the [XML reference](https://github.com/precice/precice.github.io/blob/master/_includes/xmlreference.md) using `binprecice md`
    * [ ] Look over the [Roadmap](https://www.precice.org/fundamentals-roadmap.html) and update entries.

### Marketing

* [ ] Finalize post on [Discourse](https://precice.discourse.group/)
* [ ] Write on [Matrix](https://matrix.to/#/#precice_lobby:gitter.im?web-instance[element.io]=app.gitter.im)
* [ ] CFD Online:
     * [ ] [News](https://www.cfd-online.com/Forum/news.cgi/form/0)
     * [ ] [Forum](https://www.cfd-online.com/Forums/main/261393-precice-releases.html)
* [ ] NADigest
* [ ] [NAFEMS](https://www.nafems.org/mynafems/submitnews/) (needs account, appears in [Upcoming Industry Events](https://www.nafems.org/events/industry-events/))
* [ ] Post on [Bluesky](https://bsky.app/profile/precice.org)
* [ ] Post on [Mastodon](https://fosstodon.org/@precice)
* [ ] Post in LinkedIn. Relevant places:
    * [ ] [preCICE company](https://www.linkedin.com/company/precice)
    * [ ] [OpenFOAM group](https://www.linkedin.com/groups/1920608/)
    * [ ] [HPC group](https://www.linkedin.com/groups/87791/)
    * [ ] Personal

See the [advertising a workshop](https://precice.org/precice-workshop-organizing.html#advertising) section for additional ideas.

### Misc

* [ ] Update the [PR template](https://github.com/precice/precice/blob/develop/tools/releasing/PULL_REQUEST_TEMPLATE/hotfix_pull_request_template.md)

To open a new PR with this template, use this [PR template query](https://github.com/precice/precice/compare/new?template=hotfix_pull_request_template.md)
