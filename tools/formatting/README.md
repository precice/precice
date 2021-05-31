# Formatting Tools

## check-format

This bash-script checks the format of every c(pp) & h(pp) file in the current and parent directories against the clang-format style defined in a parent `.clang-format` file.
It returns 0 if everything is formatted correctly.
Otherwise, it displays the list of files that do not match the format and returns 1.

## format-all

This bash-script applies the format of a parent `.clang-format` to every c(pp) & h(pp) file in the current and child directories. To format the complete codebase, run the script from the root directory.

This script returns 0 on success.

## format-all-dockerized

Runs `format-all` on the entire repository using the dockerimage `precice/ci-formatting:latest`.
The image contains all required tools at their required versions, so only a working docker installation is required.
