# Formatting Tools

## Moving to pre-commit

With the introduction of pre-commit, these tools became mostly obsolete.
To use, [install pre-commit](https://pre-commit.com/#install) and run `pre-commit install` at the root of the project.
You can now force the formatting on all files with `pre-commit run -a`.

## check-format

This bash-script checks the format of every c(pp) & h(pp) file in the current and parent directories against the clang-format style defined in a parent `.clang-format` file.
It returns 0 if everything is formatted correctly.
Otherwise, it displays the list of files that do not match the format and returns 1.

## cofig-formatter

This Python script formats preCICE configuration files (XML). Use as `python3 config-formatter -i precice-config.xml` to format `precice-config.xml` inline.

## format-all

This bash-script applies the format of a parent `.clang-format` to every c(pp) & h(pp) file in the current and child directories. To format the complete codebase, run the script from the root directory.

This script returns 0 on success.

## format-all-dockerized

Runs `format-all` on the entire repository using the dockerimage `precice/ci-formatting:latest`.
The image contains all required tools at their required versions, so only a working docker installation is required.
