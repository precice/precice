#! /bin/bash

set -e

USAGE="
usage: $(basename "$0") MAJOR.MINOR.PATCH

Bumps the version in all required files.
"

if (( $# != 1 )); then
    echo "$USAGE"
    exit 1
fi

VERSION=$1

echo -n "Checking version ... "
if [[ $VERSION =~ ^[0-9]+\.[0-9]+\.[0-9]+$ ]]; then
    echo "ok"
else
    echo "not matching"
    echo
    echo "The version has to match the following pattern:"
    echo "    MAJOR.MINOR.PATCH"
    echo "Example: 2.3.1"
    echo "Given version is $VERSION"
    exit 1
fi

REPO=$(git rev-parse --show-toplevel)
cd "$REPO"
echo "Root at $REPO"

echo -n "Checking location ... "
DEBLOC="tools/releasing/packaging/debian/changelog"
if [ -f "CMakeLists.txt" ] && [ -f "CHANGELOG.md" ] && [ -f "$DEBLOC" ]; then
    echo "all files found"
else
    echo "files missing"
    echo "Please execute this script in the project root directory!"
    exit 1
fi

echo "Bumping versions:"
VERSIONMAJOR=$(echo "$VERSION" | cut -d . -f 1)
VERSIONREGEX='[0-9]\+\.[0-9]\+\.[0-9]\+'
echo " - CMake"
sed "s/\(project(preCICE\s\+VERSION\s\+\)$VERSIONREGEX/\1$VERSION/" -i CMakeLists.txt
echo " - Changelog"
CL_FILES=$(find docs/changelog -type f -name "*.md")
CL_CONTENTS="## $VERSION\n"
echo "   compressing"
for cl_file in $CL_FILES; do
  # Extract the PR number from the filename
  cl_name=$(basename "$cl_file" .md)
  # Add the PR number to the end of every line in the file
  cl_entry=$(sed -s "s/$/ (https:\/\/github.com\/precice\/precice\/pull\/$cl_name)/" "$cl_file")
  if [[ -z "$cl_entry" ]]; then
    echo "WARNING $cl_file is empty"
  else
    CL_CONTENTS="$CL_CONTENTS\n$cl_entry"
  fi
  unset cl_name cl_entry
done
INSERT_LINE=4
echo -e "$(head -n${INSERT_LINE} CHANGELOG.md)\n\n$CL_CONTENTS\n$(tail +${INSERT_LINE} CHANGELOG.md)" > CHANGELOG.md
echo "   cleaning"
rm $CL_FILES
echo " - Debian Changelog"

DEBENTRY="
libprecice$VERSIONMAJOR ($VERSION) experimental; urgency=low

  * New upstream release $VERSION

 -- the preCICE developers <info@precice.org>  $(date -R)"
echo "$DEBENTRY" >> $DEBLOC

echo "done"
exit 0
