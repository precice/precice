#! /bin/bash

set -e

USAGE="
usage: $(basename $0) MAJOR.MINOR.PATCH

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

echo -n "Checking location ... "
DEBLOC="tools/packaging/debian/changelog"
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
sed "s/## develop/&\n\n## $VERSION/" -i CHANGELOG.md
echo " - Debian Changelog"

DEBENTRY="
libprecice$VERSIONMAJOR ($VERSION) experimental; urgency=low

  * New upstream release $VERSION

 -- the preCICE developers <info@precice.org>  $(date -R)"
echo "$DEBENTRY" >> $DEBLOC

echo "done"
exit 0
