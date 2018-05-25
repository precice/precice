# Debian package skeleton for Ubuntu 18.04

These files can help to produce an **EXPERIMENTAL** debian binary package,
only meant to work with Ubuntu 18.04.
Please write us about your experience in [issue #125](https://github.com/precice/precice/issues/125).

## How to produce it

* Replace the (fake) `libprecice.so.1.1.1` under `debian/usr/lib/x86_64-linux-gnu/` with the actual `libprecice.so`.
* Add in `debian/usr/include/precice` the header files decsribed in `debian/usr/include/precice/README.md` and delete the respective `README.md`.
* Compress the changelog files under `debian/usr/share/doc/libprecice1` (without timestamps) using `gzip --best -n [file]`.
* Build it with `fakeroot dpkg-deb --build debian/`. `fakeroot` is needed to assign the correct owner/permissions.
* Check it for common (policy) errors with `mv debian.deb libprecice1_1.1.1-1_amd64.deb && lintian libprecice1_1.1.1-1_amd64.deb`

Notes on building:
* The appropriate package name can be defined with this command (copy-paste from the Debian packaging manual): `objdump -p libprecice.so.1.1.1 | sed -n -e's/^[[:space:]]*SONAME[[:space:]]*//p' | sed -r -e's/([0-9])\.so\./\1-/; s/\.so(\.|$)//; y/_/-/; s/(.*)/\L&/'`
* In order to read the symbols defined in the shared library, run `readelf -d libprecice.sorecice.so`.

Some helpful readings:
* [Debian Binary Package Building HOWTO, tldp.org](http://tldp.org/HOWTO/html_single/Debian-Binary-Package-Building-HOWTO/)
* [Debian package (.deb): debian directory files, rcramer.com](http://www.rcramer.com/tech/linux/deb_debian.shtml)
* [Debian New Maintainers' Guide](https://www.debian.org/doc/manuals/maint-guide/)
* [HowToPackageForDebian](https://wiki.debian.org/HowToPackageForDebian)

## Install / Uninstall

You can install the result package with a double click (Ubuntu/Gnome Software), or with `sudo apt install ./libprecice1_1.1.1-1_amd64.deb`. This will install:

* the shared library and some verion-related links into `/usr/lib/x86_64-linux-gnu/`,
* the pkgconfig metadata into `/usr/lib/x86_64-linux-gnu/pkgconfig/`,
* the header files `Constants.hpp`, `MeshHandle.hpp`, `SolverInterface.hpp` and the headers for the C, Fortran, Fortran2003 adapters (not the Python ones, as I was not sure) into `/usr/include/precice/`.
* changelogs and a copyright file.

You can also open it with your archive manager to see the contents.

You can remove it with `sudo apt remove libprecice1`.

## Compiling and Linking adapters

In order to compile an adapter, the path `/usr/include/precice/` needs to
be discoverable by your compiler. If this does not work by default,
you can add it in your `CPLUS_INCLUDE_PATH` or run the compiler using
the flags returned by `PKG_CONFIG_ALLOW_SYSTEM_CFLAGS=1 pkg-config --cflags precice`
(this will return `-I/usr/include/precice`).
The `PKG_CONFIG_ALLOW_SYSTEM_CFLAGS=1` may be needed because pkg-config
omits the default system paths in the output, but e.g. WMake (OpenFOAM)
does not search into them.

In order to link a solver with preCICE, the path `/usr/lib/x86_64-linux-gnu/`
needs to be discoverable by your linker.

## Issues

* There are no version checks for the dependencies. This will allow the package to be installed in older Ubuntu versions, but it will fail while linking/running.
* Everything is bundled in one package. It is suggested to have a `precice` package with the shared library and a `precice-dev` package with links to it, header files, and the static library. This would be a bit more complicated.
* This package is not produced by a respective source package, which is a bit more complicated, but is the proper way to do it (and is necessary to get accepted in distribution repositories).
* Ubuntu Software reports the license as "proprietary", but I would consider this a bug of Ubuntu Software (it does the same with many other packages).
* The preCICE and preCICE package changelogs could be improved.
* There are no file hashes (`md5sums`).
* There are no post-install or pre-remove files (not needed at the moment, though).
* See the respective overriden warnings/errors in `/usr/share/lintian/overrides`.
