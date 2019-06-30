# Maintainer:  Florian Lindner <florian.lindner@xgm.de>

pkgname=precice
pkgver=1.5.0
pkgrel=1
pkgdesc="A Coupling Library for Partitioned Multi-Physics Simulations on Massively Parallel Systems"
arch=('x86_64')
url="https://www.precice.org"
license=('LGPL3')
depends=('openmpi' 'boost')
makedepends=('cmake')
optdepends=('')
provides=('precice')
source=("https://github.com/precice/precice/archive/v${pkgver}.tar.gz")
sha256sums=('a2a794becd08717e3049252134ae35692fed71966ed32e22cca796a169c16c3e')

build() {
    cd "${pkgname}-${pkgver}"
    cmake . \
	  -DCMAKE_BUILD_TYPE=Release \
	  -DCMAKE_INSTALL_PREFIX=/usr \
          -DPETSC=off \
          -DPYTHON=off
    make
}

package() {
    cd "${pkgname}-${pkgver}"
    make DESTDIR="${pkgdir}/" install
}
