# Maintainer:  Florian Lindner <florian.lindner@xgm.de>

pkgname=precice
pkgver=1.5.1
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
sha256sums=('3124642f03f181e20823da3fbaeaf85136964ea717f239509842603de9b40baa')

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
