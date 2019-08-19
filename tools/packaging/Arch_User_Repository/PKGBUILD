# Maintainer:  Florian Lindner <florian.lindner@xgm.de>

pkgname=precice
pkgver=1.4.0
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
sha256sums=('3499bfc0941fb9f004d5e32eb63d64f93e17b4057fab3ada1cde40c8311bd466')

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
