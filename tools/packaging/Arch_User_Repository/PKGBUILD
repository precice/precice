# Maintainer:  Florian Lindner <florian.lindner@xgm.de>

pkgname=precice
pkgver=1.6.1
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
sha256sums=('7d0c54faa2c69e52304f36608d93c408629868f16f3201f663a0f9b2008f0763')

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
