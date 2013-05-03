#COMSOL_INCLUDE_DIR="/import/home/student/comsol33/script/external"
#INCLUDE_DIR="/import/home/benk/include/Comsol"
#PEANO_DIR="/import/home/benk/FSI/Peano/peano/trunk/src"
#PRECICE_DIR="/import/home/benk/FSI/Peano/precice/trunk/src/precice/adapters/c"
#PEANO_LIB_DIR="/import/home/benk/FSI/Peano/peano/trunk/build/release/dim2"
#PRECICE_LIB_DIR="/import/home/benk/FSI/Peano/precice/trunk/build/release-dim2-mpi-serial"
#FSI_LIB_DIR="/import/home/benk/lib/Comsol"

COMSOL_INCLUDE_DIR="/import/home/student/comsol33/script/external"
INCLUDE_DIR="/import/home/benk/include/Comsol"
PEANO_DIR="/work/benk/Peano_build/peano/trunk/src"
#PRECICE_DIR="/work/benk/Peano_build/precice/src/precice/adapters/c"
PRECICE_DIR="/import/home/benk/FSI/Peano/precice/src/precice/adapters/c/"
PEANO_LIB_DIR="/work/benk/Peano_build/peano/build/debug/dim2/lib"
#PRECICE_LIB_DIR="/work/benk/Peano_build/precice/build/release-dim2-mpi-serial"
PRECICE_LIB_DIR="/work/benk/Peano_build/precice/build/debug-dim2-mpi-serial"
FSI_LIB_DIR="/import/home/benk/lib/Comsol"


COMP_COMAND="mpicc -c -fPIC -I$COMSOL_INCLUDE_DIR -I/usr/include -I$INCLUDE_DIR -I/usr/include/glib-2.0 -I/usr/include/glib-2.0/glib -I/usr/lib/glib-2.0/include -I$PRECICE_DIR -I$PEANO_DIR  $1.c"
# -I/usr/include/glib-2.0/glib  -I/import/home/benk/include/Comsol -I/usr/include
LINK_COMAND="mpicxx $1.o -o $1 $PRECICE_LIB_DIR/libprecice.a $PEANO_LIB_DIR/libpeano.a /usr/lib/libglib-2.0.so $FSI_LIB_DIR/libfsi_com.a $FSI_LIB_DIR/libfsi_mesh.a"
#/import/home/benk/lib/Comsol/libfsi_com.a /import/home/benk/lib/Comsol/libfsi_mesh.a /usr/lib/libglib-2.0.so
echo $COMP_COMAND
$COMP_COMAND
echo $LINK_COMAND
$LINK_COMAND
