comp=gcc #icc
#linking=gcc
linking=gfortran
COMSOL_DIR="/import/home/student/comsol33"
COMSOL_INCLUDE_DIR="/import/home/benk/include/Comsol"
LIB_DIR="/import/home/benk/lib/Comsol"
#GLIB_INCLUDE_DIR="/usr/include/glib-2.0"

other_includes="-I/usr/include -I$COMSOL_INCLUDE_DIR -I/usr/include/glib-2.0 -I/usr/include/glib-2.0/glib -I/usr/lib/glib-2.0/include"
my_include="-I$COMSOL_DIR/script/external $other_includes"

#other_lib="/import/home/benk/Work_/lib/libglib-2.0.so /import/home/benk/Work_/lib/libfsi_com.a /import/home/benk/Work_/lib/libfsi_mesh.a" 
#other_lib="/usr/lib/libglib-2.0.so /import/home/benk/Work/lib/libfsi_com.so /import/home/benk/Work/lib/libfsi_mesh.so" 
other_lib="/usr/lib/libglib-2.0.so $LIB_DIR/libfsi_com.a $LIB_DIR/libfsi_mesh.a" 
#my_lib="-L/import/ftp/installed/comsol33/lib/glnx86 -lflscriptext -L/usr/local/lib/ -lmpich $other_lib"
my_lib="-L$COMSOL_DIR/lib/glnx86 -lflscriptext -L/usr/local/lib/ -lmpich $other_lib"

comp_comand="$comp -c -fPIC $my_include $1.c"
link_comand="$linking -shared $1.o -o $1.so $my_lib"
echo $comp_comand
echo $link_comand
$comp_comand
$link_comand