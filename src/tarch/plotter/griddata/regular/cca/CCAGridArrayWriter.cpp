#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
/*
 * CCAGridArrayWriter.cpp
 *
 *  Created on: Dec 1, 2010
 *      Author: atanasoa
 */


#include "tarch/plotter/griddata/regular/cca/CCAGridArrayWriter.h"





tarch::plotter::griddata::regular::cca::CCAGridArrayWriter::CCAGridArrayWriter(
    const tarch::la::Vector<2,int>&     numberOfGridPoints,
    const tarch::la::Vector<2,double>&  domainSize,
    const tarch::la::Vector<2,double>&  origin
):
tarch::plotter::griddata::regular::CartesianGridArrayWriter(
    numberOfGridPoints,
    domainSize,
    origin
) {
}


tarch::plotter::griddata::regular::cca::CCAGridArrayWriter::CCAGridArrayWriter(
    const tarch::la::Vector<3,int>&     numberOfGridPoints,
    const tarch::la::Vector<3,double>&  domainSize,
    const tarch::la::Vector<3,double>&  origin
):
tarch::plotter::griddata::regular::CartesianGridArrayWriter(
    numberOfGridPoints,
    domainSize,
    origin
) {
}

tarch::plotter::griddata::regular::cca::CCAGridArrayWriter::~CCAGridArrayWriter() {

}

void tarch::plotter::griddata::regular::cca::CCAGridArrayWriter::writeToVertexArray( double* array, const int length) {

  if (!_vertexData.empty()) {
      int index=0;
      for (int i=0; i<static_cast<int>(_vertexData.size()); i++) {
          for (int j=0; j<tarch::la::volume(_numberOfGridPoints); j++) {
              for (int k=0; k<_vertexData[i]._recordsPerEntry; k++) {
                  array[index++]=_vertexData[i]._data[j*_vertexData[i]._recordsPerEntry+k];
              }

          }

      }
  }





}

void tarch::plotter::griddata::regular::cca::CCAGridArrayWriter::writeToCellArray( double* array, const int length) {

  if (!_cellData.empty()) {
      int index=0;
      for (int i=0; i<static_cast<int>(_cellData.size()); i++) {

          for (int j=0; j<tarch::la::volume(_numberOfGridPoints-1); j++) {
              for (int k=0; k<_cellData[i]._recordsPerEntry; k++) {
                  array[index++]= _cellData[i]._data[j*_cellData[i]._recordsPerEntry+k];
              }

          }

      }
  }

}
void tarch::plotter::griddata::regular::cca::CCAGridArrayWriter::writeToFile( const std::string& filename ) {
}

