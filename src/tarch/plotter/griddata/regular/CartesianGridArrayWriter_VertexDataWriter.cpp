#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/plotter/griddata/regular/CartesianGridArrayWriter.h"
#include <limits>

tarch::plotter::griddata::regular::CartesianGridArrayWriter::VertexDataWriter::VertexDataWriter(DataSet& dataSet, CartesianGridArrayWriter& parent):
  _dataSet( dataSet ),
  _parent(parent),
  _minValue(std::numeric_limits<double>::max()),
  _maxValue(std::numeric_limits<double>::min()) {
}


tarch::plotter::griddata::regular::CartesianGridArrayWriter::VertexDataWriter::~VertexDataWriter() {
}


void tarch::plotter::griddata::regular::CartesianGridArrayWriter::VertexDataWriter::close() {

}


void tarch::plotter::griddata::regular::CartesianGridArrayWriter::VertexDataWriter::plotVertex( int index, double value ) {
  if (value<_minValue) _minValue = value;
  if (value>_maxValue) _maxValue = value;
  _dataSet._data[ index * _dataSet._recordsPerEntry ] = value;
}


void tarch::plotter::griddata::regular::CartesianGridArrayWriter::VertexDataWriter::plotVertex( int index, const tarch::la::Vector<2,double>& value ) {
  for (int i=0; i<2; i++) {
    if (value(i)<_minValue) _minValue = value(i);
    if (value(i)>_maxValue) _maxValue = value(i);
    _dataSet._data[ index * _dataSet._recordsPerEntry + i] = value(i);
  }
}


void tarch::plotter::griddata::regular::CartesianGridArrayWriter::VertexDataWriter::plotVertex( int index, const tarch::la::Vector<3,double>& value ) {
  for (int i=0; i<3; i++) {
    if (value(i)<_minValue) _minValue = value(i);
    if (value(i)>_maxValue) _maxValue = value(i);
    _dataSet._data[ index * _dataSet._recordsPerEntry + i] = value(i);
  }
}


double tarch::plotter::griddata::regular::CartesianGridArrayWriter::VertexDataWriter::getMinValue() const {
  return _minValue;
}


double tarch::plotter::griddata::regular::CartesianGridArrayWriter::VertexDataWriter::getMaxValue() const {
  return _maxValue;
}


tarch::la::Vector<3,double> tarch::plotter::griddata::regular::CartesianGridArrayWriter::VertexDataWriter::getH() const {
  return _parent.getH();
}


int tarch::plotter::griddata::regular::CartesianGridArrayWriter::VertexDataWriter::getVertexIndex( const tarch::la::Vector<2,double>& position ) {
  return _parent.getVertexIndex( position );
}


int tarch::plotter::griddata::regular::CartesianGridArrayWriter::VertexDataWriter::getVertexIndex( const tarch::la::Vector<3,double>& position ) {
  return _parent.getVertexIndex(position);
}


bool tarch::plotter::griddata::regular::CartesianGridArrayWriter::VertexDataWriter::containsVertex( const tarch::la::Vector<3,double>& position ) const {
  return _parent.containsVertex( position );
}


bool tarch::plotter::griddata::regular::CartesianGridArrayWriter::VertexDataWriter::containsVertex( const tarch::la::Vector<2,double>& position ) const {
  return _parent.containsVertex( position );
}
