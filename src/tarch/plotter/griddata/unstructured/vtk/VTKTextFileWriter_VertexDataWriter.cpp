#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/plotter/griddata/unstructured/vtk/VTKTextFileWriter.h"

#include <limits>
#include <iomanip>

tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::VertexDataWriter::VertexDataWriter(
  const std::string& dataIdentifier, VTKTextFileWriter& writer, int recordsPerVertex
):
  _lastWriteCommandVertexNumber(-1),
  _myWriter(writer),
  _out(),
  _recordsPerVertex(recordsPerVertex),
  _minValue(std::numeric_limits<double>::max()),
  _maxValue(std::numeric_limits<double>::min()) {
  assertion(_recordsPerVertex>0);

  _out << std::setprecision(_myWriter._precision);
  if (_recordsPerVertex!=3) {
    _out << "SCALARS " << dataIdentifier << " " << _myWriter._doubleOrFloat << " " << _recordsPerVertex << std::endl
         << "LOOKUP_TABLE default" << std::endl;
  }
  else {
    _out << "VECTORS " << dataIdentifier << " " << _myWriter._doubleOrFloat << " " << std::endl;
  }

  #ifdef Asserts
  _dataIdentifier = dataIdentifier;
  #endif
}


tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::VertexDataWriter::~VertexDataWriter() {
  if (_lastWriteCommandVertexNumber>=-1) {
    close();
  }
}


void tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::VertexDataWriter::close() {
  assertionEquals2( _lastWriteCommandVertexNumber, _myWriter._numberOfVertices-1, _dataIdentifier, "perhaps not all vertices were assigned data (holds if _myWriter._numberOfVertices is bigger than local attribute)" );
  assertionMsg( _myWriter.isOpen(), "Maybe you forgot to call close() on a data writer before you destroy your writer for value " << _dataIdentifier );

  if (_lastWriteCommandVertexNumber>=-1) {
    _out << std::endl;
    _myWriter._vertexDataDescription += _out.str();
  }
  _lastWriteCommandVertexNumber = -2;
}


void tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::VertexDataWriter::plotVertex( int index, double value ) {
  assertion(_lastWriteCommandVertexNumber>=-1);
  assertion(1<=_recordsPerVertex);

  assertion2( value != std::numeric_limits<double>::infinity(), index, value);
  assertion2( value == value, index, value);  // test for not a number

  while (_lastWriteCommandVertexNumber<index-1) {
    plotVertex(_lastWriteCommandVertexNumber+1,0.0);
  }

  _lastWriteCommandVertexNumber = index;
  _out << value << " ";
  for (int i=1; i<_recordsPerVertex; i++) {
    _out << 0.0 << " ";
  }

  _out << std::endl;

  if (value<_minValue) _minValue = value;
  if (value>_maxValue) _maxValue = value;
}


void tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::VertexDataWriter::plotVertex( int index, const tarch::la::Vector<2,double>& value ) {
  assertion(_lastWriteCommandVertexNumber>=-1);
  assertion(2<=_recordsPerVertex);

  assertion1( value(0) != std::numeric_limits<double>::infinity(), value(0) );
  assertion1( value(0) == value(0), value(0) );  // test for not a number

  assertion1( value(1) != std::numeric_limits<double>::infinity(), value(1) );
  assertion1( value(1) == value(1), value(1) );  // test for not a number

  while (_lastWriteCommandVertexNumber<index-1) {
    plotVertex(_lastWriteCommandVertexNumber+1,0.0);
  }

  _lastWriteCommandVertexNumber = index;
  _out << value(0) << " ";
  _out << value(1) << " ";
  for (int i=2; i<_recordsPerVertex; i++) {
    _out << 0.0 << " ";
  }
  _out << std::endl;

  if (value(0)<_minValue) _minValue = value(0);
  if (value(0)>_maxValue) _maxValue = value(0);
  if (value(1)<_minValue) _minValue = value(1);
  if (value(1)>_maxValue) _maxValue = value(1);
}


void tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::VertexDataWriter::plotVertex( int index, const tarch::la::Vector<3,double>& value ) {
  assertion(_lastWriteCommandVertexNumber>=-1);
  assertion(3<=_recordsPerVertex);

  assertion1( value(0) != std::numeric_limits<double>::infinity(), value(0) );
  assertion1( value(0) == value(0), value(0) );  // test for not a number

  assertion1( value(1) != std::numeric_limits<double>::infinity(), value(1) );
  assertion1( value(1) == value(1), value(1) );  // test for not a number

  assertion1( value(2) != std::numeric_limits<double>::infinity(), value(2) );
  assertion1( value(2) == value(2), value(2) );  // test for not a number

  while (_lastWriteCommandVertexNumber<index-1) {
    plotVertex(_lastWriteCommandVertexNumber+1,0.0);
  }

  _lastWriteCommandVertexNumber = index;
  _out << value(0) << " ";
  _out << value(1) << " ";
  _out << value(2) << " ";
  for (int i=3; i<_recordsPerVertex; i++) {
    _out << 0.0 << " ";
  }
  _out << std::endl;

  if (value(0)<_minValue) _minValue = value(0);
  if (value(0)>_maxValue) _maxValue = value(0);
  if (value(1)<_minValue) _minValue = value(1);
  if (value(1)>_maxValue) _maxValue = value(1);
  if (value(2)<_minValue) _minValue = value(2);
  if (value(2)>_maxValue) _maxValue = value(2);
}


double tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::VertexDataWriter::getMinValue() const {
  return _minValue;
}


double tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::VertexDataWriter::getMaxValue() const {
  return _maxValue;
}
