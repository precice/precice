#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/plotter/griddata/unstructured/vtk/VTKTextFileWriter.h"

#include <iomanip>

tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::VertexWriter::VertexWriter(VTKTextFileWriter& writer):
  _currentVertexNumber(0),
  _myWriter(writer),
  _out() {
  assertion( _myWriter._numberOfVertices==0 );
  _out << std::setprecision(_myWriter._precision);
}


tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::VertexWriter::~VertexWriter() {
  if (_currentVertexNumber>=0) {
    close();
  }
}



int tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::VertexWriter::plotVertex(const tarch::la::Vector<2,double>& position) {
  assertion1( _currentVertexNumber>=0, _currentVertexNumber );

  tarch::la::Vector<3,double> p;
  p(0) = position(0);
  p(1) = position(1);
  p(2) = 0.0;

  return plotVertex(p);
}


int tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::VertexWriter::plotVertex(const tarch::la::Vector<3,double>& position) {
  assertion( _currentVertexNumber>=0 );
  _currentVertexNumber++;
  _out << position(0) << "  "
       << position(1) << "  "
       << position(2) << std::endl;
  return _currentVertexNumber-1;
}


void tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::VertexWriter::close() {
  assertion( _myWriter._numberOfVertices==0 );
  assertionMsg( _myWriter.isOpen(), "Maybe you forgot to call close() on a data writer before you destroy your writer?" );


  _myWriter._numberOfVertices  = _currentVertexNumber;
  _myWriter._vertexDescription = _out.str();
  _currentVertexNumber         = -1;
}
