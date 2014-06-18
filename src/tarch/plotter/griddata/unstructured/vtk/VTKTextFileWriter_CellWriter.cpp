#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/plotter/griddata/unstructured/vtk/VTKTextFileWriter.h"


tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::CellWriter::CellWriter(VTKTextFileWriter& writer):
  _currentCellNumber(0),
  _myWriter(writer),
  _cellListEntries(0),
  _cellOut(),
  _cellTypeOut() {
  assertion( _myWriter._numberOfCells==0 );
  assertion( _myWriter._numberOfCellEntries==0 );
}


tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::CellWriter::~CellWriter() {
  if (_currentCellNumber>=0) {
    close();
  }
}


int tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::CellWriter::plotPoint(int vertexIndex) {
  assertion( _currentCellNumber>=0 );
  assertion( _cellListEntries>=0 );

  _currentCellNumber++;
  _cellListEntries += 2;

  _cellOut << "1" << " "
           << vertexIndex << " "
           << std::endl;

  _cellTypeOut << "1" << std::endl;

  return _currentCellNumber-1;
}


int tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::CellWriter::plotHexahedron(int vertexIndex[8]) {
  assertion( _currentCellNumber>=0 );
  assertion( _cellListEntries>=0 );

  _currentCellNumber++;
  _cellListEntries += 9;

  _cellOut << "8" << " "
           << vertexIndex[0] << " "
           << vertexIndex[1] << " "
           << vertexIndex[2] << " "
           << vertexIndex[3] << " "
           << vertexIndex[4] << " "
           << vertexIndex[5] << " "
           << vertexIndex[6] << " "
           << vertexIndex[7] << " "
           << std::endl;

  _cellTypeOut << "11" << std::endl;

  return _currentCellNumber-1;
}


int tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::CellWriter::plotQuadrangle(int vertexIndex[4]) {
  assertion( _currentCellNumber>=0 );
  assertion( _cellListEntries>=0 );

  _currentCellNumber++;
  _cellListEntries += 5;

  _cellOut << "4" << " "
           << vertexIndex[0] << " "
           << vertexIndex[1] << " "
           << vertexIndex[2] << " "
           << vertexIndex[3] << " "
           << std::endl;

  _cellTypeOut << "8" << std::endl;

  return _currentCellNumber-1;
}


int tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::CellWriter::plotLine(int vertexIndex[2]) {
  assertion( _currentCellNumber>=0 );
  assertion( _cellListEntries>=0 );

  _currentCellNumber++;
  _cellListEntries += 3;

  _cellOut << "2" << " "
           << vertexIndex[0] << " "
           << vertexIndex[1] << " "
           << std::endl;

  _cellTypeOut << "3" << std::endl;

  return _currentCellNumber-1;
}


int tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::CellWriter::plotTriangle(int vertexIndex[3]) {
  assertion( _currentCellNumber>=0 );
  assertion( _cellListEntries>=0 );

  _currentCellNumber++;
  _cellListEntries += 4;

  _cellOut << "3" << " "
           << vertexIndex[0] << " "
           << vertexIndex[1] << " "
           << vertexIndex[2] << " "
           << std::endl;

  _cellTypeOut << "5" << std::endl;

  return _currentCellNumber-1;
}


void tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter::CellWriter::close() {
  assertion( _myWriter._numberOfCells==0 );
  assertion( _myWriter._numberOfCellEntries==0 );
  assertionMsg( _myWriter.isOpen(), "Maybe you forgot to call close() on a data writer before you destroy your writer?" );

  _myWriter._numberOfCells       = _currentCellNumber;
  _myWriter._numberOfCellEntries = _cellListEntries;

  _myWriter._cellDescription      = _cellOut.str();
  _myWriter._cellTypeDescription  = _cellTypeOut.str();

  _currentCellNumber = -1;
  _cellListEntries   = -1;
}
