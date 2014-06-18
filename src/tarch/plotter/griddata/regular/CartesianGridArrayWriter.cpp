#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/plotter/griddata/regular/CartesianGridArrayWriter.h"


tarch::logging::Log tarch::plotter::griddata::regular::CartesianGridArrayWriter::_log( "tarch::plotter::griddata::regular::CartesianGridArrayWriter" );


tarch::plotter::griddata::regular::CartesianGridArrayWriter::DataSet::DataSet(
  const std::string&  identifier,
  int                 recordsPerEntry,
  const tarch::la::Vector<3,int>& numberOfGridPoints
):
  _identifier( identifier ),
  _recordsPerEntry( recordsPerEntry ) {
  assertion( _recordsPerEntry > 0 );
  assertion( tarch::la::volume( numberOfGridPoints ) > 0 );
  _data = new double[ _recordsPerEntry * tarch::la::volume( numberOfGridPoints ) ];
  for (int i=0; i<_recordsPerEntry * tarch::la::volume( numberOfGridPoints ); i++) {
	_data[i] = 0.0;
  }
}


void tarch::plotter::griddata::regular::CartesianGridArrayWriter::DataSet::clear() {
  delete []_data;
  _data = 0;

}


tarch::plotter::griddata::regular::CartesianGridArrayWriter::CartesianGridArrayWriter(
  const tarch::la::Vector<2,int>&     numberOfGridPoints,
  const tarch::la::Vector<2,double>&  domainSize,
  const tarch::la::Vector<2,double>&  origin
):
  _writtenToFile( false ),
  _numberOfGridPoints(1),
  _domainSize(0.0),
  _origin(0.0),
  _vertexData(),
  _cellData() {
  for (int d=0; d<2; d++) {
    _numberOfGridPoints(d) = numberOfGridPoints(d);
    _domainSize(d)         = domainSize(d);
    _origin(d)             = origin(d);
  }
}


tarch::plotter::griddata::regular::CartesianGridArrayWriter::CartesianGridArrayWriter(
  const tarch::la::Vector<3,int>&     numberOfGridPoints,
  const tarch::la::Vector<3,double>&  domainSize,
  const tarch::la::Vector<3,double>&  origin
):
  _writtenToFile( false ),
  _numberOfGridPoints( numberOfGridPoints ),
  _domainSize( domainSize ),
  _origin( origin ),
  _vertexData(),
  _cellData() {
}


tarch::plotter::griddata::regular::CartesianGridArrayWriter::~CartesianGridArrayWriter() {
	clear();
	// if (!_writtenToFile) {
    //_log.warning( "~CartesianGridArrayWriter()", "closed writer although data has not been written to any file" );
 // }
}


bool tarch::plotter::griddata::regular::CartesianGridArrayWriter::isOpen() {
  return !_writtenToFile;
}


void tarch::plotter::griddata::regular::CartesianGridArrayWriter::clear() {
  for (
    std::vector<DataSet>::iterator p = _cellData.begin();
    p != _cellData.end();
    p++
  ) {
    p->clear();
  }

  for (
    std::vector<DataSet>::iterator p = _vertexData.begin();
    p != _vertexData.end();
    p++
  ) {
    p->clear();
  }

  _cellData.clear();
  _vertexData.clear();
}


tarch::la::Vector<3,int> tarch::plotter::griddata::regular::CartesianGridArrayWriter::getNumberOfCells() const {
  tarch::la::Vector<3,int> result;
  result = _numberOfGridPoints-1;
  if (result(2)==0) result(2)=1;
  return result;
}


tarch::plotter::griddata::Writer::CellDataWriter*
tarch::plotter::griddata::regular::CartesianGridArrayWriter::createCellDataWriter( const std::string& identifier, int recordsPerCell ) {
  _cellData.push_back(
	DataSet( identifier, recordsPerCell, getNumberOfCells() )
  );
  return new tarch::plotter::griddata::regular::CartesianGridArrayWriter::CellDataWriter(_cellData.back(),*this);
}


tarch::plotter::griddata::Writer::VertexDataWriter*
tarch::plotter::griddata::regular::CartesianGridArrayWriter::createVertexDataWriter( const std::string& identifier, int recordsPerVertex ) {
  _vertexData.push_back(
    DataSet( identifier, recordsPerVertex, _numberOfGridPoints)
  );
  return new tarch::plotter::griddata::regular::CartesianGridArrayWriter::VertexDataWriter(_vertexData.back(),*this);
}


tarch::la::Vector<3,double> tarch::plotter::griddata::regular::CartesianGridArrayWriter::getH() const {
  logTraceInWith2Arguments( "getH()", _numberOfGridPoints, _domainSize );
  tarch::la::Vector<3,double> result;
  for (int d=0; d<3; d++) {
	assertion(_numberOfGridPoints(d)>0);
	double tmp = _numberOfGridPoints(d)-1;
	if ( tarch::la::greater(tmp,0.0) ) {
      result(d) = _domainSize(d) / tmp;
	}
	else {
      result(d) = 0.0;
	}
  }
  logTraceOutWith1Argument( "getH()", result );
  return result;
}


int tarch::plotter::griddata::regular::CartesianGridArrayWriter::getVertexIndex( const tarch::la::Vector<2,double>& position ) {
  tarch::la::Vector<3,double> position3D;
  position3D(0) = position(0);
  position3D(1) = position(1);
  position3D(2) = 0.0;
  return getVertexIndex(position3D);
}


int tarch::plotter::griddata::regular::CartesianGridArrayWriter::getVertexIndex( const tarch::la::Vector<3,double>& position ) {
  logTraceInWith1Argument( "getVertexIndex(Vector)", position );
  int result(0);
  int base = 1;
  tarch::la::Vector<3,double> h;
  h = getH();
  for (int d=0; d<3; d++) {
	double translatedPointPosition = position(d)-_origin(d) +h(d)/2.0 + tarch::la::NUMERICAL_ZERO_DIFFERENCE ;
	int index = 0;
	if (tarch::la::greater(h(d),0.0)) {
	  translatedPointPosition        = translatedPointPosition/h(d);
      index = static_cast<int>(floor(translatedPointPosition));
	}
    if (index<0) {
      index = 0;
    }
    if (index>=_numberOfGridPoints(d)) {
      index = _numberOfGridPoints(d)-1;
    }
    result += index * base;
    base   *= _numberOfGridPoints(d);
  }
  logTraceOutWith1Argument( "getVertexIndex(Vector)", result );
  return result;
}


int tarch::plotter::griddata::regular::CartesianGridArrayWriter::getCellIndex( const tarch::la::Vector<2,double>& position ) {
  tarch::la::Vector<3,double> position3D;
  position3D(0) = position(0);
  position3D(1) = position(1);
  position3D(2) = 0.0;
  return getVertexIndex(position3D);
}


int tarch::plotter::griddata::regular::CartesianGridArrayWriter::getCellIndex( const tarch::la::Vector<3,double>& position ) {
  logTraceInWith1Argument( "getCellIndex(Vector)", position );
  int result(0);
  int base = 1;
  tarch::la::Vector<3,double> h;
  h = getH();
  for (int d=0; d<3; d++) {
    double translatedPointPosition = position(d)-_origin(d) +h(d)/2.0 + tarch::la::NUMERICAL_ZERO_DIFFERENCE;
    int index = 0;
    if (tarch::la::greater(h(d),0.0)) {
      translatedPointPosition        = translatedPointPosition/h(d);
      index = static_cast<int>(floor(translatedPointPosition));
    }
    if (index<0) {
      index = 0;
    }
    if (index>=_numberOfGridPoints(d)-1) {
      index = _numberOfGridPoints(d)-2;
    }
    result += index * base;
    base   *= (_numberOfGridPoints(d)-1);
  }
  logTraceOutWith1Argument( "getCellIndex(Vector)", result);
  return result;
}


bool tarch::plotter::griddata::regular::CartesianGridArrayWriter::containsVertex(
  const tarch::la::Vector<3,double>& x
) const {
  return containsCell(x,tarch::la::Vector<3,double>(0.0));
}


bool tarch::plotter::griddata::regular::CartesianGridArrayWriter::containsVertex(
  const tarch::la::Vector<2,double>& x
) const {
  tarch::la::Vector<3,double> x3D;
  x3D(0) = x(0);
  x3D(1) = x(1);
  x3D(2) = 0.0;
  return containsVertex(x3D);
}


bool tarch::plotter::griddata::regular::CartesianGridArrayWriter::containsCell(
  const tarch::la::Vector<2,double>& offset,
  const tarch::la::Vector<2,double>& boundingBox
) const {
  tarch::la::Vector<3,double> offset3D;
  tarch::la::Vector<3,double> boundingBox3D;
  offset3D(0) = offset(0);
  offset3D(1) = offset(1);
  offset3D(2) = 0.0;
  boundingBox3D(0) = boundingBox(0);
  boundingBox3D(1) = boundingBox(1);
  boundingBox3D(2) = 0.0;
  return containsCell(offset3D,boundingBox3D);
}


bool tarch::plotter::griddata::regular::CartesianGridArrayWriter::containsCell(
  const tarch::la::Vector<3,double>& offset,
  const tarch::la::Vector<3,double>& boundingBox
) const {
  tarch::la::Vector<3,double> centerOfVoxel;
  tarch::la::Vector<3,double> hHalfOfVoxel;

  hHalfOfVoxel  = boundingBox / 2.0;
  centerOfVoxel = offset + hHalfOfVoxel;

  bool result = true;
  for (int d=0; d<3; d++) {
    result &= centerOfVoxel(d) >= _origin(d)-hHalfOfVoxel(d);
    result &= centerOfVoxel(d) <= _origin(d)+_domainSize(d)+hHalfOfVoxel(d);
  }
  return result;
}
