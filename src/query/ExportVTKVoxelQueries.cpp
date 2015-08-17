// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ExportVTKVoxelQueries.hpp"
//#include "tarch/plotter/griddata/unstructured/vtk/VTKTextFileWriter.h"
#include "boost/scoped_ptr.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace query {

void ExportVTKVoxelQueries:: addQuery
(
   const utils::DynVector& voxelCenter,
   const utils::DynVector& voxelHalflengths,
   int                     containedElementsCount )
{
   _voxelCenters.push_back ( voxelCenter );
   _voxelHalflengths.push_back ( voxelHalflengths );
   _containedElementCount.push_back ( containedElementsCount );
}

void ExportVTKVoxelQueries:: exportQueries
(
  std::string filename )
{
//  using tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter;
//  using tarch::plotter::griddata::unstructured::UnstructuredGridWriter;
//  using tarch::plotter::griddata::Writer;
//  typedef boost::scoped_ptr<UnstructuredGridWriter::VertexWriter> PtrVertexWriter;
//  typedef boost::scoped_ptr<UnstructuredGridWriter::CellWriter> PtrCellWriter;
//  typedef boost::scoped_ptr<Writer::CellDataWriter> PtrCellDataWriter;

//  VTKTextFileWriter vtkWriter; // (filename + ".vtk");
//  PtrVertexWriter vertexWriter ( vtkWriter .createVertexWriter() );
//  PtrCellWriter cellWriter ( vtkWriter.createCellWriter() );
//  PtrCellDataWriter countContainedWriter (
//      vtkWriter.createCellDataWriter("ContainedVisitables", 1) );

  assertion (_voxelCenters.size() == _voxelHalflengths.size());
  assertion (_voxelCenters.size() == _containedElementCount.size());

  // for (auto & elem : _voxelCenters) {
//    int vertexIndices[utils::Def::TWO_POWER_DIM];
//    Vector result;
//    for (int iCorner=0; iCorner < utils::Def::TWO_POWER_DIM; iCorner++) {
//      Vector delin;
//      utils::delinearize(iCorner, delin);
//      Vector corner = _voxelCenters[i] + tarch::la::multiplyComponents(
//        delin, _voxelHalflengths[i], result);
//    }

//    int cellIndex = -1; // = cellWriter.getNextFreeElementNumber();
//    if ( utils::Def::DIM == 2 ) {
//      cellIndex = cellWriter->plotQuadrangle ( vertexIndices );
//    }
//    else {
//      assertion ( utils::Def::DIM == 3 );
//      cellIndex = cellWriter->plotHexahedron ( vertexIndices );
//    }
//    int containedVisitables = _containedElementCount[i];
//    assertion ( cellIndex != -1 );
//    countContainedWriter->plotCell ( cellIndex, containedVisitables );
  // }
//  vtkWriter.writeToFile ( filename + ".vtk" );
}

//void ExportVTKVoxelQueries:: exportVoxels ( std::string filename )
//{
//   tarch::plotter::VTKWriter vtkWriter (filename);
//   tarch::plotter::VTKWriter::VertexWriter vertexWriter (vtkWriter);
//   tarch::plotter::VTKWriter::ElementWriter cellWriter (vtkWriter);
//   tarch::plotter::VTKWriter::ScalarCellDataWriter countContainedWriter (vtkWriter, "ContainedVisitables");
//
//   assertion (_voxelCenters.size() == _voxelHalflengths.size());
//
//   for (size_t i=0; i < _voxelCenters.size(); i++) {
//      int vertexIndices[TWO_POWER_D];
//      for (int iCorner=0; iCorner < TWO_POWER_D; iCorner++) {
//         Vector corner = _voxelCenters[i] +
//                         utils::delinearize(iCorner) * _voxelHalflengths[i];
//         vertexIndices[iCorner] = vertexWriter.getNextFreeVertexNumber();
//         vertexWriter.plotVertex (vertexIndices[iCorner], corner);
//      }
//
//      int cellIndex = cellWriter.getNextFreeElementNumber();
//      cellWriter.plotQuadrangle (cellIndex, vertexIndices);
//   }
//   vtkWriter.plotVertices (vertexWriter);
//   vtkWriter.plotElements (cellWriter);
//}

void ExportVTKVoxelQueries:: resetQueries ()
{
  _voxelCenters.clear ();
  _voxelHalflengths.clear ();
  _containedElementCount.clear ();
}

}}
