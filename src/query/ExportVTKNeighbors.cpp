#include "ExportVTKNeighbors.hpp"
//#include "tarch/plotter/griddata/unstructured/vtk/VTKTextFileWriter.h"

namespace precice {
namespace query {

void ExportVTKNeighbors:: addNeighbors
(
  const Eigen::VectorXd& queryPoint,
  const ClosestElement&  closestNeighbor )
{
//  _neighbors += std::pair<utils::DynVector,ClosestElement>(queryPoint,closestNeighbor);
}

void ExportVTKNeighbors:: exportNeighbors
(
  const std::string& filename )
{
  //  using tarch::plotter::griddata::unstructured::vtk::VTKTextFileWriter;
  //  using tarch::plotter::griddata::unstructured::UnstructuredGridWriter;
  //  using tarch::plotter::griddata::Writer;
  //  typedef boost::scoped_ptr<UnstructuredGridWriter::VertexWriter> PtrVertexWriter;
  //  typedef boost::scoped_ptr<UnstructuredGridWriter::CellWriter> PtrCellWriter;
  //  typedef boost::scoped_ptr<Writer::CellDataWriter> PtrCellDataWriter;
  //  VTKTextFileWriter vtkWriter;
  //  PtrVertexWriter vertexWriter ( vtkWriter.createVertexWriter() );
  //  PtrCellWriter lineWriter ( vtkWriter.createCellWriter() );
  //  PtrCellDataWriter distanceWriter (
  //      vtkWriter.createCellDataWriter("DistanceToGeometry", 1) );

//  int vertexIndices[2];
//  int numLine;
//  for (size_t i=0; i < _neighbors.size(); i++) {
//    utils::DynVector searchPoint = _neighbors[i].first;
//    utils::DynVector neighborPoint = searchPoint + _neighbors[i].second.vectorToElement;
    //    vertexIndices[0] = vertexWriter->plotVertex (searchPoint);
    //    vertexIndices[1] = vertexWriter->plotVertex (neighborPoint);
    //    numLine = lineWriter->plotLine (vertexIndices);
    //    distanceWriter->plotCell ( numLine, _neighbors[i].second.distance );
//  }
  //  vtkWriter.writeToFile ( filename + ".vtk" );
}

}}
