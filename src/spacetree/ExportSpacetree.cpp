// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ExportSpacetree.hpp"
#include "io/ExportVTK.hpp"
#include <fstream>

namespace precice {
namespace spacetree {

ExportSpacetree:: ExportSpacetree
(
  const std::string& filename )
:
  _filename(filename),
  _vertexCounter(0),
  _vertices(),
  _cells(),
  _cellPositions(),
  _cellContents()
{}

void ExportSpacetree:: nodeCallback
(
  const utils::DynVector& center,
  const utils::DynVector& halflengths,
  int                     position )
{
  // do nothing
}

void ExportSpacetree:: leafCallback
(
  const utils::DynVector& center,
  const utils::DynVector& halflengths,
  int                     position,
  const mesh::Group&      content )
{
  int dim = center.size();
  int twoPowerDim = std::pow(2.0, dim);
  utils::DynVector vertex(dim);
  std::vector<int> vertexIndices(twoPowerDim);
  for (int i=0; i < twoPowerDim; i++){
    tarch::la::multiplyComponents(utils::delinearize(i,dim), halflengths, vertex);
    vertex += center;
    _vertices.push_back(vertex);
    vertexIndices[i] = _vertexCounter;
    _vertexCounter ++;
  }
  _cells.push_back(vertexIndices);
  _cellPositions.push_back(position);
  _cellContents.push_back(content.size());
}

void ExportSpacetree:: doExport
(
  Spacetree& toExport )
{
  toExport.accept(*this); // Gather spacetree information

  std::ofstream file;
  io::ExportVTK::initializeWriting(_filename, file);
  io::ExportVTK::writeHeader(file);

  // Write cell corner vertices
  file << "POINTS " << _vertices.size() << " float "<< std::endl << std::endl;
  foreach (const utils::DynVector& vertexCoords, _vertices){
    io::ExportVTK::writeVertex(vertexCoords, file);
  }
  file << std::endl;

  // Write cells
  assertion(_vertexCounter > 0);
  int dim = _vertices.front().size();
  int twoPowerDim = std::pow(2.0, dim);
  int cellSize = _cells.size();
  file << "CELLS " << cellSize << " "
       << cellSize * (twoPowerDim + 1) << std::endl << std::endl;
  foreach (const std::vector<int>& cell, _cells){
    file << twoPowerDim << " ";
    for (int i=0; i < twoPowerDim; i++){
      file << cell[i] << " ";
    }
    file << std::endl;
  }
  file << std::endl;

  // Write cell types
  file << std::endl << "CELL_TYPES " << cellSize << std::endl << std::endl;
  int cellType;
  if (dim == 2){
    cellType = 8;
  }
  else {
    assertion(dim == 3);
    cellType = 11;
  }
  for (int i=0; i < cellSize; i++){
    file << cellType << std::endl;
  }
  file << std::endl;

  // Write cell data
  file << "CELL_DATA " << cellSize << std::endl << std::endl;

  // Write cell positions
  file << "SCALARS Position(0=Undef,1=Inside,2=Outside,3=On) float 1" << std::endl
       << "LOOKUP_TABLE default" << std::endl << std::endl;
  foreach (int position, _cellPositions){
    file << position << std::endl;
  }
  file << std::endl;

  // Write cell content size
  file << "SCALARS ContentSize float 1" << std::endl
       << "LOOKUP_TABLE default" << std::endl << std::endl;
  foreach (int size, _cellContents){
    file << size << std::endl;
  }
  file << std::endl;

  file.close();
}

//void ExportSpacetree:: exportCell
//(
//  const Spacetree& spacetree )
//{
//  int dim = spacetree.getCenter().size();
//  int twoPowerDim = std::pow(2, dim);
//  if ( spacetree.isLeaf() ) {
//    std::vector<int> vertexIndices(twoPowerDim);
//    for ( int i=0; i < twoPowerDim; i++ ) {
//      utils::DynVector corner(dim);
//      tarch::la::multiplyComponents (
//          utils::delinearize(i,dim), spacetree.getHalflengths(), corner );
//      corner += spacetree.getCenter();
//      _vertices.push_back ( corner );
//      vertexIndices[i] = _vertexCounter;
//      _vertexCounter ++;
//    }
//
//    _cells.push_back ( vertexIndices );
//    _cellPositions.push_back ( spacetree.getPosition() );
//    _cellContents.push_back ( spacetree.getContent().size() );
//
////    int cellIndex = -1;
////    if ( utils::Def::DIM == 2 ) {
////      _cells.push_back (  );
////      cellIndex = _cellWriter->plotQuadrangle ( vertexIndices );
////    }
////    else {
////      assertion ( utils::Def::DIM == 3 );
//////      cellIndex = _cellWriter->plotHexahedron ( vertexIndices );
////    }
////    assertion ( cellIndex != -1 );
////    _cellPositionWriter->plotCell ( cellIndex, spacetree.getPosition() );
////    _cellContentSizeWriter->plotCell ( cellIndex, spacetree.getContent().size() );
//
////#   ifdef PRECICE_STATISTICS
////    _ambiguoutiesWriter->plotCell ( cellIndex, spacetree.getAmbiguouties() );
////    _remainingCellsVisitedWriter->plotCell ( cellIndex, spacetree.getRemainingCellsVisits() );
////    _allCellsVisitedWriter->plotCell ( cellIndex, spacetree.getAllCellsVisits() );
////    _searchDistanceCallsWriter->plotCell ( cellIndex, spacetree.getSearchDistanceCalls() );
////#   endif
//  }
//  else {
//    for ( size_t i=0; i < spacetree.getChildCount(); i++ ) {
//      exportCell ( spacetree.getChild(i) );
//    }
//  }
//}

}}
