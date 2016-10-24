#include "ExportSpacetree.hpp"
#include "io/ExportVTK.hpp"
#include <fstream>
#include <boost/filesystem.hpp>

namespace precice {
namespace spacetree {

ExportSpacetree:: ExportSpacetree
(
  const std::string& location,
  const std::string& filename )
:
  _location(location),
  _filename(filename),
  _vertexCounter(0),
  _vertices(),
  _cells(),
  _cellPositions(),
  _cellContents()
{}

void ExportSpacetree:: nodeCallback
(
  const Eigen::VectorXd& center,
  const Eigen::VectorXd& halflengths,
  int                     position )
{
  // do nothing
}

void ExportSpacetree:: leafCallback
(
  const Eigen::VectorXd& center,
  const Eigen::VectorXd& halflengths,
  int                    position,
  const mesh::Group&     content )
{
  int dim = center.size();
  int twoPowerDim = std::pow(2.0, dim);
  Eigen::VectorXd vertex(dim);
  std::vector<int> vertexIndices(twoPowerDim);
  for (int i=0; i < twoPowerDim; i++){
    vertex = utils::delinearize(i, dim).cwiseProduct(halflengths);
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

  namespace fs = boost::filesystem;
  fs::path outfile(_location);
  outfile = outfile / fs::path(_filename);
  std::ofstream outstream(outfile.string(), std::ios::trunc);

  io::ExportVTK::initializeWriting(outstream);
  io::ExportVTK::writeHeader(outstream);

  // Write cell corner vertices
  outstream << "POINTS " << _vertices.size() << " float "<< std::endl << std::endl;
  for (const Eigen::VectorXd& vertexCoords : _vertices){
    io::ExportVTK::writeVertex(vertexCoords, outstream);
  }
  outstream << std::endl;

  // Write cells
  assertion(_vertexCounter > 0);
  int dim = _vertices.front().size();
  int twoPowerDim = std::pow(2.0, dim);
  int cellSize = _cells.size();
  outstream << "CELLS " << cellSize << " "
       << cellSize * (twoPowerDim + 1) << std::endl << std::endl;
  for (const std::vector<int>& cell : _cells){
    outstream << twoPowerDim << " ";
    for (int i=0; i < twoPowerDim; i++){
      outstream << cell[i] << " ";
    }
    outstream << std::endl;
  }
  outstream << std::endl;

  // Write cell types
  outstream << std::endl << "CELL_TYPES " << cellSize << std::endl << std::endl;
  int cellType;
  if (dim == 2){
    cellType = 8;
  }
  else {
    assertion(dim == 3);
    cellType = 11;
  }
  for (int i=0; i < cellSize; i++){
    outstream << cellType << std::endl;
  }
  outstream << std::endl;

  // Write cell data
  outstream << "CELL_DATA " << cellSize << std::endl << std::endl;

  // Write cell positions
  outstream << "SCALARS Position(0=Undef,1=Inside,2=Outside,3=On) float 1" << std::endl
       << "LOOKUP_TABLE default" << std::endl << std::endl;
  for (int position : _cellPositions){
    outstream << position << std::endl;
  }
  outstream << std::endl;

  // Write cell content size
  outstream << "SCALARS ContentSize float 1" << std::endl
       << "LOOKUP_TABLE default" << std::endl << std::endl;
  for (int size : _cellContents){
    outstream << size << std::endl;
  }
  outstream << std::endl;

  outstream.close();
}

}}

