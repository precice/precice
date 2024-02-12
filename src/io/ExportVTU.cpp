#include "io/ExportVTU.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <fstream>
#include <memory>
#include <string>
#include "io/Export.hpp"
#include "logging/LogMacros.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "utils/Helpers.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::io {

std::string ExportVTU::getVTKFormat() const
{
  return "UnstructuredGrid";
}

std::string ExportVTU::getParallelExtension() const
{
  return ".pvtu";
}

std::string ExportVTU::getPieceExtension() const
{
  return ".vtu";
}

std::string ExportVTU::getPieceAttributes(const mesh::Mesh &mesh) const
{
  std::ostringstream oss;
  oss << "NumberOfPoints=\"" << mesh.nVertices() << "\" ";
  oss << "NumberOfCells=\"" << mesh.edges().size() + mesh.triangles().size() + mesh.tetrahedra().size() << "\" ";
  return oss.str();
}

void ExportVTU::writeParallelCells(std::ostream &out) const
{
  out << "      <PCells>\n";
  out << "         <PDataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\"/>\n";
  out << "         <PDataArray type=\"Int32\" Name=\"offsets\"      NumberOfComponents=\"1\"/>\n";
  out << "         <PDataArray type=\"UInt8\" Name=\"types\"        NumberOfComponents=\"1\"/>\n";
  out << "      </PCells>\n";
}

void ExportVTU::exportConnectivity(
    std::ostream &    outFile,
    const mesh::Mesh &mesh) const
{
  outFile << "         <Cells>\n";
  outFile << "            <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  outFile << "               ";
  for (const mesh::Triangle &triangle : mesh.triangles()) {
    writeTriangle(triangle, outFile);
  }
  for (const mesh::Edge &edge : mesh.edges()) {
    writeLine(edge, outFile);
  }
  for (const mesh::Tetrahedron &tetra : mesh.tetrahedra()) {
    writeTetrahedron(tetra, outFile);
  }
  outFile << '\n';
  outFile << "            </DataArray> \n";
  outFile << "            <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  outFile << "               ";
  for (size_t i = 1; i <= mesh.triangles().size(); i++) {
    outFile << 3 * i << "  ";
  }
  const auto triangleOffset = 3 * mesh.triangles().size();
  for (size_t i = 1; i <= mesh.edges().size(); i++) {
    outFile << 2 * i + triangleOffset << "  ";
  }
  const auto tetraOffset = 2 * mesh.edges().size() + triangleOffset;
  for (size_t i = 1; i <= mesh.tetrahedra().size(); i++) {
    outFile << 4 * i + tetraOffset << "  ";
  }
  outFile << '\n';
  outFile << "            </DataArray>\n";
  outFile << "            <DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  outFile << "               ";
  for (size_t i = 1; i <= mesh.triangles().size(); i++) {
    outFile << 5 << "  ";
  }
  for (size_t i = 1; i <= mesh.edges().size(); i++) {
    outFile << 3 << "  ";
  }
  for (size_t i = 1; i <= mesh.tetrahedra().size(); i++) {
    outFile << 10 << "  ";
  }
  outFile << '\n';
  outFile << "            </DataArray>\n";
  outFile << "         </Cells>\n";
}
} // namespace precice::io
