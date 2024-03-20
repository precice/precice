#include "io/ExportXML.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <memory>
#include <string>
#include "io/Export.hpp"
#include "logging/LogMacros.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Tetrahedron.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"
#include "utils/Helpers.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::io {

void ExportXML::doExport(
    const std::string &name,
    const std::string &location,
    const mesh::Mesh & mesh)
{
  PRECICE_TRACE(name, location, mesh.getName());
  processDataNamesAndDimensions(mesh);
  if (not location.empty())
    std::filesystem::create_directories(location);
  if (utils::IntraComm::isPrimary()) {
    writeParallelFile(name, location, mesh);
  }
  if (mesh.nVertices() > 0) { // only procs at the coupling interface should write output (for performance reasons)
    writeSubFile(name, location, mesh);
  }
}

void ExportXML::processDataNamesAndDimensions(const mesh::Mesh &mesh)
{
  _vectorDataNames.clear();
  _scalarDataNames.clear();
  const bool isThreeDim = (mesh.getDimensions() == 3);
  for (const mesh::PtrData &data : mesh.data()) {
    int        dataDimensions = data->getDimensions();
    const bool hasGradient    = data->hasGradient();
    PRECICE_ASSERT(dataDimensions >= 1);
    std::string dataName = data->getName();
    if (dataDimensions == 1) {
      _scalarDataNames.push_back(dataName);
      if (hasGradient) {
        _vectorDataNames.push_back(dataName + "_gradient");
      }
    } else {
      _vectorDataNames.push_back(dataName);
      if (hasGradient) {
        _vectorDataNames.push_back(dataName + "_dx");
        _vectorDataNames.push_back(dataName + "_dy");
        if (isThreeDim) {
          _vectorDataNames.push_back(dataName + "_dz");
        }
      }
    }
  }
}

void ExportXML::writeParallelFile(
    const std::string &name,
    const std::string &location,
    const mesh::Mesh & mesh) const
{
  namespace fs = std::filesystem;
  fs::path outfile(location);
  outfile = outfile / fs::path(name + getParallelExtension());
  std::ofstream outParallelFile(outfile.string(), std::ios::trunc);

  PRECICE_CHECK(outParallelFile, "{} export failed to open primary file \"{}\"", getVTKFormat(), outfile.generic_string());

  const auto formatType = getVTKFormat();
  outParallelFile << "<?xml version=\"1.0\"?>\n";
  outParallelFile << "<VTKFile type=\"P" << formatType << "\" version=\"0.1\" byte_order=\"";
  outParallelFile << (utils::isMachineBigEndian() ? "BigEndian\">" : "LittleEndian\">") << '\n';
  outParallelFile << "   <P" << formatType << " GhostLevel=\"0\">\n";

  outParallelFile << "      <PPoints>\n";
  outParallelFile << "         <PDataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"" << 3 << "\"/>\n";
  outParallelFile << "      </PPoints>\n";

  writeParallelCells(outParallelFile);

  writeParallelData(outParallelFile);

  const auto &offsets = mesh.getVertexOffsets();
  PRECICE_ASSERT(offsets.size() > 0);
  if (offsets[0] > 0) {
    outParallelFile << "      <Piece Source=\"" << name << "_" << 0 << getPieceExtension() << "\"/>\n";
  }
  for (size_t rank : utils::IntraComm::allSecondaryRanks()) {
    PRECICE_ASSERT(rank < offsets.size());
    if (offsets[rank] - offsets[rank - 1] > 0) {
      // only non-empty subfiles
      outParallelFile << "      <Piece Source=\"" << name << "_" << rank << getPieceExtension() << "\"/>\n";
    }
  }

  outParallelFile << "   </P" << formatType << ">\n";
  outParallelFile << "</VTKFile>\n";

  outParallelFile.close();
}

namespace {
std::string getPieceSuffix()
{
  if (!utils::IntraComm::isParallel()) {
    return "";
  }
  return "_" + std::to_string(utils::IntraComm::getRank());
}
} // namespace

void ExportXML::writeSubFile(
    const std::string &name,
    const std::string &location,
    const mesh::Mesh & mesh) const
{
  namespace fs = std::filesystem;
  fs::path outfile(location);
  outfile /= fs::path(name + getPieceSuffix() + getPieceExtension());
  std::ofstream outSubFile(outfile.string(), std::ios::trunc);

  PRECICE_CHECK(outSubFile, "{} export failed to open secondary file \"{}\"", getVTKFormat(), outfile.generic_string());

  const auto formatType = getVTKFormat();
  outSubFile << "<?xml version=\"1.0\"?>\n";
  outSubFile << "<VTKFile type=\"" << formatType << "\" version=\"0.1\" byte_order=\"";
  outSubFile << (utils::isMachineBigEndian() ? "BigEndian\">" : "LittleEndian\">") << '\n';

  outSubFile << "   <" << formatType << ">\n";
  outSubFile << "      <Piece " << getPieceAttributes(mesh) << "> \n";
  exportPoints(outSubFile, mesh);

  // Write Mesh
  exportConnectivity(outSubFile, mesh);

  // Write data
  exportData(outSubFile, mesh);

  outSubFile << "      </Piece>\n";
  outSubFile << "   </" << formatType << "> \n";
  outSubFile << "</VTKFile>\n";

  outSubFile.close();
}

void ExportXML::exportGradient(const mesh::PtrData data, const int spaceDim, std::ostream &outFile) const
{
  const auto &             gradients      = data->gradients();
  const int                dataDimensions = data->getDimensions();
  std::vector<std::string> suffices;
  if (dataDimensions == 1) {
    suffices = {"_gradient"};
  } else if (spaceDim == 2) {
    suffices = {"_dx", "_dy"};
  } else if (spaceDim == 3) {
    suffices = {"_dx", "_dy", "_dz"};
  }
  int counter = 0; // Counter for multicomponent
  for (const auto &suffix : suffices) {
    const std::string dataName(data->getName());
    outFile << "            <DataArray type=\"Float64\" Name=\"" << dataName << suffix << "\" NumberOfComponents=\"" << 3;
    outFile << "\" format=\"ascii\">\n";
    outFile << "               ";
    for (int i = counter; i < gradients.cols(); i += spaceDim) { // Loop over vertices
      int j = 0;
      for (; j < gradients.rows(); j++) { // Loop over components
        outFile << gradients.coeff(j, i) << " ";
      }
      if (j < 3) { // If 2D data add additional zero as third component
        outFile << "0.0"
                << " ";
      }
    }
    outFile << '\n'
            << "            </DataArray>\n";
    counter++; // Increment counter for next component
  }
}

void ExportXML::exportData(
    std::ostream &    outFile,
    const mesh::Mesh &mesh) const
{
  outFile << "         <PointData Scalars=\"Rank ";
  for (const auto &scalarDataName : _scalarDataNames) {
    outFile << scalarDataName << ' ';
  }
  outFile << "\" Vectors=\"";
  for (const auto &vectorDataName : _vectorDataNames) {
    outFile << vectorDataName << ' ';
  }
  outFile << "\">\n";

  // Export the current rank
  outFile << "            <DataArray type=\"UInt32\" Name=\"Rank\" NumberOfComponents=\"1\" format=\"ascii\">\n";
  outFile << "               ";
  const auto rank = utils::IntraComm::getRank();
  for (size_t count = 0; count < mesh.nVertices(); ++count) {
    outFile << rank << ' ';
  }
  outFile << "\n            </DataArray>\n";

  for (const mesh::PtrData &data : mesh.data()) { // Plot vertex data
    Eigen::VectorXd &values         = data->values();
    int              dataDimensions = data->getDimensions();
    std::string      dataName(data->getName());
    int              numberOfComponents = (dataDimensions == 2) ? 3 : dataDimensions;
    const bool       hasGradient        = data->hasGradient();
    outFile << "            <DataArray type=\"Float64\" Name=\"" << dataName << "\" NumberOfComponents=\"" << numberOfComponents;
    outFile << "\" format=\"ascii\">\n";
    outFile << "               ";
    if (dataDimensions > 1) {
      Eigen::VectorXd viewTemp(dataDimensions);
      for (size_t count = 0; count < mesh.nVertices(); count++) {
        size_t offset = count * dataDimensions;
        for (int i = 0; i < dataDimensions; i++) {
          viewTemp[i] = values(offset + i);
        }
        for (int i = 0; i < dataDimensions; i++) {
          outFile << viewTemp[i] << ' ';
        }
        if (dataDimensions == 2) {
          outFile << "0.0" << ' '; // 2D data needs to be 3D for vtk
        }
        outFile << ' ';
      }
    } else if (dataDimensions == 1) {
      for (size_t count = 0; count < mesh.nVertices(); count++) {
        outFile << values(count) << ' ';
      }
    }
    outFile << '\n'
            << "            </DataArray>\n";
    if (hasGradient) {
      exportGradient(data, dataDimensions, outFile);
    }
  }
  outFile << "         </PointData> \n";
}

void ExportXML::writeVertex(
    const Eigen::VectorXd &position,
    std::ostream &         outFile)
{
  outFile << "               ";
  for (int i = 0; i < position.size(); i++) {
    outFile << position(i) << "  ";
  }
  if (position.size() == 2) {
    outFile << 0.0 << "  "; // also for 2D scenario, vtk needs 3D data
  }
  outFile << '\n';
}

void ExportXML::writeTriangle(
    const mesh::Triangle &triangle,
    std::ostream &        outFile)
{
  outFile << triangle.vertex(0).getID() << "  ";
  outFile << triangle.vertex(1).getID() << "  ";
  outFile << triangle.vertex(2).getID() << "  ";
}

void ExportXML::writeTetrahedron(
    const mesh::Tetrahedron &tetra,
    std::ostream &           outFile)
{
  outFile << tetra.vertex(0).getID() << "  ";
  outFile << tetra.vertex(1).getID() << "  ";
  outFile << tetra.vertex(2).getID() << "  ";
  outFile << tetra.vertex(3).getID() << "  ";
}

void ExportXML::writeLine(
    const mesh::Edge &edge,
    std::ostream &    outFile)
{
  outFile << edge.vertex(0).getID() << "  ";
  outFile << edge.vertex(1).getID() << "  ";
}

void ExportXML::exportPoints(
    std::ostream &    outFile,
    const mesh::Mesh &mesh) const
{
  outFile << "         <Points> \n";
  outFile << "            <DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"" << 3 << "\" format=\"ascii\"> \n";
  for (const mesh::Vertex &vertex : mesh.vertices()) {
    writeVertex(vertex.getCoords(), outFile);
  }
  outFile << "            </DataArray>\n";
  outFile << "         </Points> \n\n";
}

void ExportXML::writeParallelData(std::ostream &out) const
{
  // write scalar data names
  out << "      <PPointData Scalars=\"Rank ";
  for (const auto &scalarDataName : _scalarDataNames) {
    out << scalarDataName << ' ';
  }
  // write vector data names
  out << "\" Vectors=\"";
  for (const auto &vectorDataName : _vectorDataNames) {
    out << vectorDataName << ' ';
  }
  out << "\">\n";

  out << "         <PDataArray type=\"Int32\" Name=\"Rank\" NumberOfComponents=\"1\"/>\n";

  for (const auto &scalarDataName : _scalarDataNames) {
    out << "         <PDataArray type=\"Float64\" Name=\"" << scalarDataName << "\" NumberOfComponents=\"" << 1 << "\"/>\n";
  }

  for (const auto &vectorDataName : _vectorDataNames) {
    out << "         <PDataArray type=\"Float64\" Name=\"" << vectorDataName << "\" NumberOfComponents=\"" << 3 << "\"/>\n";
  }
  out << "      </PPointData>\n";
}

} // namespace precice::io
