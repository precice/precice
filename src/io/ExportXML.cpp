#include "io/ExportXML.hpp"
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

namespace precice {
namespace io {

void ExportXML::doExport(
    const std::string &name,
    const std::string &location,
    const mesh::Mesh & mesh)
{
  PRECICE_TRACE(name, location, mesh.getName());
  processDataNamesAndDimensions(mesh);
  if (not location.empty())
    boost::filesystem::create_directories(location);
  if (utils::IntraComm::isPrimary()) {
    writeParallelFile(name, location, mesh);
  }
  if (mesh.vertices().size() > 0) { //only procs at the coupling interface should write output (for performance reasons)
    writeSubFile(name, location, mesh);
  }
}

void ExportXML::processDataNamesAndDimensions(const mesh::Mesh &mesh)
{
  _vectorDataNames.clear();
  _scalarDataNames.clear();
  for (const mesh::PtrData &data : mesh.data()) {
    int dataDimensions = data->getDimensions();
    PRECICE_ASSERT(dataDimensions >= 1);
    std::string dataName = data->getName();
    if (dataDimensions == 1) {
      _scalarDataNames.push_back(dataName);
    } else {
      _vectorDataNames.push_back(dataName);
    }
  }
}

void ExportXML::writeParallelFile(
    const std::string &name,
    const std::string &location,
    const mesh::Mesh & mesh) const
{
  namespace fs = boost::filesystem;
  fs::path outfile(location);
  outfile = outfile / fs::path(name + getParallelExtension());
  std::ofstream outPrimaryFile(outfile.string(), std::ios::trunc);

  PRECICE_CHECK(outPrimaryFile, "{} export failed to open master file \"{}\"", getVTKFormat(), outfile);

  const auto formatType = getVTKFormat();
  outPrimaryFile << "<?xml version=\"1.0\"?>\n";
  outPrimaryFile << "<VTKFile type=\"P" << formatType << "\" version=\"0.1\" byte_order=\"";
  outPrimaryFile << (utils::isMachineBigEndian() ? "BigEndian\">" : "LittleEndian\">") << '\n';
  outPrimaryFile << "   <P" << formatType << " GhostLevel=\"0\">\n";

  outPrimaryFile << "      <PPoints>\n";
  outPrimaryFile << "         <PDataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"" << 3 << "\"/>\n";
  outPrimaryFile << "      </PPoints>\n";

  writePrimaryCells(outPrimaryFile);

  writePrimaryData(outPrimaryFile);

  const auto &offsets = mesh.getVertexOffsets();
  PRECICE_ASSERT(offsets.size() > 0);
  if (offsets[0] > 0) {
    outPrimaryFile << "      <Piece Source=\"" << name << "_" << 0 << getPieceExtension() << "\"/>\n";
  }
  for (auto rank : utils::IntraComm::allSecondaries()) {
    PRECICE_ASSERT(rank < offsets.size());
    if (offsets[rank] - offsets[rank - 1] > 0) {
      //only non-empty subfiles
      outPrimaryFile << "      <Piece Source=\"" << name << "_" << rank << getPieceExtension() << "\"/>\n";
    }
  }

  outPrimaryFile << "   </P" << formatType << ">\n";
  outPrimaryFile << "</VTKFile>\n";

  outPrimaryFile.close();
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
  namespace fs = boost::filesystem;
  fs::path outfile(location);
  outfile /= fs::path(name + getPieceSuffix() + getPieceExtension());
  std::ofstream outSubFile(outfile.string(), std::ios::trunc);

  PRECICE_CHECK(outSubFile, "{} export failed to open slave file \"{}\"", getVTKFormat(), outfile);

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
  for (size_t count = 0; count < mesh.vertices().size(); ++count) {
    outFile << rank << ' ';
  }
  outFile << "\n            </DataArray>\n";

  for (const mesh::PtrData &data : mesh.data()) { // Plot vertex data
    Eigen::VectorXd &values         = data->values();
    int              dataDimensions = data->getDimensions();
    std::string      dataName(data->getName());
    int              numberOfComponents = (dataDimensions == 2) ? 3 : dataDimensions;
    outFile << "            <DataArray type=\"Float64\" Name=\"" << dataName << "\" NumberOfComponents=\"" << numberOfComponents;
    outFile << "\" format=\"ascii\">\n";
    outFile << "               ";
    if (dataDimensions > 1) {
      Eigen::VectorXd viewTemp(dataDimensions);
      for (size_t count = 0; count < mesh.vertices().size(); count++) {
        size_t offset = count * dataDimensions;
        for (int i = 0; i < dataDimensions; i++) {
          viewTemp[i] = values(offset + i);
        }
        for (int i = 0; i < dataDimensions; i++) {
          outFile << viewTemp[i] << ' ';
        }
        if (dataDimensions == 2) {
          outFile << "0.0" << ' '; //2D data needs to be 3D for vtk
        }
        outFile << ' ';
      }
    } else if (dataDimensions == 1) {
      for (size_t count = 0; count < mesh.vertices().size(); count++) {
        outFile << values(count) << ' ';
      }
    }
    outFile << '\n'
            << "            </DataArray>\n";
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
    outFile << 0.0 << "  "; //also for 2D scenario, vtk needs 3D data
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

void ExportXML::writePrimaryData(std::ostream &out) const
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

} // namespace io
} // namespace precice
