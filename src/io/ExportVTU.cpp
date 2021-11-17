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
#include "utils/MasterSlave.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace io {

void ExportVTU::doExport(
    const std::string &name,
    const std::string &location,
    mesh::Mesh &       mesh)
{
  PRECICE_TRACE(name, location, mesh.getName());
  processDataNamesAndDimensions(mesh);
  if (not location.empty())
    boost::filesystem::create_directories(location);
  if (utils::MasterSlave::isMaster()) {
    writeMasterFile(name, location, mesh);
  }
  if (mesh.vertices().size() > 0) { //only procs at the coupling interface should write output (for performance reasons)
    writeSubFile(name, location, mesh);
  }
}

void ExportVTU::processDataNamesAndDimensions(mesh::Mesh const &mesh)
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

void ExportVTU::writeMasterFile(
    const std::string &name,
    const std::string &location,
    mesh::Mesh &       mesh)
{
  namespace fs = boost::filesystem;
  fs::path outfile(location);
  outfile = outfile / fs::path(name + ".pvtu");
  std::ofstream outMasterFile(outfile.string(), std::ios::trunc);

  PRECICE_CHECK(outMasterFile, "VTU export failed to open master file \"{}\"", outfile);

  outMasterFile << "<?xml version=\"1.0\"?>\n";
  outMasterFile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"";
  outMasterFile << (utils::isMachineBigEndian() ? "BigEndian\">" : "LittleEndian\">") << '\n';
  outMasterFile << "   <PUnstructuredGrid GhostLevel=\"0\">\n";

  outMasterFile << "      <PPoints>\n";
  outMasterFile << "         <PDataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"" << 3 << "\"/>\n";
  outMasterFile << "      </PPoints>\n";

  outMasterFile << "      <PCells>\n";
  outMasterFile << "         <PDataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\"/>\n";
  outMasterFile << "         <PDataArray type=\"Int32\" Name=\"offsets\"      NumberOfComponents=\"1\"/>\n";
  outMasterFile << "         <PDataArray type=\"UInt8\" Name=\"types\"        NumberOfComponents=\"1\"/>\n";
  outMasterFile << "      </PCells>\n";

  // write scalar data names
  outMasterFile << "      <PPointData Scalars=\"Rank ";
  for (const auto &scalarDataName : _scalarDataNames) {
    outMasterFile << scalarDataName << ' ';
  }
  // write vector data names
  outMasterFile << "\" Vectors=\"";
  for (const auto &vectorDataName : _vectorDataNames) {
    outMasterFile << vectorDataName << ' ';
  }
  outMasterFile << "\">\n";

  outMasterFile << "         <PDataArray type=\"Int32\" Name=\"Rank\" NumberOfComponents=\"1\"/>\n";

  for (const auto &scalarDataName : _scalarDataNames) {
    outMasterFile << "         <PDataArray type=\"Float64\" Name=\"" << scalarDataName << "\" NumberOfComponents=\"" << 1 << "\"/>\n";
  }

  for (const auto &vectorDataName : _vectorDataNames) {
    outMasterFile << "         <PDataArray type=\"Float64\" Name=\"" << vectorDataName << "\" NumberOfComponents=\"" << 3 << "\"/>\n";
  }
  outMasterFile << "      </PPointData>\n";

  for (int i = 0; i < utils::MasterSlave::getSize(); i++) {
    if (mesh.getVertexDistribution()[i].size() > 0) { //only non-empty subfiles
      outMasterFile << "      <Piece Source=\"" << name << "_" << i << ".vtu\"/>\n";
    }
  }

  outMasterFile << "   </PUnstructuredGrid>\n";
  outMasterFile << "</VTKFile>\n";

  outMasterFile.close();
}

void ExportVTU::writeSubFile(
    const std::string &name,
    const std::string &location,
    mesh::Mesh &       mesh)
{
  int numPoints = mesh.vertices().size(); // number of vertices
  int numCells;                           // number of cells
  if (mesh.getDimensions() == 2) {
    numCells = mesh.edges().size();
  } else {
    numCells = mesh.triangles().size() + mesh.edges().size();
  }

  namespace fs = boost::filesystem;
  fs::path outfile(location);
  outfile = outfile / fs::path(name + "_" + std::to_string(utils::MasterSlave::getRank()) + ".vtu");
  std::ofstream outSubFile(outfile.string(), std::ios::trunc);

  PRECICE_CHECK(outSubFile, "VTU export failed to open slave file \"{}\"", outfile);

  outSubFile << "<?xml version=\"1.0\"?>\n";
  outSubFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"";
  outSubFile << (utils::isMachineBigEndian() ? "BigEndian\">" : "LittleEndian\">") << '\n';

  outSubFile << "   <UnstructuredGrid>\n";
  outSubFile << "      <Piece NumberOfPoints=\"" << numPoints << "\" NumberOfCells=\"" << numCells << "\"> \n";
  outSubFile << "         <Points> \n";
  outSubFile << "            <DataArray type=\"Float64\" Name=\"Position\" NumberOfComponents=\"" << 3 << "\" format=\"ascii\"> \n";
  for (const mesh::Vertex &vertex : mesh.vertices()) {
    writeVertex(vertex.getCoords(), outSubFile);
  }
  outSubFile << "            </DataArray>\n";
  outSubFile << "         </Points> \n\n";

  // Write Mesh
  exportMesh(outSubFile, mesh);

  // Write data
  exportData(outSubFile, mesh);

  outSubFile << "      </Piece>\n";
  outSubFile << "   </UnstructuredGrid> \n";
  outSubFile << "</VTKFile>\n";

  outSubFile.close();
}

void ExportVTU::exportMesh(
    std::ofstream &   outFile,
    mesh::Mesh const &mesh)
{
  if (mesh.getDimensions() == 2) { // write edges as cells
    outFile << "         <Cells>\n";
    outFile << "            <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    outFile << "               ";
    for (const mesh::Edge &edge : mesh.edges()) {
      writeLine(edge, outFile);
    }
    outFile << '\n';
    outFile << "            </DataArray> \n";
    outFile << "            <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    outFile << "               ";
    for (size_t i = 1; i <= mesh.edges().size(); i++) {
      outFile << 2 * i << "  ";
    }
    outFile << '\n';
    outFile << "            </DataArray>\n";
    outFile << "            <DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    outFile << "               ";
    for (size_t i = 1; i <= mesh.edges().size(); i++) {
      outFile << 3 << "  ";
    }
    outFile << '\n';
    outFile << "            </DataArray>\n";
    outFile << "         </Cells>\n";
  } else { // write triangles as cells

    outFile << "         <Cells>\n";
    outFile << "            <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    outFile << "               ";
    for (const mesh::Triangle &triangle : mesh.triangles()) {
      writeTriangle(triangle, outFile);
    }
    for (const mesh::Edge &edge : mesh.edges()) {
      writeLine(edge, outFile);
    }
    outFile << '\n';
    outFile << "            </DataArray> \n";
    outFile << "            <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    outFile << "               ";
    for (size_t i = 1; i <= mesh.triangles().size(); i++) {
      outFile << 3 * i << "  ";
    }
    for (size_t i = 1; i <= mesh.edges().size(); i++) {
      outFile << 2 * i << "  ";
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
    outFile << '\n';
    outFile << "            </DataArray>\n";
    outFile << "         </Cells>\n";
  }
}

void ExportVTU::exportData(
    std::ofstream &outFile,
    mesh::Mesh &   mesh)
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
  const auto rank = utils::MasterSlave::getRank();
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

void ExportVTU::writeVertex(
    const Eigen::VectorXd &position,
    std::ofstream &        outFile)
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

void ExportVTU::writeTriangle(
    const mesh::Triangle &triangle,
    std::ofstream &       outFile)
{
  outFile << triangle.vertex(0).getID() << "  ";
  outFile << triangle.vertex(1).getID() << "  ";
  outFile << triangle.vertex(2).getID() << "  ";
}

void ExportVTU::writeLine(
    const mesh::Edge &edge,
    std::ofstream &   outFile)
{
  outFile << edge.vertex(0).getID() << "  ";
  outFile << edge.vertex(1).getID() << "  ";
}

} // namespace io
} // namespace precice
