#include "ExportVTKXML.hpp"
#include <Eigen/Core>
#include <algorithm>
#include <boost/filesystem.hpp>
#include <fstream>
#include <memory>
#include <string>
#include "Constants.hpp"
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

ExportVTKXML::ExportVTKXML(
    bool writeNormals)
    : Export(),
      _writeNormals(writeNormals),
      _meshDimensions(-1)
{
}

int ExportVTKXML::getType() const
{
  return constants::exportVTKXML();
}

void ExportVTKXML::doExport(
    const std::string &name,
    const std::string &location,
    mesh::Mesh &       mesh)
{
  PRECICE_TRACE(name, location, mesh.getName());
  PRECICE_ASSERT(utils::MasterSlave::isSlave() || utils::MasterSlave::isMaster());
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

void ExportVTKXML::processDataNamesAndDimensions(mesh::Mesh const &mesh)
{
  _meshDimensions = mesh.getDimensions();
  _vectorDataNames.clear();
  _scalarDataNames.clear();
  if (_writeNormals) {
    _vectorDataNames.push_back("VertexNormals");
  }
  for (mesh::PtrData data : mesh.data()) {
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

void ExportVTKXML::writeMasterFile(
    const std::string &name,
    const std::string &location,
    mesh::Mesh &       mesh)
{
  namespace fs = boost::filesystem;
  fs::path outfile(location);
  outfile = outfile / fs::path(name + "_master.pvtu");
  std::ofstream outMasterFile(outfile.string(), std::ios::trunc);

  PRECICE_CHECK(outMasterFile, "VTKXML export failed to open master file \"" << outfile << '"');

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
  outMasterFile << "      <PPointData Scalars=\"";
  for (size_t i = 0; i < _scalarDataNames.size(); ++i) {
    outMasterFile << _scalarDataNames[i] << ' ';
  }
  // write vector data names
  outMasterFile << "\" Vectors=\"";
  for (size_t i = 0; i < _vectorDataNames.size(); ++i) {
    outMasterFile << _vectorDataNames[i] << ' ';
  }
  outMasterFile << "\">\n";

  for (size_t i = 0; i < _scalarDataNames.size(); ++i) {
    outMasterFile << "         <PDataArray type=\"Float64\" Name=\"" << _scalarDataNames[i] << "\" NumberOfComponents=\"" << 1 << "\"/>\n";
  }

  for (size_t i = 0; i < _vectorDataNames.size(); ++i) {
    outMasterFile << "         <PDataArray type=\"Float64\" Name=\"" << _vectorDataNames[i] << "\" NumberOfComponents=\"" << 3 << "\"/>\n";
  }
  outMasterFile << "      </PPointData>\n";

  for (int i = 0; i < utils::MasterSlave::getSize(); i++) {
    if (mesh.getVertexDistribution()[i].size() > 0) { //only non-empty subfiles
      outMasterFile << "      <Piece Source=\"" << name << "_r" << i << ".vtu\"/>\n";
    }
  }

  outMasterFile << "   </PUnstructuredGrid>\n";
  outMasterFile << "</VTKFile>\n";

  outMasterFile.close();
}

void ExportVTKXML::writeSubFile(
    const std::string &name,
    const std::string &location,
    mesh::Mesh &       mesh)
{
  int numPoints = mesh.vertices().size(); // number of vertices
  int numCells;                           // number of cells
  if (_meshDimensions == 2) {
    numCells = mesh.edges().size();
  } else {
    numCells = mesh.triangles().size();
  }

  namespace fs = boost::filesystem;
  fs::path outfile(location);
  outfile = outfile / fs::path(name + "_r" + std::to_string(utils::MasterSlave::getRank()) + ".vtu");
  std::ofstream outSubFile(outfile.string(), std::ios::trunc);

  PRECICE_CHECK(outSubFile, "VTKXML export failed to open slave file \"" << outfile << '"');

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

void ExportVTKXML::exportMesh(
    std::ofstream &   outFile,
    mesh::Mesh const &mesh)
{
  if (_meshDimensions == 2) { // write edges as cells
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
    outFile << '\n';
    outFile << "            </DataArray> \n";
    outFile << "            <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    outFile << "               ";
    for (size_t i = 1; i <= mesh.triangles().size(); i++) {
      outFile << 3 * i << "  ";
    }
    outFile << '\n';
    outFile << "            </DataArray>\n";
    outFile << "            <DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">\n";
    outFile << "               ";
    for (size_t i = 1; i <= mesh.triangles().size(); i++) {
      outFile << 5 << "  ";
    }
    outFile << '\n';
    outFile << "            </DataArray>\n";
    outFile << "         </Cells>\n";
  }
}

void ExportVTKXML::exportData(
    std::ofstream &outFile,
    mesh::Mesh &   mesh)
{
  outFile << "         <PointData Scalars=\"";
  for (size_t i = 0; i < _scalarDataNames.size(); i++) {
    outFile << _scalarDataNames[i] << ' ';
  }
  outFile << "\" Vectors=\"";
  for (size_t i = 0; i < _vectorDataNames.size(); i++) {
    outFile << _vectorDataNames[i] << ' ';
  }
  outFile << "\">\n";

  // Print VertexNormals
  if (_writeNormals) {
    const auto dimensions = mesh.getDimensions();
    outFile << "            <DataArray type=\"Float64\" Name=\"VertexNormals\" NumberOfComponents=\"";
    outFile << std::max(3, dimensions) << "\" format=\"ascii\">\n";
    outFile << "               ";
    for (const auto &vertex : mesh.vertices()) {
      const auto normal = vertex.getNormal();
      for (int i = 0; i < std::max(3, dimensions); i++) {
        if (i < dimensions) {
          outFile << normal[i] << ' ';
        } else {
          outFile << "0.0 ";
        }
        outFile << ' ';
      }
    }
    outFile << '\n'
            << "            </DataArray>\n";
  }
  for (mesh::PtrData data : mesh.data()) { // Plot vertex data
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

void ExportVTKXML::writeVertex(
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

void ExportVTKXML::writeTriangle(
    const mesh::Triangle &triangle,
    std::ofstream &       outFile)
{
  outFile << triangle.vertex(0).getID() << "  ";
  outFile << triangle.vertex(1).getID() << "  ";
  outFile << triangle.vertex(2).getID() << "  ";
}

void ExportVTKXML::writeLine(
    const mesh::Edge &edge,
    std::ofstream &   outFile)
{
  outFile << edge.vertex(0).getID() << "  ";
  outFile << edge.vertex(1).getID() << "  ";
}

} // namespace io
} // namespace precice
