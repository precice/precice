#include "ExportVTKXML.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Quad.hpp"
#include <Eigen/Core>
#include <string>
#include <fstream>
#include <boost/filesystem.hpp>
#include "utils/Helpers.hpp"

namespace precice {
namespace io {

logging::Logger ExportVTKXML:: _log("io::ExportVTKXML");

ExportVTKXML:: ExportVTKXML
(
  bool writeNormals )
:
  Export(),
  _writeNormals(writeNormals),
  _meshDimensions(-1)
{
}

int ExportVTKXML:: getType() const
{
  return constants::exportVTKXML();
}

void ExportVTKXML:: doExport
(
  const std::string& name,
  const std::string& location,
  mesh::Mesh&        mesh)
{
  TRACE(name, location, mesh.getName());
  assertion(utils::MasterSlave::_slaveMode || utils::MasterSlave::_masterMode);
  processDataNamesAndDimensions(mesh);
  if (utils::MasterSlave::_masterMode) {
    writeMasterFile(name, location, mesh);
  }
  if(mesh.vertices().size()>0){ //only procs at the coupling interface should write output (for performance reasons)
    writeSubFile(name, location, mesh);
  }
}

void ExportVTKXML::processDataNamesAndDimensions
(
  mesh::Mesh& mesh)
{
  _meshDimensions = mesh.getDimensions();
  _vectorDataNames.clear();
  _scalarDataNames.clear();
  if (_writeNormals) {
    _vectorDataNames.push_back("VertexNormals ");
  }
  for (mesh::PtrData data : mesh.data()) {
    int dataDimensions = data->getDimensions();
    assertion(dataDimensions>=1);
    std::string dataName = data->getName();
    if ( dataDimensions == 1) {
      _scalarDataNames.push_back(dataName);
    } else {
      _vectorDataNames.push_back(dataName);
    }
  }
}

void ExportVTKXML::writeMasterFile
(
  const std::string& name,
  const std::string& location,
  mesh::Mesh&        mesh)
{
  namespace fs = boost::filesystem;
  fs::path outfile(location);
  outfile = outfile / fs::path(name + "_master.pvtu");
  std::ofstream outMasterFile(outfile.string(), std::ios::trunc);

  CHECK(outMasterFile, "Could not open master file \"" << outfile.c_str() << "\" for VTKXML export!");

  outMasterFile << "<?xml version=\"1.0\"?>" << std::endl;
  outMasterFile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"";
  outMasterFile << (utils::isMachineBigEndian() ? "BigEndian\">" : "LittleEndian\">")  << std::endl;
  outMasterFile << "   <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;

  outMasterFile << "      <PPoints>" << std::endl;
  outMasterFile << "         <PDataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"" << 3 << "\"/>" << std::endl;
  outMasterFile << "      </PPoints>" << std::endl;

  outMasterFile << "      <PCells>" << std::endl;
  outMasterFile << "         <PDataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\"/>" << std::endl;
  outMasterFile << "         <PDataArray type=\"Int32\" Name=\"offsets\"      NumberOfComponents=\"1\"/>" << std::endl;
  outMasterFile << "         <PDataArray type=\"UInt8\" Name=\"types\"        NumberOfComponents=\"1\"/>" << std::endl;
  outMasterFile << "      </PCells>" << std::endl;

  // write scalar data names
  outMasterFile << "      <PPointData Scalars=\"";
  for (size_t i = 0; i < _scalarDataNames.size(); ++i) {
    outMasterFile << _scalarDataNames[i] << " ";
  }
  // write vector data names
  outMasterFile << "\" Vectors=\"";
  for (size_t i = 0; i < _vectorDataNames.size(); ++i) {
    outMasterFile << _vectorDataNames[i] << " ";
  }
  outMasterFile << "\">" << std::endl;

  for (size_t i = 0; i < _scalarDataNames.size(); ++i) {
    outMasterFile << "         <PDataArray type=\"Float32\" Name=\""<< _scalarDataNames[i] << "\" NumberOfComponents=\"" << 1 << "\"/>" << std::endl;
  }

  for (size_t i = 0; i < _vectorDataNames.size(); ++i) {
    outMasterFile << "         <PDataArray type=\"Float32\" Name=\""<< _vectorDataNames[i] << "\" NumberOfComponents=\"" << 3 << "\"/>" << std::endl;
  }
  outMasterFile << "      </PPointData>" << std::endl;

  for (int i = 0; i < utils::MasterSlave::_size; i++) {
    if(mesh.getVertexDistribution()[i].size()>0){ //only non-empty subfiles
      outMasterFile << "      <Piece Source=\"" << name << "_r" << i << ".vtu\"/>" << std::endl;
    }
  }

  outMasterFile << "   </PUnstructuredGrid>" << std::endl;
  outMasterFile << "</VTKFile>" << std::endl;

  outMasterFile.close();
}

void ExportVTKXML::writeSubFile
(
  const std::string& name,
  const std::string& location,
  mesh::Mesh&        mesh)
{
  int numPoints = mesh.vertices().size(); // number of vertices
  int numCells; // number of cells
  if (_meshDimensions == 2) {
    numCells = mesh.edges().size();
  } else {
    numCells = mesh.triangles().size() + mesh.quads().size();
  }

  namespace fs = boost::filesystem;
  fs::path outfile(location);
  outfile = outfile / fs::path(name + "_r" + std::to_string(utils::MasterSlave::_rank) + ".vtu");
  std::ofstream outSubFile(outfile.string(), std::ios::trunc);

  CHECK(outSubFile, "Could not open slave file \"" << outfile.c_str() << "\" for VTKXML export!");

  outSubFile << "<?xml version=\"1.0\"?>" << std::endl;
  outSubFile << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"";
  outSubFile << (utils::isMachineBigEndian() ? "BigEndian\">" : "LittleEndian\">")  << std::endl;

  outSubFile << "   <UnstructuredGrid>" << std::endl;
  outSubFile << "      <Piece NumberOfPoints=\"" << numPoints << "\" NumberOfCells=\"" << numCells << "\"> " << std::endl;
  outSubFile << "         <Points> " << std::endl;
  outSubFile << "            <DataArray type=\"Float32\" Name=\"Position\" NumberOfComponents=\"" << 3 << "\" format=\"ascii\"> " << std::endl;
  for (mesh::Vertex& vertex : mesh.vertices()) {
    writeVertex(vertex.getCoords(), outSubFile);
  }
  outSubFile << "            </DataArray>" << std::endl;
  outSubFile << "         </Points> " << std::endl << std::endl;

  // Write Mesh
  exportMesh(outSubFile, mesh);

  // Write data
  exportData(outSubFile, mesh);

  outSubFile << "      </Piece>" << std::endl;
  outSubFile << "   </UnstructuredGrid> " << std::endl;
  outSubFile << "</VTKFile>" << std::endl;

  outSubFile.close();
}

void ExportVTKXML::exportMesh
(
  std::ofstream& outFile,
  mesh::Mesh&    mesh)
{
  if (_meshDimensions == 2) { // write edges as cells
    outFile << "         <Cells>" << std::endl;
    outFile << "            <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
    outFile << "               ";
    for (mesh::Edge & edge : mesh.edges()) {
      writeLine(edge, outFile);
    }
    outFile << std::endl;
    outFile << "            </DataArray> " << std::endl;
    outFile << "            <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
    outFile << "               ";
    for (size_t i = 1; i <= mesh.edges().size(); i++) {
      outFile << 2*i << "  ";
    }
    outFile << std::endl;
    outFile << "            </DataArray>" << std::endl;
    outFile << "            <DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
    outFile << "               ";
    for (size_t i = 1; i <= mesh.edges().size(); i++) {
      outFile << 3 << "  ";
    }
    outFile << std::endl;
    outFile << "            </DataArray>" << std::endl;
    outFile << "         </Cells>" << std::endl;
  } else { // write triangles and quads as cells

    outFile << "         <Cells>" << std::endl;
    outFile << "            <DataArray type=\"Int32\" Name=\"connectivity\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
    outFile << "               ";
    for (mesh::Triangle& triangle : mesh.triangles()) {
      writeTriangle(triangle, outFile);
    }
    for (mesh::Quad& quad : mesh.quads()) {
      writeQuadrangle(quad, outFile);
    }
    outFile << std::endl;
    outFile << "            </DataArray> " << std::endl;
    outFile << "            <DataArray type=\"Int32\" Name=\"offsets\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
    outFile << "               ";
    for (size_t i = 1; i <= mesh.triangles().size(); i++) {
      outFile << 3*i << "  ";
    }
    for (size_t i = 1; i <= mesh.quads().size(); i++) {
      outFile << 4*i << "  ";
    }
    outFile << std::endl;
    outFile << "            </DataArray>" << std::endl;
    outFile << "            <DataArray type=\"UInt8\"  Name=\"types\" NumberOfComponents=\"1\" format=\"ascii\">" << std::endl;
    outFile << "               ";
    for (size_t i = 1; i <= mesh.triangles().size(); i++) {
      outFile << 5 << "  ";
    }
    for (size_t i = 1; i <= mesh.quads().size(); i++) {
      outFile << 9 << "  ";
    }
    outFile << std::endl;
    outFile << "            </DataArray>" << std::endl;
    outFile << "         </Cells>" << std::endl;
  }
}

void ExportVTKXML:: exportData
(
  std::ofstream& outFile,
  mesh::Mesh&    mesh)
{
  outFile << "         <PointData Scalars=\"";
  for (size_t i = 0; i < _scalarDataNames.size(); i++) {
    outFile << _scalarDataNames[i] << " ";
  }
  outFile << "\" Vectors=\"";
  for (size_t i = 0; i < _vectorDataNames.size(); i++) {
    outFile << _vectorDataNames[i] << " ";
  }
  outFile << "\">" << std::endl;

  for (mesh::PtrData data : mesh.data()) { // Plot vertex data
    Eigen::VectorXd& values = data->values();
    int dataDimensions = data->getDimensions();
    std::string dataName(data->getName());
    int numberOfComponents = (dataDimensions==2) ? 3 : dataDimensions;
    outFile << "            <DataArray type=\"Float32\" Name=\"" << dataName << "\" NumberOfComponents=\"" << numberOfComponents;
    outFile << "\" format=\"ascii\">" << std::endl;
    outFile << "               ";
    if(dataDimensions > 1) {
      Eigen::VectorXd viewTemp(dataDimensions);
      for (size_t count = 0; count < mesh.vertices().size(); count++) {
        size_t offset = count * dataDimensions;
        for(int i=0; i < dataDimensions; i++){
          viewTemp[i] = values(offset + i);
        }
        for(int i = 0; i < dataDimensions; i++){
          outFile << viewTemp[i] << " ";
        }
        if(dataDimensions == 2){
          outFile << "0.0" << " "; //2D data needs to be 3D for vtk
        }
        outFile << " ";
      }
    }
    else if(dataDimensions == 1) {
      for (size_t count = 0; count < mesh.vertices().size(); count++) {
        outFile << values(count) << " ";
      }
    }
    outFile << std::endl << "            </DataArray>" << std::endl;
  }
  outFile << "         </PointData> " << std::endl;
}

void ExportVTKXML::writeVertex
(
  const Eigen::VectorXd& position,
  std::ofstream&         outFile)
{
  outFile << "               ";
  for (int i = 0; i < position.size(); i++){
    outFile << position(i) << "  ";
  }
  if (position.size() ==2)
  {
    outFile << 0.0 << "  ";  //also for 2D scenario, vtk needs 3D data
  }
  outFile << std::endl;
}


void ExportVTKXML:: writeTriangle
(
  mesh::Triangle& triangle,
  std::ofstream&  outFile)
{
  outFile << triangle.vertex(0).getID() << "  ";
  outFile << triangle.vertex(1).getID() << "  ";
  outFile << triangle.vertex(2).getID() << "  ";

}

void ExportVTKXML:: writeQuadrangle
(
  mesh::Quad&    quad,
  std::ofstream& outFile)
{
  outFile << quad.vertex(0).getID() << "  ";
  outFile << quad.vertex(1).getID() << "  ";
  outFile << quad.vertex(2).getID() << "  ";
  outFile << quad.vertex(3).getID() << "  ";
}

void ExportVTKXML:: writeLine
(
  mesh::Edge& edge,
  std::ofstream& outFile)
{
  outFile << edge.vertex(0).getID() << "  ";
  outFile << edge.vertex(1).getID() << "  ";
}

}} // namespace precice, io
