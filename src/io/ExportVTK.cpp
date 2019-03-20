#include "ExportVTK.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Quad.hpp"
#include "Constants.hpp"
#include <Eigen/Core>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <boost/filesystem.hpp>

namespace precice {
namespace io {


ExportVTK:: ExportVTK
(
  bool writeNormals )
:
  Export(),
  _writeNormals(writeNormals)
{}

int ExportVTK:: getType() const
{
  return constants::exportVTK();
}

void ExportVTK:: doExport
(
  const std::string& name,
  const std::string& location,
  mesh::Mesh&        mesh)
{
  TRACE(name, location, mesh.getName());
  assertion(name != std::string(""));

  namespace fs = boost::filesystem;
  fs::path outfile(location);
  outfile = outfile / fs::path(name + ".vtk");
  std::ofstream outstream(outfile.string(), std::ios::trunc);
  CHECK(outstream, "Could not open file \"" << outfile.c_str() << "\" for VTK export!");

  initializeWriting(outstream);
  writeHeader(outstream);
  exportMesh(outstream, mesh);
  exportData(outstream, mesh);
  outstream.close();
}

void ExportVTK::exportMesh
(
  std::ofstream& outFile,
  mesh::Mesh&    mesh)
{
  TRACE(mesh.getName());

  // Plot vertices
  outFile << "POINTS " << mesh.vertices().size() << " float \n\n";
  for (mesh::Vertex& vertex : mesh.vertices()) {
    writeVertex(vertex.getCoords(), outFile);
  }
  outFile << '\n';


  // Plot edges
  if(mesh.getDimensions() == 2) {
    outFile << "CELLS " << mesh.edges().size() << ' ' << mesh.edges().size() * 3
            << '\n' << '\n';
    for (mesh::Edge & edge : mesh.edges()) {
      int internalIndices[2];
      internalIndices[0] = edge.vertex(0).getID();
      internalIndices[1] = edge.vertex(1).getID();
      writeLine(internalIndices, outFile);
    }
    outFile << "\nCELL_TYPES " << mesh.edges().size() << "\n\n";
    for(size_t i = 0; i < mesh.edges().size(); ++i) {
      outFile << "3\n";
    }
  }

  // Plot triangles
  if(mesh.getDimensions() == 3) {
    size_t sizeTriangles = mesh.triangles().size();
    size_t sizeQuads = mesh.quads().size();
    outFile << "CELLS " << sizeTriangles + sizeQuads << ' '
            << sizeTriangles * 4 + sizeQuads * 5 << "\n\n";
    for (mesh::Triangle& triangle : mesh.triangles()) {
      int internalIndices[3];
      internalIndices[0] = triangle.vertex(0).getID();
      internalIndices[1] = triangle.vertex(1).getID();
      internalIndices[2] = triangle.vertex(2).getID();
      writeTriangle(internalIndices, outFile);
    }
    for (mesh::Quad& quad : mesh.quads()) {
      int internalIndices[4];
      internalIndices[0] = quad.vertex(0).getID();
      internalIndices[1] = quad.vertex(1).getID();
      internalIndices[2] = quad.vertex(2).getID();
      internalIndices[3] = quad.vertex(3).getID();
      writeQuadrangle(internalIndices, outFile);
    }

    outFile << "\nCELL_TYPES " << sizeTriangles + sizeQuads << "\n\n";
    for(size_t i=0; i < sizeTriangles; i++){
      outFile << "5\n";
    }
    for(size_t i=0; i < sizeQuads; i++){
      outFile << "9\n";
    }


    // OLD
//    outFile << "CELLS " << mesh.triangles().size() << ' '
//            << mesh.triangles().size() * 4  << "\n\n";
//    for (mesh::Triangle & triangle : mesh.triangles()){
//      int internalIndices[3];
//      internalIndices[0] = triangle.vertex(0).getID();
//      internalIndices[1] = triangle.vertex(1).getID();
//      internalIndices[2] = triangle.vertex(2).getID();
//      writeTriangle(internalIndices, outFile);
//    }
//    outFile << "\nCELL_TYPES " << mesh.triangles().size() << "\n\n";
//    for(size_t i=0; i < mesh.triangles().size(); i++){
//      outFile << "5\n";
//    }

  }

  outFile << '\n';
}

void ExportVTK:: exportData
(
  std::ofstream& outFile,
  mesh::Mesh&    mesh)
{
  outFile << "POINT_DATA " << mesh.vertices().size() << '\n';
  outFile << '\n';

  if(_writeNormals) { // Plot vertex normals
    outFile << "VECTORS VertexNormals float\n";
    outFile << '\n';
    for (mesh::Vertex& vertex : mesh.vertices()) {
      int i = 0;
      for(; i < mesh.getDimensions(); i++){
        outFile << vertex.getNormal()[i] << ' ';
      }
      if(i < 3){
        outFile << '0';
      }
      outFile << '\n';
    }
    outFile << '\n';

    // Plot edge normals
//    if(_plotNormals) {
//      VTKTextFileWriter vtkWriterEdgeNormals; //(fileName + "-edgenormals.vtk");
//      PtrVertexWriter normalsOriginWriter (vtkWriterEdgeNormals.createVertexWriter());
//      PtrVertexDataWriter normalsWriter (
//          vtkWriterEdgeNormals.createVertexDataWriter("EdgeNormals", utils::Def::DIM));
//      for (mesh::Edge & edge : mesh.edges()) {
//  //      int vertexID = normalsOriginWriter.getNextFreeVertexNumber ();
//        int vertexID = normalsOriginWriter->plotVertex(edge.getCenter());
//        normalsWriter->plotVertex(vertexID, edge.getNormal());
//      }
//      vtkWriterEdgeNormals.writeToFile(fileName + "-edgenormals.vtk");
//  //    vtkWriterEdgeNormals.plotVertices(normalsOriginWriter);
//  //    vtkWriterEdgeNormals.plotPointData(normalsWriter);
//    }
  }

  for (mesh::PtrData data : mesh.data()) { // Plot vertex data
    Eigen::VectorXd& values = data->values();
    if(data->getDimensions() > 1) {
      Eigen::VectorXd viewTemp(data->getDimensions());
      outFile << "VECTORS " << data->getName() << " float\n";
      for (mesh::Vertex& vertex : mesh.vertices()) {
        int offset = vertex.getID() * data->getDimensions();
        for(int i=0; i < data->getDimensions(); i++){
          viewTemp[i] = values(offset + i);
        }
        int i = 0;
        for(; i < data->getDimensions(); i++){
          outFile << viewTemp[i] << ' ';
        }
        if(i < 3){
          outFile << '0';
        }
        outFile << '\n';
      }
      outFile << '\n';
    }
    else if(data->getDimensions() == 1) {
      outFile << "SCALARS " << data->getName() << " float\n";
      outFile << "LOOKUP_TABLE default\n";
      for (mesh::Vertex& vertex : mesh.vertices()) {
        outFile << values(vertex.getID()) << '\n';
      }
      outFile << '\n';
    }
  }
}

void ExportVTK:: initializeWriting
(
  std::ofstream&     filestream)
{
  //size_t pos = fullFilename.rfind(".vtk");
  //if ((pos == std::string::npos) || (pos != fullFilename.size()-4)){
  //  fullFilename += ".vtk";
  //}
  filestream.setf(std::ios::showpoint);
  filestream.setf(std::ios::scientific);
  filestream << std::setprecision(16);
}

void ExportVTK:: writeHeader
(
  std::ostream& outFile)
{
  outFile << "# vtk DataFile Version 2.0\n\n"
          << "ASCII\n\n"
          << "DATASET UNSTRUCTURED_GRID\n\n";
}

void ExportVTK:: writeVertex
(
  const Eigen::VectorXd&  position,
  std::ostream&           outFile)
{
  if(position.size() == 2) {
    outFile << position(0) << "  " << position(1) << "  " << 0.0 << '\n';
  }
  else {
    assertion(position.size() == 3);
    outFile << position(0) << "  " << position(1) << "  " << position(2) << '\n';
  }
}


void ExportVTK:: writeTriangle
(
  int           vertexIndices[3],
  std::ostream& outFile)
{
  outFile << 3 << ' ';
  for(int i=0; i < 3; i++) {
    outFile << vertexIndices[i] << ' ';
  }
  outFile << '\n';
}

void ExportVTK:: writeQuadrangle
(
  int           vertexIndices[4],
  std::ostream& outFile)
{
  outFile << 4 << ' ';
  for(int i=0; i < 4; i++) {
    outFile << vertexIndices[i] << ' ';
  }
  outFile << '\n';
}

void ExportVTK:: writeLine
(
  int           vertexIndices[2],
  std::ostream& outFile)
{
  outFile << 2 << ' ';
  for(int i=0; i<2; i++) {
    outFile << vertexIndices[i] << ' ';
  }
  outFile << '\n';
}

}} // namespace precice, io

