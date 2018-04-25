#include "ExportVRML.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/PropertyContainer.hpp"
#include "utils/String.hpp"
#include <iostream>
#include <Eigen/Core>
#include <fstream>
#include <map>
#include <boost/filesystem.hpp>
#include "utils/Helpers.hpp"

namespace precice {
namespace io {

logging::Logger ExportVRML:: _log ("io::ExportVRML");

ExportVRML:: ExportVRML
(
  bool plotNormals )
:
  Export()
{}

int ExportVRML:: getType() const
{
  return constants::exportVRML();
}

void ExportVRML:: doExport
(
  const std::string& name,
  const std::string& location,
  mesh::Mesh&        mesh )
{
  namespace fs = boost::filesystem;
  fs::path outfile(location);
  outfile = outfile / fs::path(name);
  std::ofstream outstream(outfile.string(), std::ios::trunc);
  CHECK(outstream, "Could not open file \"" << outfile.c_str() << "\" for VTK export!");

  writeHeader(outstream);
  writeMesh(outstream, mesh);
  outstream << "}" << std::endl;
  outstream.close();
}

void ExportVRML:: doExportCheckpoint
(
  const std::string& filename,
  mesh::Mesh&        mesh )
{
  //currently checkpoints do not have a dedicated location, but write straight to the main directory
  std::ofstream outFile;
  std::string fullFilename(filename);
  openFile(outFile, utils::checkAppendExtension(fullFilename, std::string(".wrl")));
  writeHeader(outFile);
  writeMesh(outFile, mesh);
  writeVertexData(outFile, mesh);
  writePropertyContainer(outFile, mesh);
  outFile << "}" << std::endl;
  outFile.close();
}

void ExportVRML:: openFile
(
  std::ofstream&     outFile,
  const std::string& filename ) const
{
  outFile.open ( filename.c_str() );
  CHECK(outFile, "Could not open file \"" << filename << "\" for VRML export!" );
  outFile.setf(std::ios::showpoint);
  outFile.setf(std::ios::scientific);
  outFile << std::setprecision(16);
}

void ExportVRML:: writeHeader
(
  std::ofstream& outFile ) const
{
  outFile  << "#VRML V1.0 ascii"  << std::endl
           << "Separator {"       << std::endl
           << "   Coordinate3 {"  << std::endl
           << "      point ["     << std::endl;
}

void ExportVRML:: writeMesh
(
  std::ofstream& outFile,
  mesh::Mesh&    mesh ) const
{
  std::map<int,int> vertexIndices; // first: precice id, second: vrml file id
  int dimensions = mesh.getDimensions();
  assertion ( (dimensions == 2) || (dimensions == 3), dimensions );

  // Export vertices
  int vertexIndex = 0;
  for (mesh::Vertex& vertex : mesh.vertices()) {
    if ( dimensions == 2 ) {
      outFile << "         " << vertex.getCoords()[0]
              << " " << vertex.getCoords()[1]
              << " " << "0," << std::endl;
    }
    else {
      outFile << "         " << vertex.getCoords()[0]
              << " "      << vertex.getCoords()[1]
              << " "      << vertex.getCoords()[2]
              << ","      << std::endl;
    }
    assertion ( vertexIndices.count(vertex.getID()) == 0 );
    vertexIndices[vertex.getID()] = vertexIndex;
    vertexIndex ++;
  }

  outFile << " ]"                   << std::endl
          << "   }"                 << std::endl
          << std::endl;

  // Export faces (edges or triangles)
  if ( dimensions == 2 ) {
    outFile << "   IndexedLineSet  {" << std::endl
      << "      coordIndex  [ " << std::endl;
    for (mesh::Edge& edge : mesh.edges()) {
      assertion ( vertexIndices.count(edge.vertex(0).getID()) > 0 );
      assertion ( vertexIndices.count(edge.vertex(1).getID()) > 0 );
      outFile << "         "
              << vertexIndices.find(edge.vertex(0).getID())->second
              << ", "
              << vertexIndices.find(edge.vertex(1).getID())->second
              << ", -1,"
              << std::endl;
    }
  }
  else {
    assertion ( dimensions == 3 );
    outFile << "   IndexedFaceSet  {" << std::endl
            << "      coordIndex  [ " << std::endl;
    for (mesh::Triangle& triangle : mesh.triangles()) {
      outFile << "         ";
      for ( int i=0; i < 3; i++ ) {
        int vertexIndex = triangle.vertex(i).getID();
        assertion ( utils::contained(vertexIndex,vertexIndices) );
        outFile << vertexIndices[vertexIndex] << ", ";
      }
      outFile << "-1," << std::endl;
    }

//   int indexEdge0Vertex0 = vertexIndices[ mesh.getTriangle ( i ).edge ( 0 ).vertex (0 ).getID() ];
//   int indexEdge0Vertex1 = vertexIndices[ mesh.getTriangle ( i ).edge ( 0 ).vertex (1 ).getID() ];
//
//   int indexEdge1Vertex0 = vertexIndices [ mesh.getTriangle ( i ).edge ( 1 ).vertex (0 ).getID() ];
//   int indexEdge1Vertex1 = vertexIndices [ mesh.getTriangle ( i ).edge ( 1 ).vertex (1 ).getID() ];
//
//   // ensure clockweise
//   if ( indexEdge0Vertex0 == indexEdge1Vertex0 ) {
//      outFile << "   " << indexEdge0Vertex1
//              << ", " << indexEdge0Vertex0
//              << ", " << indexEdge1Vertex1
//              << ", -1," << endl;
//   }
//   else if ( indexEdge0Vertex0 == indexEdge1Vertex1 ) {
//      outFile << "   " << indexEdge0Vertex1
//              << ", " << indexEdge0Vertex0
//              << ", " << indexEdge1Vertex0
//              << ", -1," << endl;
//   }
//   else if ( indexEdge0Vertex1 == indexEdge1Vertex0 ) {
//      outFile << "   " << indexEdge0Vertex0
//              << ", " << indexEdge0Vertex1
//              << ", " << indexEdge1Vertex1
//              << ", -1," << endl;
//   }
//   else if ( indexEdge0Vertex1 == indexEdge1Vertex1 ) {
//      outFile << "   " << indexEdge0Vertex0
//              << ", " << indexEdge0Vertex1
//              << ", " << indexEdge1Vertex0
//              << ", -1," << endl;
//   }
   }
   outFile << "      ]" << std::endl
           << "   }"    << std::endl;
}

void ExportVRML:: writeVertexData
(
  std::ofstream& outFile,
   mesh::Mesh&   mesh ) const
{
  for (mesh::PtrData data : mesh.data()) {
    outFile << "   Info {"                         << std::endl
            << "      string \"preCICE data\""     << std::endl
            << "      fields [ SFString dataname SFString datatype "
            << "MFFloat datavalues ]"              << std::endl
            << "      dataname \"" << data->getName() << "\"" << std::endl
            << "      datadimensions " << data->getDimensions() << std::endl
            << "      datavalues ["                << std::endl;

    const Eigen::VectorXd& values = data->values();
    int dims = data->getDimensions();
    for ( int i=0; i < values.size(); i+=dims ) {
      outFile << "         ";
      for ( int dim=0; dim < dims; dim++ ) {
        outFile << values(i+dim);
        if ( dim + 1 < dims ) {
          outFile << " ";
        }
      }
      outFile << "," << std::endl;
    }

//      for ( mesh::Vertex & vertex : mesh.vertices() ) {
//         assertion ( vertex.hasProperty(data->getID()) );
//         if ( data->getType() == mesh::Data::TYPE_VECTOR ) {
//            using utils::Vector;
//            Vector value ( vertex.getProperty<Vector>(data->getID()) );
//            if ( utils::Def::DIM == 2 ) {
//               outFile << "         " << value[0] << " "
//                       << value[1] << "," << std::endl;
//            }
//            else {
//               assertion ( utils::Def::DIM == 3 );
//               outFile << "         " << value[0] << " "
//                       << value[1] << " " << value[2] << "," << std::endl;
//            }
//         }
//         else {
//            preciceCheck ( false, "writeVertexData()",
//                           "Only vector data is supported!" );
//         }
//      }

      outFile << "      ]" << std::endl
              << "   }"    << std::endl;

   }
}

void ExportVRML:: writePropertyContainer
(
  std::ofstream& outFile,
  mesh::Mesh&    mesh ) const
{
  using mesh::PropertyContainer;
  int dimensions = mesh.getDimensions();
  for (PropertyContainer & container : mesh.propertyContainers()) {
    int id = container.getProperty<int> ( container.INDEX_GEOMETRY_ID );
    std::string idName ( "" );
    for (auto &nameID : mesh.getNameIDPairs()) {
      if ( nameID.second == id ) {
        idName = nameID.first;
        break;
      }
    }
    assertion ( idName != std::string("") );

    outFile << "   Info {"                                << std::endl
      << "      string \"preCICE property container\""    << std::endl
      << "      fields [ SFString idname MFInt32 faces ]" << std::endl
      << "      idname \"" << idName << "\"" << std::endl
      << "      faces [ "                    << std::endl;

//    using mesh::Vertex;
//    int vertexIndex = 0;
//    for ( Vertex & vertex : mesh.vertices() ) {
//      for ( int i=0; i < vertex.getParentCount(); i++ ) {
//        if ( &vertex.getParent(i) == &container  ) {
//          outFile << "         " << vertexIndex << "," << std::endl;
//        }
//      }
//      vertexIndex ++;
//    }
//    outFile << "      ]" << std::endl;

    int index = 0;
    if (dimensions == 2){
      for (mesh::Edge& edge : mesh.edges()) {
        for (int i=0; i < edge.getParentCount(); i++){
          if (&edge.getParent(i) == &container ){
            outFile << "         " << index << "," << std::endl;
          }
        }
        index++;
      }
    }
    else {
      for (mesh::Triangle& triangle : mesh.triangles()) {
        for (int i=0; i < triangle.getParentCount(); i++){
          if (&triangle.getParent(i) == &container ){
            outFile << "         " << index << "," << std::endl;
          }
        }
        index++;
      }

    }
    outFile << "      ]" << std::endl;
    outFile << "   }"    << std::endl;
  }
}

}} // namespace precice, io


