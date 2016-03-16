// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ImportVRML.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/PropertyContainer.hpp"
#include "utils/Dimensions.hpp"
#include "utils/String.hpp"
#include "tarch/la/WrappedVector.h"
#ifndef PRECICE_NO_SPIRIT2
#include <boost/spirit/include/qi_parse.hpp>
#include <boost/spirit/include/support_multi_pass.hpp>
namespace spirit = boost::spirit;
#endif // not PRECICE_NO_SPIRIT2
#include <fstream>
#include <vector>
#include <list>

namespace precice {
namespace io {

logging::Logger ImportVRML:: _log ( "precice::io::ImportVRML" );

ImportVRML:: ImportVRML
(
  const std::string& location )
:
  Import(location),
  _createMesh(true)
{
# ifdef PRECICE_NO_SPIRIT2
  preciceError("ImportVRML()",
               "VRML import can only be used with Boost.Spirit V2.0!");
# endif
}

void ImportVRML:: doImport
(
  const std::string& name,
  mesh::Mesh&        mesh )
{
  doImport(name, mesh, false);
}

void ImportVRML:: doImportCheckpoint
(
  const std::string& name,
  mesh::Mesh&        mesh,
  bool               createMesh)
{
  _createMesh = createMesh;
  doImport(name, mesh, true);
  _createMesh = true;
}

void ImportVRML:: doImport
(
  const std::string& name,
  mesh::Mesh&        mesh,
  bool               isCheckpoint )
{
# ifndef PRECICE_NO_SPIRIT2
  using namespace tarch::la;
  int dimensions = mesh.getDimensions();
  assertion1((dimensions == 2) || (dimensions == 3), dimensions);

  // Create input filestream to VRML file
  std::string filename(getLocation() + name);
  std::ifstream in(utils::checkAppendExtension(filename, ".wrl").c_str());
  preciceCheck(in.is_open(), "doImport()",
               "Could not open input file " << filename << "!");

  // Wrap input file stream into a multi pass iterator (requirement)
  typedef std::istreambuf_iterator<char> FileIter;
  typedef spirit::multi_pass<FileIter> MultiPassFileIter;
  MultiPassFileIter first = spirit::make_default_multi_pass(FileIter(in));
  MultiPassFileIter last = spirit::make_default_multi_pass(FileIter());

  // Parse VRML file and check validity
  impl::VRML10Parser<MultiPassFileIter> vrmlParser(dimensions);
  bool success = spirit::qi::phrase_parse (
      first, last, vrmlParser, spirit::qi::space );
  preciceCheck(success && (first == last), "doImport()",
               "Parsing of file " << filename << " failed! Left over: "
               << std::endl << std::string(first, last));

  // Construct vertex coordinates from parsed information
  if(_createMesh){
    std::vector<mesh::Vertex*> vertices;
    if (dimensions == 2){
      for (size_t i=0; i < vrmlParser.coordinates.size(); i+=2){
        assertion2(i + 1 < vrmlParser.coordinates.size(),
                   i + 1, vrmlParser.coordinates.size());
        vertices += &mesh.createVertex(wrap<2,double>(&vrmlParser.coordinates[i]));
      }
      // Construct edge indices from parsed data.
      // The parsed data has the form: i0, ..., in, -1, i0, ..., im, -1, ...., -1
      //assertion(vrmlParser.indices.size() > 2);
      std::vector<tarch::la::Vector<2,int> > indices;
      for (size_t i=0; i < vrmlParser.indices.size(); i++){
        if (vrmlParser.indices[i+1] == -1){
          i++;
        }
        else {
          assertion(vertices[vrmlParser.indices[i]] != nullptr);
          assertion(vertices[vrmlParser.indices[i+1]] != nullptr);
          mesh.createEdge(*vertices[vrmlParser.indices[i]],
                          *vertices[vrmlParser.indices[i+1]]);
        }
      }
    }
    else { // 3D
      // Create vertices
      for (size_t i=0; i < vrmlParser.coordinates.size(); i+=3){
        assertion2(i + 2 < vrmlParser.coordinates.size(),
                   i + 2, vrmlParser.coordinates.size());
        vertices += &mesh.createVertex(wrap<3,double>(&vrmlParser.coordinates[i]));
      }

      // Construct triangle indices from parsed data.
      std::vector<std::list<mesh::Edge*> > adjacencyList(vertices.size());
      for (size_t i=0; i < vrmlParser.indices.size(); i += 3){
        assertion(i + 2 < vrmlParser.indices.size());
        mesh::Vertex* v0 = vertices[vrmlParser.indices[i]];
        mesh::Vertex* v1 = vertices[vrmlParser.indices[i+1]];
        mesh::Vertex* v2 = vertices[vrmlParser.indices[i+2]];
        mesh::Edge& edge0 = getEdge(*v0, *v1, mesh, adjacencyList);
        mesh::Edge& edge1 = getEdge(*v1, *v2, mesh, adjacencyList);
        mesh::Edge& edge2 = getEdge(*v2, *v0, mesh, adjacencyList);
        mesh.createTriangle(edge0, edge1, edge2);
      }
    }
  }
  else{
    //for provided meshes, no vertices need to be created,
    //but we check if the once given by the user coincide with the once read from file

    preciceCheck(vrmlParser.coordinates.size()/dimensions == mesh.vertices().size(),
        "doImport()",
        "For the mesh " << mesh.getName() << ", " << mesh.vertices().size()
        << " vertices were set, while " << vrmlParser.coordinates.size()/dimensions
        << " vertices are read from file for restart.");
    if (dimensions == 2){
      for (size_t i=0; i < mesh.vertices().size(); i++){
        preciceCheck((mesh.vertices()[i].getCoords()[0] == vrmlParser.coordinates[i*2])
             && (mesh.vertices()[i].getCoords()[1] == vrmlParser.coordinates[i*2+1])                                                                ,
            "doImport()","For mesh " << mesh.getName() << " the vertices that were set"
            << " do not coincide with those read from file.");
      }
    }
    else { // 3D
      for (size_t i=0; i < mesh.vertices().size(); i++){
        preciceCheck((mesh.vertices()[i].getCoords()[0] == vrmlParser.coordinates[i*3])
             && (mesh.vertices()[i].getCoords()[1] == vrmlParser.coordinates[i*3+1])
             && (mesh.vertices()[i].getCoords()[2] == vrmlParser.coordinates[i*3+2])   ,
            "doImport()","For mesh " << mesh.getName() << " the vertices that were set"
            << " do not coincide with those read from file.");
      }
    }
  }

   // Stop here, if no checkpoint data is to be read
   if (not isCheckpoint){
      return;
   }

   // Set data values of vertices
   if (mesh.data().size() != vrmlParser.data.size()){
      preciceError("doImport()", "Number of vertex data sets configured does "
                   << "not match with that in VRML file!");
   }
   mesh.allocateDataValues();
   for (size_t iData=0; iData < vrmlParser.data.size(); iData++){
      mesh::PtrData meshData = mesh.data()[iData];
      const impl::VRML10Parser<MultiPassFileIter>::Data&
          vrmlData = vrmlParser.data[iData];
      if (meshData->getName() != vrmlData.name){
         preciceError("doImport()", "Name of configured data set " << iData
                      << " \"" << meshData->getName()
                      << "\" does not match with that in VRML file "
                      << "\"" << vrmlData.name << "\"!");
      }
      if ( meshData->getDimensions() != vrmlData.dimensions ) {
         preciceError ( "doImport()", "Dimension of configured data set " << iData
                         << " \"" << meshData->getDimensions()
                         << "\" does not match with that in VRML file "
                         << "\"" << vrmlData.dimensions << "\"!" );
      }
      Eigen::VectorXd& values = meshData->values();
      preciceCheck(values.size() == (int)vrmlData.values.size(),
                   "doImport()", "Number of data values from VRML file ("
                   << vrmlData.values.size()
                   << ") does not fit to number of expected data values ("
                   << values.size() << ") for data \""
                   << meshData->getName() << "\"!");
      for ( size_t i=0; i < vrmlData.values.size(); i++ ) {
         values(i) = vrmlData.values[i];
      }
   }

   // Add property containers
   preciceCheck ( vrmlParser.propertyContainers.size()
                  == mesh.propertyContainers().size(), "doImport()",
                  "Number of property containers in VRML file does not equals "
                  << "that of mesh from configuration!" );
   for (auto & elem : vrmlParser.propertyContainers) {
      std::string subIDName ( elem.subIDName );
      mesh::PropertyContainer& cont = mesh.getPropertyContainer(subIDName);
      if (dimensions == 2){
        mesh::Mesh::EdgeContainer& edges = mesh.edges();
        size_t size = elem.faces.size();
        for (size_t iFace=0; iFace < size; iFace++){
          int index = elem.faces[iFace];
          mesh::Edge& edge = edges[index];
          edge.addParent(cont);
          addParentIfMissing(edge.vertex(0), cont);
          addParentIfMissing(edge.vertex(1), cont);
        }
      }
      else {
        assertionMsg(dimensions == 3, dimensions);
        mesh::Mesh::TriangleContainer& triangles = mesh.triangles();
        size_t size = elem.faces.size();
        for (size_t iFace=0; iFace < size; iFace++){
          int index = elem.faces[iFace];
          mesh::Triangle& triangle = triangles[index];
          triangle.addParent(cont);
          addParentIfMissing(triangle.edge(0), cont);
          addParentIfMissing(triangle.edge(1), cont);
          addParentIfMissing(triangle.edge(2), cont);
          addParentIfMissing(triangle.vertex(0), cont);
          addParentIfMissing(triangle.vertex(1), cont);
          addParentIfMissing(triangle.vertex(2), cont);
        }
      }
   }
#  endif // not PRECICE_NO_SPIRIT2
}

mesh::Edge& ImportVRML:: getEdge
(
  mesh::Vertex& vertexOne,
  mesh::Vertex& vertexTwo,
  mesh::Mesh&   mesh,
  std::vector<std::list<mesh::Edge*> >& adjacencyList )
{
  // Edge might be created already
  std::list<mesh::Edge*>& adjEdgesOne = adjacencyList[vertexOne.getID()];
  for (mesh::Edge* edge : adjEdgesOne){
    assertion(edge != nullptr);
    if (   (edge->vertex(0).getID() == vertexTwo.getID())
        || (edge->vertex(1).getID() == vertexTwo.getID()))
    {
      return *edge;
    }
  }
  std::list<mesh::Edge*>& adjEdgesTwo = adjacencyList[vertexTwo.getID()];
  for (mesh::Edge* edge : adjEdgesTwo){
    assertion(edge != nullptr);
    if (   (edge->vertex(0).getID() == vertexOne.getID())
        || (edge->vertex(1).getID() == vertexOne.getID()))
    {
      return *edge;
    }
  }
  // Edge does not exist yet, create and add it to adjacency list
  mesh::Edge& edge = mesh.createEdge(vertexOne, vertexTwo);
  adjacencyList[vertexOne.getID()].push_front(&edge);
  adjacencyList[vertexTwo.getID()].push_front(&edge);
  return edge;
}

void ImportVRML:: addParentIfMissing
(
  mesh::PropertyContainer& object,
  mesh::PropertyContainer& parent )
{
  for (int i=0; i < object.getParentCount(); i++){
    if (&object.getParent(i) == &parent){
      return;
    }
  }
  object.addParent(parent);
}


}} // namespace precice, io

