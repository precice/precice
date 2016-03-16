// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IO_IMPORTVRML_HPP_
#define PRECICE_IO_IMPORTVRML_HPP_

#include "Import.hpp"
#include "utils/Dimensions.hpp"
#include "impl/VRML10Parser.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Helpers.hpp"
#include "logging/Logger.hpp"
#include <map>

namespace precice {
  namespace mesh {
    class Mesh;
    class Edge;
    class Vertex;
    class PropertyContainer;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace io {

/**
 * @brief Imports a geometry from a vrml file.
 */
class ImportVRML : public Import
{
public:

   /**
    * @brief Constructor.
    */
   ImportVRML ( const std::string& location );


   /**
    * @brief Destructor.
    */
   virtual ~ImportVRML() {}

   /**
    * @brief Imports the geometry from an vrml file into a Mesh object.
    *
    * @param mesh [IN/OUT] The imported elements are added to this mesh.
    */
   virtual void doImport (
      const std::string& name,
      mesh::Mesh&        mesh );

   void doImportCheckpoint (
      const std::string& name,
      mesh::Mesh&        mesh,
      bool               createMesh);

private:

   // @brief Logging device.
   static logging::Logger _log;

   bool _createMesh;

   void doImport (
      const std::string& name,
      mesh::Mesh&        mesh,
      bool                isCheckpoint );

   mesh::Edge& getEdge (
     mesh::Vertex& vertexOne,
     mesh::Vertex& vertexTwo,
     mesh::Mesh&   mesh,
     std::vector<std::list<mesh::Edge*> >& adjacencyList );

   void addParentIfMissing(
     mesh::PropertyContainer& object,
     mesh::PropertyContainer& parent );

//   int getIndexExistingEdge (
//      const std::vector<mesh::Edge*>& edges,
//      const mesh::Vertex*             vertex0,
//      const mesh::Vertex*             vertex1 ) const;

};

}} // namespace precice, io

#endif /* PRECICE_IO_IMPORTVRML_HPP_ */
