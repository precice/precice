// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef GEOMETRY_CONFIG_GEOMETRYCONFIGURATION_HPP_
#define GEOMETRY_CONFIG_GEOMETRYCONFIGURATION_HPP_

#include "geometry/Geometry.hpp"
#include "geometry/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "com/Communication.hpp"
#include "utils/xml/XMLTag.hpp"
#include "tarch/logging/Log.h"
#include <vector>
#include <set>
#include "boost/tuple/tuple.hpp"
#include <string>

namespace precice {
  namespace mesh {
    class DataConfiguration;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace geometry {

class GeometryConfiguration : public utils::XMLTag::Listener
{
public:

   // @brief Name of the xml-tag corresponding to the GeometryConfiguration.
   //static const std::string& getTag();

   /**
    * @brief Constructor
    */
   GeometryConfiguration (
     utils::XMLTag&             parent,
     mesh::PtrMeshConfiguration meshConfig );

   /**
    * @brief Sets the spatial dimensions of the geometries to be configured.
    */
   void setDimensions ( int dimensions );

   /**
    * Parse a TAG-tag. The parser reads the tag, validates the context, ensures
    * the corresponding compiler switch is set or not (depending on tag) and
    * returns. The argument may not be 0. If the either the validation or the
    * compiler switch check fails, all successing isValid() calls fail.
    */
   //bool parseSubtag ( utils::XMLTag::XMLReader* xmlReader );

   /**
    * @returns Whether all compiler switches checked using the parseSubtag()
    *          routine were valid compared to the values specified in the
    *          configuration file.
    */
   //bool isValid() const;

   std::vector<PtrGeometry>& geometries()
   {
      return _geometries;
   }

   /**
    * @brief Returns the geometry of the given mesh, or a NULL initialized Ptr.
    */
   PtrGeometry getGeometry ( const std::string& meshName );

   /**
    * @brief Callback from automated xml-configuration for opening tag.
    */
   virtual void xmlTagCallback ( utils::XMLTag& callingTag);

   /**
    * @brief Callback from automated xml-configuration for closing tag.
    */
   virtual void xmlEndTagCallback ( utils::XMLTag& callingTag );

   /**
    * @brief For manual configuration, as done in test-cases
    */
   void addGeometry (
      PtrGeometry geometry,
      const std::string&  meshName );

private:

   // @brief Log device
   static tarch::logging::Log _log;

   // @brief XML tag corresponding to class
   const std::string TAG;
   const std::string TAG_LENGTH;
   const std::string TAG_DISCRETIZATION_WIDTH;
   const std::string TAG_OFFSET;
   const std::string TAG_PORES;
   const std::string TAG_RADIUS;
   const std::string TAG_DEFORMATION;
   const std::string TAG_FILENAME;
   const std::string TAG_FILETYPE;
   const std::string TAG_PROVIDER;
   const std::string TAG_RECEIVER;

   const std::string ATTR_MESH;
   const std::string ATTR_TYPE;
   const std::string ATTR_VALUE;
   const std::string ATTR_NAME;

   const std::string VALUE_BUILTIN_CUBOID;
   const std::string VALUE_BUILTIN_SPHERE;
   const std::string VALUE_BUILTIN_BUBBLE;
   const std::string VALUE_BUILTIN_DRATCHET;
   const std::string VALUE_BUILTIN_FACE;
   const std::string VALUE_IMPORT;
   const std::string VALUE_NONE;
   const std::string VALUE_AUTO;

   mesh::PtrMeshConfiguration _meshConfig;

   int _dimensions;

   /**
    * @brief Stores already read data, to be used when creating a geometry
    */
   struct ReadData
   {
      ReadData (int dimensions)
      :
         noReadData(true),
         type(""),
         mesh(""),
         discretizationWidth(0.0),
         offset(utils::DynVector(dimensions, 0.0)),
         radius(0.0),
         deformation(0.0),
         scalarLength(0.0),
         length(utils::DynVector(dimensions, 0.0)),
         pores(0.0),
         filename(""),
         filetype("")
      {}

      ReadData& operator= (const ReadData& toCopy)
      {
        noReadData = toCopy.noReadData;
        type = toCopy.type;
        mesh = toCopy.mesh;
        discretizationWidth = toCopy.discretizationWidth;
        offset.clear();
        offset.append(toCopy.offset);
        radius = toCopy.radius;
        deformation = toCopy.deformation;
        scalarLength = toCopy.scalarLength;
        length.clear();
        length.append(toCopy.offset);
        pores = toCopy.pores;
        filename = toCopy.filename;
        filetype = toCopy.filetype;
        return *this;
      }

      bool noReadData; // To note, whether any data has been read yet
      std::string type; // Type of geometry
      std::string mesh; // Name of mesh to build with geometry.
      double discretizationWidth;
      utils::DynVector offset; // Translational displacement of the geometry
      double radius; // Radius for sphere or drift ratchet or bubble
      double deformation; // Deformation for bubble
      double scalarLength; // Length of drift ratchet
      utils::DynVector length; // Sidelengths of cuboid
      double pores; // Number of pores for drift ratchet
      std::string filename; // Filename of geometry to import
      std::string filetype;
   };

   ReadData _readData;

   // @brief Geometries configured (besides custom geometries) and meshes.
   std::vector<PtrGeometry> _geometries;

   std::vector<std::string> _meshNames;

   // @brief Holds true if the filter configuration is valid
   //bool _isValid;

   /**
    * @brief Adds geometry type specific XML attributes.
    */
   void addTypeSpecificAttributes (
      const std::string& type,
      utils::XMLTag&     tag );

   bool addCuboid ();

   bool addDriftRatchet ();

   void addSphere ();

   void addBubble ();

   void addImportGeometry ();

   void checkMeshName ( const std::string& meshName );
};

}} // namespace precice, geoemtry

#endif /*GEOMETRY_CONFIG_GEOMETRYCONFIGURATION_HPP_*/
