// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IO_IMPL_VRML10PARSER_HPP_
#define PRECICE_IO_IMPL_VRML10PARSER_HPP_

#ifndef PRECICE_NO_SPIRIT2

#include "mesh/PropertyContainer.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Helpers.hpp"
#include "tarch/logging/Log.h"
#include <boost/spirit/include/qi.hpp>
#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/include/phoenix_stl.hpp>
#include <boost/bind.hpp>
#include <iostream>
#include <fstream>
#include <vector>

namespace precice {
namespace io {
namespace impl {

namespace qi = boost::spirit::qi;
namespace phoenix = boost::phoenix;

#ifdef Debug
#define PRECICE_PARSER_DEBUG(message) >> qi::eps[boost::bind(&VRML10Parser::debug, this, message)]
#else
#define PRECICE_PARSER_DEBUG(message)
#endif

/**
 * @brief Defines grammer for parsing VRML 1.0 files.
 */
template< typename ITERATOR_T >
struct VRML10Parser : public qi::grammar<ITERATOR_T, qi::space_type>
{
   struct Data {
      std::string name;
      int dimensions;
      std::vector<double> values;
   };

   struct PropertyContainer {
     std::string subIDName;
     //std::vector<int> vertices;
     std::vector<int> faces; // Can be edge (2D) or triangle (3D) indices
   };

   // @brief Logging device.
   static tarch::logging::Log _log;

   // @brief The entry point for parsing.
   qi::rule<ITERATOR_T, qi::space_type> start;

   // @brief One valid vrml tag.
   qi::rule<ITERATOR_T, qi::space_type> ruleVRMLTag;

   // @brief Hierarchical tag for nesting further vrml tags.
   qi::rule<ITERATOR_T, qi::space_type> ruleSeparator;

   // @brief Hierarchical tag for nesting further vrml tags.
   qi::rule<ITERATOR_T, qi::space_type> ruleGroup;

   // @brief Hierarchical tag for nesting further vrml tags.
   qi::rule<ITERATOR_T, qi::space_type> ruleSwitch;

   // @brief
   qi::rule<ITERATOR_T, qi::space_type> ruleCoordinateSet;

   // @brief
   qi::rule<ITERATOR_T, qi::space_type> ruleCoordinates;

   // @brief
   qi::rule<ITERATOR_T, qi::space_type> ruleCoordinate;

   // @brief
   qi::rule<ITERATOR_T, qi::space_type> ruleIndexSet;

   // @brief
   qi::rule<ITERATOR_T, qi::space_type> ruleIndex;

   // @brief
   qi::rule<ITERATOR_T, qi::space_type> ruleData;

   qi::rule<ITERATOR_T, qi::space_type> rulePropertyContainer;

   // @brief
   qi::rule<ITERATOR_T, std::string(), qi::space_type> ruleQuotedString;

   // @brief
   qi::rule<ITERATOR_T, qi::space_type> ruleComment;

   // @brief
   qi::rule<ITERATOR_T, qi::space_type> ruleSkip;

   // @brief
   qi::rule<ITERATOR_T, qi::space_type> ruleDEF;

   // @brief
   qi::rule<ITERATOR_T, qi::space_type> ruleShapeHints;

   // @brief
   qi::rule<ITERATOR_T, qi::space_type> ruleMaterial;

   // @brief Coordinate components of vertices (V0C0, V0C1, V1C0, V1C2, ...)
   std::vector<double> coordinates;

   // @brief Indices of vertices of edges/triangles (T0I0, T0I1, T0I2, T1I0, ...)
   std::vector<int> indices;

   // @brief Data of vertices.
   std::vector<Data> data;

   std::vector<PropertyContainer> propertyContainers;

   // @brief Dimension of mesh/data to be parsed.
   int dimensions;

   /**
    * @brief Constructor.
    */
   VRML10Parser ( int dimensions );

   /**
    * @brief Adds name of a new vertex data set.
    */
   void addDataName ( const std::string & name );

   /**
    * @brief Adds type of a new vertex data set.
    */
   void addDataDimension ( int dimensions );

   /**
    * @brief Adds one value of a vertex data set.
    */
   void addDataValue ( double value );

   /**
    * @brief Sets the geoemtry ID of the property container.
    */
   void setPropertyContainerSubIDName ( const std::string& subIDName );

   /**
    * @brief Adds a vertex ID belonging to the property container.
    */
   void addPropertyContainerFace ( int faceID );

   /**
    * @brief Checks the validity of the parsed vrml file.
    */
   void checkInputValidity ();

   /**
    * @brief For debug output during the parsing process.
    */
   void debug ( const std::string& message );
};


// --------------------------------------------------------- HEADER DEFINITIONS


template< typename ITERATOR_T >
tarch::logging::Log VRML10Parser<ITERATOR_T>:: _log ( "precice::io::impl::VRML10Parser" );

template< typename ITERATOR_T >
VRML10Parser<ITERATOR_T>:: VRML10Parser
(
  int dimensions )
:
   VRML10Parser::base_type(start),
   dimensions ( dimensions )
{
   start =
      qi::lit("#VRML V1.0 ascii")
      >> *( ruleVRMLTag )
      >> qi::eps[boost::bind(&VRML10Parser::checkInputValidity, this)];

   ruleVRMLTag =
        ruleSeparator
      | ruleGroup
      | ruleSwitch
      | ruleCoordinateSet
      | ruleIndexSet
      | ruleData
      | rulePropertyContainer
      | ruleSkip;

   ruleSeparator =
      qi::lit("Separator")
      PRECICE_PARSER_DEBUG("Applying ruleSeparator")
      >> qi::char_('{')
      >> *ruleVRMLTag
      >> qi::char_('}')
      PRECICE_PARSER_DEBUG("Leaving ruleSeparator");

   ruleGroup =
      qi::lit("Group")
      PRECICE_PARSER_DEBUG("Applying ruleGroup")
      >> qi::char_('{')
      >> *ruleVRMLTag
      >> qi::char_('}')
      PRECICE_PARSER_DEBUG("Leaving ruleGroup");

   ruleSwitch =
      qi::lit("Switch")
      PRECICE_PARSER_DEBUG("Applying ruleSwitch")
      >> qi::char_('{')
      >> *ruleVRMLTag
      >> qi::char_('}')
      PRECICE_PARSER_DEBUG("Leaving ruleSwitch");

   ruleCoordinateSet =
      "Coordinate3"
      PRECICE_PARSER_DEBUG("Applying ruleCoordinateSet")
      >> qi::char_('{')
      >> "point"
      >> qi::char_('[')
      >> *(ruleCoordinates || ',')
      >> qi::char_(']')
      >> qi::char_('}');

   ruleData =
      "Info"
      >> qi::char_('{')
      >> qi::lit("string")
      >> qi::lit("\"preCICE data\"")
      PRECICE_PARSER_DEBUG("Applying ruleData")
      >> qi::lit("fields [ SFString dataname SFString datatype MFFloat datavalues ]")
      >> qi::eps[phoenix::push_back(phoenix::ref(data), Data())]  // Create new data entry
      >> qi::lit("dataname")
      >> ruleQuotedString[ boost::bind(&VRML10Parser::addDataName, this, ::_1) ]
      >> qi::lit("datadimensions")
      >> qi::int_[ boost::bind(&VRML10Parser::addDataDimension, this, ::_1) ]
      >> qi::lit("datavalues")
      >> qi::char_('[')
      >> *(
         qi::double_[ boost::bind(&VRML10Parser::addDataValue, this, ::_1) ]
         || qi::char_(',')
         )
      >> qi::char_(']')
      >> qi::char_('}')
      PRECICE_PARSER_DEBUG("Leaving ruleData");

   rulePropertyContainer =
      "Info"
      >> qi::char_("{")
      >> qi::lit("string")
      >> qi::lit("\"preCICE property container\"")
      PRECICE_PARSER_DEBUG("Applying rulePropertyContainer")
      >> qi::lit("fields [ SFString idname MFInt32 faces ]")
      >> qi::eps[phoenix::push_back(
                 phoenix::ref(propertyContainers), PropertyContainer())]  // Create new data entry
      >> qi::lit("idname")
      >> ruleQuotedString[ boost::bind(&VRML10Parser::setPropertyContainerSubIDName, this, ::_1) ]
      >> qi::lit("faces")
      >> qi::char_('[')
      >> *(
         qi::int_[ boost::bind(&VRML10Parser::addPropertyContainerFace, this, ::_1) ]
         || qi::char_(',')
         )
      >> qi::char_(']')
      >> qi::char_('}')
      PRECICE_PARSER_DEBUG("Leaving rulePropertyContainer");


   if ( dimensions == 2 ) {
      ruleCoordinates =
         ruleCoordinate
         >> ruleCoordinate
         >> qi::double_;

      ruleIndexSet =
         "IndexedLineSet"
         PRECICE_PARSER_DEBUG("Applying ruleIndexSet (2D)")
         >> qi::char_('{')
         >> "coordIndex"
         >> qi::char_('[')
         >> *( ruleIndex || ',' )
         >> qi::char_(']')
         >> qi::char_('}')
         PRECICE_PARSER_DEBUG("Leaving ruleIndexSet (2D)");
   }
   else {
      assertion ( dimensions == 3 );
      ruleCoordinates =
         ruleCoordinate
         >> ruleCoordinate
         >> ruleCoordinate;

      ruleIndexSet =
         "IndexedFaceSet"
         PRECICE_PARSER_DEBUG("Applying ruleIndexSet (3D)")
         >> qi::char_('{')
         >> "coordIndex"
         >> qi::char_('[')
         >> *(
            ruleIndex >> ','
            >> ruleIndex >> ','
            >> ruleIndex >> ','
            >> qi::lit("-1") || ','
            )
         >> qi::char_(']')
         >> qi::char_('}')
         PRECICE_PARSER_DEBUG("Leaving ruleIndexSet (3D)");
   }

   ruleCoordinate =
      qi::double_[phoenix::push_back(phoenix::ref(coordinates), qi::_1)];

   ruleIndex =
      qi::int_[phoenix::push_back(phoenix::ref(indices), qi::_1)];

   ruleQuotedString =
      qi::eps
      >> qi::lexeme[
         '"'
         >> +(qi::char_ - '"')[qi::_val += qi::_1]
         >> '"'
         ];

   ruleSkip =
        ruleComment
      | ruleDEF
      | ruleShapeHints
      | ruleMaterial;

   ruleComment =
      qi::char_('#')
      PRECICE_PARSER_DEBUG("Applying ruleComment")
      >> qi::lexeme[ *(qi::char_ - qi::eol) >> qi::eol ]
      PRECICE_PARSER_DEBUG("Leaving ruleComment");

   ruleDEF =
      qi::lit("DEF")
      PRECICE_PARSER_DEBUG("Applying ruleDEF")
      >> qi::lexeme[ *(qi::char_ - qi::eol) >> qi::eol ]
      PRECICE_PARSER_DEBUG("Leaving ruleDEF");

   ruleShapeHints =
      qi::lit("ShapeHints")
      PRECICE_PARSER_DEBUG("Applying ruleShapeHints")
      >> qi::char_('{')
      >> *(qi::char_ - qi::char_('}'))
      >> qi::char_('}')
      PRECICE_PARSER_DEBUG("Leaving ruleShapeHints");

   ruleMaterial =
      qi::lit("Material")
      PRECICE_PARSER_DEBUG("Applying ruleMaterial")
      >> qi::char_('{')
      >> *(qi::char_ - qi::char_('}'))
      >> qi::char_('}')
      PRECICE_PARSER_DEBUG("Leaving ruleMaterial");
}

template< typename ITERATOR_T >
void VRML10Parser<ITERATOR_T>:: addDataName ( const std::string & name )
{
  data.back().name = name;
}

template< typename ITERATOR_T >
void VRML10Parser<ITERATOR_T>:: addDataDimension ( int dimensions )
{
  data.back().dimensions = dimensions;
}

template< typename ITERATOR_T >
void VRML10Parser<ITERATOR_T>:: addDataValue ( double value )
{
  data.back().values.push_back ( value );
}

template< typename ITERATOR_T >
void VRML10Parser<ITERATOR_T>:: setPropertyContainerSubIDName ( const std::string & name )
{
   propertyContainers.back().subIDName = name;
}

template< typename ITERATOR_T >
void VRML10Parser<ITERATOR_T>:: addPropertyContainerFace ( int faceID )
{
   propertyContainers.back().faces.push_back ( faceID );
}

template< typename ITERATOR_T >
void VRML10Parser<ITERATOR_T>:: debug ( const std::string & message )
{
  std::string preciceMethodName ("debug()");
  preciceDebug ( message );
}

template< typename ITERATOR_T >
void VRML10Parser<ITERATOR_T>:: checkInputValidity ()
{
   bool valid = true;
   if ( dimensions == 2 ) {
      valid &= (coordinates.size() % 2) == 0 ? true : false;
   }
   else {
      assertion ( dimensions == 3 );
      valid &= (coordinates.size() % 3) == 0 ? true : false;
   }
}

}}} // namespace precice, io, impl

#endif // not PRECICE_NO_SPIRIT2

#endif // PRECICE_IO_IMPL_VRML10PARSER_HPP_


