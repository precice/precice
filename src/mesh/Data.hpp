// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_MESH_DATA_HPP_
#define PRECICE_MESH_DATA_HPP_

#include "SharedPointer.hpp"
#include "utils/Globals.hpp"
#include "utils/Dimensions.hpp"
#include "utils/PointerVector.hpp"
#include "tarch/logging/Log.h"
#include "tarch/la/DynamicVector.h"
#include <string>

namespace precice {
  namespace mesh {
    class Mesh;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace mesh {

/**
 * @brief Describes a set of data values belonging to the vertices of a mesh.
 */
class Data
{
public:

  // @brief Possible types of data values.
//  enum DataType {
//    TYPE_UNDEFINED,
//    TYPE_DOUBLE,
//    TYPE_VECTOR
//  };

  // @brief Name of an undefined data type.
  //static const std::string TYPE_NAME_UNDEFINED;
  // @brief Name of a double data type.
  //static const std::string TYPE_NAME_DOUBLE;
  // @brief Name of a vector data type.
  //static const std::string TYPE_NAME_VECTOR;

  /**
   * @brief Returns the number of created (and still existing) Data objects.
   *
   * Used to give Data objects unique IDs.
   */
  static size_t getDataCount ();

  /**
   * @brief Sets the data counter to zero.
   *
   * Used in test cases where multiple scenarios with data are run.
   */
  static void resetDataCount ();

  /**
   * @brief Do not use this constructor! Only there for compatibility with std::map.
   */
  Data ();

  /**
   * @brief Constructor.
   */
  Data (
    const std::string& name,
    int                id,
    int                dimension );

  /**
   * @brief Destructor, decrements data count.
   */
  ~Data();

  /**
   * @brief Returns a reference to the data values.
   */
  utils::DynVector& values ();

  /**
   * @rief Returns a const reference to the data values.
   */
  const utils::DynVector& values () const;

  /**
   * @brief Returns the name of the data set.
   */
  const std::string& getName () const;

  /**
   * @brief Returns the ID of the data set (supposed to be unique).
   */
  int getID () const;

  /**
   * @brief Returns the type constant of the data set.
   */
//  DataType getType () const;

  /**
   * @brief Returns the dimension (i.e., number of components) of one data value.
   */
  int getDimensions () const;

  /**
   * @brief Returns the name constant of the type of the data set.
   */
  //const std::string& getTypeName () const;


private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  // @brief Counter for existing Data objects.
  static size_t _dataCount;

  // @brief
  utils::DynVector _values;

  // @brief Name of the data set.
  std::string _name;

  // @brief ID of the data set (supposed to be unique).
  int _id;

//  // @brief Type of data (scalar or vector).
//  DataType _type;

  // @brief Dimensionality of one data value.
  int _dimensions;

};

}} // namespace precice, mesh

#endif /* PRECICE_MESH_DATA_HPP_ */
