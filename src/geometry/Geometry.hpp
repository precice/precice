// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_GEOMETRY_GEOMETRY_HPP_
#define PRECICE_GEOMETRY_GEOMETRY_HPP_

#include "utils/Dimensions.hpp"
#include "tarch/logging/Log.h"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"

namespace precice {
   namespace mesh {
      class Mesh;
   }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace geometry {

/**
 * @brief Abstract base class for all geometries.
 *
 * A geometry is there to create geometrical objects, based on the data
 * structure defined in component container. It uses the Mesh class to setup
 * a specific kind of geometry. In addition, every geometry is enanced by a
 * unique name and id.
 */
class Geometry
{
public:

  /**
   * @brief Constructor.
   *
   * @param name [IN] Unique name for the geometry.
   * @param isVolumeEnclosed [IN] If true, the volume of the geometry is inside
   *                              of its (closed) surface.
   * @param offset [IN] Constant offset of all vertex coordinates.
   */
  Geometry ( const utils::DynVector& offset );

  /**
   * @brief Destructor.
   */
  virtual ~Geometry() {}

  /**
   * @brief Creates the geometry into the given seed Mesh.
   */
  void create ( mesh::Mesh& seed );

  /**
   * @brief Sets an offset from zero for the geometry.
   */
  void setOffset ( const utils::DynVector& offset )
  {
    _offset = offset;
  }

  /**
   * @brief Returns the offset of the geometry from zero.
   */
  const utils::DynVector& getOffset () const
  {
    return _offset;
  }


protected:

  /**
   * @brief Is called from create() to actually create the geometry.
   */
  virtual void specializedCreate ( mesh::Mesh& seed ) =0;

  /**
   * @brief By default, mesh data is allocated after specializedCreate.
   *
   * If a subclass of Geometry already allocated the mesh data and
   * assigned values, this method has to be overwritten to not to allocate and
   * overwrite the data values.
   */
  virtual void allocateDataValues ( mesh::Mesh& mesh );

private:

  static tarch::logging::Log _log;

  // @brief Offset of reference point of geometry from zero point
  utils::DynVector _offset;

};

}} // namespace precice, geometry

#endif /* PRECICE_GEOMETRY_GEOMETRY_HPP_ */
