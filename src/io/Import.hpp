// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IO_IMPORT_HPP_
#define PRECICE_IO_IMPORT_HPP_

#include "utils/Dimensions.hpp"
#include <string>
#include <map>

namespace precice {
  namespace mesh {
    class Mesh;
  }
}

// ------------------------------------------------------------ CLASS DEFINTION

namespace precice {
namespace io {

/**
 * @brief Abstract base class for import classes.
 *
 * Defines a unique interface to import geometries.
 */
class Import
{
public:

  Import ( const std::string& location );

  virtual ~Import() {}

  /**
   * @brief Does the import.
   *
   * @param mesh [IN/OUT] The importet elements are added to the mesh.
   */
  virtual void doImport (
    const std::string& name,
    mesh::Mesh&        mesh ) =0;

protected:

  const std::string& getLocation () const
  {
    return _location;
  }

private:

  std::string _location;
};

}} // namespace precice, io

#endif /* PRECICE_IO_IMPORT_HPP_ */
