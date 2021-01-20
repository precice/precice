#pragma once

#include <Eigen/Core>
#include "mesh/Data.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace cplscheme {

struct CouplingData { // @todo: should be a class from a design standpoint. See https://github.com/precice/precice/pull/865#discussion_r495825098
  using DataMatrix = Eigen::MatrixXd;

  /// Returns a reference to the data values.
  Eigen::VectorXd &values()
  {
    PRECICE_ASSERT(data != nullptr);
    return data->values();
  }

  /// Returns a const reference to the data values.
  const Eigen::VectorXd &values() const
  {
    PRECICE_ASSERT(data != nullptr);
    return data->values();
  }

  /// Data values of previous iteration (1st col) and previous time windows.
  DataMatrix oldValues;

  mesh::PtrData data;

  mesh::PtrMesh mesh;

  ///  True, if the data values if this CouplingData requires to be initialized by a participant.
  bool requiresInitialization;

  int getDimensions()
  {
    PRECICE_ASSERT(data != nullptr);
    return data->getDimensions();
  }

  /**
   * @brief Default constructor, not to be used!
   *
   * Necessary when compiler creates template code for std::map::operator[].
   */
  CouplingData()
  {
    PRECICE_ASSERT(false);
  }

  CouplingData(
      mesh::PtrData data,
      mesh::PtrMesh mesh,
      bool          requiresInitialization)
      : data(data),
        mesh(mesh),
        requiresInitialization(requiresInitialization)
  {
    PRECICE_ASSERT(data != nullptr);
    PRECICE_ASSERT(mesh != nullptr);
    PRECICE_ASSERT(mesh.use_count() > 0);
  }
};

} // namespace cplscheme
} // namespace precice
