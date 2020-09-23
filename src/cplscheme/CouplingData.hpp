#pragma once

#include <Eigen/Core>
#include "mesh/Data.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace cplscheme {

struct CouplingData {
  using DataMatrix = Eigen::MatrixXd;

  [[deprecated("Use &(CouplingData::dataValues()). Even better: Don't use a pointer.")]]  // @todo: remove before merge, since it is C++14
  /// Data values of current iteration.
  Eigen::VectorXd *values;

  /// Returns a reference to the data values.
  Eigen::VectorXd &dataValues()
  {
    return data->values();
  }

  /// Returns a const reference to the data values.
  const Eigen::VectorXd &dataValues() const
  {
    return data->values();
  }

  /// Data values of previous iteration (1st col) and previous time windows.
  DataMatrix oldValues;

  mesh::PtrMesh mesh;

  mesh::PtrData data;

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
      : values(&(data->values())),
        mesh(mesh),
        requiresInitialization(requiresInitialization)
  {
    PRECICE_ASSERT(values != NULL);
    PRECICE_ASSERT(mesh.use_count() > 0);
  }
};

} // namespace cplscheme
} // namespace precice
