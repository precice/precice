#pragma once

#include <Eigen/Core>
#include "mesh/Data.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace cplscheme {

struct CouplingData {
  using DataMatrix = Eigen::MatrixXd;

  /// Data values of current iteration.
  Eigen::VectorXd *values;

  /// Data values of previous iteration (1st col) and previous time windows.
  DataMatrix oldValues;

  mesh::PtrMesh mesh;

  ///  True, if the data values if this CouplingData requires to be initialized by a participant.
  bool requiresInitialization;

  /// dimension of one data value (scalar=1, or vectorial=interface-dimension)
  int dimension;

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
      Eigen::VectorXd *values,
      mesh::PtrMesh    mesh,
      bool             requiresInitialization,
      int              dimension)
      : values(values),
        mesh(mesh),
        requiresInitialization(requiresInitialization),
        dimension(dimension)
  {
    PRECICE_ASSERT(values != NULL);
    PRECICE_ASSERT(mesh.use_count() > 0);
  }
};

} // namespace cplscheme
} // namespace precice
