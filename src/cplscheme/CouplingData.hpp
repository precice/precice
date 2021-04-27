#pragma once

#include <Eigen/Core>
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace cplscheme {

class CouplingData {
public:
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

  int getDimensions() const
  {
    PRECICE_ASSERT(data != nullptr);
    return data->getDimensions();
  }

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

  int getMeshID()
  {
    return mesh->getID();
  }

  std::vector<int> getVertexOffsets()
  {
    return mesh->getVertexOffsets();
  }

  /// Data values of previous iteration (1st col) and previous time windows.
  DataMatrix oldValues;

  ///  True, if the data values if this CouplingData requires to be initialized by a participant.
  const bool requiresInitialization;

private:
  /**
   * @brief Default constructor, not to be used!
   *
   * Necessary when compiler creates template code for std::map::operator[].
   */
  CouplingData()
      : requiresInitialization(false)
  {
    PRECICE_ASSERT(false);
  }

  mesh::PtrData data;

  mesh::PtrMesh mesh;
};

} // namespace cplscheme
} // namespace precice
