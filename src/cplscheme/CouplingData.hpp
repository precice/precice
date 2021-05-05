#pragma once

#include <Eigen/Core>
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "time/Waveform.hpp"
#include "utils/EigenHelperFunctions.hpp"
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
    _previousIteration = data->values(); // initialize previous iteration with current data.
  }

  int getDimensions() const
  {
    PRECICE_ASSERT(data != nullptr);
    return data->getDimensions();
  }

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

  void storeIteration()
  {
    _previousIteration = this->values();
  }

  const Eigen::VectorXd previousIteration() const
  {
    return _previousIteration;
  }

  int getMeshID()
  {
    return mesh->getID();
  }

  std::vector<int> getVertexOffsets()
  {
    return mesh->getVertexOffsets();
  }

  ///  True, if the data values of this CouplingData require to be initialized by this participant.
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

  /// Data values of previous iteration.
  Eigen::VectorXd _previousIteration;

  mesh::PtrData data;

  mesh::PtrMesh mesh;
};

} // namespace cplscheme
} // namespace precice
