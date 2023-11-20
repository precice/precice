#pragma once

#include "mapping/Mapping.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace action {

/**
 * @brief Abstract base class for configurable actions on data and/or meshes.
 *
 * Actions are executed on call of precice::Participant::initialize() and precice::Participant::advance().
 * They can change meshes and in particular data values.
 */
class Action {
public:
  /// Defines the time and place of application of the action.
  enum Timing {
    WRITE_MAPPING_POST, // At the end of a time window, after write mapping (if existent) and before advancing cpl scheme
    READ_MAPPING_POST   // At the end of a time window, after read mapping (if existent)
  };

  Action(
      Timing                            timing,
      const mesh::PtrMesh &             mesh,
      mapping::Mapping::MeshRequirement requirement)
      : _timing(timing),
        _mesh(mesh),
        _meshRequirement(requirement)
  {
  }

  Action(
      Timing               timing,
      const mesh::PtrMesh &mesh)
      : _timing(timing),
        _mesh(mesh)
  {
  }

  Action &operator=(Action &&) = delete;

  /// Destructor, empty.
  virtual ~Action() {}

  /**
   * @brief Performs the action, to be overwritten by subclasses.
   *
   * @param[in] time the current total simulation time.
   */
  virtual void performAction(double time) = 0;

  /// Returns the timing of the action.
  Timing getTiming() const
  {
    return _timing;
  }

  /// Returns the mesh carrying the data used in the action.
  const mesh::PtrMesh &getMesh() const
  {
    return _mesh;
  }

  /// Returns the mesh requirement of this action
  mapping::Mapping::MeshRequirement getMeshRequirement() const
  {
    return _meshRequirement;
  }

private:
  /// Determines when the action will be executed.
  Timing _timing;

  /// Mesh carrying the data used in the action.
  mesh::PtrMesh _mesh;

  /// The mesh requirements for the mesh
  mapping::Mapping::MeshRequirement _meshRequirement = mapping::Mapping::MeshRequirement::UNDEFINED;
};

} // namespace action
} // namespace precice
