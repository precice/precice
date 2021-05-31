#pragma once

#include "mapping/Mapping.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace action {

/**
 * @brief Abstract base class for configurable actions on data and/or meshes.
 *
 * Actions are executed on call of precice::SolverInterface::initialize(),
 * precice::SolverInterface::initializeData(), and precice::SolverInterface::advance(). They can change meshes and in particular
 * data values.
 */
class Action {
public:
  /// Defines the time and place of application of the action.
  enum Timing {
    ON_TIME_WINDOW_COMPLETE_POST, // On advancing to next time window, after adv. cpl scheme
    WRITE_MAPPING_PRIOR,          // Everytime, before write mapping
    WRITE_MAPPING_POST,           // Everytime, after write mapping and before advancing cpl scheme
    READ_MAPPING_PRIOR,           // Everytime, after advancing cpl scheme and before read mapping
    READ_MAPPING_POST             // Everytime, after read mapping
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
    * @param[in] timeStepSize Length of last time step computed.
    * @param[in] computedTimeWindowPart Sum of all time steps within current time window, i.e. part that is already computed.
    * @param[in] timeWindowSize Current time window size.
    */
  virtual void performAction(
      double time,
      double timeStepSize,
      double computedTimeWindowPart,
      double timeWindowSize) = 0;

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
