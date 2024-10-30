#pragma once

#include <Eigen/Core>
#include <iosfwd>

#include "mapping/MappingDataCache.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "precice/span.hpp"

namespace precice::mapping {

/**
 * @brief Abstract base class for mapping of data from one mesh to another.
 */
class Mapping {
public:
  /**
   * @brief Specifies additional constraints for a mapping.
   *
   * A consistent mapping retains mean values. When mapping displacements, e.g.
   * rigid body motions are retained. A conservative mapping retains the sum of
   * the values. The scaled-consistent-surface/volume mappings first map the values consistently,
   * then scales the mapped such that the integrals on both meshes are equal.
   * Integrals are either done on surfaces or volumes depending on the mode.
   * - Continuous fields such as displacements or temperatures should use consistent maps.
   * - Quantities whose sum is preserved such as forces should use conservative maps.
   * - Continuous fields whose integral matters, such as pressure or heat fluxes should be consistent or scaled-consistent.
   */
  enum Constraint {
    CONSISTENT,
    CONSERVATIVE,
    SCALED_CONSISTENT_SURFACE,
    SCALED_CONSISTENT_VOLUME
  };

  /**
   * @brief Specifies requirements for the input and output meshes of a mapping.
   *
   * Different mapping types have different requirements on the meshes involved
   * in the mapping, while the input and output mesh holding the data to map can
   * have different requirements. FULL requires a mesh consisting of vertices
   * connected by edges and faces. VERTEX requires a mesh consisting of vertices
   * only.
   */
  enum class MeshRequirement {
    UNDEFINED = 0,
    /// Vertices only.
    VERTEX = 1,
    /// Full mesh.
    FULL = 2
  };

  /**
   * @brief Specifies whether the mapping requires an initial guess
   *
   * Iterative mappings use an additional initial guess to perform the mapping.
   * When calling the version of map with the initialGuess, derived classes of Mapping can
   * access and update this initial guess using \ref initialGuess().
   * Note that the size of the initial guess is controlled by the Mapping.
   *
   * The first initial guess is expected to be an empty VectorXd.
   */
  enum class InitialGuessRequirement : bool {
    Required = true,
    None     = false
  };

  /// Constructor, takes mapping constraint.
  Mapping(Constraint constraint, int dimensions, bool requiresGradientData, InitialGuessRequirement initialGuessRequirement);

  Mapping &operator=(Mapping &&) = delete;

  /// Destructor, empty.
  virtual ~Mapping() = default;

  /**
   * @brief Sets input and output meshes carrying data to be mapped.
   *
   * @param[in] input Mesh with known data values to be mapped.
   * @param[in] output Mesh with unknown data values to be computed from input.
   */
  void setMeshes(
      const mesh::PtrMesh &input,
      const mesh::PtrMesh &output);

  const mesh::PtrMesh &getInputMesh() const;

  const mesh::PtrMesh &getOutputMesh() const;

  /// Returns the constraint (consistent/conservative) of the mapping.
  Constraint getConstraint() const;

  /// Returns the requirement on the input mesh.
  MeshRequirement getInputRequirement() const;

  /// Returns the requirement on the output mesh.
  MeshRequirement getOutputRequirement() const;

  /// Computes the mapping coefficients from the in- and output mesh.
  virtual void computeMapping() = 0;

  /**
   * @brief Returns true, if the mapping has been computed.
   *
   * After a call to clear(), a computed mapping is removed and false returned.
   */
  bool hasComputedMapping() const;

  /// Checks whether the mapping has the given constraint or not
  virtual bool hasConstraint(const Constraint &constraint) const;

  /// Returns true if mapping is a form of scaled consistent mapping
  bool isScaledConsistent() const;

  /// Return true if the mapping requires an initial guess
  bool requiresInitialGuess() const;

  bool isIndirectMapping() const;

  /// Return the provided initial guess of a mapping using an initialGuess
  const Eigen::VectorXd &initialGuess() const;

  Eigen::VectorXd &initialGuess();

  /// True if initialGuess().size() == 0
  bool hasInitialGuess() const;

  /// Removes a computed mapping.
  virtual void clear() = 0;

  /// @deprecated
  void map(int inputDataID, int outputDataID);
  /// @deprecated
  void map(int inputDataID, int outputDataID, Eigen::VectorXd &initialGuess);

  /**
   * @brief Maps an input \ref Sample to output data from input mesh to output mesh.
   *
   * Derived classed must implement the mapping functionality using mapConsistent() and mapConservative()
   *
   * @param[in] input sample to map
   * @param[out] output result data
   *
   * @pre \ref hasComputedMapping() == true
   * @pre \ref requiresInitialGuess() == false
   *
   * @post output contains the mapped data
   */
  void map(const time::Sample &input, Eigen::VectorXd &output);

  /**
   * @brief Maps an input \ref Sample to output data from input mesh to output mesh, given an initialGuess.
   *
   * The initialGuess must be either an empty VectorXd, or the result of a previous invocation of \ref map.
   *
   * Derived classed must implement the mapping functionality using mapConsistent() and mapConservative()
   *
   * @warning this is the map version which requires an initial guess
   *
   * @param[in] input sample to map
   * @param[out] output result data
   * @param[inout] initialGuess guess to use during map, may be an empty VectorXd
   *
   * @pre initialGuess is either an empty VectorXd or the result of a previous invocation of \ref map
   * @pre \ref hasComputedMapping() == true
   * @pre \ref requiresInitialGuess() == true
   *
   * @post output contains the mapped data
   * @post \ref initialGuess() contains the initial guess for the next call to \ref map
   */
  void map(const time::Sample &input, Eigen::VectorXd &output, Eigen::VectorXd &initialGuess);

  /// Method used by partition. Tags vertices that could be owned by this rank.
  virtual void tagMeshFirstRound() = 0;

  /// Method used by partition. Tags vertices that can be filtered out.
  virtual void tagMeshSecondRound() = 0;

  /**
   * @brief Scales the consistently mapped output data such that the surface integral
   * of the values on input mesh and output mesh are equal
   *
   *
   * @pre Input and output mesh should have full connectivity information.
   */
  virtual void scaleConsistentMapping(const Eigen::VectorXd &input, Eigen::VectorXd &output, Constraint type) const;

  /// Returns whether the mapping requires gradient data
  bool requiresGradientData() const;

  /// Returns the name of the mapping method for logging purpose
  virtual std::string getName() const = 0;

  virtual void writeConservativeAt(::precice::span<const double> coordinates, Eigen::Map<const Eigen::MatrixXd> &source, Eigen::Map<Eigen::MatrixXd> &target);

  // @todo consider making this a private method in the RBFMapping/PUM mapping class
  virtual void updateMappingDataCache(MappingDataCache &cache, Eigen::VectorXd &in);

  // For now only for read-consistent
  virtual void evaluateMappingDataCacheAt(::precice::span<const double> coordinates, const MappingDataCache &cache, ::precice::span<double> values);

protected:
  /// Returns pointer to input mesh.
  mesh::PtrMesh input() const;

  /// Returns pointer to output mesh.
  mesh::PtrMesh output() const;

  /// Sets the mesh requirement for the input mesh.
  void setInputRequirement(MeshRequirement requirement);

  /// Sets the mesh requirement for the output mesh.
  void setOutputRequirement(MeshRequirement requirement);

  int getDimensions() const;

  /// Flag to indicate whether computeMapping() has been called.
  bool _hasComputedMapping = false;

  /// Flag if gradient data is required for the mapping
  bool _requiresGradientData;

  /**
   * @brief Maps data using a conservative constraint
   *
   * @param[in] input Sample to map data from
   * @param[in] output Values to map to
   *
   * If requiresInitialGuess(), then the initial guess is available via initialGuess().
   * Provide a new initial guess by overwriting it.
   * The mapping has full control over its size.
   *
   * @see For mappings requiring an initialGuess: initialGuess() hasInitialGuess()
   */
  virtual void mapConservative(const time::Sample &input, Eigen::VectorXd &output) = 0;
  /**
   * @brief Maps data using a consistent constraint
   *
   * @param[in] input Sample to map data from
   * @param[in] output Values to map to
   *
   * If requiresInitialGuess(), then the initial guess is available via initialGuess().
   * Provide a new initial guess by overwriting it.
   * The mapping has full control over its size.
   *
   * @see For mappings requiring an initialGuess: initialGuess() hasInitialGuess()
   */
  virtual void mapConsistent(const time::Sample &input, Eigen::VectorXd &output) = 0;

private:
  /// Determines whether mapping is consistent or conservative.
  Constraint _constraint;

  /// Requirement on input mesh.
  MeshRequirement _inputRequirement;

  /// Requirement on output mesh.
  MeshRequirement _outputRequirement;

  /// Pointer to input mesh.
  mesh::PtrMesh _input;

  /// Pointer to output mesh.
  mesh::PtrMesh _output;

  int _dimensions;

  /// The InitialGuessRequirement of the Mapping
  InitialGuessRequirement _initialGuessRequirement;

  /// Pointer to the initialGuess set and unset by \ref map.
  Eigen::VectorXd *_initialGuess = nullptr;
};

/** Defines an ordering for MeshRequirement in terms of specificality
* @param[in] lhs the left-hand side of the binary operator
* @param[in] rhs the right-hand side of the binary operator
*/
bool operator<(Mapping::MeshRequirement lhs, Mapping::MeshRequirement rhs);

/** Defines the output operation to streams
* @param[in,out] out stream to output to.
* @param[in] val the value to output.
*/
std::ostream &operator<<(std::ostream &out, Mapping::MeshRequirement val);

} // namespace precice::mapping
