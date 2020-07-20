#pragma once

#include <iosfwd>
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace mapping {

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
   * the values. Values integrated over some area should be mapped conservative,
   * while area independent values such as pressure or stresses should be mapped
   * consistent.
   */
  enum Constraint {
    CONSISTENT,
    CONSERVATIVE
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

  /// Constructor, takes mapping constraint.
  Mapping(Constraint constraint, int dimensions);

  Mapping &operator=(Mapping &&) = delete;

  /// Destructor, empty.
  virtual ~Mapping() {}

  /**
   * @brief Sets input and output meshes carrying data to be mapped.
   *
   * @param[in] input Mesh with known data values to be mapped.
   * @param[in] output Mesh with unknwon data values to be computed from input.
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
  virtual bool hasComputedMapping() const = 0;

  /// Removes a computed mapping.
  virtual void clear() = 0;

  /**
   * @brief Maps input data to output data from input mesh to output mesh.
   *
   * Pre-conditions:
   * - hasComputedMapping() returns true
   *
   * Post-conditions:
   * - output values are computed from input values
   */
  virtual void map(
      int inputDataID,
      int outputDataID) = 0;

  /// Method used by partition. Tags vertices that could be owned by this rank.
  virtual void tagMeshFirstRound() = 0;

  /// Method used by partition. Tags vertices that can be filtered out.
  virtual void tagMeshSecondRound() = 0;

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

private:
  /// Determines wether mapping is consistent or conservative.
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

} // namespace mapping
} // namespace precice
