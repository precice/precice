#pragma once

#include "spacetree/Spacetree.hpp"
#include "spacetree/impl/OctreeCell.hpp"

namespace precice {
namespace spacetree {

/**
 * @brief Quad- (2D), Octree (3D) that is built completely on initialization.
 */
class StaticOctree : public Spacetree
{
public:

  /**
   * @brief Constructor.
   *
   * @param[in] center Center of the root cell.
   * @param[in] halflength  Half sidelength of the root cell.
   * @param[in] refinementLimit  Fully refined cell is guaranteed to have smaller
   *        or equal sidelengths.
   */
  StaticOctree (
    const Eigen::VectorXd& center,
    double halflength,
    double refinementLimit );

  /**
   * @brief Destructor, empty.
   */
  virtual ~StaticOctree() {}

  /**
   * @brief Adds a mesh to the spacetree.
   *
   * The spacetree is rebuilt on the next search call.
   */
  virtual void addMesh ( const mesh::PtrMesh& content );

  /**
   * @brief Returns the meshes added to the spacetree.
   */
  virtual const std::vector<mesh::PtrMesh>& meshes() const;

  /**
   * @brief (Re)Initializes the spacetree structure.
   *
   * Is called automatically, after meshChanged was called, and a search method
   * is envoked.
   */
  virtual void initialize();

  /**
   * @brief Sets flag, which leads to a rebuild of the tree when issuing search.
   */
  virtual void meshChanged ( mesh::Mesh& mesh );

  /**
   * @brief Searches for the position (inside/outside/on) of a point in space.
   */
  virtual int searchPosition ( const Eigen::VectorXd& point );

  /**
   * @brief Searches for the distance to the next mesh.
   */
  virtual void searchDistance ( query::FindClosest& findClosest );

  /**
   * @brief Searches for the content and position of a hexahedron.
   */
  virtual int searchContent ( query::FindVoxelContent& find );

  /**
   * @brief Visitor entry to visit all cells of the spacetree.
   */
  virtual void accept ( Visitor& visitor );

  /**
   * @brief Clears the spacetree.
   */
  virtual void clear();

private:

  static logging::Logger _log;

  // @brief Contained meshes.
  std::vector<mesh::PtrMesh> _meshes;

  // @brief Center of root cell.
  Eigen::VectorXd _center;

  // @brief Half sidelengths of root cell.
  double _halflength;

  // @brief Guaranteed max sidelengths of fully refined cell.
  double _refinementLimit;

  // @brief Root cell data structure.
  impl::OctreeCell _rootCell;

  // @brief Flag to signal that a contained mesh has changed.
  bool _meshChanged;
};

}} // namespace precice, spacetree
