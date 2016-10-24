#pragma once

#include "spacetree/Spacetree.hpp"
#include "spacetree/impl/OctreeCell.hpp"
#include <vector>

namespace precice {
namespace spacetree {

/**
 * @brief Quad- (2D), DynamicOctree (3D).
 */
class DynamicOctree : public Spacetree
{
public:

  DynamicOctree (
    const Eigen::VectorXd& center,
    double halflength,
    double refinementLimit );

  virtual ~DynamicOctree() {}

  virtual void addMesh (  const mesh::PtrMesh& content );

  virtual const std::vector<mesh::PtrMesh>& meshes() const;

  virtual void initialize();

  virtual void meshChanged ( mesh::Mesh& mesh );

  virtual int searchPosition ( const Eigen::VectorXd& point );

  /**
   * @brief Spacetree accelerated distance search.
   */
  virtual void searchDistance ( query::FindClosest& findClosest );

  /**
   * @brief Spacetree accelerated voxel content search.
   */
  virtual int searchContent ( query::FindVoxelContent& find );

  virtual void accept ( Visitor& visitor );

  virtual void clear();

private:

  static logging::Logger _log;

  // @brief Contained meshes.
  std::vector<mesh::PtrMesh> _meshes;

  Eigen::VectorXd _center;

  double _halflength;

  double _refinementLimit;

  impl::OctreeCell _rootCell;

  bool _meshChanged;
};

}} // namespace precice, spacetree
