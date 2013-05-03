// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_NEWSPACETREE_SPACETREE_HPP_
#define PRECICE_NEWSPACETREE_SPACETREE_HPP_

#include "spacetree/SharedPointer.hpp"
#include "mesh/Mesh.hpp"
#include "utils/Dimensions.hpp"
#include "boost/noncopyable.hpp"

namespace precice {
  namespace query {
    class FindClosest;
    //class FindClosestVertex;
    class FindVoxelContent;
  }
  namespace spacetree {
    //class ExportSpacetree;
  }
}

// ---------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace spacetree {

/**
 * @brief Abstract base class for space trees.
 *
 * A spacetree is used to accellerate a search for container objects representing
 * a 2D or 3D geometry. It does so by recursively subdividing the objects forming
 * the geometry into smaller groups. The subdivision is such that every group is
 * located in a predefined spatial neighborhood which is known. A search can be
 * tracked down to a small part of involved objects then, and the computational
 * effort is reduced.
 *
 * There are several partitioning schemes possible: the standard scheme subdivides
 * every spatial dimension into two equally sized parts, which results in four
 * childs in 2D and eight in 3D. This scheme is implemented in the class
 * StandardSpacetree. In future, a subdivision into three parts for every
 * dimension (nonal-tree) or a non-equal partition (kd-tree) seems to be interesting.
 */
class Spacetree : public mesh::Mesh::MeshListener, private boost::noncopyable
{
public:

  /**
   * @brief Interface for spacetree visitors, used to export spacetree.
   */
  struct Visitor
  {
    /**
     * @brief Callback from spacetree node.
     */
    virtual void nodeCallback (
      const utils::DynVector& center,
      const utils::DynVector& halflengths,
      int                     position ) = 0;

    /**
     * @brief Callback from spacetree leaf.
     */
    virtual void leafCallback (
      const utils::DynVector& center,
      const utils::DynVector& halflengths,
      int                     position,
      const mesh::Group&      content ) = 0;
  };

  static int positionUndefined() { return 0; }
  static int positionInsideOfGeometry() { return 1; }
  static int positionOutsideOfGeometry() { return 2; }
  static int positionOnGeometry() { return 3; }
  static int minElementsToRefineCell;

  /**
   * @brief Destructor.
   */
  virtual ~Spacetree() {}

  /**
   * @brief Adds a mesh to be contained in the spacetree.
   */
  virtual void addMesh ( const mesh::PtrMesh& content ) =0;

  virtual const std::vector<mesh::PtrMesh>& meshes() const =0;

  /**
   * @brief Initializes the spacetree.
   *
   * This method can be called only once after the creation of the spacetree,
   * or after clear() has been called.
   */
  virtual void initialize() =0;

  /**
   * @brief Callback method for mesh object that changed.
   */
  virtual void meshChanged ( mesh::Mesh& mesh ) =0;

  /**
   * @brief Spacetree accelerated position search.
   *
   * @return Position (see positionXYZ() methods).
   */
  virtual int searchPosition ( const utils::DynVector& point ) =0;

  /**
   * @brief Wrapper for searchPosition(utils::DynVector).
   */
  int searchPosition ( const utils::Vector3D& point )
  {
    utils::DynVector pointCopy(point);
    return searchPosition(pointCopy);
  }

  /**
   * @brief Wrapper for searchPosition(utils::DynVector).
   */
  int searchPosition ( const utils::Vector2D& point )
  {
    utils::DynVector pointCopy(point);
    return searchPosition(pointCopy);
  }

  /**
   * @brief Spacetree accelerated distance search.
   */
  virtual void searchDistance ( query::FindClosest& findClosest ) =0;

  /**
   * @brief Spacetree accelerated voxel content search.
   */
  virtual int searchContent ( query::FindVoxelContent& find ) =0;

  /**
   * @brief Traverses each cell of the spacetree
   */
  virtual void accept ( Visitor& visitor ) =0;

  /**
   * @brief Removes all cells from the spacetree.
   */
  virtual void clear() =0;
};

}} // namespace precice, spacetree

#endif // PRECICE_NEWSPACETREE_SPACETREE_HPP_
