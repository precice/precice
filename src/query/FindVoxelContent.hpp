#ifndef PRECICE_QUERY_FINDVOXELCONTAINEDVISITOR_HPP_
#define PRECICE_QUERY_FINDVOXELCONTAINEDVISITOR_HPP_

#include "mesh/Group.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/logging/Log.h"
#include "utils/Helpers.hpp"
#include <vector>

namespace precice {
namespace query {

/**
 * @brief Find visitables that are contained in a given voxel.
 *
 */
class FindVoxelContent
{
public:

  /**
   * @brief Deciders to include/exclude elements touching the voxels bounds.
   */
  enum BoundaryInclusion {
    INCLUDE_BOUNDARY,
    EXCLUDE_BOUNDARY
  };

  /**
   * @brief Constructor.
   *
   * @param voxelCenter     [IN] Coordinates of the voxel center
   * @param halflengths     [IN] Half of the sidelengths of the voxel
   * @param includeTouching [IN] Decides whether touching objects are included
   */
  template<typename VECTORA_T, typename VECTORB_T>
  FindVoxelContent(
    const VECTORA_T&  voxelCenter,
    const VECTORB_T&  halflengths,
    BoundaryInclusion boundaryInclusion );

  /**
   * @brief Performs the find operation on the content of the container.
   *
   * In two dimensions, the content edges need to have valid normals. In three
   * dimensions, the triangles need to have valid normals.
   */
  template<typename CONTAINER_T>
  bool operator()( CONTAINER_T& container );

  /**
   * @brief Returns voxel center coordinates
   */
  const utils::DynVector& getVoxelCenter() const;

  /**
   * @brief Returns halflength of voxel
   */
  const utils::DynVector& getVoxelHalflengths() const;

  /**
   * @brief Returns true, if touching objects are included in content
   */
  BoundaryInclusion getBoundaryInclusion() const;

  /**
   * @brief Returns a Group with all contained visitables
   */
  mesh::Group& content();

  void clear();

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  // @brief Center point of the voxel
  utils::DynVector _voxelCenter;

  // @brief Half the sidelengths of the voxel
  utils::DynVector _voxelHalflengths;

  // @brief Determines, whether objects that touch the voxel are contained
  BoundaryInclusion _boundaryInclusion;

  // @brief Spatial dimensions of the voxel (2D or 3D).
  int _dimensions;

  // @brief Contained and partially contained visitables of the voxel
  mesh::Group _content;

  void checkVertex( mesh::Vertex& vertex );

  void checkEdge( mesh::Edge& edge );

  void checkTriangle( mesh::Triangle& triangle );

  /**
   * @brief Returns true, if a plane square and a segment intersect.
   */
  bool computeIntersection(
    const utils::Vector3D&  squareCenter,
    const utils::DynVector& halflengths,
    int                     squareNormalDirection,
    const utils::Vector3D&  firstPointSegment,
    const utils::Vector3D&  secondPointSegment,
    bool                    countTouchingAsIntersection ) const;

  bool computeIntersection(
    const mesh::Triangle&  triangle,
    const utils::Vector3D& firstPointSegment,
    const utils::Vector3D& secondPointSegment,
    bool                   countTouchingAsIntersection );
};

// ---------------------------------------------------------- HEAER DEFINITIONS

template<typename VECTORA_T, typename VECTORB_T>
FindVoxelContent:: FindVoxelContent
(
  const VECTORA_T&  voxelCenter,
  const VECTORB_T&  halflengths,
  BoundaryInclusion boundaryInclusion )
:
  _voxelCenter ( voxelCenter ),
  _voxelHalflengths ( halflengths ),
  _boundaryInclusion ( boundaryInclusion ),
  _dimensions ( voxelCenter.size() ),
  _content ()
{
  assertion ( voxelCenter.size() == halflengths.size(),
               voxelCenter.size(), halflengths.size() );
  assertion ( (_dimensions == 2) || (_dimensions == 3), _dimensions );
}

template< typename CONTAINER_T >
bool FindVoxelContent:: operator()
(
  CONTAINER_T& container )
{
  for ( mesh::Vertex& vertex : container.vertices() ){
    checkVertex ( vertex );
  }
  for ( mesh::Edge& edge : container.edges() ){
    checkEdge ( edge );
  }
  for ( mesh::Triangle& triangle : container.triangles() ){
    checkTriangle ( triangle );
  }
  return not _content.empty ();
}

}} // namespace precice, query

#endif /* PRECICE_QUERY_FINDVOXELCONTAINEDVISITOR_HPP_ */


