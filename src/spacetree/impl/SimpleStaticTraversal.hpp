#pragma once

#include "spacetree/Spacetree.hpp"
#include "spacetree/impl/Environment.hpp"
#include "mesh/Merge.hpp"
#include "query/FindClosest.hpp"
#include "query/FindVoxelContent.hpp"
#include "utils/PointerVector.hpp"
// #include "utils/Dimensions.hpp"
#include "logging/Logger.hpp"
// #include "utils/Helpers.hpp"
#include "math/math.hpp"
#include <list>
// #include <memory>

namespace precice {
namespace spacetree {
namespace impl {

template<typename CELL_T>
class StaticTraversal
{
public:

  void refineAll (
    CELL_T&                 cell,
    const Eigen::VectorXd& cellCenter,
    const Eigen::VectorXd& cellHalflengths,
    double                  refinementLimit );

  int searchPosition (
    CELL_T&                 cell,
    const Eigen::VectorXd& searchPoint,
    const Eigen::VectorXd& cellCenter,
    const Eigen::VectorXd& cellHalflengths ) const;

  bool searchDistance (
    CELL_T&                 cell,
    query::FindClosest&     findClosest,
    const Eigen::VectorXd& cellCenter,
    const Eigen::VectorXd& cellHalflengths ) const;

  int searchContent (
    CELL_T&                  cell,
    query::FindVoxelContent& findContent,
    const Eigen::VectorXd&  cellCenter,
    const Eigen::VectorXd&  cellHalflengths ) const;

private:

  struct Root
  {
    CELL_T& cell;
    Eigen::VectorXd center;
    Eigen::VectorXd halflengths;
    Root (
        CELL_T& cell,
        const Eigen::VectorXd& center,
        const Eigen::VectorXd& halflengths )
    :
      cell(cell),
      center(center),
      halflengths(halflengths)
    {}
  };

  struct SearchPositionResult
  {
    int position;
    bool ambiguous;
    std::list<CELL_T*> uncachedCells;
    std::list<Eigen::VectorXd> uncachedCellCenters;
    query::FindClosest findClosest;

    SearchPositionResult ( const Eigen::VectorXd& searchPoint )
    : position(Spacetree::positionUndefined()), ambiguous(false), uncachedCells(),
      uncachedCellCenters(), findClosest(searchPoint) {}
  };

  struct SearchContentResult
  {
    int position;
    std::list<CELL_T*> uncachedCells;
    std::list<Eigen::VectorXd> uncachedCellCenters;

    SearchContentResult()
    : position(Spacetree::positionUndefined()), uncachedCells(),
      uncachedCellCenters() {}
  };

  static logging::Logger _log;

  void setCellPositions (
    CELL_T&                 cell,
    const Eigen::VectorXd& cellCenter,
    const Eigen::VectorXd& cellHalflengths,
    const Root&             root );

  void refineAllInternal (
    CELL_T&                 cell,
    const Eigen::VectorXd& cellCenter,
    const Eigen::VectorXd& cellHalflengths,
    double                  refinementLimit );

  std::shared_ptr<SearchPositionResult> searchPositionInternal (
    CELL_T&                 cell,
    const Eigen::VectorXd& searchPoint,
    const Eigen::VectorXd& cellCenter,
    const Eigen::VectorXd& cellHalflengths ) const;

  std::shared_ptr<SearchContentResult> searchContentInternal (
    CELL_T&                  cell,
    query::FindVoxelContent& findContent,
    const Eigen::VectorXd&  cellCenter,
    const Eigen::VectorXd&  cellHalflengths ) const;

  double distanceToBoundary (
    const Eigen::VectorXd& cellCenter,
    const Eigen::VectorXd& cellHalflengths,
    const Eigen::VectorXd& searchPoint ) const;

  bool isCovered (
    const Eigen::VectorXd& cellCenter,
    const Eigen::VectorXd& cellHalflengths,
    const Eigen::VectorXd& voxelCenter,
    const Eigen::VectorXd& voxelHalflengths ) const;

  bool isOverlapped (
    const Eigen::VectorXd& cellCenter,
    const Eigen::VectorXd& cellHalflengths,
    const Eigen::VectorXd& voxelCenter,
    const Eigen::VectorXd& voxelHalflengths ) const;

  template<typename VISITOR_T>
  void visitAllCells (
    CELL_T&    cell,
    VISITOR_T& visitor ) const;

  template<typename VISITOR_T>
  void visitRemainingCells (
    const CELL_T& toExclude,
    CELL_T&       cell,
    VISITOR_T&    visitor ) const;
};

// ----------------------------------------------------- HEADER IMPLEMENTATIONS

template<typename CELL_T>
logging::Logger StaticTraversal<CELL_T>::
  _log("spacetree::StaticTraversal");

template<typename CELL_T>
void StaticTraversal<CELL_T>:: refineAll
(
  CELL_T&                 cell,
  const Eigen::VectorXd& cellCenter,
  const Eigen::VectorXd& cellHalflengths,
  double                  refinementLimit )
{
  preciceTrace ( "refineAll()", cellCenter, cellHalflengths, refinementLimit );

//  // The environment gives information on the position of the cells surrounding
//  // a current cell of consideration.
//  utils::DynVector outsidePoint(cellHalflengths);
//  outsidePoint *= 2.0;
//  outsidePoint += cellCenter;
//  query::FindClosest findClosest(outsidePoint);
//  findClosest(cell.content());
//  assertion(not tarch::la::equals(findClosest.getClosest().distance, 0.0));
//  int env = findClosest.getClosest().distance > 0
//            ? Spacetree::positionOutsideOfGeometry()
//            : Spacetree::positionInsideOfGeometry();

  refineAllInternal(cell, cellCenter, cellHalflengths, refinementLimit);
  Root root(cell, cellCenter, cellHalflengths);
  setCellPositions(cell, cellCenter, cellHalflengths, root);
//  cell.setPosition(env);
//  SetCellPositionsResult result;
//  utils::DynVector childCenter(cellCenter.size());
//  utils::DynVector childHalflengths(cellCenter.size());
//  setCellPositions ( cell, cellCenter, cellHalflengths, result );
//  while ( not result.cells.empty() ){
//    CELL_T& leftCell = *result.cells.front();
//    assertion(not leftCell.isLeaf());
//    const utils::DynVector& center = result.cellCenters.front();
//    const utils::DynVector& h = result.cellHalflengths.front();
//    for (int i=0; i < leftCell.getChildCount(); i++){
//      if (leftCell.child(i).getPosition() == Spacetree::positionUndefined()){
//        leftCell.getChildData(i, center, h, childCenter, childHalflengths);
//        query::FindClosest findDistance ( childCenter );
//        searchDistance ( cell, findDistance, cellCenter, cellHalflengths );
//        assertion(not tarch::la::equals(findDistance.getEuclidianDistance(), 0.0));
//        int pos = findDistance.getClosest().distance > 0
//                  ? Spacetree::positionOutsideOfGeometry()
//                  : Spacetree::positionInsideOfGeometry();
//        leftCell.setPosition(pos);
//        break;
//      }
//    }
//    assertion(leftCell.getPosition() != Spacetree::positionOnGeometry());
//    setCellPositions ( leftCell, center, h, result );
//    result.cells.pop_front();
//    result.cellCenters.pop_front();
//    result.cellHalflengths.pop_front();
//  }
}

template<typename CELL_T>
void StaticTraversal<CELL_T>:: refineAllInternal
(
  CELL_T&                 cell,
  const Eigen::VectorXd& cellCenter,
  const Eigen::VectorXd& cellHalflengths,
  double                  refinementLimit )
{
  preciceTrace ( "refineAllInternal()", cellCenter, cellHalflengths,
                  refinementLimit, cell.isLeaf() );
  if ( cell.isLeaf() && cell.needsRefinement(cellHalflengths, refinementLimit)){
    DEBUG ( "  Refine leaf" );
    assertion ( not cell.content().empty() );
    assertion ( cell.getPosition() == Spacetree::positionOnGeometry() );
    cell.refine(cellCenter, cellHalflengths);
  }

  if ( not cell.isLeaf() ){
    DEBUG ( "  Node" );
    assertion ( cell.getPosition() != Spacetree::positionUndefined() );
    Eigen::VectorXd childCenter(cellCenter.size());
    Eigen::VectorXd childHalflengths(cellCenter.size());
    for ( int i=0; i < cell.getChildCount(); i++ ){
      cell.getChildData(i, cellCenter, cellHalflengths, childCenter, childHalflengths);
      CELL_T& childCell = cell.child(i);
      refineAllInternal ( childCell, childCenter, childHalflengths, refinementLimit );
    }
  }
}

template<typename CELL_T>
void StaticTraversal<CELL_T>:: setCellPositions
(
  CELL_T&                 cell,
  const Eigen::VectorXd& cellCenter,
  const Eigen::VectorXd& cellHalflengths,
  const Root&             root )
{
  preciceTrace("setCellPositions()", cellCenter, cellHalflengths, root.center, root.halflengths);
  assertion(not cell.isLeaf());
  int position = Spacetree::positionUndefined();
  Eigen::VectorXd childCenter(cellCenter.size());
  Eigen::VectorXd childH(cellCenter.size());
  for (int i=0; i < cell.getChildCount(); i++){
    CELL_T& child = cell.child(i);
    if (child.getPosition() == Spacetree::positionUndefined()){
      assertion(child.isLeaf());
      if (position == Spacetree::positionUndefined()){
        cell.getChildData(i, cellCenter, cellHalflengths, childCenter, childH);
        query::FindClosest findClosest(childCenter);
        searchDistance(root.cell, findClosest, root.center, root.halflengths);
        assertion(findClosest.hasFound());
        assertion(math::greater(findClosest.getEuclidianDistance(), 0.0));
        int position = findClosest.getClosest().distance > 0
                   ? Spacetree::positionOutsideOfGeometry()
                   : Spacetree::positionInsideOfGeometry();
        DEBUG("  Computed position " << position << " for child " << i
                     << " at " << childCenter);
        child.setPosition(position);
      }
      else {
        DEBUG("  Setting known position " << position << " for child " << i);
        assertion((position == Spacetree::positionOutsideOfGeometry())
                  || (position == Spacetree::positionInsideOfGeometry()));
        child.setPosition(position);
      }
    }
    else if (not child.isLeaf()){
      assertion(child.getPosition() == Spacetree::positionOnGeometry());
      cell.getChildData(i, cellCenter, cellHalflengths, childCenter, childH);
      DEBUG("  Recursing to child " << i);
      setCellPositions(child, childCenter, childH, root);
    }
  }
}

template<typename CELL_T>
int StaticTraversal<CELL_T>:: searchPosition
(
  CELL_T&                 cell,
  const Eigen::VectorXd& searchPoint,
  const Eigen::VectorXd& cellCenter,
  const Eigen::VectorXd& cellHalflengths ) const
{
  preciceTrace ( "searchPosition()", searchPoint, cellCenter, cellHalflengths );
  typedef std::shared_ptr<SearchPositionResult> PtrResult;
  PtrResult result = searchPositionInternal ( cell, searchPoint, cellCenter,
                                              cellHalflengths );
  assertion ( result.use_count() > 0 );
  assertion ( result->position != Spacetree::positionUndefined() );
  assertion ( result->uncachedCells.empty() );

  return result->position;
}

template<typename CELL_T>
bool StaticTraversal<CELL_T>:: searchDistance
(
  CELL_T&                 cell,
  query::FindClosest&     findClosest,
  const Eigen::VectorXd& cellCenter,
  const Eigen::VectorXd& cellHalflengths ) const
{
  preciceTrace ( "searchDistance()", cellCenter, cellHalflengths, findClosest.getSearchPoint() );
  if ( cell.isLeaf() ){
    DEBUG ( "  Leaf" );
    findClosest ( cell.content() );
  }

  if ( not cell.isLeaf() ){
    DEBUG ( "  Node" );
    int childIndex = cell.getChildIndex (findClosest.getSearchPoint(), cellCenter,
                                         cellHalflengths);
    Eigen::VectorXd newCenter(cellCenter);
    Eigen::VectorXd newHalflengths(cellHalflengths);
    cell.getChildData (childIndex, cellCenter, cellHalflengths, newCenter,
                       newHalflengths);
    CELL_T& subtree = cell.child ( childIndex );
    bool ambiguous = searchDistance(subtree, findClosest, newCenter,
                                    newHalflengths );
    if ( ambiguous ){
      visitRemainingCells ( subtree, cell, findClosest );
    }
    //return false; // Found non-ambiguous
  }

  if ( findClosest.hasFound() ){
    double distance = distanceToBoundary(cellCenter, cellHalflengths,
                                         findClosest.getSearchPoint());
    bool isAmbiguous = math::greater ( findClosest.getEuclidianDistance(), distance );
    DEBUG ( "  hasfound, distance = " << findClosest.getEuclidianDistance()
                   << ", return " << isAmbiguous );
    return isAmbiguous;
  }
  DEBUG ( "  return true" );
  return true;
}

template<typename CELL_T>
int StaticTraversal<CELL_T>:: searchContent
(
  CELL_T&                  cell,
  query::FindVoxelContent& findContent,
  const Eigen::VectorXd&  cellCenter,
  const Eigen::VectorXd&  cellHalflengths ) const
{
  preciceTrace ( "searchContent()", findContent.getVoxelCenter(),
                  findContent.getVoxelHalflengths() );
  std::shared_ptr<SearchContentResult> result = searchContentInternal (
      cell, findContent, cellCenter, cellHalflengths );
  assertion ( result.use_count() > 0 );
  assertion ( result->uncachedCells.empty() );
  if ( result->position == Spacetree::positionUndefined() ){
    // Some/all searched cells had content, but not in the search voxel
    if ( result->position == Spacetree::positionUndefined() ){
      DEBUG ( "Computing position of search voxel" );
      query::FindClosest findDistance ( findContent.getVoxelCenter() );
      searchDistance ( cell, findDistance, cellCenter, cellHalflengths );
      double distance = findDistance.getClosest().distance;
      assertion ( not math::equals(distance, 0.0) );
      result->position = distance > 0 ? Spacetree::positionOutsideOfGeometry()
                                      : Spacetree::positionInsideOfGeometry();
    }
  }

  assertion ( result->position != Spacetree::positionUndefined() );
  DEBUG ( "return content().size() = " << findContent.content().size() );
  return result->position;
}

template<typename CELL_T>
std::shared_ptr<typename StaticTraversal<CELL_T>::SearchPositionResult>
StaticTraversal<CELL_T>:: searchPositionInternal
(
  CELL_T&                 cell,
  const Eigen::VectorXd& searchPoint,
  const Eigen::VectorXd& cellCenter,
  const Eigen::VectorXd& cellHalflengths ) const
{
  preciceTrace ( "searchPositionInternal()", searchPoint, cellCenter,
                  cellHalflengths );
  using namespace math;
  std::shared_ptr<SearchPositionResult> data;
  double distance = 0.0;
  if ( cell.isLeaf() ){
    DEBUG ( "  Leaf" );
    assertion ( data.use_count() == 0 );
    data = std::shared_ptr<SearchPositionResult> ( new SearchPositionResult(searchPoint) );
    if ( cell.getPosition() == Spacetree::positionOnGeometry() ){
      DEBUG ( "    Has content" );
      //query::FindClosest findClosest ( searchPoint );
      data->findClosest ( cell.content() );
      if ( data->findClosest.hasFound() ){
        DEBUG ( "    Found elements" );
        distance = data->findClosest.getClosest().distance;
        data->position = Spacetree::positionOnGeometry(); // May by altered later
      }
    }
    else {
      DEBUG ( "    Is empty" );
      assertion ( cell.content().empty() );
      assertion ( cell.getPosition() != Spacetree::positionUndefined() );
      data->position = cell.getPosition();
    }
  }

  if ( not cell.isLeaf() ) { // could be a leaf on entrance to searchPositionInternal
    DEBUG ( "  Node" );
    //assertion ( cell.getPosition() == Spacetree::positionOnGeometry() );
    int childIndex = cell.getChildIndex(searchPoint, cellCenter, cellHalflengths);
    Eigen::VectorXd newCenter(cellCenter);
    Eigen::VectorXd newHalflengths(cellHalflengths);
    cell.getChildData (childIndex, cellCenter, cellHalflengths, newCenter,
                       newHalflengths);
    CELL_T& subtree = cell.child(childIndex);
    assertion ( data.use_count() == 0 );
    data = searchPositionInternal ( subtree, searchPoint, newCenter, newHalflengths );
    assertion ( data->position != Spacetree::positionUndefined() );
    if ( data->ambiguous ){
      DEBUG ( "    Did not find elements or ambiguous, visit others" );
      //query::FindClosest findClosest ( searchPoint );
      visitRemainingCells(subtree, cell, data->findClosest);
      if ( data->findClosest.hasFound() ){
        DEBUG ( "    Found elements in others" );
        distance = data->findClosest.getClosest().distance;
        data->position = Spacetree::positionOnGeometry();
        data->ambiguous = false;
      }
    }
  }

  // Set inside/outside and check for ambiguities
  if ( not equals(distance, 0.0) ){
    DEBUG( "  Checking for ambiguities of found objects");
    assertion ( data.use_count() > 0 );
    if ( greater(distance, 0.0) ){
      data->position = Spacetree::positionOutsideOfGeometry();
    }
    else if ( greater(0.0, distance) ){
      data->position = Spacetree::positionInsideOfGeometry();
    }
    DEBUG ( "  found pos = " << data->position );
    double distanceToBound =  distanceToBoundary(cellCenter, cellHalflengths,
                                                 searchPoint);
    if ( greater(std::abs(distance), distanceToBound) ){
      DEBUG ( "  is ambigious" );
      data->ambiguous = true;
    }
  }
  DEBUG ( "  return position = " << data->position << ", ambiguous = "
                 << data->ambiguous );
  return data;
}

template<typename CELL_T>
std::shared_ptr<typename StaticTraversal<CELL_T>::SearchContentResult>
StaticTraversal<CELL_T>:: searchContentInternal
(
  CELL_T&                  cell,
  query::FindVoxelContent& findContent,
  const Eigen::VectorXd&  cellCenter,
  const Eigen::VectorXd&  cellHalflengths ) const
{
  preciceTrace ( "searchContentInternal()", cellCenter, cellHalflengths,
                  findContent.getVoxelCenter(), findContent.getVoxelHalflengths() );
  if ( cell.isLeaf() ){
    DEBUG ( "Leaf..." );
    std::shared_ptr<SearchContentResult> data ( new SearchContentResult() );
    assertion ( cell.getPosition() != Spacetree::positionUndefined() );
    if ( (cell.getPosition() == Spacetree::positionOutsideOfGeometry())
         || (cell.getPosition() == Spacetree::positionInsideOfGeometry()) )
    {
      DEBUG ( "empty, position = " << cell.getPosition() );
      data->position = cell.getPosition();
    }
    else {
      DEBUG ( "content size = " << cell.content().size() );
      assertion ( cell.getPosition() == Spacetree::positionOnGeometry() );
      assertion ( not cell.content().empty() );
      bool set = isCovered ( cellCenter, cellHalflengths,
                             findContent.getVoxelCenter(),
                             findContent.getVoxelHalflengths() );
      set &= findContent.getBoundaryInclusion() == query::FindVoxelContent::INCLUDE_BOUNDARY;
      if ( set ){
        DEBUG ( "Is covered by voxel, add cell content to find content" );
        findContent.content().add ( cell.content() );
        data->position = Spacetree::positionOnGeometry();
      }
      else {
        DEBUG ( "Isn't covered by voxel, apply find content" );
        findContent ( cell.content() );
        if ( not findContent.content().empty() ){
          DEBUG("Some content is contained in voxel");
          data->position = Spacetree::positionOnGeometry();
        }
      }
    }
    DEBUG ( "return size = " << findContent.content().size()
                   << ", pos = " << data->position );
    return data;
  }

  if ( not cell.isLeaf() ) {
    DEBUG ( "Node..." );
    int searchCount = 0;
    std::shared_ptr<SearchContentResult> data ( new SearchContentResult() );
    Eigen::VectorXd childCenter(cellCenter.size());
    Eigen::VectorXd childHalflengths(cellCenter.size());
    for ( int i=0; i < cell.getChildCount(); i++ ){
      cell.getChildData(i, cellCenter, cellHalflengths, childCenter, childHalflengths);
      if ( isOverlapped(childCenter, childHalflengths,
           findContent.getVoxelCenter(), findContent.getVoxelHalflengths()) )
      {
        std::shared_ptr<SearchContentResult> tempData ( new SearchContentResult() );
        searchCount++;
        CELL_T& childCell = cell.child(i);
        tempData = searchContentInternal ( childCell, findContent, childCenter,
                                           childHalflengths );
        if ( (tempData->position != Spacetree::positionUndefined())
             && (data->position == Spacetree::positionUndefined()) )
        {
          data->position = tempData->position;
        }
        else if ( tempData->position == Spacetree::positionOnGeometry() ){
          data->position = Spacetree::positionOnGeometry();
        }
//#       ifdef Asserts
//        else if ( (tempData->position != Spacetree::positionUndefined())
//                  && (data->position != Spacetree::positionOnGeometry()) )
//        {
//          // If the position is already inside or outside, it cannot be the
//          // opposite at the same time due to the shared boundary of the
//          // spacetree cells.
//          assertion ( tempData->position == data->position );
//        }
//#       endif // Asserts
      }
    }
//    if ( searchCount > 1 ) {
//      DEBUG ( "Merging found content of size = " << findContent.content().size() );
//      mesh::Merge mergeContent;
//      mergeContent ( findContent.content() );
//      findContent.content() = mergeContent.content();
//      DEBUG ( "Merged size = " << findContent.content().size() );
//    }
    return data;
  }
  preciceError("findContentInternal()", "Reached invalid code location!");
}

template<typename CELL_T>
double StaticTraversal<CELL_T>:: distanceToBoundary
(
  const Eigen::VectorXd& cellCenter,
  const Eigen::VectorXd& cellHalflengths,
  const Eigen::VectorXd& searchPoint ) const
{
  Eigen::VectorXd centerDiff ( cellCenter );
  centerDiff -= searchPoint;
  centerDiff = centerDiff.cwiseAbs();
  double maxDistanceToCenter = centerDiff.maxCoeff();
  assertion ( math::greaterEquals(cellHalflengths[0], maxDistanceToCenter) );
  return cellHalflengths[0] - maxDistanceToCenter;
}

template<typename CELL_T>
bool StaticTraversal<CELL_T>:: isCovered
(
  const Eigen::VectorXd& cellCenter,
  const Eigen::VectorXd& cellHalflengths,
  const Eigen::VectorXd& voxelCenter,
  const Eigen::VectorXd& voxelHalflengths ) const
{
  preciceTrace ( "isCovered()", cellCenter, cellHalflengths, voxelCenter,
                  voxelHalflengths );
  Eigen::VectorXd coverage ( cellCenter );
  coverage -= voxelCenter;
  coverage = coverage.cwiseAbs();
  coverage += cellHalflengths;
  DEBUG ( "return " << not math::oneGreater(coverage, voxelHalflengths) );
  return not math::oneGreater(coverage, voxelHalflengths);
}

template<typename CELL_T>
bool StaticTraversal<CELL_T>:: isOverlapped
(
  const Eigen::VectorXd& cellCenter,
  const Eigen::VectorXd& cellHalflengths,
  const Eigen::VectorXd& voxelCenter,
  const Eigen::VectorXd& voxelHalflengths ) const
{
  preciceTrace ( "isOverlapped()", cellCenter, cellHalflengths, voxelCenter,
                  voxelHalflengths );
  Eigen::VectorXd overlap = cellCenter;
  overlap -= voxelCenter;
  overlap = overlap.cwiseAbs();
  overlap -= cellHalflengths;
  DEBUG ( "return = " << math::allGreater(voxelHalflengths, overlap) );
  return math::allGreater(voxelHalflengths, overlap);
}

template<typename CELL_T>
template<typename VISITOR_T>
void StaticTraversal<CELL_T>:: visitAllCells
(
  CELL_T&    cell,
  VISITOR_T& visitor ) const
{
  //preciceTrace ( "visitAllCells()" );
  if ( cell.isLeaf() ){
    //DEBUG ( "  Applying visitor to leaf cell content with size = "
    //               << cell.content().size() );
    visitor(cell.content());
  }
  else {
    for ( int i=0; i < cell.getChildCount(); i++ ){
      visitAllCells ( cell.child(i), visitor );
    }
  }
}

template<typename CELL_T>
template<typename VISITOR_T>
void StaticTraversal<CELL_T>:: visitRemainingCells
(
  const CELL_T& toExclude,
  CELL_T&       cell,
  VISITOR_T&    visitor ) const
{
  assertion ( not cell.isLeaf() );
  for ( int i=0; i < cell.getChildCount(); i++ ){
    if ( &toExclude != &cell.child(i) ) {
      visitAllCells ( cell.child(i), visitor );
    }
  }
}

}}} // namespace precice, spacetree, impl

