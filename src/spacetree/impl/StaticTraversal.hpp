// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_NEWSPACETREE_STATICTRAVERSAL_HPP_
#define PRECICE_NEWSPACETREE_STATICTRAVERSAL_HPP_

#include "spacetree/Spacetree.hpp"
#include "spacetree/impl/Environment.hpp"
#include "mesh/Merge.hpp"
#include "query/FindClosest.hpp"
#include "query/FindVoxelContent.hpp"
#include "utils/PointerVector.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/logging/Log.h"
#include "utils/Helpers.hpp"
#include <list>
#include <memory>

namespace precice {
namespace spacetree {
namespace impl {

template<typename CELL_T>
class StaticTraversal
{
public:

  void refineAll (
    CELL_T&                 cell,
    const utils::DynVector& cellCenter,
    const utils::DynVector& cellHalflengths,
    double                  refinementLimit,
    Environment&            environment );

  int searchPosition (
    CELL_T&                 cell,
    const utils::DynVector& searchPoint,
    const utils::DynVector& cellCenter,
    const utils::DynVector& cellHalflengths );

  bool searchDistance (
    CELL_T&                 cell,
    query::FindClosest&     findClosest,
    const utils::DynVector& cellCenter,
    const utils::DynVector& cellHalflengths );

  int searchContent (
    CELL_T&                  cell,
    query::FindVoxelContent& findContent,
    const utils::DynVector&  cellCenter,
    const utils::DynVector&  cellHalflengths );

private:

  struct RefineAllResult
  {
    int position;
    std::list<CELL_T*> cells;
    std::list<utils::DynVector> cellCenters;
    std::list<utils::DynVector> cellHalflengths;
    std::list<Environment> cellEnvironments;

    RefineAllResult()
    : position(Spacetree::positionUndefined()), cells(),
      cellCenters(), cellHalflengths(), cellEnvironments() {}
  };

  struct SearchPositionResult
  {
    int position;
    bool ambiguous;
    std::list<CELL_T*> uncachedCells;
    std::list<utils::DynVector> uncachedCellCenters;
    query::FindClosest findClosest;

    SearchPositionResult ( const utils::DynVector& searchPoint )
    : position(Spacetree::positionUndefined()), ambiguous(false), uncachedCells(),
      uncachedCellCenters(), findClosest(searchPoint) {}
  };

  struct SearchContentResult
  {
    int position;
    std::list<CELL_T*> uncachedCells;
    std::list<utils::DynVector> uncachedCellCenters;

    SearchContentResult()
    : position(Spacetree::positionUndefined()), uncachedCells(),
      uncachedCellCenters() {}
  };

  static tarch::logging::Log _log;

  void refineUndefinedCells (
    RefineAllResult&        refineAllResult,
    CELL_T&                 cell,
    const utils::DynVector& cellCenter,
    const utils::DynVector& cellHalflengths,
    double                  refinementLimit );

  void refineAllInternal (
    CELL_T&                 cell,
    const utils::DynVector& cellCenter,
    const utils::DynVector& cellHalflengths,
    double                  refinementLimit,
    Environment&            environment,
    RefineAllResult&        result );

  std::shared_ptr<SearchPositionResult> searchPositionInternal (
    CELL_T&                 cell,
    const utils::DynVector& searchPoint,
    const utils::DynVector& cellCenter,
    const utils::DynVector& cellHalflengths );

  std::shared_ptr<SearchContentResult> searchContentInternal (
    CELL_T&                  cell,
    query::FindVoxelContent& findContent,
    const utils::DynVector&  cellCenter,
    const utils::DynVector&  cellHalflengths );

  double distanceToBoundary (
    const utils::DynVector& cellCenter,
    const utils::DynVector& cellHalflengths,
    const utils::DynVector& searchPoint ) const;

  bool isCovered (
    const utils::DynVector& cellCenter,
    const utils::DynVector& cellHalflengths,
    const utils::DynVector& voxelCenter,
    const utils::DynVector& voxelHalflengths ) const;

  bool isOverlapped (
    const utils::DynVector& cellCenter,
    const utils::DynVector& cellHalflengths,
    const utils::DynVector& voxelCenter,
    const utils::DynVector& voxelHalflengths ) const;

  template<typename VISITOR_T>
  void visitAllCells (
    CELL_T&    cell,
    VISITOR_T& visitor );

  template<typename VISITOR_T>
  void visitRemainingCells (
    const CELL_T& toExclude,
    CELL_T&       cell,
    VISITOR_T&    visitor );
};

// ----------------------------------------------------- HEADER IMPLEMENTATIONS

template<typename CELL_T>
tarch::logging::Log StaticTraversal<CELL_T>::
  _log("precice::spacetree::impl::StaticTraversal");

template<typename CELL_T>
void StaticTraversal<CELL_T>:: refineAll
(
  CELL_T&                 cell,
  const utils::DynVector& cellCenter,
  const utils::DynVector& cellHalflengths,
  double                  refinementLimit,
  Environment&            env )
{
  preciceTrace3("refineAll()", cellCenter, cellHalflengths, refinementLimit);
  // The environment gives information on the position of the cells surrounding
  // a current cell of consideration. Since a mixture of in- and out-cells is
  // not possible, the 0th component of the environment vector is used to
  // indicate, whether there are only on-geometry, out-geometry (and on), or
  // in-geometry (and on) cells surrounding.
  utils::DynVector outsidePoint(cellHalflengths);
  outsidePoint *= 2.0;
  outsidePoint += cellCenter;
  query::FindClosest findClosest(outsidePoint);
  findClosest(cell.content());
  assertion(not tarch::la::equals(findClosest.getClosest().distance, 0.0));
  int pos = findClosest.getClosest().distance > 0
            ? Spacetree::positionOutsideOfGeometry()
            : Spacetree::positionInsideOfGeometry();
  env.setAllNeighborCellPositions(pos);
  env.computePosition();
  RefineAllResult result;
  refineAllInternal(cell, cellCenter, cellHalflengths, refinementLimit, env, result);
  refineUndefinedCells(result, cell, cellCenter, cellHalflengths, refinementLimit);
}

template<typename CELL_T>
void StaticTraversal<CELL_T>:: refineAllInternal
(
  CELL_T&                 cell,
  const utils::DynVector& cellCenter,
  const utils::DynVector& cellHalflengths,
  double                  refinementLimit,
  Environment&            env,
  RefineAllResult&        result )
{
  preciceTrace4 ( "refineAllInternal()", cellCenter, cellHalflengths,
                  refinementLimit, env.getNeighborCellPositions() );
  using namespace tarch::la;
  bool environmentIncomplete = false;
  if (cell.isLeaf()){
    preciceDebug("  Leaf");
    if (cell.needsRefinement(cellHalflengths, refinementLimit)){
      preciceDebug("    Needs refinement");
      assertion(not cell.content().empty());
      assertion(cell.getPosition() == Spacetree::positionOnGeometry());
      cell.refine(cellCenter, cellHalflengths);

      Environment oldEnvironment(env);
      for (int i=0; i < cell.getChildCount(); i++){
        CELL_T& childCell = cell.child(i);
        if (childCell.getPosition() == Spacetree::positionUndefined()){
          // Modify environment positions
          const DynamicVector<int>& cellIndices = env.getNeighborCellIndices(i);
          const DynamicVector<int>& sideIndices = env.getNeighborSideIndices(i);
          assertion2(cellIndices.size() == sideIndices.size(),
                     cellIndices.size(), sideIndices.size());
          for (int j=0; j < (int)cellIndices.size(); j++){
            env.setNeighborCellPosition(
                sideIndices[j], cell.child(cellIndices[j]).getPosition());
          }
          env.computePosition();
          assertion(env.getPosition() != Spacetree::positionUndefined());
          if (env.getPosition() != Spacetree::positionOnGeometry()){
            preciceDebug("    Derive cell position " << env.getPosition()
                         << " from environment = " << env.getNeighborCellPositions());
            // If some of the surrounding cells are either outside or inside,
            // the new empty cell has to be also outside or inside respectively.
            childCell.setPosition(env.getPosition());
          }
          else {
            preciceDebug("    Environment incomplete to derive position");
            environmentIncomplete = true;
          }
          env = oldEnvironment;
        }
      }
    }
  }

  if ( environmentIncomplete ){
    preciceDebug ( "  Incomplete environment, storing cell" );
    result.cells.push_back(&cell);
    result.cellCenters.push_back(cellCenter);
    result.cellHalflengths.push_back(cellHalflengths);
    result.cellEnvironments.push_back(env);
  }
  else if ( not cell.isLeaf() ){
    preciceDebug ( "  Node" );
    assertion ( cell.getPosition() != Spacetree::positionUndefined() );
    utils::DynVector childCenter(cellCenter.size());
    utils::DynVector childHalflengths(cellCenter.size());

    Environment oldEnvironment(env);
    for ( int i=0; i < cell.getChildCount(); i++ ){
      cell.getChildData(i, cellCenter, cellHalflengths, childCenter, childHalflengths);
      CELL_T& childCell = cell.child(i);
      // Modify environment positions
      const DynamicVector<int>& cellIndices = env.getNeighborCellIndices(i);
      const DynamicVector<int>& sideIndices = env.getNeighborSideIndices(i);
      assertion2 ( cellIndices.size() == sideIndices.size(),
                   cellIndices.size(), sideIndices.size() );
      for (int j=0; j < (int)cellIndices.size(); j++){
        env.setNeighborCellPosition(
            sideIndices[j], cell.child(cellIndices[j]).getPosition());
      }
      env.computePosition();
      refineAllInternal ( childCell, childCenter, childHalflengths, refinementLimit,
                          env, result );
      env = oldEnvironment;
    }
  }
}

template<typename CELL_T>
void StaticTraversal<CELL_T>:: refineUndefinedCells
(
  RefineAllResult&        refineAllResult,
  CELL_T&                 cell,
  const utils::DynVector& cellCenter,
  const utils::DynVector& cellHalflengths,
  double                  refinementLimit )
{
  preciceTrace3("refineUndefinedCells()", cellCenter, cellHalflengths, refinementLimit);
  typename std::list<CELL_T*>::iterator cellIter;
  std::list<utils::DynVector>::iterator centerIter;
  std::list<utils::DynVector>::iterator hIter;
  std::list<Environment>::iterator envIter;
  cellIter = refineAllResult.cells.begin();
  centerIter = refineAllResult.cellCenters.begin();
  hIter = refineAllResult.cellHalflengths.begin();
  envIter = refineAllResult.cellEnvironments.begin();
  while ( cellIter != refineAllResult.cells.end() ){
    assertion(centerIter != refineAllResult.cellCenters.end());
    assertion(hIter != refineAllResult.cellHalflengths.end());
    assertion(envIter != refineAllResult.cellEnvironments.end());
    preciceDebug("  Compute child positions of cell with center = " << *centerIter
                 << ", h = " << *hIter);
    bool knowPosition = false;
    int pos = Spacetree::positionUndefined();
    for ( int i=0; i < (*cellIter)->getChildCount(); i++ ){
      preciceDebug("    Child number " << i);
      CELL_T& child = (*cellIter)->child(i);
      if (child.getPosition() == Spacetree::positionUndefined()){
        if (knowPosition){
          preciceDebug("    Know position already, position = " << pos);
          child.setPosition(pos);
        }
        else {
          preciceDebug("    Compute position by findDistance");
          query::FindClosest findDistance ( *centerIter );
          searchDistance ( cell, findDistance, cellCenter, cellHalflengths );
          assertion(not tarch::la::equals(findDistance.getEuclidianDistance(), 0.0));
          pos = findDistance.getClosest().distance > 0
                ? Spacetree::positionOutsideOfGeometry()
                : Spacetree::positionInsideOfGeometry();
          preciceDebug("    Set computed position = " << pos);
          child.setPosition(pos);
          knowPosition = true;
        }
      }
    }
    RefineAllResult result;
    preciceDebug("  Go on computing subcell positions");
    refineAllInternal ( **cellIter, *centerIter, *hIter, refinementLimit, *envIter, result );
    refineUndefinedCells ( result, cell, cellCenter, cellHalflengths, refinementLimit );
    cellIter++;
    centerIter++;
    hIter++;
    envIter++;
  }
}

template<typename CELL_T>
int StaticTraversal<CELL_T>:: searchPosition
(
  CELL_T&                 cell,
  const utils::DynVector& searchPoint,
  const utils::DynVector& cellCenter,
  const utils::DynVector& cellHalflengths )
{
  preciceTrace3 ( "searchPosition()", searchPoint, cellCenter, cellHalflengths );
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
  const utils::DynVector& cellCenter,
  const utils::DynVector& cellHalflengths )
{
  preciceTrace2 ( "searchDistance()", cellCenter, cellHalflengths );
  if ( cell.isLeaf() ){
    preciceDebug ( "  Leaf" );
    findClosest ( cell.content() );
  }

  if ( not cell.isLeaf() ){
    preciceDebug ( "  Node" );
    int childIndex = cell.getChildIndex (findClosest.getSearchPoint(), cellCenter,
                                         cellHalflengths);
    utils::DynVector newCenter(cellCenter);
    utils::DynVector newHalflengths(cellHalflengths);
    cell.getChildData (childIndex, cellCenter, cellHalflengths, newCenter,
                       newHalflengths);
    CELL_T& subtree = cell.child ( childIndex );
    bool ambiguous = searchDistance(subtree, findClosest, newCenter,
                                    newHalflengths );
    if ( ambiguous ){
      visitRemainingCells ( subtree, cell, findClosest );
    }
  }

  if ( findClosest.hasFound() ){
    double distance = distanceToBoundary(cellCenter, cellHalflengths,
                                         findClosest.getSearchPoint());
    using tarch::la::greater;
    bool isAmbiguous = greater ( findClosest.getEuclidianDistance(), distance );
    preciceDebug ( "  hasfound, return " << isAmbiguous );
    return isAmbiguous;
  }
  preciceDebug ( "  return true" );
  return true;
}

template<typename CELL_T>
int StaticTraversal<CELL_T>:: searchContent
(
  CELL_T&                  cell,
  query::FindVoxelContent& findContent,
  const utils::DynVector&  cellCenter,
  const utils::DynVector&  cellHalflengths )
{
  preciceTrace2 ( "searchContent()", findContent.getVoxelCenter(),
                  findContent.getVoxelHalflengths() );
  std::shared_ptr<SearchContentResult> result = searchContentInternal (
      cell, findContent, cellCenter, cellHalflengths );
  assertion ( result.use_count() > 0 );
  assertion ( result->uncachedCells.empty() );
  if ( result->position == Spacetree::positionUndefined() ){
    // Some/all searched cells had content, but not in the search voxel
    if ( result->position == Spacetree::positionUndefined() ){
      preciceDebug ( "Computing position of search voxel" );
      query::FindClosest findDistance ( findContent.getVoxelCenter() );
      searchDistance ( cell, findDistance, cellCenter, cellHalflengths );
      double distance = findDistance.getClosest().distance;
      assertion ( not tarch::la::equals(distance, 0.0) );
      result->position = distance > 0 ? Spacetree::positionOutsideOfGeometry()
                                      : Spacetree::positionInsideOfGeometry();
    }
  }

  assertion ( result->position != Spacetree::positionUndefined() );
  preciceDebug ( "return content().size() = " << findContent.content().size() );
  return result->position;
}

template<typename CELL_T>
std::shared_ptr<typename StaticTraversal<CELL_T>::SearchPositionResult>
StaticTraversal<CELL_T>:: searchPositionInternal
(
  CELL_T&                 cell,
  const utils::DynVector& searchPoint,
  const utils::DynVector& cellCenter,
  const utils::DynVector& cellHalflengths )
{
  preciceTrace3 ( "searchPositionInternal()", searchPoint, cellCenter,
                  cellHalflengths );
  using namespace tarch::la;
  std::shared_ptr<SearchPositionResult> data;
  double distance = 0.0;
  if ( cell.isLeaf() ){
    preciceDebug ( "  Leaf" );
    assertion ( data.use_count() == 0 );
    data = std::shared_ptr<SearchPositionResult> (
        new SearchPositionResult(searchPoint) );
    if ( cell.getPosition() == Spacetree::positionOnGeometry() ){
      preciceDebug ( "    Has content" );
      //query::FindClosest findClosest ( searchPoint );
      data->findClosest ( cell.content() );
      if ( data->findClosest.hasFound() ){
        preciceDebug ( "    Found elements" );
        distance = data->findClosest.getClosest().distance;
        data->position = Spacetree::positionOnGeometry(); // May by altered later
      }
    }
    else {
      preciceDebug ( "    Is empty" );
      assertion ( cell.content().empty() );
      assertion ( cell.getPosition() != Spacetree::positionUndefined() );
      data->position = cell.getPosition();
    }
  }

  if ( not cell.isLeaf() ) { // could be a leaf on entrance to searchPositionInternal
    preciceDebug ( "  Node" );
    assertion ( cell.getPosition() == Spacetree::positionOnGeometry() );
    int childIndex = cell.getChildIndex(searchPoint, cellCenter, cellHalflengths);
    utils::DynVector newCenter(cellCenter);
    utils::DynVector newHalflengths(cellHalflengths);
    cell.getChildData (childIndex, cellCenter, cellHalflengths, newCenter,
                       newHalflengths);
    CELL_T& subtree = cell.child(childIndex);
    assertion ( data.use_count() == 0 );
    data = searchPositionInternal ( subtree, searchPoint, newCenter, newHalflengths );
    if ( (data->position == Spacetree::positionUndefined()) || data->ambiguous ){
      preciceDebug ( "    Did not find elements or ambiguous, visit others" );
      visitRemainingCells(subtree, cell, data->findClosest);
      if ( data->findClosest.hasFound() ){
        preciceDebug ( "    Found elements in others" );
        distance = data->findClosest.getClosest().distance;
        data->position = Spacetree::positionOnGeometry();
        data->ambiguous = false;
      }
    }
  }

  // Set inside/outside and check for ambiguities
  if ( not equals(distance, 0.0) ){
    preciceDebug( "  Checking for ambiguities of found objects");
    assertion ( data.use_count() > 0 );
    if ( greater(distance, 0.0) ){
      data->position = Spacetree::positionOutsideOfGeometry();
    }
    else if ( tarch::la::greater(0.0, distance) ){
      data->position = Spacetree::positionInsideOfGeometry();
    }
    preciceDebug ( "  found pos = " << data->position );
    double distanceToBound =  distanceToBoundary(cellCenter, cellHalflengths,
                                                 searchPoint);
    if ( greater(std::abs(distance), distanceToBound) ){
      preciceDebug ( "  is ambigious" );
      data->ambiguous = true;
    }
  }
  preciceDebug ( "  return position = " << data->position << ", ambiguous = "
                 << data->ambiguous );
  return data;
}

template<typename CELL_T>
std::shared_ptr<typename StaticTraversal<CELL_T>::SearchContentResult>
StaticTraversal<CELL_T>:: searchContentInternal
(
  CELL_T&                  cell,
  query::FindVoxelContent& findContent,
  const utils::DynVector&  cellCenter,
  const utils::DynVector&  cellHalflengths )
{
  preciceTrace4 ( "searchContentInternal()", cellCenter, cellHalflengths,
                  findContent.getVoxelCenter(), findContent.getVoxelHalflengths() );
  if ( cell.isLeaf() ){
    preciceDebug ( "Leaf..." );
    std::shared_ptr<SearchContentResult> data ( new SearchContentResult() );
    assertion ( cell.getPosition() != Spacetree::positionUndefined() );
    if ( (cell.getPosition() == Spacetree::positionOutsideOfGeometry())
         || (cell.getPosition() == Spacetree::positionInsideOfGeometry()) )
    {
      preciceDebug ( "empty, position = " << cell.getPosition() );
      data->position = cell.getPosition();
    }
    else {
      preciceDebug ( "content size = " << cell.content().size() );
      assertion ( cell.getPosition() == Spacetree::positionOnGeometry() );
      assertion ( not cell.content().empty() );
      bool set = isCovered ( cellCenter, cellHalflengths,
                             findContent.getVoxelCenter(),
                             findContent.getVoxelHalflengths() );
      set &= findContent.getBoundaryInclusion() == query::FindVoxelContent::INCLUDE_BOUNDARY;
      if ( set ){
        preciceDebug ( "Is covered by voxel, add cell content to find content" );
        findContent.content().add ( cell.content() );
        data->position = Spacetree::positionOnGeometry();
      }
      else {
        preciceDebug ( "Isn't covered by voxel, apply find content" );
        findContent ( cell.content() );
        if ( not findContent.content().empty() ){
          preciceDebug("Some content is contained in voxel");
          data->position = Spacetree::positionOnGeometry();
        }
      }
    }
    preciceDebug ( "return size = " << findContent.content().size()
                   << ", pos = " << data->position );
    return data;
  }

  if ( not cell.isLeaf() ) {
    preciceDebug ( "Node..." );
    int searchCount = 0;
    std::shared_ptr<SearchContentResult> data ( new SearchContentResult() );
    utils::DynVector childCenter(cellCenter.size());
    utils::DynVector childHalflengths(cellCenter.size());
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
//          assertion ( tempData->position == data->position );
//        }
//#       endif // Asserts
      }
    }
    return data;
  }
  preciceError("findContentInternal()", "Reached invalid code location!");
}

template<typename CELL_T>
double StaticTraversal<CELL_T>:: distanceToBoundary
(
  const utils::DynVector& cellCenter,
  const utils::DynVector& cellHalflengths,
  const utils::DynVector& searchPoint ) const
{
  utils::DynVector centerDiff ( cellCenter );
  centerDiff -= searchPoint;
  tarch::la::abs(centerDiff, centerDiff);
  double maxDistanceToCenter = tarch::la::max ( centerDiff );
  assertion ( tarch::la::greaterEquals(cellHalflengths[0], maxDistanceToCenter) );
  return cellHalflengths[0] - maxDistanceToCenter;
}

template<typename CELL_T>
bool StaticTraversal<CELL_T>:: isCovered
(
  const utils::DynVector& cellCenter,
  const utils::DynVector& cellHalflengths,
  const utils::DynVector& voxelCenter,
  const utils::DynVector& voxelHalflengths ) const
{
  preciceTrace4 ( "isCovered()", cellCenter, cellHalflengths, voxelCenter,
                  voxelHalflengths );
  utils::DynVector coverage ( cellCenter );
  coverage -= voxelCenter;
  tarch::la::abs ( coverage, coverage );
  coverage += cellHalflengths;
  preciceDebug ( "return " << not tarch::la::oneGreater(coverage, voxelHalflengths) );
  return not tarch::la::oneGreater(coverage, voxelHalflengths);
}

template<typename CELL_T>
bool StaticTraversal<CELL_T>:: isOverlapped
(
  const utils::DynVector& cellCenter,
  const utils::DynVector& cellHalflengths,
  const utils::DynVector& voxelCenter,
  const utils::DynVector& voxelHalflengths ) const
{
  preciceTrace4 ( "isOverlapped()", cellCenter, cellHalflengths, voxelCenter,
                  voxelHalflengths );
  utils::DynVector overlap = cellCenter;
  overlap -= voxelCenter;
  tarch::la::abs ( overlap, overlap );
  overlap -= cellHalflengths;
  preciceDebug ( "return = " << tarch::la::allGreater(voxelHalflengths, overlap) );
  return tarch::la::allGreater(voxelHalflengths, overlap);
}

template<typename CELL_T>
template<typename VISITOR_T>
void StaticTraversal<CELL_T>:: visitAllCells
(
  CELL_T&    cell,
  VISITOR_T& visitor )
{
  if ( cell.isLeaf() ){
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
  VISITOR_T&    visitor )
{
  assertion ( not cell.isLeaf() );
  for ( int i=0; i < cell.getChildCount(); i++ ){
    if ( &toExclude != &cell.child(i) ) {
      visitAllCells ( cell.child(i), visitor );
    }
  }
}

}}} // namespace precice, spacetree, impl

#endif /* PRECICE_NEWSPACETREE_STATICTRAVERSAL_HPP_ */
