// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_NEWSPACETREE_DYNAMICTRAVERSAL_HPP_
#define PRECICE_NEWSPACETREE_DYNAMICTRAVERSAL_HPP_

#include "spacetree/Spacetree.hpp"
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
class DynamicTraversal
{
public:

  int searchPosition (
    CELL_T&                 cell,
    const utils::DynVector& searchPoint,
    const utils::DynVector& cellCenter,
    const utils::DynVector& cellHalflengths,
    double                  refinementLimit );

  bool searchDistance (
    CELL_T&                 cell,
    query::FindClosest&     findClosest,
    const utils::DynVector& cellCenter,
    const utils::DynVector& cellHalflengths,
    double                  refinementLimit );

  int searchContent (
    CELL_T&                  cell,
    query::FindVoxelContent& findContent,
    const utils::DynVector&  cellCenter,
    const utils::DynVector&  cellHalflengths,
    double                   refinementLimit );

private:

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

  std::shared_ptr<SearchPositionResult> searchPositionInternal (
    CELL_T&                 cell,
    const utils::DynVector& searchPoint,
    const utils::DynVector& cellCenter,
    const utils::DynVector& cellHalflengths,
    double                  refinementLimit );

  std::shared_ptr<SearchContentResult> searchContentInternal (
    CELL_T&                  cell,
    query::FindVoxelContent& findContent,
    const utils::DynVector&  cellCenter,
    const utils::DynVector&  cellHalflengths,
    double                   refinementLimit );

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
tarch::logging::Log DynamicTraversal<CELL_T>::
  _log("precice::spacetree::impl::DynamicTraversal");

template<typename CELL_T>
int DynamicTraversal<CELL_T>:: searchPosition
(
  CELL_T&                 cell,
  const utils::DynVector& searchPoint,
  const utils::DynVector& cellCenter,
  const utils::DynVector& cellHalflengths,
  double                  refinementLimit )
{
  preciceTrace4("searchPosition()", searchPoint, cellCenter, cellHalflengths,
                refinementLimit);
  typedef std::shared_ptr<SearchPositionResult> PtrResult;
  PtrResult result = searchPositionInternal(cell, searchPoint, cellCenter,
                                            cellHalflengths, refinementLimit);
  assertion(result.use_count() > 0);
  assertion(result->position != Spacetree::positionUndefined());

  // Compute positions of uncached empty cells
  typename std::list<CELL_T*>::iterator cellIter;
  std::list<utils::DynVector>::iterator centerIter;
  cellIter = result->uncachedCells.begin();
  centerIter = result->uncachedCellCenters.begin();
  while (cellIter != result->uncachedCells.end()){
    assertion(centerIter != result->uncachedCellCenters.end());
    preciceDebug("Computing position of cell at center " << *centerIter);
    assertion((*cellIter)->getPosition() == Spacetree::positionUndefined());
    assertion((*cellIter)->content().empty());
    PtrResult tempResult = searchPositionInternal(cell, *centerIter, cellCenter,
                                                  cellHalflengths, refinementLimit);
    assertion(tempResult->position != Spacetree::positionUndefined());
    assertion(tempResult->position != Spacetree::positionOnGeometry());
    (*cellIter)->setPosition(tempResult->position);
    cellIter++;
    centerIter++;
  }
  preciceDebug("Return position = " << result->position);
  return result->position;
}

template<typename CELL_T>
bool DynamicTraversal<CELL_T>:: searchDistance
(
  CELL_T&                 cell,
  query::FindClosest&     findClosest,
  const utils::DynVector& cellCenter,
  const utils::DynVector& cellHalflengths,
  double                  refinementLimit )
{
  preciceTrace3 ( "searchDistance()", cellCenter, cellHalflengths, refinementLimit );
  if ( cell.isLeaf() ){
    preciceDebug ( "  Leaf" );
    if ( cell.needsRefinement(cellHalflengths, refinementLimit) ){
      preciceDebug ( "    Needs refinement" );
      assertion ( not cell.content().empty() );
      cell.refine(cellCenter, cellHalflengths);
    }
    else {
      preciceDebug ( "    Needs no refinement, apply findClosest" );
      findClosest ( cell.content() );
    }
  }

  if ( not cell.isLeaf() ) {
    preciceDebug ( "  Node" );
    int childIndex = cell.getChildIndex (findClosest.getSearchPoint(), cellCenter,
                                         cellHalflengths);
    utils::DynVector newCenter(cellCenter);
    utils::DynVector newHalflengths(cellHalflengths);
    cell.getChildData (childIndex, cellCenter, cellHalflengths, newCenter,
                       newHalflengths);
    CELL_T& subtree = cell.child ( childIndex );
    bool ambiguous = searchDistance(subtree, findClosest, newCenter,
                                    newHalflengths, refinementLimit);
    if ( ambiguous ){
      visitRemainingCells ( subtree, cell, findClosest );
    }
  }

  if ( findClosest.hasFound() ) {
    double distance = distanceToBoundary(cellCenter, cellHalflengths,
                                         findClosest.getSearchPoint());
    using tarch::la::greater;
    bool isAmbiguous = greater ( findClosest.getEuclidianDistance(), distance );
    preciceDebug ( "  hasfound, return ambiguous = " << isAmbiguous );
    return isAmbiguous;
  }
  preciceDebug ( "  return ambiguous or not found" );
  return true;
}

template<typename CELL_T>
int DynamicTraversal<CELL_T>:: searchContent
(
  CELL_T&                  cell,
  query::FindVoxelContent& findContent,
  const utils::DynVector&  cellCenter,
  const utils::DynVector&  cellHalflengths,
  double                   refinementLimit )
{
  preciceTrace2 ( "searchContent()", findContent.getVoxelCenter(),
                  findContent.getVoxelHalflengths() );
  std::shared_ptr<SearchContentResult> result = searchContentInternal (
      cell, findContent, cellCenter, cellHalflengths, refinementLimit );
  assertion ( result.use_count() > 0 );
  if ( result->position == Spacetree::positionUndefined() ){
    preciceDebug("Position undefined after searchContentInternal()");
    // Compute positions of empty cells
    typename std::list<CELL_T*>::iterator cellIter;
    std::list<utils::DynVector>::iterator centerIter;
    cellIter = result->uncachedCells.begin();
    centerIter = result->uncachedCellCenters.begin();
    while ( cellIter != result->uncachedCells.end() ){
      preciceDebug ( "Computing position of uncached cell at center " << *centerIter );
      assertion ( (*cellIter)->getPosition() == Spacetree::positionUndefined() );
      query::FindClosest findDistance ( *centerIter );
      searchDistance ( cell, findDistance, cellCenter, cellHalflengths, refinementLimit );
      double distance = findDistance.getClosest().distance;
      int pos = distance > 0 ? Spacetree::positionOutsideOfGeometry()
                             : Spacetree::positionInsideOfGeometry();
      (*cellIter)->setPosition ( pos );
      if ( result->position == Spacetree::positionUndefined() ){
        // Set only once, since all positions have to coincide for the voxel
        preciceDebug ( "Set result position to " << pos );
        result->position = pos;
      }
#     ifdef Asserts
      else {
        // If the voxel overlaps with cells that are inside and outside, that
        // would mean the voxel is on the geometry.
        assertion ( result->position == pos );
      }
#     endif // Asserts
      cellIter++;
      centerIter++;
    }

    // Some/all searched cells had content, but not in the search voxel
    if ( result->position == Spacetree::positionUndefined() ){
      preciceDebug ( "Computing position of search voxel" );
      query::FindClosest findDistance ( findContent.getVoxelCenter() );
      searchDistance ( cell, findDistance, cellCenter, cellHalflengths, refinementLimit );
      double distance = findDistance.getClosest().distance;
      assertion ( not tarch::la::equals(distance, 0.0) );
      result->position = distance > 0 ? Spacetree::positionOutsideOfGeometry()
                                      : Spacetree::positionInsideOfGeometry();
    }
  }

  assertion ( result->position != Spacetree::positionUndefined() );
  preciceDebug("return content().size() = " << findContent.content().size()
               << ", position = " << result->position);
  return result->position;
}

template<typename CELL_T>
std::shared_ptr<typename DynamicTraversal<CELL_T>::SearchPositionResult>
DynamicTraversal<CELL_T>:: searchPositionInternal
(
  CELL_T&                 cell,
  const utils::DynVector& searchPoint,
  const utils::DynVector& cellCenter,
  const utils::DynVector& cellHalflengths,
  double                  refinementLimit )
{
  preciceTrace4("searchPositionInternal()", searchPoint, cellCenter,
                cellHalflengths, refinementLimit);
  using namespace tarch::la;
  std::shared_ptr<SearchPositionResult> data;
  double distance = 0.0;
  if (cell.isLeaf()){
    preciceDebug("  Leaf");
    if (cell.needsRefinement(cellHalflengths, refinementLimit)){
      preciceDebug("    Needs refinement");
      assertion(not cell.content().empty());
      cell.refine(cellCenter, cellHalflengths);
    }
    else {
      preciceDebug("    Needs no refinement");
      assertion(data.use_count() == 0);
      data = std::shared_ptr<SearchPositionResult> (
             new SearchPositionResult(searchPoint) );
      if (cell.getPosition() == Spacetree::positionOnGeometry()){
        preciceDebug("    Has content");
        //query::FindClosest findClosest ( searchPoint );
        data->findClosest(cell.content());
        if (data->findClosest.hasFound()){
          preciceDebug("    Found elements");
          distance = data->findClosest.getClosest().distance;
          data->position = Spacetree::positionOnGeometry(); // May by altered later
        }
      }
      else {
        preciceDebug("    Is empty");
        assertion(cell.content().empty());
        data->position = cell.getPosition();
        if (data->position == Spacetree::positionUndefined()){
          preciceDebug("    Has undefined position");
          data->uncachedCells.push_back(&cell);
          data->uncachedCellCenters.push_back(cellCenter);
        }
      }
    }
  }

  if (not cell.isLeaf()){ // could be a leaf on entrance to searchPositionInternal
    preciceDebug("  Node");
    assertion(cell.getPosition() == Spacetree::positionOnGeometry());
    int childIndex = cell.getChildIndex(searchPoint, cellCenter, cellHalflengths);
    utils::DynVector newCenter(cellCenter);
    utils::DynVector newHalflengths(cellHalflengths);
    cell.getChildData(childIndex, cellCenter, cellHalflengths, newCenter,
                      newHalflengths);
    CELL_T& subtree = cell.child(childIndex);
    assertion(data.use_count() == 0);
    data = searchPositionInternal(subtree, searchPoint, newCenter, newHalflengths,
                                  refinementLimit);
    if ((data->position == Spacetree::positionUndefined()) || data->ambiguous){
      preciceDebug("    Did not find elements or ambiguous, visit others");
      visitRemainingCells(subtree, cell, data->findClosest);
      if (data->findClosest.hasFound()){
        preciceDebug("    Found elements in others");
        distance = data->findClosest.getClosest().distance;
        data->position = Spacetree::positionOnGeometry();
        data->ambiguous = false;
      }
    }
  }

  // Set inside/outside and check for ambiguities
  if (not equals(distance, 0.0)){
    preciceDebug("  Checking for ambiguities of found objects");
    assertion(data.use_count() > 0);
    if (greater(distance, 0.0)){
      data->position = Spacetree::positionOutsideOfGeometry();
    }
    else if (tarch::la::greater(0.0, distance)){
      data->position = Spacetree::positionInsideOfGeometry();
    }
    preciceDebug("  found pos = " << data->position);
    double distanceToBound =  distanceToBoundary(cellCenter, cellHalflengths,
                                                 searchPoint);
    if (greater(std::abs(distance), distanceToBound)){
      preciceDebug("  is ambigious");
      data->ambiguous = true;
    }
  }
  preciceDebug("  return position = " << data->position << ", ambiguous = "
               << data->ambiguous);
  return data;
}

template<typename CELL_T>
std::shared_ptr<typename DynamicTraversal<CELL_T>::SearchContentResult>
DynamicTraversal<CELL_T>:: searchContentInternal
(
  CELL_T&                  cell,
  query::FindVoxelContent& findContent,
  const utils::DynVector&  cellCenter,
  const utils::DynVector&  cellHalflengths,
  double                   refinementLimit )
{
  preciceTrace4 ( "searchContentInternal()", cellCenter, cellHalflengths,
                  findContent.getVoxelCenter(), findContent.getVoxelHalflengths() );
  if ( cell.isLeaf() ){
    preciceDebug ( "Leaf..." );
    if ( cell.needsRefinement(cellHalflengths, refinementLimit) ){
      preciceDebug ( "Needs refinement..." );
      cell.refine(cellCenter, cellHalflengths);
    }
    else {
      preciceDebug ( "Don't needs refinement..." );
      std::shared_ptr<SearchContentResult> data ( new SearchContentResult() );
      if ( (cell.getPosition() == Spacetree::positionOutsideOfGeometry())
           || (cell.getPosition() == Spacetree::positionInsideOfGeometry()) )
      {
        preciceDebug ( "empty, position = " << cell.getPosition() );
        data->position = cell.getPosition();
      }
      else if ( cell.getPosition() == Spacetree::positionUndefined() ) {
        preciceDebug ( "empty, undefined position" );
        data->uncachedCells.push_back(&cell);
        data->uncachedCellCenters.push_back(cellCenter);
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
          preciceDebug ( "Is covered by voxel..." );
          findContent.content().add ( cell.content() );
          data->position = Spacetree::positionOnGeometry();
        }
        else {
          preciceDebug ( "Isn't covered by voxel..." );
          findContent ( cell.content() );
          if ( not findContent.content().empty() ) {
            data->position = Spacetree::positionOnGeometry();
          }
        }
      }
      preciceDebug ( "return size = " << findContent.content().size()
                     << ", pos = " << data->position );
      return data;
    }
  }

  if ( not cell.isLeaf() ){
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
                                           childHalflengths, refinementLimit );
        if ( (tempData->position != Spacetree::positionUndefined())
             && (data->position == Spacetree::positionUndefined()) )
        {
          data->position = tempData->position;
        }
        else if ( tempData->position == Spacetree::positionOnGeometry() ){
          data->position = Spacetree::positionOnGeometry();
        }
#       ifdef Asserts
        else if ( (tempData->position != Spacetree::positionUndefined())
                  && (data->position != Spacetree::positionOnGeometry()) )
        {
          assertion ( tempData->position == data->position );
        }
#       endif // Asserts
      }
    }
    return data;
  }
  preciceError("findContentInternal()", "Reached invalid code location!");
//  assertion ( false );
//
//  std::shared_ptr<SearchContentResult> neverUsedSharedPtr( new SearchContentResult() );
//  return neverUsedSharedPtr;
}

template<typename CELL_T>
double DynamicTraversal<CELL_T>:: distanceToBoundary
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
bool DynamicTraversal<CELL_T>:: isCovered
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
bool DynamicTraversal<CELL_T>:: isOverlapped
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
void DynamicTraversal<CELL_T>:: visitAllCells
(
  CELL_T&    cell,
  VISITOR_T& visitor )
{
  //preciceTrace ( "visitAllCells()" );
  if ( cell.isLeaf() ){
    //preciceDebug ( "  Applying visitor to leaf cell content with size = "
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
void DynamicTraversal<CELL_T>:: visitRemainingCells
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

#endif /* PRECICE_NEWSPACETREE_DYNAMICTRAVERSAL_HPP_ */
