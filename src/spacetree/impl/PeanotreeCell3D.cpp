// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "PeanotreeCell3D.hpp"

namespace precice {
namespace spacetree {
namespace impl {

tarch::logging::Log PeanotreeCell3D:: _log("precice::spacetree::impl::PeanotreeCell3D");

const double PeanotreeCell3D::_oneThird(1.0/3.0);

PeanotreeCell3D:: PeanotreeCell3D()
:
  _content(new mesh::Group()),
  _position(Spacetree::positionUndefined()),
  _childs()
{
  assertion(false, "Class implementation not yet finished!!");
}

PeanotreeCell3D:: ~PeanotreeCell3D()
{
  if (_content != nullptr){
    delete _content;
  }
  _childs.deleteElements();
}

bool PeanotreeCell3D:: needsRefinement
(
  const utils::DynVector& cellHalflengths,
  double                  refinementLimit )
{
  assertion(_content != nullptr);
  // All halflengths are assumed to be equal
  if (tarch::la::smaller(cellHalflengths[0], refinementLimit)){
    return false;
  }
  if (_content->triangles().size() <= 1){
    return false;
  }
  return true;
}

void PeanotreeCell3D:: refine
(
  const utils::DynVector& cellCenter,
  const utils::DynVector& cellHalflengths )
{
  tpreciceTrace2 ( "refine()", cellCenter, cellHalflengths );
  assertion ( _content != nullptr );
  assertion ( _childs.size() == 0, _childs.size() );
  int dim = cellCenter.size();
  int threePowerD = dim==2 ? 9 : 27;  // only 2D and 3D supported
  //_subtrees = new RegularSpacetree*[twoPowerDim];
  utils::DynVector newCenter(dim);
  utils::DynVector newHalflengths(dim, 1.0 / 3.0 * cellHalflengths[0]);
  for ( int i=0; i < threePowerD; i++ ){
    tpreciceDebug ( "Creating child " << i );
    newCenter = PEANO_FINE_CELL_CENTER_POSITIONS_3D[i];
    newCenter *= newHalflengths[0];
    newCenter += cellCenter;
    _childs.push_back(new PeanotreeCell3D());
    query::FindVoxelContent findVoxel (
      newCenter, newHalflengths, query::FindVoxelContent::INCLUDE_BOUNDARY );
    findVoxel (*_content);
    _childs[i]._content->add(findVoxel.content());
    if ( not _childs[i]._content->empty() ){
      _childs[i]._position = Spacetree::positionOnGeometry();
    }
  }
  delete _content;
  _content = nullptr; // Important to recognize the cell as node
}

int PeanotreeCell3D:: getChildIndex
(
  const utils::DynVector& searchPoint,
  const utils::DynVector& cellCenter,
  const utils::DynVector& cellHalflengths )
{
  assertion ( cellCenter.size() == cellHalflengths.size(),
               cellCenter.size(), cellHalflengths.size());
  assertion ( cellCenter.size() == searchPoint.size(),
               cellCenter.size(), searchPoint.size());
  utils::DynVector halfspaces(searchPoint.size());
  for ( int i=0; i < searchPoint.size(); i++ ) {
    halfspaces[i] = searchPoint[i] > cellCenter[i] ? 1.0 : -1.0;
  }
  return utils::linearize(halfspaces);
}

void PeanotreeCell3D:: getChildData
(
  int                     childIndex,
  const utils::DynVector& cellCenter,
  const utils::DynVector& cellHalflengths,
  utils::DynVector&       childCenter,
  utils::DynVector&       childHalflengths )
{
  tpreciceTrace3 ( "getChildData()", childIndex, cellCenter, cellHalflengths );
  assign(childHalflengths) = 0.5 * cellHalflengths[0];
  childCenter = utils::delinearize(childIndex, childCenter.size());
  childCenter *= childHalflengths[0];
  childCenter += cellCenter;
}

void PeanotreeCell3D:: accept
(
  Spacetree::Visitor& visitor,
  const utils::DynVector& center,
  const utils::DynVector& halflengths )
{
  if ( isLeaf() ){
    visitor.leafCallback(center, halflengths, _position, content());
  }
  else {
    visitor.nodeCallback(center, halflengths, _position);
    utils::DynVector childCenter(center.size());
    utils::DynVector childHalflengths(center.size());
    for (int i=0; i < (int)_childs.size(); i++){
      getChildData(i, center, halflengths, childCenter, childHalflengths);
      _childs[i].accept(visitor, childCenter, childHalflengths);
    }
  }
}

void PeanotreeCell3D:: clear()
{
  if (_content != nullptr){
    _content->clear();
  }
  else {
    _content = new mesh::Group();
  }
  _childs.deleteElements();
  _childs.clear();
  _position = Spacetree::positionUndefined();
  assertion ( _childs.size() == 0, _childs.size() );
  assertion ( _content->empty() );
}

}}} // namespace precice, spacetree, impl
