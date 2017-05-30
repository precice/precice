#include "OctreeCell.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"
#include "math/math.hpp"

namespace precice {
namespace spacetree {
namespace impl {

logging::Logger OctreeCell:: _log("spacetree::impl::OctreeCell");

OctreeCell:: OctreeCell()
:
  _content(new mesh::Group()),
  _position(Spacetree::positionUndefined()),
  _childs()
{}

OctreeCell:: ~OctreeCell()
{
  if (_content != nullptr){
    delete _content;
  }
  _childs.deleteElements();
}

bool OctreeCell:: needsRefinement
(
  const Eigen::VectorXd& cellHalflengths,
  double                 refinementLimit )
{
  assertion(_content != nullptr);
  // All halflengths are assumed to be equal
  if (math::smaller(cellHalflengths[0], refinementLimit)){
    return false;
  }
  else if ((int)_content->size() < Spacetree::minElementsToRefineCell){
    return false;
  }
  return true;
}

void OctreeCell:: refine
(
  const Eigen::VectorXd& cellCenter,
  const Eigen::VectorXd& cellHalflengths )
{
  TRACE(cellCenter, cellHalflengths);
  assertion(_content != nullptr);
  assertion(_childs.size() == 0, _childs.size());
  int dim = cellCenter.size();
  int twoPowerDim = std::pow(2.0,dim);
  //_subtrees = new RegularSpacetree*[twoPowerDim];
  Eigen::VectorXd newCenter(dim);
  Eigen::VectorXd newHalflengths = Eigen::VectorXd::Constant(dim, 0.5*cellHalflengths[0]);
  for (int i=0; i < twoPowerDim; i++){
    //DEBUG ( "Creating child " << i );
    newCenter = utils::delinearize(i, dim);
    newCenter *= newHalflengths[0];
    newCenter += cellCenter;
    _childs.push_back(new OctreeCell());
    query::FindVoxelContent findVoxel (
      newCenter, newHalflengths, query::FindVoxelContent::INCLUDE_BOUNDARY );
    findVoxel(*_content);
    _childs[i]._content->add(findVoxel.content());
    if (not _childs[i]._content->empty()){
      DEBUG("  Refined child cell with center " << newCenter
                   << " is on geometry");
//      if (_childs[i]._content->size() == 1){
//        INFO("  Cell at x=" << newCenter << ", h=" << newHalflengths
//                     << ", vertices=" << _childs[i]._content->vertices().size()
//                     << ", edges=" << _childs[i]._content->edges().size()
//                     << ", triangles=" << _childs[i]._content->triangles().size());
//      }
//      if(_childs[i]._content->vertices().size() == 0
//         && _childs[i]._content->edges().size() == 1
//         && _childs[i]._content->triangles().size() == 1)
//      {
//        INFO("Edge from " << _childs[i]._content->edges()[0].vertex(0).getCoords()
//                     << " to " << _childs[i]._content->edges()[0].vertex(1).getCoords());
//        INFO("Triangle at " << _childs[i]._content->triangles()[0].vertex(0).getCoords()
//                     << ", " << _childs[i]._content->triangles()[0].vertex(1).getCoords()
//                     << ", " << _childs[i]._content->triangles()[0].vertex(2).getCoords());
//      }
      _childs[i]._position = Spacetree::positionOnGeometry();
    }
  }
  delete _content;
  _content = nullptr; // Important to recognize the cell as node
}

int OctreeCell:: getChildIndex
(
  const Eigen::VectorXd& searchPoint,
  const Eigen::VectorXd& cellCenter,
  const Eigen::VectorXd& cellHalflengths )
{
  assertion(cellCenter.size() == cellHalflengths.size(),
            cellCenter.size(), cellHalflengths.size());
  assertion(cellCenter.size() == searchPoint.size(),
            cellCenter.size(), searchPoint.size());
  Eigen::VectorXd halfspaces(searchPoint.size());
  for (int i=0; i < searchPoint.size(); i++){
    halfspaces[i] = searchPoint[i] > cellCenter[i] ? 1.0 : -1.0;
  }
  return utils::linearize(halfspaces);
}

void OctreeCell:: getChildData
(
  int                    childIndex,
  const Eigen::VectorXd& cellCenter,
  const Eigen::VectorXd& cellHalflengths,
  Eigen::VectorXd&       childCenter,
  Eigen::VectorXd&       childHalflengths )
{
  TRACE(childIndex, cellCenter, cellHalflengths);
  childHalflengths.setConstant(0.5 * cellHalflengths[0]);
  childCenter = utils::delinearize(childIndex, childCenter.size());
  childCenter *= childHalflengths[0];
  childCenter += cellCenter;
  DEBUG("  Computed center = " << childCenter << ", h = " << childHalflengths);
}

void OctreeCell:: accept
(
  Spacetree::Visitor& visitor,
  const Eigen::VectorXd& center,
  const Eigen::VectorXd& halflengths )
{
  if (isLeaf()){
    visitor.leafCallback(center, halflengths, _position, content());
  }
  else {
    visitor.nodeCallback(center, halflengths, _position);
    Eigen::VectorXd childCenter(center.size());
    Eigen::VectorXd childHalflengths(center.size());
    for (int i=0; i < (int)_childs.size(); i++){
      getChildData(i, center, halflengths, childCenter, childHalflengths);
      _childs[i].accept(visitor, childCenter, childHalflengths);
    }
  }
}

void OctreeCell:: clear()
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
  assertion(_childs.size() == 0, _childs.size());
  assertion(_content->empty());
}

}}} // namespace precice, spacetree, impl
