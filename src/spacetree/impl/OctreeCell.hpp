#pragma once

#include "spacetree/Spacetree.hpp"
#include "mesh/Group.hpp"
#include "logging/Logger.hpp"
#include "query/FindVoxelContent.hpp"

namespace precice {
  namespace spacetree {
    namespace tests {
      class SpacetreeTestScenarios;
    }
  }
}

// ------------------------------------------------------------ CLASS DEFINITION

namespace precice {
namespace spacetree {
namespace impl {

class OctreeCell
{
public:

  OctreeCell();

  ~OctreeCell();

  bool isLeaf() const
  {
    return _content != NULL;
  }

  mesh::Group& content()
  {
    assertion(_content != NULL);
    return *_content;
  }

  bool needsRefinement (
    const Eigen::VectorXd& cellHalflengths,
    double                  refinementLimit );

  int getPosition() const
  {
    return _position;
  }

  void setPosition ( int position )
  {
    _position = position;
  }

  void refine (
    const Eigen::VectorXd& cellCenter,
    const Eigen::VectorXd& cellHalflengths );

  int getChildCount()
  {
    return _childs.size();
  }

  OctreeCell& child ( int index )
  {
    assertion ( (index >= 0) && (index < (int)_childs.size()), index, _childs.size() );
    return _childs[index];
  }

  int getChildIndex (
    const Eigen::VectorXd& searchPoint,
    const Eigen::VectorXd& cellCenter,
    const Eigen::VectorXd& cellHalflengths );

  /**
   * @brief Computes and sets child center and halflengths from parent cell.
   */
  void getChildData (
    int                    childIndex,
    const Eigen::VectorXd& cellCenter,
    const Eigen::VectorXd& cellHalflengths,
    Eigen::VectorXd&       childCenter,
    Eigen::VectorXd&       childHalflengths );

  void accept (
    Spacetree::Visitor& visitor,
    const Eigen::VectorXd& center,
    const Eigen::VectorXd& halflengths );

  void clear();

private:

  static logging::Logger _log;

  mesh::Group* _content;

  int _position;

  utils::ptr_vector<OctreeCell> _childs;
};

}}} // namespace precice, spacetree, impl

