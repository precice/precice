#ifndef PRECICE_CPLSCHEME_COUPLINGDATA_HPP_
#define PRECICE_CPLSCHEME_COUPLINGDATA_HPP_

#include "mesh/SharedPointer.hpp"
#include "tarch/la/DynamicColumnMatrix.h"
#include "utils/Helpers.hpp"
#include "mesh/Data.hpp"
#include "Eigen/Dense"
#include <vector>

namespace precice {
namespace cplscheme {

struct CouplingData
{
  typedef Eigen::MatrixXd DataMatrix;

  /// @brief Data values of current iteration.
  Eigen::VectorXd* values;

  /// @brief Data values of previous iteration (1st col) and previous timesteps.
  DataMatrix oldValues;

  mesh::PtrMesh mesh;

  /// @brief True, if the data values are initialized by a participant.
  bool initialize;

  /// @ dimension of one data value (scalar=1, or vectorial=interface-dimension)
  int dimension;

  /**
   * @brief Default constructor, not to be used!
   *
   * Necessary when compiler creates template code for std::map::operator[].
   */
  CouplingData ()
    {
      assertion ( false );
    }

  /**
   * @brief Constructor.
   */
  CouplingData (
    Eigen::VectorXd*  values,
    mesh::PtrMesh     mesh,
    bool              initialize,
    int               dimension)
    :
    values ( values ),
    oldValues (),
    mesh(mesh),
    initialize ( initialize ),
    dimension(dimension)
    {
      assertion ( values != NULL );
      assertion ( mesh.use_count()>0);
    }
};

}} // namespace precice, cplscheme

#endif /* PRECICE_CPLSCHEME_COUPLINGDATA_HPP_ */
