#include "Dimensions.hpp"
#include "assertion.hpp"

namespace precice {
namespace utils {

const Eigen::VectorXd DELINEARIZE_2D[8] =
    {
        Eigen::VectorXd(Eigen::Vector2d(-1.0, -1.0)),
        Eigen::VectorXd(Eigen::Vector2d(1.0, -1.0)),
        Eigen::VectorXd(Eigen::Vector2d(-1.0, 1.0)),
        Eigen::VectorXd(Eigen::Vector2d(1.0, 1.0))};

const Eigen::VectorXd DELINEARIZE_3D[8] =
    {
        Eigen::VectorXd(Eigen::Vector3d(-1.0, -1.0, -1.0)),
        Eigen::VectorXd(Eigen::Vector3d(1.0, -1.0, -1.0)),
        Eigen::VectorXd(Eigen::Vector3d(-1.0, 1.0, -1.0)),
        Eigen::VectorXd(Eigen::Vector3d(1.0, 1.0, -1.0)),
        Eigen::VectorXd(Eigen::Vector3d(-1.0, -1.0, 1.0)),
        Eigen::VectorXd(Eigen::Vector3d(1.0, -1.0, 1.0)),
        Eigen::VectorXd(Eigen::Vector3d(-1.0, 1.0, 1.0)),
        Eigen::VectorXd(Eigen::Vector3d(1.0, 1.0, 1.0))};

const Eigen::VectorXd &delinearize(
    int toDelinearize,
    int dimensions)
{
  if (dimensions == 2) {
    PRECICE_ASSERT((toDelinearize >= 0) && (toDelinearize < 4), toDelinearize);
    return DELINEARIZE_2D[toDelinearize];
  } else {
    PRECICE_ASSERT(dimensions == 3, dimensions);
    PRECICE_ASSERT((toDelinearize >= 0) && (toDelinearize < 8), toDelinearize);
    return DELINEARIZE_3D[toDelinearize];
  }
}

const int IndexMaps<2>::CUBOID_EDGE_VERTICES[4][2] =
    {
        {0, 2}, // edge 0
        {1, 3}, // edge 1
        {0, 1}, // edge 2
        {2, 3}  // edge 3
};

const int IndexMaps<3>::CUBOID_FACE_VERTICES[6][4] =
    {
        {0, 2, 4, 6},
        {1, 3, 5, 7},
        {0, 1, 4, 5},
        {2, 3, 6, 7},
        {0, 1, 2, 3},
        {4, 5, 6, 7}};

const int IndexMaps<3>::CUBOID_FACE_EDGES[6][4] =
    {
        {4, 6, 8, 10},
        {5, 7, 9, 11},
        {0, 2, 8, 9},
        {1, 3, 10, 11},
        {0, 1, 4, 5},
        {2, 3, 6, 7}};

const int IndexMaps<3>::CUBOID_EDGE_VERTICES[12][2] =
    {
        /*  0 */ {0, 1},
        /*  1 */ {2, 3},
        /*  2 */ {4, 5},
        /*  3 */ {6, 7},
        /*  4 */ {0, 2},
        /*  5 */ {1, 3},
        /*  6 */ {4, 6},
        /*  7 */ {5, 7},
        /*  8 */ {0, 4},
        /*  9 */ {1, 5},
        /* 10 */ {2, 6},
        /* 11 */ {3, 7}};

} // namespace utils
} // namespace precice
