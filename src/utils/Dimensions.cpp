// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "Dimensions.hpp"
//#include "utils/Globals.hpp"

namespace precice {
namespace utils {


//void initializeZero ( Vector & toInitialize )
//{
//   toInitialize = Vector(0.0);
//}

std::string getTypeName(const DynVector& var)
{
  return std::string("vector");
}

std::string getTypeName(const Vector2D& var)
{
  return std::string("2d-vector");
}

std::string getTypeName(const Vector3D& var)
{
  return std::string("3d-vector");
}

const utils::DynVector DELINEARIZE_2D[4] =
{
 utils::DynVector(utils::Vector2D(-1.0, -1.0)),
 utils::DynVector(utils::Vector2D( 1.0, -1.0)),
 utils::DynVector(utils::Vector2D(-1.0,  1.0)),
 utils::DynVector(utils::Vector2D( 1.0,  1.0))
};

const utils::DynVector DELINEARIZE_3D[8] =
{
 utils::DynVector(utils::Vector3D(-1.0, -1.0, -1.0)),
 utils::DynVector(utils::Vector3D( 1.0, -1.0, -1.0)),
 utils::DynVector(utils::Vector3D(-1.0,  1.0, -1.0)),
 utils::DynVector(utils::Vector3D( 1.0,  1.0, -1.0)),
 utils::DynVector(utils::Vector3D(-1.0, -1.0,  1.0)),
 utils::DynVector(utils::Vector3D( 1.0, -1.0,  1.0)),
 utils::DynVector(utils::Vector3D(-1.0,  1.0,  1.0)),
 utils::DynVector(utils::Vector3D( 1.0,  1.0,  1.0))
};

const utils::DynVector& delinearize
(
  int toDelinearize,
  int dimensions )
{
  if ( dimensions == 2 ){
    assertion1 ( (toDelinearize >= 0) && (toDelinearize < 4), toDelinearize );
    return DELINEARIZE_2D[toDelinearize];
  }
  else {
    assertion1 ( dimensions == 3, dimensions );
    assertion1 ( (toDelinearize >= 0) && (toDelinearize < 8), toDelinearize );
    return DELINEARIZE_3D[toDelinearize];
  }
}

//void delinearize
//(
//  int              index,
//  utils::Vector2D& result )
//{
//   assertion (index < 4);
//   double dblIndex = (double)index;
//   for (int dim = 1; dim >= 0 ; dim--) {
//      if ((dblIndex / std::pow(2.0, dim)) >= 1.0) {
//         result(dim) = 1.0;
//         dblIndex -= std::pow(2.0, dim);
//         assertion (dblIndex >= 0);
//      }
//      else {
//         result(dim) = -1.0;
//      }
//   }
//}
//
//void delinearize
//(
//  int              index,
//  utils::Vector3D& result )
//{
//   assertion (index < 8);
//   double dblIndex = (double)index;
//   for (int dim = 2; dim >= 0 ; dim--) {
//      if ((dblIndex / std::pow(2.0, dim)) >= 1.0) {
//         result(dim) = 1.0;
//         dblIndex -= std::pow(2.0, dim);
//         assertion (dblIndex >= 0);
//      }
//      else {
//        result(dim) = -1.0;
//      }
//   }
//}

//tarch::la::Vector<Def::DIM-1, double> delinearizeDimensionReduced
//(
//  int index )
//{
//   assertion (index < Def::TWO_POWER_DIM_M1);
//   double dblIndex = (double)index;
//   tarch::la::Vector<Def::DIM-1, double> delinearized;
//   for (int dim = Def::DIM - 2; dim >= 0 ; dim--) {
//      if ((dblIndex / std::pow(2.0, dim)) >= 1.0) {
//         delinearized(dim) = 1.0;
//         dblIndex -= std::pow(2.0, dim);
//         assertion (dblIndex >= 0);
//      }
//      else {
//         delinearized(dim) = -1.0;
//      }
//   }
//   return delinearized;
//}

//int linearize
//(
//  const utils::Vector2D& dimensions )
//{
//   int index = 0;
//   for (int dim=0; dim < 2; dim++) {
//      if (dimensions(dim) > 0.0) {
//         index += (int) std::pow(2.0, dim);
//      }
//   }
//   return index;
//}
//
//int linearize
//(
//  const utils::Vector3D& dimensions )
//{
//   int index = 0;
//   for (int dim=0; dim < 3; dim++) {
//      if (dimensions(dim) > 0.0) {
//         index += (int) std::pow(2.0, dim);
//      }
//   }
//   return index;
//}

//void getHyperfaceCornerIndices
//(
//  int faceIndex,
//  int cornerIndices[CompilePower<2,Def::DIM-1>::VALUE])
//{
//   assertion (faceIndex >= 0);
//   assertion (faceIndex < 2 * Def::DIM);
//
//   double dblIndex = faceIndex;
//   int faceNormalDimension = -1;
//   double coordinateSide = 0.0;
//   for (int dim = Def::DIM - 1; dim >= 0 ; dim--) {
//      if ((dblIndex / 2.0) >= 1.0) {
//         dblIndex -= 2.0;
//      }
//      else {
//         faceNormalDimension = Def::DIM - dim - 1;
//         coordinateSide = dblIndex < 0.5 ? -1.0 : 1.0;
//         break;
//      }
//   }
//   assertion ( faceNormalDimension >= 0 );
//   assertion ( faceNormalDimension < Def::DIM );
//
//   for (int i=0; i < CompilePower<2,Def::DIM-1>::VALUE; i++) {
//      tarch::la::Vector<Def::DIM-1, double> remainingDimensions =
//         delinearizeDimensionReduced ( i );
//      Vector vertexLocation;
//      int remainingDim = 0;
//      for (int dim=0; dim < Def::DIM; dim++) {
//         if (dim == faceNormalDimension) {
//            vertexLocation(dim) = coordinateSide;
//         }
//         else {
//            vertexLocation(dim) = remainingDimensions(remainingDim);
//            remainingDim++;
//         }
//      }
//      cornerIndices[i] = linearize (vertexLocation);
//   }
//}

const int IndexMaps<2>:: CUBOID_EDGE_VERTICES[4][2] =
{
   {0, 2}, // edge 0
   {1, 3}, // edge 1
   {0, 1}, // edge 2
   {2, 3}  // edge 3
};

const int IndexMaps<3>:: CUBOID_FACE_VERTICES[6][4] =
{
   {0, 2, 4, 6},
   {1, 3, 5, 7},
   {0, 1, 4, 5},
   {2, 3, 6, 7},
   {0, 1, 2, 3},
   {4, 5, 6 ,7}
};

const int IndexMaps<3>:: CUBOID_FACE_EDGES[6][4] =
{
   {4,  6,  8, 10},
   {5,  7,  9, 11},
   {0,  2,  8,  9},
   {1,  3, 10, 11},
   {0,  1,  4,  5},
   {2,  3,  6,  7}
};

const int IndexMaps<3>:: CUBOID_EDGE_VERTICES[12][2] =
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
   /* 11 */ {3, 7}
};

//const utils::DynVector& DELINEARIZE_3D[8] =
//{
// utils::DynVector(utils::Vector2D()),
// utils::DynVector(utils::Vector2D()),
// utils::DynVector(utils::Vector2D()),
// utils::DynVector(utils::Vector2D())
//};

}} // namespace precice, utils
