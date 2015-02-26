// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "BroydenPostProcessing.hpp"
#include "cplscheme/CouplingData.hpp"
#include "utils/Globals.hpp"
#include "tarch/la/GramSchmidt.h"
#include "tarch/la/LUDecomposition.h"
#include "tarch/la/MatrixVectorOperations.h"
#include "tarch/la/TransposedMatrix.h"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/la/Scalar.h"
#include "io/TXTWriter.hpp"
#include "io/TXTReader.hpp"

#include <time.h>
#include <sstream>
#include <fstream>
//#include "utils/NumericalCompare.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

// tarch::logging::Log BroydenPostProcessing::
 //      _log("precice::cplscheme::impl::BroydenPostProcessing");

      
BroydenPostProcessing:: BroydenPostProcessing
(
  double initialRelaxation,
  int    maxIterationsUsed,
  int    timestepsReused,
  double singularityLimit,
  std::vector<int> dataIDs,
  std::map<int,double> scalings)
:
  BaseQNPostProcessing(initialRelaxation, maxIterationsUsed, timestepsReused,
		       singularityLimit, dataIDs, scalings),
//  _secondaryOldXTildes(),
  _invJacobian(),
  k(0),
  t(0)
//  _matrixWriter("jacobain.m")
{}



void BroydenPostProcessing:: initialize
(
  DataMap& cplData )
{
  // do common QN post processing initialization
  BaseQNPostProcessing::initialize(cplData);
  
  double init = 0.0;
  size_t entries= _residuals.size();
  
  _invJacobian = Matrix(entries, entries, init);
}



void BroydenPostProcessing::computeUnderrelaxationSecondaryData
(
  DataMap& cplData)
{
    //Store x_tildes for secondary data
  //  foreach (int id, _secondaryDataIDs){
  //    assertion2(_secondaryOldXTildes[id].size() == cplData[id]->values->size(),
  //               _secondaryOldXTildes[id].size(), cplData[id]->values->size());
  //    _secondaryOldXTildes[id] = *(cplData[id]->values);
  //  }

    // Perform underrelaxation with initial relaxation factor for secondary data
    foreach (int id, _secondaryDataIDs){
      PtrCouplingData data = cplData[id];
      DataValues& values = *(data->values);
      values *= _initialRelaxation;                   // new * omg
      DataValues& secResiduals = _secondaryResiduals[id];
      secResiduals = data->oldValues.column(0);    // old
      secResiduals *= 1.0 - _initialRelaxation;       // (1-omg) * old
      values += secResiduals;                      // (1-omg) * old + new * omg
    }
}




void BroydenPostProcessing::updateDifferenceMatrices
(
  DataMap& cplData)
{
  using namespace tarch::la;

//   // Compute residuals of secondary data
//   foreach (int id, _secondaryDataIDs){
//     DataValues& secResiduals = _secondaryResiduals[id];
//     PtrCouplingData data = cplData[id];
//     assertion2(secResiduals.size() == data->values->size(),
//                secResiduals.size(), data->values->size());
//     secResiduals = *(data->values);
//     secResiduals -= data->oldValues.column(0);
//   }

  /*
   * ATTETION: changed the condition from _firstIteration && _firstTimeStep
   * to the following: 
   * underrelaxation has to be done, if the scheme has converged without even
   * entering post processing. In this case the V, W matrices would still be empty.
   * This case happended in the open foam example beamInCrossFlow.
   */ 
  if(_firstIteration && (_firstTimeStep ||  (_matrixCols.size() < 2))){
    k++;
    // Perform underrelaxation with initial relaxation factor for secondary data
//     foreach (int id, _secondaryDataIDs){
//       PtrCouplingData data = cplData[id];
//       DataValues& values = *(data->values);
//       values *= _initialRelaxation;                   // new * omg
//       DataValues& secResiduals = _secondaryResiduals[id];
//       secResiduals = data->oldValues.column(0);    // old
//       secResiduals *= 1.0 - _initialRelaxation;       // (1-omg) * old
//       values += secResiduals;                      // (1-omg) * old + new * omg
//     }
  }
  else {
    if (not _firstIteration){
      k++;
    }
  }
  
  // call the base method for common update of V, W matrices
  BaseQNPostProcessing::updateDifferenceMatrices(cplData);
}



void BroydenPostProcessing::performPPSecondaryData
(
  DataMap& cplData)
{}


void BroydenPostProcessing::computeQNUpdate
    (PostProcessing::DataMap& cplData, DataValues& xUpdate)
{
  preciceTrace("computeQNUpdate()");
  using namespace tarch::la;
  
    // ------------- update inverse Jacobian -----------
    // J_inv = J_inv_n + (w- J_inv_n*v)*v^T/|v|_l2
    // ----------------------------------------- -------
  
  DataValues v = _matrixV.column(0);
  DataValues w = _matrixW.column(0);
  Matrix JUpdate(_invJacobian.rows(),_invJacobian.cols(), 0.);

  preciceDebug("took latest column of V,W");

  double dotproductV = v*v;
  DataValues tmp = _invJacobian * v;    // J_inv*v
  tmp = w - tmp;                        // (w-J_inv*v)
  tmp = tmp/dotproductV;                // (w-J_inv*v)/|v|_l2
  preciceDebug("did step (W-J_inv*v)/|v|");
  //preciceDebug("v.size = "<<v.size());
  //preciceDebug("tmp.size() = "<<tmp.size());
  Matrix tmp_mat(tmp.size(),1);
  Matrix vT_mat(1,v.size());
  //preciceDebug("tmp_mat.size = "<<tmp_mat.rows()<<", "<<tmp_mat.cols());
  //preciceDebug("vT_mat.size = "<<vT_mat.rows()<<", "<<vT_mat.cols());
  for(int i = 0; i < v.size(); i++)      // transform vectors in matrices
  {
    tmp_mat(i,0) = tmp(i);
    vT_mat(0,i) = v(i);
  }
  preciceDebug("converted vectors into matrices");
  assertion2(tmp_mat.cols() == vT_mat.rows(), tmp_mat.cols(), vT_mat.rows());
  multiply(tmp_mat, vT_mat, JUpdate);   // (w-J_inv*v)/|v| * v^T
  preciceDebug("multiplied (w-J_inv*v)/|v| * v^T");
  assertion2(_invJacobian.rows() == JUpdate.rows(), _invJacobian.rows(), JUpdate.rows());
  _invJacobian = _invJacobian + JUpdate;
  
  DataValues negRes(_residuals);
  negRes *= -1.;
 
  // solve delta_x = - J_inv*residuals
  multiply(_invJacobian, negRes, xUpdate);   
}



void BroydenPostProcessing:: specializedIterationsConverged
(
   DataMap & cplData)
{
  
//   // ---- DEBUG --------------------------
//   
//   using namespace tarch::la;
//   // compute frobenius norm of difference between Jacobian matrix from current
//   // time step and Jcobian from old time step
//   DataMatrix jacobianDiff(_invJacobian.rows(), _invJacobian.cols(), 0.0);
//   jacobianDiff = _oldInvJacobian;
//   jacobianDiff *= -1.;
//   jacobianDiff = _invJacobian + jacobianDiff;
//   double frob = frobeniusNorm(jacobianDiff); 
//   f<<t<<"  "<<frob<<"\n";
//   if(t >= 100) f.close();
//   // -------------------------------------

   //----------------------------DEBUG------------------------------------
//   std::stringstream sk; sk <<k;
//   std::stringstream st; st <<t;
//   std::string jfile("j_"+st.str()+".m");
//   _invJacobian.printm(jfile.c_str());
  
  // --------------------------------------------------------------------
  
  k = 0;
  t++;
  // store inverse Jacobian
//  _matrixWriter.write(_invJacobian);
  //_oldInvJacobian = _invJacobian;
}

}}} // namespace precice, cplscheme, impl
