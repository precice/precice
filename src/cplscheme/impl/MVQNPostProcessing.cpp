// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "MVQNPostProcessing.hpp"
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
#include "Eigen/Dense"

#include <time.h>
#include <sstream>
#include <fstream>
//#include "utils/NumericalCompare.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

// tarch::logging::Log MVQNPostProcessing::
//       _log("precice::cplscheme::impl::MVQNPostProcessing");

      
MVQNPostProcessing:: MVQNPostProcessing
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
  _oldInvJacobian(),
  k(0),
  t(0)
{}



void MVQNPostProcessing:: initialize
(
  DataMap& cplData )
{
  // do common QN post processing initialization
  BaseQNPostProcessing::initialize(cplData);
  
  double init = 0.0;
  size_t entries= _residuals.size();
  
  _invJacobian = Matrix(entries, entries, init);
  _oldInvJacobian = Matrix(entries, entries, init);
}



void MVQNPostProcessing::computeUnderrelaxationSecondaryData
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




void MVQNPostProcessing::updateDifferenceMatrices
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



void MVQNPostProcessing::computeQNUpdate
    (PostProcessing::DataMap& cplData, DataValues& xUpdate)
{
  preciceTrace("computeQNUpdate()");
  using namespace tarch::la;

    // ------------- update inverse Jacobian -----------
    // J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
    // ----------------------------------------- -------

    preciceDebug("   Compute Newton factors");
    //computeNewtonFactorsLUDecomposition(cplData, xUpdate);
    
    DataValues xUpdate2(xUpdate.size(),0.0);
    
    // computes xUpdate using updatedQR decompositon (does not modify _invJacobian)
    computeNewtonFactorsUpdatedQRDecomposition(cplData, xUpdate2);
    // computes xUpdate using modifiedGramSchmidt QR-dec
    computeNewtonFactorsQRDecomposition(cplData,xUpdate);

    // validation
    for(int i = 0; i<xUpdate.size(); i++)
    {
      if(!tarch::la::equals(xUpdate(i), xUpdate2(i), 1e-12))
      {
	std::cerr<<"xUpdates were not the same for standart and updatedQR decomposition:\n"<<"standart:\n   "<<xUpdate<<"updated:\n   "<<xUpdate2<<std::endl;
	_infostream<<"validation failed for xUpdate.\n   - modifiedGramSchmidt:\n"<<xUpdate<<"\n  - updatedQR:\n"<<xUpdate2<<"\n"<<std::flush;
      }
    }
}


void MVQNPostProcessing::computeNewtonFactorsUpdatedQRDecomposition
(PostProcessing::DataMap& cplData, DataValues& xUpdate)
{
  preciceTrace("computeNewtonFactorsQRDecomposition()");
  using namespace tarch::la;
 
  // ------------- update inverse Jacobian -----------
  // J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
  // ----------------------------------------- -------

  DataMatrix v;
  bool linearDependence = true;
  
  while (linearDependence)
  {
    linearDependence = false;
    v.clear();
    
    Matrix __R(_matrixV.cols(), _matrixV.cols(), 0.0);
    auto r = _qrV.matrixR();
    for(int i = 0; i<r.rows(); i++)
      for(int j = 0; j<r.cols(); j++)
      {
	__R(i,j) = r(i,j);
      }

    if (_matrixV.cols() > 1){
      for (int i=0; i < _matrixV.cols(); i++){
	if (__R(i,i) < _singularityLimit){
	  preciceDebug("   Removing linear dependent column " << i);
	  _infostream<<"(updatedQR) removing linear dependent column "<<i<<"  time step: "<<t<<" iteration: "<<k<<"\n"<<std::flush;
	  linearDependence = true;
	  removeMatrixColumn(i);
	}
      }
    }
    if(not linearDependence)
    {
      Matrix __Q(_matrixV.rows(), _matrixV.cols(), 0.0);
      
      DataValues __ytmpVec(_matrixV.cols(), 0.0);
      DataValues __matrixQRow;
      auto q = _qrV.matrixQ();
      for(int i = 0; i<q.rows(); i++)
	for(int j = 0; j<q.cols(); j++)
	{
	  __Q(i,j) = q(i,j);
	}
	
      for(int i = 0; i < __Q.rows(); i++)
      {
	for(int j=0; j < __Q.cols(); j++){
	  __matrixQRow.append(__Q(i,j));
	}
	
	backSubstitution(__R, __matrixQRow, __ytmpVec);
	v.append(__ytmpVec);  
	__matrixQRow.clear();
      }
    }
  }

  // tmpMatrix = J_inv_n*V
  Matrix tmpMatrix(_matrixV.rows(), _matrixV.cols(), 0.0);
  assertion2(_oldInvJacobian.cols() == _matrixV.rows(), _oldInvJacobian.cols(), _matrixV.rows());
  multiply(_oldInvJacobian, _matrixV, tmpMatrix);

  // tmpMatrix = (W-J_inv_n*V)
  tmpMatrix *= -1.;
  tmpMatrix = tmpMatrix + _matrixW;
  
  // invJacobian = (W - J_inv_n*V)*(V^T*V)^-1*V^T
  assertion2(tmpMatrix.cols() == v.rows(), tmpMatrix.cols(), v.rows());

  Matrix tmp_invJacobian(_invJacobian.rows(), _invJacobian.cols(), 0.0);
  multiply(tmpMatrix, v, tmp_invJacobian);
  tmp_invJacobian = tmp_invJacobian + _oldInvJacobian;

  DataValues negRes(_residuals);
  negRes *= -1.;
  // solve delta_x = - J_inv*residuals
  multiply(tmp_invJacobian, negRes, xUpdate);  
}


void MVQNPostProcessing::computeNewtonFactorsQRDecomposition
(PostProcessing::DataMap& cplData, DataValues& xUpdate)
{
  preciceTrace("computeNewtonFactorsQRDecomposition()");
  using namespace tarch::la;
 
  // ------------- update inverse Jacobian -----------
  // J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
  // ----------------------------------------- -------
  
  DataMatrix v;
  bool linearDependence = true;
  while (linearDependence){
    linearDependence = false;
    v.clear();
    
    DataMatrix Vcopy(_matrixV);
    DataMatrix Q(Vcopy.rows(), Vcopy.cols(), 0.0);
    DataMatrix R(Vcopy.cols(), Vcopy.cols(), 0.0);
 
    modifiedGramSchmidt(Vcopy, Q, R);
    
    if (_matrixV.cols() > 1){
      for (int i=0; i < _matrixV.cols(); i++){
	if (R(i,i) < _singularityLimit){
	  preciceDebug("   Removing linear dependent column " << i);
	  _infostream<<"(modifiedGramSchmidt) removing linear dependent column "<<i<<"  time step: "<<t<<" iteration: "<<k<<"\n"<<std::flush;
	  linearDependence = true;
	  removeMatrixColumn(i);
	}
      }
    }
    if(not linearDependence)
    {
      DataValues ytmpVec(_matrixV.cols(), 0.0);
      DataValues _matrixQRow;
      for(int i = 0; i < Q.rows(); i++)
      {
	for(int j=0; j < Q.cols(); j++){
	_matrixQRow.append(Q(i,j));
        }
	backSubstitution(R, _matrixQRow, ytmpVec);
	v.append(ytmpVec);  
        _matrixQRow.clear();
      }
    }
  }

  // tmpMatrix = J_inv_n*V
  Matrix tmpMatrix(_matrixV.rows(), _matrixV.cols(), 0.0);
  assertion2(_oldInvJacobian.cols() == _matrixV.rows(), _oldInvJacobian.cols(), _matrixV.rows());
  multiply(_oldInvJacobian, _matrixV, tmpMatrix);

  // tmpMatrix = (W-J_inv_n*V)
  tmpMatrix *= -1.;
  tmpMatrix = tmpMatrix + _matrixW;
  
  // invJacobian = (W - J_inv_n*V)*(V^T*V)^-1*V^T
  assertion2(tmpMatrix.cols() == v.rows(), tmpMatrix.cols(), v.rows());
  multiply(tmpMatrix, v, _invJacobian);
  _invJacobian = _invJacobian + _oldInvJacobian;
  
  DataValues negRes(_residuals);
  negRes *= -1.;
  
  // solve delta_x = - J_inv*residuals
  multiply(_invJacobian, negRes, xUpdate); 
}


void MVQNPostProcessing::computeNewtonFactorsLUDecomposition
(PostProcessing::DataMap& cplData, DataValues& xUpdate)
{
  preciceTrace("computeNewtonFactorsLUDecomposition()");
  using namespace tarch::la;
  
  // ------------- update inverse Jacobian -----------
  // J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
  // ----------------------------------------- -------

  //----------------------------DEBUG------------------------------------
  /*
  std::string Vfile("matrixV"+st.str()+"_"+sk.str()+".m");
  std::string Wfile("matrixW"+st.str()+"_"+sk.str()+".m");
  _matrixV.printm(Vfile.c_str());
  _matrixW.printm(Wfile.c_str());
  */
  // --------------------------------------------------------------------
  
  DataMatrix VTVLU(_matrixV.cols(), _matrixV.cols(), 0.0);
  DataMatrix v;
  multiply(transpose(_matrixV), _matrixV, VTVLU);  // VTV = V^T*V
  
  
  //----------------------------DEBUG------------------------------------
  /*
  std::string VTVfile("VTV"+st.str()+"_"+sk.str()+".m");
  VTVLU.printm(VTVfile.c_str());
  */
  // --------------------------------------------------------------------
  
  
  DataValues pivots(_matrixV.cols(), 0.0);
  lu(VTVLU,pivots);
  
  //----------------------------DEBUG------------------------------------
  /*
  std::string VTVLUfile("VTVLU"+st.str()+"_"+sk.str()+".m");
  VTVLU.printm(VTVLUfile.c_str());
  */
  // --------------------------------------------------------------------
  
  
  DataValues ytmpVec(_matrixV.cols(), 0.0);
  DataValues xtmpVec(_matrixV.cols(), 0.0);
  DataValues _matrixVRow;
  for(int i = 0; i < _matrixV.rows(); i++)
  {
    for(int j=0; j < _matrixV.cols(); j++){
      _matrixVRow.append(_matrixV(i,j));
    }
    
    //----------------------------DEBUG------------------------------------
    /*
    std::stringstream si; si <<i;
    std::string mrowfile("mrow"+st.str()+"_"+sk.str()+"_"+si.str()+".m");
    _matrixVRow.printm(mrowfile.c_str());
    */
    // --------------------------------------------------------------------
    
    
    // account for pivoting in lu-decomposition
    assertion2(_matrixVRow.size() == pivots.size(), _matrixVRow.size(), pivots.size());
    for ( int i=0; i < _matrixVRow.size(); i++ ){
      double temp = _matrixVRow[i];
      _matrixVRow[i] = _matrixVRow[pivots[i]];
      _matrixVRow[pivots[i]] = temp;
    }
    forwardSubstitution(VTVLU, _matrixVRow, ytmpVec);
    
    //----------------------------DEBUG------------------------------------
    /*
    std::string ytmpfile("ytmpVec"+st.str()+"_"+sk.str()+"_"+si.str()+".m");
    ytmpVec.printm(ytmpfile.c_str());
    */
    // --------------------------------------------------------------------
    
    backSubstitution(VTVLU, ytmpVec, xtmpVec);
    
    //----------------------------DEBUG------------------------------------
    /*
    std::string xtmpfile("xtmpVec"+st.str()+"_"+sk.str()+"_"+si.str()+".m");
    xtmpVec.printm(xtmpfile.c_str());
    */
    // --------------------------------------------------------------------
    
    v.append(xtmpVec);  
    _matrixVRow.clear();
  }
  
  
  //----------------------------DEBUG------------------------------------
  /*
  std::string vfile("v"+st.str()+"_"+sk.str()+".m");
  v.printm(vfile.c_str());
  */
  // --------------------------------------------------------------------
  
    
  // tmpMatrix = J_inv_n*V
  DataMatrix tmpMatrix(_matrixV.rows(), _matrixV.cols(), 0.0);
  assertion2(_oldInvJacobian.cols() == _matrixV.rows(), _oldInvJacobian.cols(), _matrixV.rows());
  multiply(_oldInvJacobian, _matrixV, tmpMatrix);
  // tmpMatrix = (W-J_inv_n*V)
  tmpMatrix *= -1.;
  tmpMatrix = _matrixW + tmpMatrix;
  // invJacobian = (W - J_inv_n*V)*(V^T*V)^-1*V^T
  assertion2(tmpMatrix.cols() == v.rows(), tmpMatrix.cols(), v.rows());
  multiply(tmpMatrix, v, _invJacobian);
  _invJacobian = _invJacobian + _oldInvJacobian;
  
  DataValues negRes(_residuals);
  negRes *= -1.;
  
  // solve delta_x = - J_inv*residuals
  multiply(_invJacobian, negRes, xUpdate); 
}



void MVQNPostProcessing:: specializedIterationsConverged
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
  _oldInvJacobian = _invJacobian;
}

}}} // namespace precice, cplscheme, impl