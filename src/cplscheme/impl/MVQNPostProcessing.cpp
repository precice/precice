#ifndef PRECICE_NO_MPI

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
#include "utils/MasterSlave.hpp"
#include "utils/EventTimings.hpp"
#include "tarch/la/Scalar.h"
#include "io/TXTWriter.hpp"
#include "io/TXTReader.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "com/Communication.hpp"
#include "Eigen/Dense"

#include <time.h>
#include <sstream>
#include <fstream>
#include <cstring>

#include <thread>
#include <chrono>
//#include "utils/NumericalCompare.hpp"

using precice::utils::Event;


namespace precice {
namespace cplscheme {
namespace impl {

 //tarch::logging::Log MVQNPostProcessing::_log("precice::cplscheme::impl::MVQNPostProcessing");

      
MVQNPostProcessing:: MVQNPostProcessing
(
  double initialRelaxation,
  bool forceInitialRelaxation,
  int    maxIterationsUsed,
  int    timestepsReused,
  int 	 filter,
  double singularityLimit,
  std::vector<int> dataIDs,
  PtrPreconditioner preconditioner)
:
  BaseQNPostProcessing(initialRelaxation, forceInitialRelaxation, maxIterationsUsed, timestepsReused,
		       filter, singularityLimit, dataIDs, preconditioner),
//  _secondaryOldXTildes(),
  _invJacobian(),
  _oldInvJacobian(),
  _cyclicCommLeft(nullptr),
  _cyclicCommRight(nullptr),
  _parMatrixOps()
{}

MVQNPostProcessing::~MVQNPostProcessing()
{
	//if(utils::MasterSlave::_masterMode ||utils::MasterSlave::_slaveMode){ // not possible because of tests, MasterSlave is deactivated when PP is killed

	// close and shut down cyclic communication connections
	if(_cyclicCommRight != nullptr || _cyclicCommLeft != nullptr){
		if((utils::MasterSlave::_rank % 2) == 0)
		{
		  _cyclicCommLeft->closeConnection();
		  _cyclicCommRight->closeConnection();
		}else{
		  _cyclicCommRight->closeConnection();
		  _cyclicCommLeft->closeConnection();
		}
		_cyclicCommRight = nullptr;
		_cyclicCommLeft = nullptr;
	}
}


void MVQNPostProcessing:: initialize
(
  DataMap& cplData )
{
  preciceTrace(__func__);
  Event e("MVQNPostProcessing::initialize()", true, true); // time measurement, barrier
  // do common QN post processing initialization
  BaseQNPostProcessing::initialize(cplData);
  

  if(utils::MasterSlave::_masterMode ||utils::MasterSlave::_slaveMode){
		/*
		 * TODO: FIXME: This is a temporary and hacky realization of the cyclic commmunication between slaves
		 * 				Therefore the requesterName and accessorName are not given (cf solverInterfaceImpl).
		 * 				The master-slave communication should be modified such that direct communication between
		 * 				slaves is possible (via MPIDirect)
		 */


		 _cyclicCommLeft = com::Communication::SharedPointer(new com::MPIPortsCommunication("."));
		 _cyclicCommRight = com::Communication::SharedPointer(new com::MPIPortsCommunication("."));

		 // initialize cyclic communication between successive slaves
		int prevProc = (utils::MasterSlave::_rank-1 < 0) ? utils::MasterSlave::_size-1 : utils::MasterSlave::_rank-1;
		if((utils::MasterSlave::_rank % 2) == 0)
		{
		  _cyclicCommLeft->acceptConnection("cyclicComm-" + std::to_string(prevProc), "", 0, 1 );
		  _cyclicCommRight->requestConnection("cyclicComm-" +  std::to_string(utils::MasterSlave::_rank), "", 0, 1 );
		}else{
		  _cyclicCommRight->requestConnection("cyclicComm-" +  std::to_string(utils::MasterSlave::_rank), "", 0, 1 );
		  _cyclicCommLeft->acceptConnection("cyclicComm-" + std::to_string(prevProc), "", 0, 1 );
		}
  }

  // initialize parallel matrix-matrix operation module
  _parMatrixOps.initialize(_cyclicCommLeft, _cyclicCommRight);

  int entries = _residuals.size();
  int global_n = 0;

	if (not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode) {
		global_n = entries;
	}else{
		global_n = _dimOffsets.back();
	}
  
  _invJacobian = Eigen::MatrixXd::Zero(global_n, entries);
  _oldInvJacobian = Eigen::MatrixXd::Zero(global_n, entries);
  _Wtil = Eigen::MatrixXd::Zero(entries, 0);

  _preconditioner->triggerGlobalWeights(global_n);
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
    for (int id: _secondaryDataIDs){
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
  preciceTrace(__func__);

  // call the base method for common update of V, W matrices
  BaseQNPostProcessing::updateDifferenceMatrices(cplData);

  if (_firstIteration && (_firstTimeStep || _forceInitialRelaxation)) {
    // do nothing: constant relaxation
  } else {
    if (not _firstIteration) {
      // Update matrix _Wtil = (W - J_prev*V) with newest information

      Eigen::VectorXd v = _matrixV.column(0);
      Eigen::VectorXd w = _matrixW.column(0);


      bool columnLimitReached = _Wtil.cols() == _maxIterationsUsed;
      bool overdetermined = _Wtil.cols() <= getLSSystemRows();

      Eigen::VectorXd wtil = Eigen::VectorXd::Zero(_matrixV.rows());

      std::cout<<"update W_til ... "<<std::endl;

      // compute J_prev * V(0) := wtil the new column in _Wtil of dimension: (n x n) * (n x 1) = (n x 1),
      //                                        parallel: (n_global x n_local) * (n_local x 1) = (n_local x 1)
      _parMatrixOps.multiply(_oldInvJacobian, v, wtil, _dimOffsets, getLSSystemRows(), getLSSystemRows(), 1, false);

      wtil = w - wtil;

      if (not columnLimitReached && overdetermined) {

        appendFront(_Wtil, wtil);   std::cout<<"      append front "<<std::endl;
      }else {
        shiftSetFirst(_Wtil, wtil); std::cout<<"      shift set first"<<std::endl;
      }
    }
  }
}



void MVQNPostProcessing::computeQNUpdate
    (PostProcessing::DataMap& cplData, DataValues& xUpdate)
{
  preciceTrace(__func__);
  Event e(__func__, true, true); // time measurement, barrier
  preciceDebug("Compute Newton factors ");

    /**      --- update inverse Jacobian ---
     *
     * J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
     */

  Event ePrecond_1("precond J (1)", true, true); // time measurement, barrier
  _preconditioner->apply(_Wtil);
  _preconditioner->apply(_oldInvJacobian,false);
  _preconditioner->revert(_oldInvJacobian,true);
  ePrecond_1.stop();
  computeNewtonFactors(cplData, xUpdate);
  Event ePrecond_2("precond J (2)", true, true); // time measurement, barrier
  _preconditioner->revert(_Wtil);
  _preconditioner->revert(_oldInvJacobian,false);
  _preconditioner->apply(_oldInvJacobian,true);
  ePrecond_2.stop();
}


void MVQNPostProcessing::buildWtil()
{
  preciceTrace(__func__);
  /**
   * assumes that V, W, J_prev are already preconditioned,
   */
  Eigen::MatrixXd V(_matrixV.rows(), _matrixV.cols());
  Eigen::MatrixXd W(_matrixW.rows(), _matrixW.cols());

  // convert tarch matrices to Eigen matrices
  for (int i = 0; i < V.rows(); i++)
    for (int j = 0; j < V.cols(); j++) {
      V(i, j) = _matrixV(i, j);
    }
  for (int i = 0; i < W.rows(); i++)
    for (int j = 0; j < W.cols(); j++) {
      W(i, j) = _matrixW(i, j);
    }

  std::cout<<"build Wtil  ... "<<std::endl;

  Event e_WtilV("compute W_til = (W - J_prev*V)", true, true); // time measurement, barrier
  assertion2(_matrixV.rows() == _qrV.rows(), _matrixV.rows(), _qrV.rows());  assertion2(getLSSystemCols() == _qrV.cols(), getLSSystemCols(), _qrV.cols());

  // TODO: transpose V efficiently using blocking in parallel
  //       such that multiplication is cache efficient
  _Wtil = Eigen::MatrixXd::Zero(_qrV.rows(), _qrV.cols());

  // multiply J_prev * V = W_til of dimension: (n x n) * (n x m) = (n x m),
  //                                    parallel:  (n_global x n_local) * (n_local x m) = (n_local x m)
  _parMatrixOps.multiply(_oldInvJacobian, V, _Wtil, _dimOffsets, getLSSystemRows(), getLSSystemRows(), getLSSystemCols(), false);

  // W_til = (W-J_inv_n*V) = (W-V_tilde)
  _Wtil *= -1.;
  _Wtil = _Wtil + W;

  _resetLS = false;
  e_WtilV.stop();
}


void MVQNPostProcessing::buildJacobian()
{
  preciceTrace(__func__);
  /**      --- compute inverse Jacobian ---
  *
  * J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
  *
  * assumes that J_prev, V, W are already preconditioned
  */


  Eigen::MatrixXd Z(_qrV.cols(), _qrV.rows());
  Eigen::MatrixXd V(_matrixV.rows(), _matrixV.cols());
  Eigen::MatrixXd W(_matrixW.rows(), _matrixW.cols());

  // convert tarch matrices to Eigen matrices
  for (int i = 0; i < V.rows(); i++)
    for (int j = 0; j < V.cols(); j++) {
      V(i, j) = _matrixV(i, j);
    }
  for (int i = 0; i < W.rows(); i++)
    for (int j = 0; j < W.cols(); j++) {
      W(i, j) = _matrixW(i, j);
    }

  auto Q = _qrV.matrixQ();
  auto R = _qrV.matrixR();

  Event e(__func__, true, true); // time measurement, barrier
  Eigen::VectorXd yVec(_qrV.cols());

  // assertions for the case of processors with no vertices
  if(!_hasNodesOnInterface){
      assertion2(_qrV.cols() == getLSSystemCols(), _qrV.cols(), getLSSystemCols());
      assertion1(_qrV.rows() == 0, _qrV.rows());
      assertion1(Q.size() == 0, Q.size());
  }

  std::cout<<"build Jacobian ... "<<std::endl;

  /**
   *  (1) computation of matrix Z = (V^TV)^-1 * V^T as solution to the equation
   *      R*z = Q^T(i) for all columns i,  via back substitution.
   */
  Event e_qr("solve Z = (V^TV)^-1V^T via QR", true, true); // time measurement, barrier
  for (int i = 0; i < Q.rows(); i++) {
    Eigen::VectorXd Qrow = Q.row(i);
    yVec = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheLeft>(Qrow);
    Z.col(i) = yVec;
  }
  e_qr.stop();


  /**
  *  (2) Multiply J_prev * V =: W_tilde
  */
  assertion2(_matrixV.rows() == _qrV.rows(), _matrixV.rows(), _qrV.rows());  assertion2(getLSSystemCols() == _qrV.cols(), getLSSystemCols(), _qrV.cols());
  if(_resetLS)
    buildWtil();

  /**
  *  (3) compute invJacobian = W_til*Z
  *
  *  where Z = (V^T*V)^-1*V^T via QR-dec and back-substitution
  *  and W_til = (W - J_inv_n*V)
  *
  *  dimension: (n x n) * (n x m) = (n x m),
  *  parallel:  (n_global x n_local) * (n_local x m) = (n_local x m)
  */
  Event e_WtilZ("compute J = W_til*Z", true, true); // time measurement, barrier

  _parMatrixOps.multiply(_Wtil, Z, _invJacobian, _dimOffsets, getLSSystemRows(), getLSSystemCols(), getLSSystemRows());
  e_WtilZ.stop();

  // update Jacobian
  _invJacobian = _invJacobian + _oldInvJacobian;
}

void MVQNPostProcessing::computeNewtonFactors
(PostProcessing::DataMap& cplData, DataValues& xUpdate)
{
  preciceTrace(__func__);
  /**      --- update inverse Jacobian ---
  *
  * J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
  */


  Eigen::MatrixXd Z(_qrV.cols(), _qrV.rows());
  Eigen::MatrixXd V(_matrixV.rows(), _matrixV.cols());
  Eigen::MatrixXd W(_matrixW.rows(), _matrixW.cols());

  // convert tarch matrices to Eigen matrices
  for (int i = 0; i < V.rows(); i++)
    for (int j = 0; j < V.cols(); j++) {
      V(i, j) = _matrixV(i, j);
    }
  for (int i = 0; i < W.rows(); i++)
    for (int j = 0; j < W.cols(); j++) {
      W(i, j) = _matrixW(i, j);
    }

  auto Q = _qrV.matrixQ();
  auto R = _qrV.matrixR();

  Event e(__func__, true, true); // time measurement, barrier
  Eigen::VectorXd yVec(_qrV.cols());

  // assertions for the case of processors with no vertices
  if(!_hasNodesOnInterface){
      assertion2(_qrV.cols() == getLSSystemCols(), _qrV.cols(), getLSSystemCols());
      assertion1(_qrV.rows() == 0, _qrV.rows());
      assertion1(Q.size() == 0, Q.size());
  }

  /**
   *  (1) computation of matrix Z = (V^TV)^-1 * V^T as solution to the equation
   *      R*z = Q^T(i) for all columns i,  via back substitution.
   */
  Event e_qr("solve Z = (V^TV)^-1V^T via QR", true, true); // time measurement, barrier
  for (int i = 0; i < Q.rows(); i++) {
    Eigen::VectorXd Qrow = Q.row(i);
    yVec = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheLeft>(Qrow);
    Z.col(i) = yVec;
  }
  e_qr.stop();


  /**
  *  (2) Multiply J_prev * V =: W_tilde
  */
  assertion2(_matrixV.rows() == _qrV.rows(), _matrixV.rows(), _qrV.rows());  assertion2(getLSSystemCols() == _qrV.cols(), getLSSystemCols(), _qrV.cols());

  // rebuild matrix Wtil if V changes completely.
  if(_resetLS)
    buildWtil();

  /**
  *  (3) compute r_til = Z*(-residual)
  *
  *  where Z = (V^T*V)^-1*V^T via QR-dec and back-substitution
  *
  *  dimension: (m x n) * (n x 1) = (m x 1),
  *  parallel:  (m x n_local) * (n x 1) = (m x 1)
  */
  Eigen::VectorXd negRes(_residuals.size());
  for(int i = 0; i < negRes.size(); i++)
    negRes(i) = - _residuals(i);

  Eigen::VectorXd r_til_loc = Eigen::VectorXd::Zero(getLSSystemCols());
  Eigen::VectorXd r_til = Eigen::VectorXd::Zero(getLSSystemCols());

  std::cout<<"          compute Z*(-res) "<<std::endl;

  r_til_loc.noalias() = Z * negRes;

  std::this_thread::sleep_for (std::chrono::seconds(1* (1+utils::MasterSlave::_rank)));
  std::cout<<"r_til_loc.size() on proc "<<utils::MasterSlave::_rank<<": "<<r_til_loc.size()<<std::endl;
  std::cout<<"r_til.size() on proc "<<utils::MasterSlave::_rank<<": "<<r_til.size()<<std::endl;

  // if serial computation on single processor, i.e, no master-slave mode
  if( not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode){
    r_til = r_til_loc;
  }else{
    utils::MasterSlave::allreduceSum(r_til_loc.data(), r_til.data(), r_til_loc.size());
  }

  //std::this_thread::sleep_for (std::chrono::seconds(1*utils::MasterSlave::_rank));
  //std::cout<<"r_til:loc: proc"<<utils::MasterSlave::_rank<<": "<<r_til_loc<<std::endl;
  //std::cout<<"r_til: "<<r_til<<std::endl;


  /**
   * (4) compute _Wtil * r_til
   *
   * dimension: (n x m) * (m x 1) = (n x 1),
   * parallel:  (n_local x m) * (m x 1) = (n_local x 1)
   *
   * Note: r_til is not distributed but locally stored on each proc (dimension m x 1)
   */
  std::cout<<"          compute W_til * tmp "<<std::endl;
  Eigen::VectorXd xUptmp(_residuals.size());
  xUptmp = _Wtil * r_til;                      // local product, result is naturally distributed.

  /**
   *  (5) xUp = J_prev * (-res) + Wtil*Z*(-res)
   */
  std::cout<<"          compute J_prev*(-res) "<<std::endl;
  Eigen::VectorXd xUp(_residuals.size());
  _parMatrixOps.multiply(_oldInvJacobian, negRes, xUp, _dimOffsets, getLSSystemRows(), getLSSystemRows(), 1, false);

  xUp = xUp + xUptmp;

  for(int i = 0; i < xUp.size(); i++)
    xUpdate(i) = xUp(i);

  // pending deletion: delete Wtil
  if (_firstIteration && _timestepsReused == 0 && not _forceInitialRelaxation) {
    //_Wtil.conservativeResize(0,0);
    _resetLS = true;
  }
}



void MVQNPostProcessing::computeNewtonFactorsUpdatedQRDecomposition
(PostProcessing::DataMap& cplData, DataValues& xUpdate)
{
	preciceTrace("computeNewtonFactorsQRDecomposition()");

	/**      --- update inverse Jacobian ---
	*
	* J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
	*/
  

	Eigen::MatrixXd Z(_qrV.cols(), _qrV.rows());
	Eigen::MatrixXd V(_matrixV.rows(), _matrixV.cols());
	Eigen::MatrixXd W(_matrixW.rows(), _matrixW.cols());

	// convert tarch matrices to Eigen matrices
	for (int i = 0; i < V.rows(); i++)
		for (int j = 0; j < V.cols(); j++) {
			V(i, j) = _matrixV(i, j);
		}
	for (int i = 0; i < W.rows(); i++)
		for (int j = 0; j < W.cols(); j++) {
			W(i, j) = _matrixW(i, j);
		}

	auto Q = _qrV.matrixQ();
	auto R = _qrV.matrixR();

	Event e(__func__, true, true); // time measurement, barrier
	Eigen::VectorXd yVec(_qrV.cols());

	// assertions for the case of processors with no vertices
	if(!_hasNodesOnInterface){
			assertion2(_qrV.cols() == getLSSystemCols(), _qrV.cols(), getLSSystemCols());
			assertion1(_qrV.rows() == 0, _qrV.rows());
			assertion1(Q.size() == 0, Q.size());
	}

	/**
	 *  (1) computation of matrix Z = (V^TV)^-1 * V^T as solution to the equation
	 *      R*z = Q^T(i) for all columns i,  via back substitution.
	 */
	Event e_qr("solve Z = (V^TV)^-1V^T via QR", true, true); // time measurement, barrier
	for (int i = 0; i < Q.rows(); i++) {
		Eigen::VectorXd Qrow = Q.row(i);
		yVec = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheLeft>(Qrow);
		Z.col(i) = yVec;
	}
	e_qr.stop();


	/**
	*  (2) Multiply J_prev * V =: W_tilde
	*/
	Event e_WtilV("compute W_til = (W - J_prev*V)", true, true); // time measurement, barrier
	assertion2(_matrixV.rows() == _qrV.rows(), _matrixV.rows(), _qrV.rows());  assertion2(getLSSystemCols() == _qrV.cols(), getLSSystemCols(), _qrV.cols());

	// TODO: transpose V efficiently using blocking in parallel
	//       such that multiplication is cache efficient
	Eigen::MatrixXd W_til = Eigen::MatrixXd::Zero(_qrV.rows(), _qrV.cols());

	// multiply J_prev * V = W_til of dimension: (n x n) * (n x m) = (n x m),
	//                                    parallel:  (n_global x n_local) * (n_local x m) = (n_local x m)
	_parMatrixOps.multiply(_oldInvJacobian, V, W_til, _dimOffsets, getLSSystemRows(), getLSSystemRows(), getLSSystemCols(), false);


	// W_til = (W-J_inv_n*V) = (W-V_tilde)
	W_til *= -1.;
	W_til = W_til + W;

	e_WtilV.stop();

	/**
	*  (3) compute invJacobian = W_til*Z
	*
	*  where Z = (V^T*V)^-1*V^T via QR-dec and back-substitution
	*  and W_til = (W - J_inv_n*V)
	*
	*  dimension: (n x n) * (n x m) = (n x m),
	*  parallel:  (n_global x n_local) * (n_local x m) = (n_local x m)
	*/
	Event e_WtilZ("compute J = W_til*Z", true, true); // time measurement, barrier

	_parMatrixOps.multiply(W_til, Z, _invJacobian, _dimOffsets, getLSSystemRows(), getLSSystemCols(), getLSSystemRows());
	e_WtilV.stop();

	// update Jacobian
	_invJacobian = _invJacobian + _oldInvJacobian;

	/**
 	 *  (4) solve delta_x = - J_inv * res
	 */
	Eigen::VectorXd res_tilde(_residuals.size());
  Eigen::VectorXd xUp(_residuals.size());
  for(int i = 0; i < res_tilde.size(); i++)
    res_tilde(i) = _residuals(i);

	res_tilde *= -1.;

	Event e_up("compute update = J*(-res)", true, true); // time measurement, barrier
	// multiply J_inv * (-res) = x_Update of dimension: (n x n) * (n x 1) = (n x 1),
	//                                        parallel: (n_global x n_local) * (n_local x 1) = (n_local x 1)
	_parMatrixOps.multiply(_invJacobian, res_tilde, xUp, _dimOffsets, getLSSystemRows(), getLSSystemRows(), 1, false);
  e_up.stop();

	for(int i = 0; i < xUp.size(); i++)
	  xUpdate(i) = xUp(i);
}


void MVQNPostProcessing:: specializedIterationsConverged
(
   DataMap & cplData)
{

  Event ePrecond_1("precond J (1)", true, true); // time measurement, barrier
  _preconditioner->apply(_residuals);
  _preconditioner->apply(_matrixV);
  _preconditioner->apply(_matrixW);

  _preconditioner->apply(_Wtil);
  _preconditioner->apply(_oldInvJacobian,false);
  _preconditioner->revert(_oldInvJacobian,true);

  if(_preconditioner->requireNewQR()){
    if(not (_filter==PostProcessing::QR2FILTER)){ //for QR2 filter, there is no need to do this twice
      _qrV.reset(_matrixV, getLSSystemRows());
    }
    _preconditioner->newQRfulfilled();
  }
  BaseQNPostProcessing::applyFilter();  // apply the configured filter to the LS system
  ePrecond_1.stop();

  // compute explicit representation of Jacobian
  buildJacobian();
  // store inverse Jacobian from last time step
  _oldInvJacobian = _invJacobian;

  Event ePrecond_2("precond J (2)", true, true); // time measurement, barrier
  _preconditioner->revert(_matrixW);
  _preconditioner->revert(_matrixV);
  _preconditioner->revert(_residuals);

  _preconditioner->revert(_oldInvJacobian,false);
  _preconditioner->apply(_oldInvJacobian,true);
  _preconditioner->revert(_Wtil);
  ePrecond_2.stop();


  // delete columns from matrix _Wtil according to reused time steps
  if (_matrixCols.front() == 0) { // Did only one iteration
    _matrixCols.pop_front();
  }

  if (_timestepsReused == 0) {
    if (_forceInitialRelaxation)
    {  // reset _Wtil if initial relaxation is enforced
      _Wtil.conservativeResize(0, 0);
    }
    //else: pending deletion
  }
  else if ((int) _matrixCols.size() > _timestepsReused) {
    int toRemove = _matrixCols.back();
    assertion1(toRemove > 0, toRemove);  assertion2(getLSSystemCols() > toRemove, getLSSystemCols(), toRemove);
    for (int i = 0; i < toRemove; i++) {
      removeColumnFromMatrix(_Wtil, _Wtil.cols() - 1);
    }
  }
}


void MVQNPostProcessing:: removeMatrixColumn
(
  int columnIndex)
{
  assertion(_matrixV.cols() > 1); assertion(_Wtil.cols() > 1);

  // remove column from matrix _Wtil
  removeColumnFromMatrix(_Wtil, columnIndex);

  BaseQNPostProcessing::removeMatrixColumn(columnIndex);
}

// ================== move this to helper class/module ================================
void MVQNPostProcessing::shiftSetFirst
(
    Eigen::MatrixXd& A, Eigen::VectorXd& v)
{
  assertion2(v.size() == A.rows(), v.size(), A.rows());
  int n = A.rows(), m = A.cols();
  //A.bottomRightCorner(n, m - 1) = A.topLeftCorner(n, m - 1);
  for(auto i = A.cols()-1; i > 0; i--)
        A.col(i) = A.col(i-1);
  A.col(0) = v;
}

void MVQNPostProcessing::appendFront
(
    Eigen::MatrixXd& A, Eigen::VectorXd& v)
{
  int n = A.rows(), m = A.cols();
  if (n <= 0 && m <= 0) {
    A = v;
  } else {
    assertion2(v.size() == n, v.size(), A.rows());
    A.conservativeResize(n, m + 1);
    //A.topRightCorner(n, m) = A.topLeftCorner(n, m); // bad error, reason unknown!
    for(auto i = A.cols()-1; i > 0; i--)
      A.col(i) = A.col(i-1);
    A.col(0) = v;
  }
}

void MVQNPostProcessing::removeColumnFromMatrix
(
    Eigen::MatrixXd& A, int col)
{
  assertion2(col < A.cols() && col >= 0, col, A.cols())
  for (int j = col; j < A.cols() - 1; j++)
    A.col(j) = A.col(j + 1);

  A.conservativeResize(A.rows(), A.cols() - 1);
}
// ================== move this to helper class/module ================================

}}} // namespace precice, cplscheme, impl

#endif
