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
//#include "utils/NumericalCompare.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

// tarch::logging::Log MVQNPostProcessing::
//       _log("precice::cplscheme::impl::MVQNPostProcessing");

      
MVQNPostProcessing:: MVQNPostProcessing
(
  double initialRelaxation,
  bool forceInitialRelaxation,
  int    maxIterationsUsed,
  int    timestepsReused,
  int 	 filter,
  double singularityLimit,
  std::vector<int> dataIDs,
  std::map<int,double> scalings)
:
  BaseQNPostProcessing(initialRelaxation, forceInitialRelaxation, maxIterationsUsed, timestepsReused,
		       filter, singularityLimit, dataIDs, scalings),
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





/**
 * This is a no-op at the moment, add implementation to handle secondary data
 */
void MVQNPostProcessing::updateDifferenceMatrices
(
  DataMap& cplData)
{
  using namespace tarch::la;

  // call the base method for common update of V, W matrices
  BaseQNPostProcessing::updateDifferenceMatrices(cplData);
}



void MVQNPostProcessing::computeQNUpdate
    (PostProcessing::DataMap& cplData, DataValues& xUpdate)
{
  preciceTrace("computeQNUpdate()");
  preciceDebug("Compute Newton factors ");

    /**      --- update inverse Jacobian ---
     *
     * J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
     */
    computeNewtonFactorsUpdatedQRDecomposition(cplData, xUpdate);
}



void MVQNPostProcessing::computeNewtonFactorsUpdatedQRDecomposition
(PostProcessing::DataMap& cplData, DataValues& xUpdate)
{
	preciceTrace("computeNewtonFactorsQRDecomposition()");

	/**      --- update inverse Jacobian ---
	*
	* J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
	*/
  
	// do: filtering of least-squares system to maintain good conditioning
	 if(_filter == PostProcessing::NOFILTER){
	    	// do nothing
	}else{
		// do: filtering of least-squares system to maintain good conditioning
		std::vector<int> delIndices(0);
		_qrV.applyFilter(_singularityLimit, delIndices, _matrixV);
		// start with largest index (as V,W matrices are shrinked and shifted
		for(int i = delIndices.size()-1; i >= 0; i--){
			preciceDebug("   Removing linear dependent column " << delIndices[i]);
			removeMatrixColumn(delIndices[i]);
			std::stringstream ss;
			ss << "(updatedQR) removing linear dependent column "<< delIndices[i] << "  time step: " << tSteps
			   << " iteration: " << its<< "\n" << std::endl;
			preciceDebug(ss.str()); writeInfo(ss.str()); std::cout<<ss.str()<<std::endl;
		}
		assertion2(_matrixV.cols() == _qrV.cols(), _matrixV.cols(), _qrV.cols());
	}
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

	Eigen::VectorXd yVec(_qrV.cols());

	// assertions for the case of processors with no vertices
	if(!_hasNodesOnInterface){
			assertion2(_qrV.cols() == getLSSystemCols(), _qrV.cols(), getLSSystemCols());
			assertion1(_qrV.rows() == 0, _qrV.rows());
			assertion1(Q.size() == 0, Q.size());
	}

	for (int i = 0; i < Q.rows(); i++) {
		Eigen::VectorXd Qrow = Q.row(i);
		yVec = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheLeft>(Qrow);
		Z.col(i) = yVec;
	}


	/**
	*  (1) Multiply J_prev * V =: V_tilde
	*/
	assertion2(_matrixV.rows() == _qrV.rows(), _matrixV.rows(), _qrV.rows());
	assertion2(getLSSystemCols() == _qrV.cols(), getLSSystemCols(), _qrV.cols());

	// TODO: transpose V efficiently using blocking in parallel
	//       such that multiplication is cache efficient
	Eigen::MatrixXd W_til = Eigen::MatrixXd::Zero(_qrV.rows(), _qrV.cols());

	// multiply J_prev * V = W_til of dimension: (n x n) * (n x m) = (n x m),
	//                                    parallel:  (n_global x n_local) * (n_local x m) = (n_local x m)
	_parMatrixOps.multiply(_oldInvJacobian, V, W_til, _dimOffsets, getLSSystemRows(), getLSSystemRows(), getLSSystemCols(), false);


	// W_til = (W-J_inv_n*V) = (W-V_tilde)
	W_til *= -1.;
	W_til = W_til + W;


	/**
	*  (2) compute invJacobian = W_til*Z
	*
	*  where Z = (V^T*V)^-1*V^T via QR-dec and back-substitution
	*  and W_til = (W - J_inv_n*V)
	*
	*  dimension: (n x n) * (n x m) = (n x m),
	*  parallel:  (n_global x n_local) * (n_local x m) = (n_local x m)
	*/

	_parMatrixOps.multiply(W_til, Z, _invJacobian, _dimOffsets, getLSSystemRows(), getLSSystemCols(), getLSSystemRows());

	// update Jacobian
	_invJacobian = _invJacobian + _oldInvJacobian;

	/**
 	 *  (3) solve delta_x = - J_inv * res
	 */
	Eigen::VectorXd res_tilde(_residuals.size());
	Eigen::VectorXd xUp(_residuals.size());
	for(int i = 0; i < res_tilde.size(); i++)
	  res_tilde(i) = _residuals(i) - _designSpecification(i);

	res_tilde *= -1.;

	// multiply J_inv * (-res) = x_Update of dimension: (n x n) * (n x 1) = (n x 1),
	//                                        parallel:  (n_global x n_local) * (n_local x 1) = (n_local x 1)
	_parMatrixOps.multiply(_invJacobian, res_tilde, xUp, _dimOffsets, getLSSystemRows(), getLSSystemRows(), 1);

	for(int i = 0; i < xUp.size(); i++)
	  xUpdate(i) = xUp(i);

}


//void MVQNPostProcessing::computeNewtonFactorsUpdatedQRDecomposition
//(PostProcessing::DataMap& cplData, DataValues& xUpdate)
//{
//  preciceTrace("computeNewtonFactorsQRDecomposition()");
//  using namespace tarch::la;
//
//  /**      --- update inverse Jacobian ---
//   *
//   * J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
//   */
//  Eigen::MatrixXd Z;
//  Eigen::MatrixXd V;
//  Eigen::MatrixXd W;
//  bool linearDependence = true;
//
//  while (linearDependence) {
//		linearDependence = false;
//
//		auto R =  _qrV.matrixR();
//		if (getLSSystemCols() > 1) {
//			for (int i = 0; i < _qrV.cols(); i++) {
//				//if (std::fabs(R(i, i)) < _singularityLimit) {
//				if (std::fabs(R(i, i)) < 0.0) {
//					std::stringstream ss;
//					ss << "(updatedQR) removing linear dependent column "<< i << "  time step: " << tSteps
//					   << " iteration: " << its<< "\n" << std::endl;
//					preciceDebug(ss.str()); writeInfo(ss.str()); std::cout<<ss.str()<<std::endl;
//
//					linearDependence = true;
//					removeMatrixColumn(i);
//				}
//			}
//		}
//		if (not linearDependence) {
//			//Matrix __Q(_qrV.rows(), _qrV.cols(), 0.0);
//			//Z = Matrix(_qrV.cols(), _qrV.rows(), 0.0);
//			Z.resize(_qrV.cols(), _qrV.rows());
//			V.resize(_matrixV.rows(), _matrixV.cols());
//			W.resize(_matrixW.rows(), _matrixW.cols());
//
//			for (int i = 0; i < V.rows(); i++)
//				for (int j = 0; j < V.cols(); j++) {
//					V(i, j) = _matrixV(i, j);
//				}
//			for (int i = 0; i < V.rows(); i++)
//				for (int j = 0; j < V.cols(); j++) {
//					W(i, j) = _matrixW(i, j);
//				}
//
//			auto Q = _qrV.matrixQ();
//			R = _qrV.matrixR();
//
//			Eigen::VectorXd yVec(_qrV.cols());
//
//			// assertions for the case of processors with no vertices
//			if(!_hasNodesOnInterface){
//					assertion2(_qrV.cols() == getLSSystemCols(), _qrV.cols(), getLSSystemCols());
//					assertion1(_qrV.rows() == 0, _qrV.rows());
//					assertion1(Q.size() == 0, Q.size());
//			}
//
//			for (int i = 0; i < Q.rows(); i++) {
//				Eigen::VectorXd Qrow = Q.row(i);
//				yVec = R.triangularView<Eigen::Upper>().solve<Eigen::OnTheLeft>(Qrow);
//				Z.col(i) = yVec;
//			}
//		}
//	}
//
//
//  /**
//   *  (1) Multiply J_prev * V =: V_tilde
//   */
//  assertion2(_matrixV.rows() == _qrV.rows(), _matrixV.rows(), _qrV.rows());
//  assertion2(getLSSystemCols() == _qrV.cols(), getLSSystemCols(), _qrV.cols());
//
//  // TODO: transpose V efficiently using blocking in parallel
//  //       such that multiplication is cache efficient
//  Eigen::MatrixXd W_til = Eigen::MatrixXd::Zero(_qrV.rows(), _qrV.cols());
//
//  // multiply J_prev * V = W_til of dimension: (n x n) * (n x m) = (n x m),
//  //                                    parallel:  (n_global x n_local) * (n_local x m) = (n_local x m)
//  _parMatrixOps.multiply(_oldInvJacobian, V, W_til, _dimOffsets, getLSSystemRows(), getLSSystemRows(), getLSSystemCols(), false);
//
//
//  // W_til = (W-J_inv_n*V) = (W-V_tilde)
//  W_til *= -1.;
//  W_til = W_til + W;
//
//
//  /**
//   *  (2) compute invJacobian = W_til*Z
//   *
//   *  where Z = (V^T*V)^-1*V^T via QR-dec and back-substitution
//   *  and W_til = (W - J_inv_n*V)
//   *
//   *  dimension: (n x n) * (n x m) = (n x m),
//   *  parallel:  (n_global x n_local) * (n_local x m) = (n_local x m)
//   */
//
//  _parMatrixOps.multiply(W_til, Z, _invJacobian, _dimOffsets, getLSSystemRows(), getLSSystemCols(), getLSSystemRows());
//
//  // update Jacobian
//  _invJacobian = _invJacobian + _oldInvJacobian;
//
//  /**
//   *  (3) solve delta_x = - J_inv * res
//   */
//
//  Eigen::VectorXd res(_residuals.size());
//  Eigen::VectorXd xUp(_residuals.size());
//  for(int i = 0; i < res.size(); i++)
//	  res(i) = _residuals(i);
//
//  res *= -1.;
//
//  // multiply J_inv * (-res) = x_Update of dimension: (n x n) * (n x 1) = (n x 1),
//  //                                        parallel:  (n_global x n_local) * (n_local x 1) = (n_local x 1)
//  _parMatrixOps.multiply(_invJacobian, res, xUp, _dimOffsets, getLSSystemRows(), getLSSystemRows(), 1);
//
//  for(int i = 0; i < xUp.size(); i++)
//	  xUpdate(i) = xUp(i);
//
//}

void MVQNPostProcessing:: specializedIterationsConverged
(
   DataMap & cplData)
{
  // store inverse Jacobian from last time step
  _oldInvJacobian = _invJacobian;
}

}}} // namespace precice, cplscheme, impl

#endif
