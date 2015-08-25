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


  double init = 0.0;
  int entries = _residuals.size();
  int global_n = 0;

	if (not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode) {
		global_n = entries;
	}else{
		global_n = _dimOffsets.back();
	}
  
  _invJacobian = Matrix(global_n, entries, init);
  _oldInvJacobian = Matrix(global_n, entries, init);
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
  using namespace tarch::la;

    // ------------- update inverse Jacobian -----------
    // J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
    // ----------------------------------------- -------
    preciceDebug("   Compute Newton factors ");
    computeNewtonFactorsUpdatedQRDecomposition(cplData, xUpdate);
}


void MVQNPostProcessing::computeNewtonFactorsUpdatedQRDecomposition
(PostProcessing::DataMap& cplData, DataValues& xUpdate)
{
  preciceTrace("computeNewtonFactorsQRDecomposition()");
  using namespace tarch::la;
 
  // ------------- update inverse Jacobian -----------
  // J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
  // ----------------------------------------- -------

  Matrix Z;
  bool linearDependence = true;
  
  while (linearDependence) {
		linearDependence = false;

		Matrix __R(_qrV.cols(), _qrV.cols(), 0.0);
		auto r = _qrV.matrixR();
		for (int i = 0; i < r.rows(); i++)
			for (int j = 0; j < r.cols(); j++) {
				__R(i, j) = r(i, j);
			}

		if (getLSSystemCols() > 1) {
			for (int i = 0; i < __R.rows(); i++) {
				//if (std::fabs(__R(i, i)) < _singularityLimit) {
				if (std::fabs(__R(i, i)) < 0.0) {
					std::stringstream ss;
					ss << "(updatedQR) removing linear dependent column "<< i << "  time step: " << tSteps
					   << " iteration: " << its<< "\n" << std::endl;
					preciceDebug(ss.str()); writeInfo(ss.str()); std::cout<<ss.str()<<std::endl;

					linearDependence = true;
					removeMatrixColumn(i);
				}
			}
		}
		if (not linearDependence) {
			Matrix __Q(_qrV.rows(), _qrV.cols(), 0.0);
			Z = Matrix(_qrV.cols(), _qrV.rows(), 0.0);

			DataValues __ytmpVec(_qrV.cols(), 0.0);
			DataValues __matrixQRow;
			auto q = _qrV.matrixQ();
			for (int i = 0; i < q.rows(); i++)
				for (int j = 0; j < q.cols(); j++) {
					__Q(i, j) = q(i, j);
				}

			// assertions for the case of processors with no vertices
			if(!_hasNodesOnInterface){
					assertion2(_qrV.cols() == getLSSystemCols(), _qrV.cols(), getLSSystemCols());
					assertion1(_qrV.rows() == 0, _qrV.rows());
					assertion1(__Q.size() == 0, __Q.size());
			}

			r = _qrV.matrixR();
			for (int i = 0; i < r.rows(); i++)
				for (int j = 0; j < r.cols(); j++) {
					__R(i, j) = r(i, j);
				}

			for (int i = 0; i < __Q.rows(); i++) {
				for (int j = 0; j < __Q.cols(); j++) {
					__matrixQRow.append(__Q(i, j));
				}


				backSubstitution(__R, __matrixQRow, __ytmpVec);
				for(int p = 0; p < __ytmpVec.size(); p++)
					Z(p,i) = __ytmpVec(p);
				__matrixQRow.clear();
			}
		}
	}


  /**
   *  (1) Multiply J_prev * V =: V_tilde
   */
  assertion2(_matrixV.rows() == _qrV.rows(), _matrixV.rows(), _qrV.rows());
  assertion2(getLSSystemCols() == _qrV.cols(), getLSSystemCols(), _qrV.cols());

  // TODO: transpose V efficiently using blocking in parallel
  //       such that multiplication is cache efficient
  Matrix W_til(_qrV.rows(), _qrV.cols(), 0.0);

  // multiply J_prev * V = W_til of dimension: (n x n) * (n x m) = (n x m),
  //                                    parallel:  (n_global x n_local) * (n_local x m) = (n_local x m)
  _parMatrixOps.multiply(_oldInvJacobian, _matrixV, W_til, _dimOffsets, getLSSystemRows(), getLSSystemRows(), getLSSystemCols());


  // W_til = (W-J_inv_n*V) = (W-V_tilde)
  W_til *= -1.;
  W_til = W_til + _matrixW;

  
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

  DataValues res(_residuals);
  res *= -1.;

  // multiply J_inv * (-res) = x_Update of dimension: (n x n) * (n x 1) = (n x 1),
  //                                        parallel:  (n_global x n_local) * (n_local x 1) = (n_local x 1)
  _parMatrixOps.multiply(_invJacobian, res, xUpdate, _dimOffsets, getLSSystemRows(), getLSSystemRows());

}

void MVQNPostProcessing:: specializedIterationsConverged
(
   DataMap & cplData)
{
  // store inverse Jacobian from last time step
  _oldInvJacobian = _invJacobian;
}

}}} // namespace precice, cplscheme, impl

#endif
