// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "IQNILSPostProcessing.hpp"
#include "cplscheme/CouplingData.hpp"
#include "utils/Globals.hpp"
#include "tarch/la/GramSchmidt.h"
#include "tarch/la/MatrixVectorOperations.h"
#include "tarch/la/TransposedMatrix.h"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/la/Scalar.h"
#include "io/TXTWriter.hpp"
#include "io/TXTReader.hpp"
#include "utils/MasterSlave.hpp"
#include "QRFactorization.hpp"
#include "Eigen/Dense"

#include "tarch/tests/TestMacros.h"

#include <time.h>

//#include "utils/NumericalCompare.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

// tarch::logging::Log IQNILSPostProcessing::
//       _log("precice::cplscheme::impl::IQNILSPostProcessing");

IQNILSPostProcessing:: IQNILSPostProcessing
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
  _secondaryOldXTildes(),
  _secondaryMatricesW()
{
}

void IQNILSPostProcessing:: initialize
(
  DataMap& cplData )
{
  // do common QN post processing initialization
  BaseQNPostProcessing::initialize(cplData);

  double init = 0.0;
  // Fetch secondary data IDs, to be relaxed with same coefficients from IQN-ILS
  foreach (DataMap::value_type& pair, cplData){
    if (not utils::contained(pair.first, _dataIDs)){
      int secondaryEntries = pair.second->values->size();
      _secondaryOldXTildes[pair.first].append(DataValues(secondaryEntries, init));
    }
  }
  
}


void IQNILSPostProcessing::updateDifferenceMatrices
(
  DataMap& cplData)
{
  // Compute residuals of secondary data
  foreach (int id, _secondaryDataIDs){
    DataValues& secResiduals = _secondaryResiduals[id];
    PtrCouplingData data = cplData[id];
    assertion2(secResiduals.size() == data->values->size(),
               secResiduals.size(), data->values->size());
    secResiduals = *(data->values);
    secResiduals -= data->oldValues.column(0);
  }

  /*
   * ATTETION: changed the condition from _firstIteration && _firstTimeStep
   * to the following: 
   * underrelaxation has to be done, if the scheme has converged without even
   * entering post processing. In this case the V, W matrices would still be empty.
   * This case happended in the open foam example beamInCrossFlow.
   */ 
  if(_firstIteration && (_firstTimeStep ||  (_matrixCols.size() < 2))){
//     // Store x_tildes for secondary data
//     foreach (int id, _secondaryDataIDs){
//       assertion2(_secondaryOldXTildes[id].size() == cplData[id]->values->size(),
//                  _secondaryOldXTildes[id].size(), cplData[id]->values->size());
//       _secondaryOldXTildes[id] = *(cplData[id]->values);
//     }
// 
//     // Perform underrelaxation with initial relaxation factor for secondary data
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
      bool columnLimitReached = _matrixV.cols() == _maxIterationsUsed;
      bool overdetermined = _matrixV.cols() <= _matrixV.rows();
      if (not columnLimitReached && overdetermined){
        
	// Append column for secondary W matrices
        foreach (int id, _secondaryDataIDs){
          _secondaryMatricesW[id].appendFront(_secondaryResiduals[id]);
        }
      }
      else {
        // Shift column for secondary W matrices
        foreach (int id, _secondaryDataIDs){
          _secondaryMatricesW[id].shiftSetFirst(_secondaryResiduals[id]);
        }
      }

      // Compute delta_x_tilde for secondary data
      foreach (int id, _secondaryDataIDs){
        DataMatrix& secW = _secondaryMatricesW[id];
        assertion2(secW.column(0).size() == cplData[id]->values->size(),
                   secW.column(0).size(), cplData[id]->values->size());
        secW.column(0) = *(cplData[id]->values);
        secW.column(0) -= _secondaryOldXTildes[id];
      }
    }

    // Store x_tildes for secondary data
    foreach (int id, _secondaryDataIDs){ 
      assertion2(_secondaryOldXTildes[id].size() == cplData[id]->values->size(),
                 _secondaryOldXTildes[id].size(), cplData[id]->values->size());
      _secondaryOldXTildes[id] = *(cplData[id]->values);
    }
  }
  
  
  // call the base method for common update of V, W matrices
  BaseQNPostProcessing::updateDifferenceMatrices(cplData);
}


void IQNILSPostProcessing::computeUnderrelaxationSecondaryData
(
  DataMap& cplData)
{
    //Store x_tildes for secondary data
    foreach (int id, _secondaryDataIDs){
      assertion2(_secondaryOldXTildes[id].size() == cplData[id]->values->size(),
                 _secondaryOldXTildes[id].size(), cplData[id]->values->size());
      _secondaryOldXTildes[id] = *(cplData[id]->values);
    }

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


void IQNILSPostProcessing::computeQNUpdate
(PostProcessing::DataMap& cplData, DataValues& xUpdate)
{
  preciceTrace("computeQNUpdate()");
    using namespace tarch::la;
    
    // Calculate QR decomposition of matrix V and solve Rc = -Qr
    DataValues c;
    bool linearDependence = true;
    while (linearDependence){
      preciceDebug("   Compute Newton factors");
      linearDependence = false;
      
      DataMatrix Vcopy(_matrixV);
      DataMatrix Q(Vcopy.rows(), Vcopy.cols(), 0.0);
      DataMatrix R(Vcopy.cols(), Vcopy.cols(), 0.0);
      modifiedGramSchmidt(Vcopy, Q, R);
      
      if (_matrixV.cols() > 1){
        for (int i=0; i < _matrixV.cols(); i++){
          if (R(i,i) < _singularityLimit){
            preciceDebug("   Removing linear dependent column " << i);
	    _infostream<<"removing linear dependent column "<<i<<"\n"<<std::flush;
            linearDependence = true;
            removeMatrixColumn(i);
          }
        }
      }
      if (not linearDependence){
        
	preciceDebug("   Apply Newton factors");
        DataValues b(Q.cols(), 0.0);
        multiply(transpose(Q), _residuals, b); // = Qr
        b *= -1.0; // = -Qr
        assertion1(c.size() == 0, c.size());
        c.append(b.size(), 0.0);
	
        backSubstitution(R, b, c);
	
	multiply(_matrixW, c, xUpdate); // = Wc
	
	// ---------- QN factors with updatedQR -----------
	Matrix __Qt(_matrixV.cols(), _matrixV.rows(), 0.0);
        Matrix __R(_matrixV.cols(), _matrixV.cols(), 0.0);
	
	auto q = _qrV.matrixQ();
	for(int i = 0; i<q.rows(); i++)
	  for(int j = 0; j<q.cols(); j++)
	  {
	    __Qt(j,i) = q(i,j);
	  }
	auto r = _qrV.matrixR();
	for(int i = 0; i<r.rows(); i++)
	  for(int j = 0; j<r.cols(); j++)
	  {
	    __R(i,j) = r(i,j);
	  }
	  
	DataValues __b(__Qt.rows(), 0.0);
        multiply(__Qt, _residuals, __b); 
        __b *= -1.0; // = -Qr
        DataValues __c(__b.size(), 0.0);
        backSubstitution(__R, __b, __c);
	
	// validation
	DataValues update(_residuals.size(), 0.0);
	multiply(_matrixW, __c, update); 
	for(int i = 0; i<xUpdate.size(); i++)
	{
	  if(!tarch::la::equals(xUpdate(i), update(i), 1e-12))
	  {
	    std::cerr<<"xUpdates were not the same for standart and updatedQR decomposition:\n"<<"standart:\n   "<<xUpdate<<"updated:\n   "<<update<<std::endl;
	    _infostream<<"validation failed for xUpdate.\n   - modifiedGramSchmidt:\n"<<xUpdate<<"\n  - updatedQR:\n"<<update<<"\n"<<std::flush;
	  }
	}
	// ------------------------------------------------
	
        preciceDebug("c = " << c);
      }
    }
    
    // Perform QN relaxation for secondary data
    foreach (int id, _secondaryDataIDs){
      PtrCouplingData data = cplData[id];
      DataValues& values = *(data->values);
      assertion2(_secondaryMatricesW[id].cols() == c.size(),
                 _secondaryMatricesW[id].cols(), c.size());
      multiply(_secondaryMatricesW[id], c, values);
      assertion2(values.size() == data->oldValues.column(0).size(),
                 values.size(), data->oldValues.column(0).size());
      values += data->oldValues.column(0);
      assertion2(values.size() == _secondaryResiduals[id].size(),
                 values.size(), _secondaryResiduals[id].size());
      values += _secondaryResiduals[id];
    }
}


void IQNILSPostProcessing:: specializedIterationsConverged
(
   DataMap & cplData)
{
  
  if (_matrixCols.front() == 0){ // Did only one iteration
    _matrixCols.pop_front(); 
  }
  
  if (_timestepsReused == 0){
    foreach (int id, _secondaryDataIDs){
      _secondaryMatricesW[id].clear();
    }
  }
  else if ((int)_matrixCols.size() > _timestepsReused){
    int toRemove = _matrixCols.back();
    foreach (int id, _secondaryDataIDs){
      DataMatrix& secW = _secondaryMatricesW[id];
      assertion3(secW.cols() > toRemove, secW, toRemove, id);
      for (int i=0; i < toRemove; i++){
        secW.remove(secW.cols() - 1);
      }
    }
  }
  
}


void IQNILSPostProcessing:: removeMatrixColumn
(
  int columnIndex)
{
  assertion(_matrixV.cols() > 1);
  // remove column from secondary Data Matrix W
 foreach (int id, _secondaryDataIDs){
    _secondaryMatricesW[id].remove(columnIndex);
  }
  
  BaseQNPostProcessing::removeMatrixColumn(columnIndex);
}

}}} // namespace precice, cplscheme, impl
