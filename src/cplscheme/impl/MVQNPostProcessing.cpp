// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "MVQNPostProcessing.hpp"
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
//#include "utils/NumericalCompare.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

tarch::logging::Log MVQNPostProcessing::
      _log("precice::cplscheme::impl::MVQNPostProcessing");

      
// ------------ UNTOUCHED --------------------------------------------------------------
// except secondary data, output in case of error
// -------------------------------------------------------------------------------------
MVQNPostProcessing:: MVQNPostProcessing
(
  double initialRelaxation,
  int    maxIterationsUsed,
  int    timestepsReused,
  double singularityLimit,
  std::vector<int> dataIDs,
  std::map<int,double> scalings)
:
  PostProcessing(),
  _initialRelaxation(initialRelaxation),
  _maxIterationsUsed(maxIterationsUsed),
  _timestepsReused(timestepsReused),
  _singularityLimit(singularityLimit),
  _dataIDs(dataIDs),
//  _secondaryDataIDs(),
  _scalings(scalings),
  _firstIteration(true),
  _oldXTilde(),
//  _secondaryOldXTildes(),
  _residuals(),
//  _secondaryResiduals(),
  _scaledValues(),
  _scaledOldValues(),
  _oldResiduals(),
  _matrixV(),
  _matrixW(),
  _invJacobian(),
  _oldInvJacobian(),
//  _secondaryMatricesW(),
  _matrixCols()
{
   preciceCheck((_initialRelaxation > 0.0) && (_initialRelaxation <= 1.0),
                "MVQNPostProcessing()",
                "Initial relaxation factor for MVQN post-processing has to "
                << "be larger than zero and smaller or equal than one!");
   preciceCheck(_maxIterationsUsed > 0, "MVQNPostProcessing()",
                "Maximal iterations used for MVQN post-processing has to "
                << "be larger than zero!");
   preciceCheck(_timestepsReused >= 0, "MVQNPostProcessing()",
                "Number of old timesteps to be reused for MVQN "
                 << "post-processing has to be >= 0!");
   preciceCheck(tarch::la::greater(_singularityLimit, 0.0),
                "MVQNPostProcessing()", "Singularity limit for MVQN "
                << "post-processing has to be larger than numerical zero ("
                << tarch::la::NUMERICAL_ZERO_DIFFERENCE << ")!");
}


// ------------ UNTOUCHED --------------------------------------------------------------
// except secondary data
// -------------------------------------------------------------------------------------
void MVQNPostProcessing:: initialize
(
  DataMap& cplData )
{
  preciceTrace1("initialize()", cplData.size());
  size_t entries=0;

  for(size_t i=0;i<_dataIDs.size();i++){
    preciceCheck(utils::contained(_dataIDs[i], cplData), "initialize()",
                 "Data with ID " << _dataIDs[i] << " is not contained in data "
                 "given at initialization!");
    entries += cplData[_dataIDs[i]]->values->size();
  }



  assertion(entries > 0);
  double init = 0.0;
  assertion(_oldXTilde.size() == 0);
  assertion(_oldResiduals.size() == 0);
  _oldXTilde.append(DataValues(entries, init));
  _oldResiduals.append(DataValues(entries, init));
  _residuals.append(DataValues(entries, init));
  _scaledValues.append(DataValues(entries, init));
  _scaledOldValues.append(DataValues(entries, init));
  _invJacobian.append(DataMatrix(entries, entries, init));
  _oldInvJacobian.append(DataMatrix(entries, entries, init));
  _matrixCols.push_front(0);
  _firstIteration = true;

  // Fetch secondary data IDs, to be relaxed with same coefficients from MVQN
//  foreach (DataMap::value_type& pair, cplData){
//    if (not utils::contained(pair.first, _dataIDs)){
//      _secondaryDataIDs.push_back(pair.first);
//      int secondaryEntries = pair.second->values->size();
//      _secondaryOldXTildes[pair.first].append(DataValues(secondaryEntries, init));
//      _secondaryResiduals[pair.first].append(DataValues(secondaryEntries, init));
//    }
//  }

  // Append old value columns, if not done outside of post-processing already
  foreach (DataMap::value_type& pair, cplData){
    int cols = pair.second->oldValues.cols();
    if (cols < 1){ // Add only, if not already done
      assertion1(pair.second->values->size() > 0, pair.first);
      pair.second->oldValues.append(
        CouplingData::DataMatrix(pair.second->values->size(), 1, 0.0));
    }
  }
}

void MVQNPostProcessing:: performPostProcessing
(
  DataMap& cplData)
{
  preciceTrace2("performPostProcessing()", _dataIDs.size(), cplData.size());
  using namespace tarch::la;
  assertion2(_dataIDs.size() == _scalings.size(), _dataIDs.size(), _scalings.size());
  assertion2(_oldResiduals.size() == _oldXTilde.size(),
             _oldResiduals.size(), _oldXTilde.size());
  assertion2(_scaledValues.size() == _oldXTilde.size(),
             _scaledValues.size(), _oldXTilde.size());
  assertion2(_scaledOldValues.size() == _oldXTilde.size(),
             _scaledOldValues.size(), _oldXTilde.size());
  assertion2(_residuals.size() == _oldXTilde.size(),
             _residuals.size(), _oldXTilde.size());

  int offset = 0;
  foreach (int id, _dataIDs){
    double factor = _scalings[id];
    preciceDebug("Scaling Factor " << factor << " for id: " << id);
    int size = cplData[id]->values->size();
    DataValues& values = *cplData[id]->values;
    DataValues& oldValues = cplData[id]->oldValues.column(0);
    for (int i=0; i < size; i++){
      _scaledValues[i+offset] = values[i]/factor;
      _scaledOldValues[i+offset] = oldValues[i]/factor;
    }
    offset += size;
  }

  // Compute current residual: vertex-data - oldData
  _residuals = _scaledValues;
  _residuals -= _scaledOldValues;

  // Compute residuals of secondary data
//  foreach (int id, _secondaryDataIDs){
//    DataValues& secResiduals = _secondaryResiduals[id];
//    PtrCouplingData data = cplData[id];
//    assertion2(secResiduals.size() == data->values->size(),
//               secResiduals.size(), data->values->size());
//    secResiduals = *(data->values);
//    secResiduals -= data->oldValues.column(0);
//  }

  if (_firstIteration && (_matrixCols.size() < 2)){
    preciceDebug("   Performing underrelaxation");
    _oldInvJacobian = _invJacobian; // store inverse Jacobian
    _oldXTilde = _scaledValues; // Store x tilde
    _oldResiduals = _residuals; // Store current residual
    // Perform underrelaxation with residual: x_new = x_old + omega * res
    _residuals *= _initialRelaxation;
    _residuals += _scaledOldValues;
    _scaledValues = _residuals;

    // Store x_tildes for secondary data
//    foreach (int id, _secondaryDataIDs){
//      assertion2(_secondaryOldXTildes[id].size() == cplData[id]->values->size(),
//                 _secondaryOldXTildes[id].size(), cplData[id]->values->size());
//      _secondaryOldXTildes[id] = *(cplData[id]->values);
//    }

    // Perform underrelaxation with initial relaxation factor for rest of data.
//    foreach (int id, _secondaryDataIDs){
//      PtrCouplingData data = cplData[id];
//      DataValues& values = *(data->values);
//      values *= _initialRelaxation;                   // new * omg
//      DataValues& secResiduals = _secondaryResiduals[id];
//      secResiduals = data->oldValues.column(0);    // old
//      secResiduals *= 1.0 - _initialRelaxation;       // (1-omg) * old
//      values += secResiduals;                      // (1-omg) * old + new * omg
//    }
    
  } else {
    preciceDebug("   Performing QN step");

    if (not _firstIteration){ 
      // Update matrices V, W with newest information
      assertion2(_matrixV.cols() == _matrixW.cols(), _matrixV.cols(), _matrixW.cols());
      assertion2(_matrixV.cols() <= _maxIterationsUsed,
                 _matrixV.cols(), _maxIterationsUsed);
      bool columnLimitReached = _matrixV.cols() == _maxIterationsUsed;
      bool overdetermined = _matrixV.cols() <= _matrixV.rows();
      
      if (not columnLimitReached && overdetermined){
        _matrixV.appendFront(_residuals); // Will be modified to delta_r
        _matrixW.appendFront(_residuals); // Will be overwritten by delta_x_tilde
        // Append column for secondary W matrices
//        foreach (int id, _secondaryDataIDs){
//          _secondaryMatricesW[id].appendFront(_secondaryResiduals[id]);
//        }
        _matrixCols.front()++;
      }
      else {
        _matrixV.shiftSetFirst(_residuals); // Will be modified to delta_r
        _matrixW.shiftSetFirst(_residuals); // Will be overwritten by delta_x_tilde
        // Shift column for secondary W matrices
//        foreach (int id, _secondaryDataIDs){
//          _secondaryMatricesW[id].shiftSetFirst(_secondaryResiduals[id]);
//        }
        _matrixCols.front()++;
        _matrixCols.back()--;
        if (_matrixCols.back() == 0){
          _matrixCols.pop_back();
        }
      }
      
      // Compute delta_residual = residual - residual_old
      _matrixV.column(0) -= _oldResiduals;
      // Compute delta_x_tilde = x_tilde - x_tilde_old
      _matrixW.column(0) = _scaledValues;
      _matrixW.column(0) -= _oldXTilde;
      
      // Compute delta_x_tilde for secondary data
//      foreach (int id, _secondaryDataIDs){
//        DataMatrix& secW = _secondaryMatricesW[id];
//        assertion2(secW.column(0).size() == cplData[id]->values->size(),
//                   secW.column(0).size(), cplData[id]->values->size());
//        secW.column(0) = *(cplData[id]->values);
//        secW.column(0) -= _secondaryOldXTildes[id];
//      }
    }

    _oldResiduals = _residuals;   // Store residuals
    _oldXTilde = _scaledValues;   // Store x_tilde
//    foreach (int id, _secondaryDataIDs){ // Store x_tildes for secondary data
//      assertion2(_secondaryOldXTildes[id].size() == cplData[id]->values->size(),
//                 _secondaryOldXTildes[id].size(), cplData[id]->values->size());
//      _secondaryOldXTildes[id] = *(cplData[id]->values);
//    }



    // ------------- update inverse Jacobian -----------
    // J_inv = J_inv_n + (W - J_inv_n*V)*(V^T*V)^-1*V^T
    // ----------------------------------------- -------

    DataValues xUpdate(_residuals.size(), 0.0);

    preciceDebug("   Compute Newton factors");
    DataMatrix VTV(_matrixV.cols(), _matrixV.cols(), 0.0);
    DataMatrix Q(_matrixV.cols(), _matrixV.cols(), 0.0);
    DataMatrix R(_matrixV.cols(), _matrixV.cols(), 0.0);
    DataMatrix v;
    multiply(transpose(_matrixV), _matrixV, VTV);  // VTV = V^T*V
    modifiedGramSchmidt(VTV, Q, R);
    
    preciceDebug("   Apply Newton factors");
    DataMatrix QTVT(_matrixV.cols(), _matrixV.rows(), 0.0);
    multiply(transpose(Q), transpose(_matrixV), QTVT); // QTVT = Q^T*V^T
    DataValues tmp(_matrixV.cols(), 0.0);
    
    // v = (V^T*V)^-1*V^T
    for(int i = 0; i < QTVT.cols(); i++)
    {
      backSubstitution(R, QTVT.column(i), tmp);
      v.append(tmp);  
    }
    // JCopy = J_inv_n*V
    DataMatrix JCopy(_invJacobian);
    multiply(_oldInvJacobian, _matrixV, JCopy);
    // JCopy = (W-J_inv_n*V)
    JCopy *= -1.;
    JCopy = JCopy + _matrixW;
    // invJacobian = (W - J_inv_n*V)*(V^T*V)^-1*V^T
    multiply(JCopy, v, _invJacobian);
    _invJacobian = _invJacobian + _oldInvJacobian;
    
    tmp = _residuals;
    tmp *= -1.;
    
    // solve delta_x = - J_inv*residuals
    multiply(_invJacobian, tmp, xUpdate); 
  
    //preciceDebug("performPostProcessing()", "   oldValues = " << oldValues);
    _scaledValues = _scaledOldValues;  // = x^k
    _scaledValues += xUpdate;        // = x^k + Wc
    _scaledValues += _residuals; // = x^k + Wc + r^k

    // Perform QN relaxation for secondary data
//    foreach (int id, _secondaryDataIDs){
//      PtrCouplingData data = cplData[id];
//      DataValues& values = *(data->values);
//      assertion2(_secondaryMatricesW[id].cols() == c.size(),
//                 _secondaryMatricesW[id].cols(), c.size());
//      multiply(_secondaryMatricesW[id], c, values);
//      assertion2(values.size() == data->oldValues.column(0).size(),
//                 values.size(), data->oldValues.column(0).size());
//      values += data->oldValues.column(0);
//      assertion2(values.size() == _secondaryResiduals[id].size(),
//                 values.size(), _secondaryResiduals[id].size());
//      values += _secondaryResiduals[id];
//    }
  }

  // Undo scaling of values and overwrite originals
  offset = 0;
  foreach(int id, _dataIDs){
    double factor = _scalings[id];
    int size = cplData[id]->values->size();
    preciceDebug("Copying values back, size: " << size);
    utils::DynVector& valuesPart = *(cplData[id]->values);
    utils::DynVector& oldValuesPart = cplData[id]->oldValues.column(0);
    for(int i=0; i < size; i++){
      valuesPart[i] = _scaledValues[i+offset]*factor;
      oldValuesPart[i] = _scaledOldValues[i+offset]*factor;
    }
    offset += size;
  }
  _firstIteration = false;
}


// ------------ UNTOUCHED --------------------------------------------------------------
// except secondary data
// -------------------------------------------------------------------------------------
void MVQNPostProcessing:: iterationsConverged
(
   DataMap & cplData)
{
  preciceTrace("iterationsConverged()");
# ifdef Debug
  std::ostringstream stream;
  stream << "Matrix column counters: ";
  foreach (int cols, _matrixCols){
    stream << cols << ", ";
  }
  preciceDebug(stream.str());
# endif // Debug

  if (_matrixCols.front() == 0){ // Did only one iteration
    _matrixCols.pop_front();
  }

  if (_timestepsReused == 0){
    //precicePrint("Removing all columns from V, W");
    _matrixV.clear();
    _matrixW.clear();
//    foreach (int id, _secondaryDataIDs){
//      _secondaryMatricesW[id].clear();
//    }
    _matrixCols.clear();
  }
  else if ((int)_matrixCols.size() > _timestepsReused){
    int toRemove = _matrixCols.back();
    assertion1(toRemove > 0, toRemove);
    preciceDebug("Removing " << toRemove << " cols from IQN-ILS matrices with "
                 << _matrixV.cols() << " cols");
    assertion2(_matrixV.cols() == _matrixW.cols(), _matrixV.cols(), _matrixW.cols());
    assertion2(_matrixV.cols() > toRemove, _matrixV.cols(), toRemove);
    for (int i=0; i < toRemove; i++){
      _matrixV.remove(_matrixV.cols() - 1);
      _matrixW.remove(_matrixW.cols() - 1);
    }
//    foreach (int id, _secondaryDataIDs){
//      DataMatrix& secW = _secondaryMatricesW[id];
//      assertion3(secW.cols() > toRemove, secW, toRemove, id);
//      for (int i=0; i < toRemove; i++){
//        secW.remove(secW.cols() - 1);
//      }
//    }
    _matrixCols.pop_back();
  }
  _matrixCols.push_front(0);
  _firstIteration = true;
}

// ------------ UNTOUCHED --------------------------------------------------------------
void MVQNPostProcessing:: exportState(io::TXTWriter& writer)
{
//  tarch::la::Vector<1,int> colSize(_matrixCols.size());
//  writer.write(colSize);
//  writer.write(_matrixCols);
//  writer.write(_oldXTilde);
//  writer.write(_oldResiduals);
//  writer.write(_matrixV);
//  writer.write(_matrixW);
}

// ------------ UNTOUCHED --------------------------------------------------------------
void MVQNPostProcessing:: importState(io::TXTReader& reader)
{
//  tarch::la::Vector<1,int> colSize;
//  reader.read(colSize);
//  _matrixCols.resize(colSize[0]);
//  reader.read(_oldXTilde);
//  reader.read(_oldResiduals);
//  for (int i=0; i<colSize[0]; i++ ){
//    // Using _oldXTilde to have a vector with appropriate size.
//    // Values are overwritten afterwards in file read.
//    _matrixV.append(_oldXTilde);
//    _matrixW.append(_oldXTilde);
//  }
//  reader.read(_matrixV);
//  reader.read(_matrixW);
//  if (colSize[0] > 1){
//    _firstIteration = false;
//  }
}

// ------------ UNTOUCHED --------------------------------------------------------------
void MVQNPostProcessing:: removeMatrixColumn
(
  int columnIndex)
{
  preciceTrace2("removeMatrixColumn()", columnIndex, _matrixV.cols());

  // Remove matrix columns
  assertion(_matrixV.cols() > 1);
  _matrixV.remove(columnIndex);
  _matrixW.remove(columnIndex);
//  foreach (int id, _secondaryDataIDs){
//    _secondaryMatricesW[id].remove(columnIndex);
//  }

  // Reduce column count
  std::deque<int>::iterator iter = _matrixCols.begin();
  int cols = 0;
  while(iter != _matrixCols.end()) {
    cols += *iter;
    if(cols > columnIndex) {
      assertion(*iter > 0);
      *iter -= 1;
      if(*iter == 0) {
        _matrixCols.erase(iter);
      }
      break;
    }
    iter++;
  }
}

}}} // namespace precice, cplscheme, impl
