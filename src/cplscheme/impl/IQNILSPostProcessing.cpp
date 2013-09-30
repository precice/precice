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
//#include "utils/NumericalCompare.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

tarch::logging::Log IQNILSPostProcessing::
      _log("precice::cplscheme::impl::IQNILSPostProcessing");

IQNILSPostProcessing:: IQNILSPostProcessing
(
  double initialRelaxation,
  int    maxIterationsUsed,
  int    timestepsReused,
  double singularityLimit,
  std::vector<int>    dataIDs)
:
  PostProcessing(),
  _initialRelaxation(initialRelaxation),
  _maxIterationsUsed(maxIterationsUsed),
  _timestepsReused(timestepsReused),
  _singularityLimit(singularityLimit),
  _dataIDs(dataIDs),
  _firstIteration(true),
  _oldXTilde(),
  _oldResiduals(),
  _matrixV(),
  _matrixW(),
  _matrixCols()
{
   preciceCheck((_initialRelaxation > 0.0) && (_initialRelaxation <= 1.0),
                "IQNILSPostProcessing()",
                "Initial relaxation factor for IQN-ILS post-processing has to "
                << "be larger than zero and smaller or equal than one!");
   preciceCheck(_maxIterationsUsed > 0, "IQNILSPostProcessing()",
                "Maximal iterations used for IQN-ILS post-processing has to "
                << "be larger than zero!");
   preciceCheck(_timestepsReused >= 0, "IQNILSPostProcessing()",
                "Number of old timesteps to be reused for IQN-ILS "
                 << "post-processing has to be >= 0!");
   preciceCheck(tarch::la::greater(_singularityLimit, 0.0),
                "IQNILSPostProcessing()", "Singularity limit for IQN-ILS "
                << "post-processing has to be larger than numerical zero ("
                << tarch::la::NUMERICAL_ZERO_DIFFERENCE << ")!");
}

void IQNILSPostProcessing:: initialize
(
   DataMap& cplData )
{
   preciceCheck(utils::contained(*_dataIDs.begin(), cplData), "initialize()",
                "Data with ID " << *_dataIDs.begin() << " is not contained in data "
                "given at initialization!");

   size_t entries=0;
   if(_dataIDs.size()==1){
     entries = cplData[_dataIDs.at(0)]->values->size();
   }
   else{
     assertion(_dataIDs.size()==2);
     entries = cplData[_dataIDs.at(0)]->values->size() +
         cplData[_dataIDs.at(1)]->values->size();
   }

   assertion(entries > 0);
   double init = 0.0;
   assertion(_oldXTilde.size() == 0);
   assertion(_oldResiduals.size() == 0);
   _oldXTilde.append(DataValues(entries, init));
   _oldResiduals.append(DataValues(entries, init));
   _matrixCols.push_front(0);
   _firstIteration = true;

   // Append column for old values if not done by coupling scheme yet
   foreach (DataMap::value_type& pair, cplData){
     int cols = pair.second->oldValues.cols();
     if (cols < 1){
       assertion1(pair.second->values->size() > 0, pair.first);
       pair.second->oldValues.append(CouplingData::DataMatrix(
         pair.second->values->size(), 1, 0.0));
     }
   }
}

void IQNILSPostProcessing:: performPostProcessing
(
   DataMap& cplData)
{
  preciceTrace("performPostProcessing()");
  //there is already a preciceCheck in cplscheme->initialize
  assertion1(_dataIDs.size()<=2 && _dataIDs.size()>=1,
                      "The number of coupling data vectors should be 1 or 2");
  using namespace tarch::la;
  assertion1(utils::contained(*_dataIDs.begin(), cplData), *_dataIDs.begin());
  assertion2(_oldResiduals.size() == _oldXTilde.size(),
             _oldResiduals.size(), _oldXTilde.size());

  DataValues values;
  DataValues oldValues;
  preciceDebug("dataId size " << _dataIDs.size());
  preciceDebug("cplData size " << cplData.size());
  foreach (int id, _dataIDs){
    values.append(*(cplData[id]->values));
    oldValues.append(cplData[id]->oldValues.column(0));
  }

  //preciceDebug("Untouched values = " << values);
  //preciceDebug("Old values = " << oldValues);

  // Compute current residual: vertex-data - oldData
  DataValues residuals(values);
  residuals -= oldValues;

  if (_firstIteration && (_matrixCols.size() < 2)){
    preciceDebug("   Performing underrelaxation");
    _oldXTilde = values; // Store x tilde
    _oldResiduals = residuals; // Store current residual
    // Perform underrelaxation with residual: x_new = x_old + omega * res
    residuals *= _initialRelaxation;
    residuals += oldValues;
    values = residuals;
    //      precicePrint("   Performing constant relaxation with omg = " << _initialRelaxation);

    // Perform underrelaxation with initial relaxation factor for rest of data.
    foreach (DataMap::value_type& pair, cplData){
      if (pair.first != *_dataIDs.begin() && pair.first != *(_dataIDs.end()-1)){
        //            precicePrint("   More data ...");
        DataValues & values = *pair.second->values;
        values *= _initialRelaxation;                   // new * omg
        DataValues & residuals = pair.second->oldValues.column(0);  // old
        residuals *= 1.0 - _initialRelaxation;          // (1-omg) * old
        values += residuals;                   // (1-omg) * old + new * omg
      }
    }
  }
  else {
    preciceDebug("   Performing QN step");

    if (not _firstIteration){ // Update matrices V, W with newest information
      assertion2(_matrixV.cols() == _matrixW.cols(), _matrixV.cols(), _matrixW.cols());
      assertion2(_matrixV.cols() <= _maxIterationsUsed,
                 _matrixV.cols(), _maxIterationsUsed);
      bool columnLimitReached = _matrixV.cols() == _maxIterationsUsed;
      bool overdetermined = _matrixV.cols() <= _matrixV.rows();
      if (not columnLimitReached && overdetermined){
        //preciceDebug("performPostProcessing()",
        //               "   Adding column to V,W. Cols = " << _matrixV.cols());
        _matrixV.appendFront(residuals); // Will be modified to delta_r
        _matrixW.appendFront(residuals); // Will be overwritten by delta_x_tilde
        _matrixCols.front() ++;
        //precicePrint("   Added column to V,W. V.col(0)= " << _matrixV.column(0));
      }
      else {
        //preciceDebug("performPostProcessing()", "   Shift inserting new column to V,W");
        _matrixV.shiftSetFirst(residuals); // Will be modified to delta_r
        _matrixW.shiftSetFirst(residuals); // Will be overwritten by delta_x_tilde
        _matrixCols.front() ++;
        _matrixCols.back() --;
        if (_matrixCols.back() == 0){
          _matrixCols.pop_back ();
        }
      }
      // Compute delta_residual = residual - residual_old
      _matrixV.column(0) -= _oldResiduals;

      // Compute delta_x_tilde = x_tilde - x_tilde_old
      _matrixW.column(0) = values;
      _matrixW.column(0) -= _oldXTilde;
    }

    _oldResiduals = residuals;   // Store residuals
    _oldXTilde = values;         // Store x_tilde

    DataValues Wc(residuals.size(), 0.0);

    // Calculate QR decomposition of matrix V and solve Rc = -Qr
    //precicePrint("V.cols() = " << _matrixV.cols());
    //precicePrint("V.rows() = " << _matrixV.rows());
    bool linearDependence = true;
    while (linearDependence){
      preciceDebug("   Compute Newton factors");
      linearDependence = false;
      DataMatrix Vcopy(_matrixV);
      //precicePrint("Vcopy.cols() = " << Vcopy.cols());
      //precicePrint("Vcopy.rows() = " << Vcopy.rows());
      //precicePrint("Vcopy = " << Vcopy);
      DataMatrix Q(Vcopy.rows(), Vcopy.cols(), 0.0);
      DataMatrix R(Vcopy.cols(), Vcopy.cols(), 0.0);
      //precicePrint("Q.cols() = " << Q.cols() << ", Q.rows() = " << Q.rows());
      modifiedGramSchmidt(Vcopy, Q, R);
      if (_matrixV.cols() > 1){
        for (int i=0; i < _matrixV.cols(); i++){
          if (R(i,i) < _singularityLimit){
            preciceDebug("   Removing linear dependent column " << i);
            linearDependence = true;
            removeMatrixColumn(i);
          }
        }
      }
      if (not linearDependence){
        preciceDebug("   Apply Newton factors");
        //preciceDebug("performPostProcessing()", "   Q = " << Q);
        DataValues b(Q.cols(), 0.0);
        //precicePrint("Q.cols() = " << Q.cols() << ", Q.rows() = " << Q.rows());
        //precicePrint("Q_T.cols() = " << transpose(Q).cols() << ", Q_T.rows() = " << transpose(Q).rows());
        multiply(transpose(Q), residuals, b); // = Qr
        b *= -1.0; // = -Qr
        DataValues c(b.size(), 0.0);
        //preciceDebug("performPostProcessing()", "   R = " << R);
        //preciceDebug("performPostProcessing()", "   b = " << b);
        backSubstitution(R, b, c);
        multiply(_matrixW, c, Wc); // = Wc
        //preciceDebug("performPostProcessing()", "   _matrixW = " << _matrixW);
        preciceDebug("c = " << c);
      }
    }
    //preciceDebug("performPostProcessing()", "   oldValues = " << oldValues);
    values = oldValues;  // = x^k
    values += Wc;        // = x^k + Wc
    values += residuals; // = x^k + Wc + r^k

    foreach (DataMap::value_type & pair, cplData){
      if (pair.first != *_dataIDs.begin() && pair.first != *(_dataIDs.end()-1)){
        DataValues & residuals = *pair.second->values;                 // = x_tilde
        residuals -= pair.second->oldValues.column(0); // = x_tilde - x^k = r^k
        // Perform update with Wc
        residuals += Wc;                                 // = Wc + r^k
        residuals += pair.second->oldValues.column(0); // = x^k + Wc + r^k
        *pair.second->values = residuals;
      }
    }
    //preciceDebug("performPostprocessing()", "Postprocessed values = " << values);
  }


  // Set all values back from copies to original
  int offset = 0;
  foreach(int id, _dataIDs){
    int size = cplData[id]->values->size();
    preciceDebug("Copying values back, size: " << size);
    utils::DynVector& valuesPart = *(cplData[id]->values);
    utils::DynVector& oldValuesPart = cplData[id]->oldValues.column(0);
    for(int i=0; i<size; i++){
      //preciceDebug("Copying values back, values, id: " << id <<" i: " << i);
      valuesPart[i] = values[i + offset];
      oldValuesPart[i] = oldValues[i + offset];
    }
    offset += size;
  }
  //preciceDebug("copied values back, values size: " << cplData[*_dataIDs.begin()]->values->size());
  //preciceDebug("copied values back, oldValues size: " << cplData[*_dataIDs.begin()]->oldValues.column(0).size());
  _firstIteration = false;
}

void IQNILSPostProcessing:: iterationsConverged
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
    _matrixCols.pop_back();
  }
  _matrixCols.push_front(0);
  _firstIteration = true;
}

void IQNILSPostProcessing:: exportState(io::TXTWriter& writer)
{
  tarch::la::Vector<1,int> colSize(_matrixCols.size());
  writer.write(colSize);
  writer.write(_matrixCols);
  writer.write(_oldXTilde);
  writer.write(_oldResiduals);
  writer.write(_matrixV);
  writer.write(_matrixW);
}

void IQNILSPostProcessing:: importState(io::TXTReader& reader)
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

void IQNILSPostProcessing:: removeMatrixColumn
(
   int columnIndex)
{
  preciceTrace2("removeMatrixColumn()", columnIndex, _matrixV.cols());

   // Remove matrix columns
   assertion(_matrixV.cols() > 1);
   _matrixV.remove(columnIndex);
   _matrixW.remove(columnIndex);

   // Reduce column count
   std::deque<int>::iterator iter = _matrixCols.begin();
   int cols = 0;
   while(iter != _matrixCols.end()) {
      cols += *iter;
//      precicePrint("   cols = " << cols << ", iter = " << *iter);
      if(cols > columnIndex) {
//         precicePrint("   Decreasing iter");
         assertion(*iter > 0);
         *iter -= 1;
         if(*iter == 0) {
            _matrixCols.erase(iter);
         }
         break;
      }
      iter ++;
   }
}

}}} // namespace precice, cplscheme, impl
