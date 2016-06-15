// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#pragma once

#include "utils/Dimensions.hpp"
#include "utils/Helpers.hpp"
#include "utils/Globals.hpp"
#include "tarch/la/DynamicColumnMatrix.h"
#include <Eigen/Dense>
#include "../SharedPointer.hpp"
#include <vector>
#include "utils/MasterSlave.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

/**
 * @brief Interface for preconditioner variants that can be applied to quasi-Newton post-processing schemes.
 *
 * Preconditioning (formerly also known as scaling) improves the balance between several sub-vectors for parallel or multi coupling.
 *
 * apply() applies the weighting, i.e. transforms from physical values to balanced values.
 * revert() reverts the weighting, i.e. transforms from balanced values back to physical values.
 * update() updates the preconditioner, after every FSI iteration (though some variants might only be updated after a complete timestep)
 */
class Preconditioner
{
public:

  //see post-processing definitions
  typedef tarch::la::DynamicVector<double> DataValues;
  typedef std::map<int,PtrCouplingData> DataMap;
  typedef tarch::la::DynamicColumnMatrix<double> DataMatrix;
  typedef Eigen::MatrixXd EigenMatrix;

  Preconditioner(
      std::vector<int> dimensions,
      int maxNonConstTimesteps)
  :
    _weights(),
    _invWeights(),
    _dimensions(dimensions),
    _sizeOfSubVector(-1),
    _maxNonConstTimesteps(maxNonConstTimesteps),
    _nbNonConstTimesteps(0),
    _requireNewQR(false),
    _freezed(false)
  {}


  /**
   * @brief Destructor, empty.
   */
  virtual ~Preconditioner() {}


  /**
   * @brief initialize the preconditioner
   * @param size of the pp system (e.g. rows of V)
   */
  virtual void initialize(int N){
    preciceTrace1("initialize()", N);

    assertion(_weights.size()==0);
    // cannot do this already in the constructor as the size is unknown at that point
    _weights.resize(N, 1.0);
    _invWeights.resize(N, 1.0);

    int numberOfParts = 0;
    for (int dim : _dimensions){
      numberOfParts += dim;
    }
    _sizeOfSubVector = N / numberOfParts;
    assertion(numberOfParts * _sizeOfSubVector == N);
  }

  /**
   * @brief To transform physical values to balanced values. Matrix version
   */
  void apply(DataMatrix& M){
    preciceTrace("apply()");

    assertion(M.column(0).size()==(int)_weights.size());

    // scale matrix M
    for(int i=0; i<M.cols(); i++){
      for(int j=0; j<M.column(0).size(); j++){
        M(j,i) *= _weights[j];
      }
    }
  }

  /**
   * @brief Apply preconditioner to matrix
   * @param transpose: false = from left, true = from right
   */
  void apply(EigenMatrix& M, bool transpose){
    preciceTrace(__func__);
    if(transpose){
      assertion(M.cols()==(int)_weights.size(), M.cols(), _weights.size());
      for(int i=0; i<M.cols(); i++){
        for(int j=0; j<M.rows(); j++){
          M(j,i) *= _weights[i];
        }
      }
    }
    else{
      assertion(M.rows()==(int)_weights.size(), M.rows(), (int)_weights.size());
      for(int i=0; i<M.cols(); i++){
        for(int j=0; j<M.rows(); j++){
          M(j,i) *= _weights[j];
        }
      }
    }
  }

  /**
   * @brief Apply inverse preconditioner to matrix
   * @param transpose: false = from left, true = from right
   */
  void revert(EigenMatrix& M, bool transpose){
    preciceTrace(__func__);
    //assertion(_needsGlobalWeights);
    if (transpose) {
      assertion(M.cols()==(int)_invWeights.size());
      for (int i = 0; i < M.cols(); i++) {
        for (int j = 0; j < M.rows(); j++) {
          M(j, i) *= _invWeights[i];
        }
      }
    }
    else {
      assertion(M.rows()==(int)_invWeights.size(), M.rows(), (int)_invWeights.size());
      for (int i = 0; i < M.cols(); i++) {
        for (int j = 0; j < M.rows(); j++) {
          M(j, i) *= _invWeights[j];
        }
      }
    }
  }



  /**
   * @brief To transform physical values to balanced values. Vector version
   */
  void apply(DataValues& v){
    preciceTrace("apply()");
    assertion(v.size()==(int)_weights.size());

    // scale residual
    for(int j=0; j<v.size(); j++){
      v[j] *= _weights[j];
    }
  }

  /**
   * @brief To transform balanced values back to physical values. Matrix version
   */
  void revert(DataMatrix& M){
    preciceTrace("revert()");
    assertion(M.column(0).size()==(int)_weights.size());

    // scale matrix M
    for(int i=0; i<M.cols(); i++){
      for(int j=0; j<M.column(0).size(); j++){
        M(j,i) *= _invWeights[j];
      }
    }
  }

  /**
   * @brief To transform balanced values back to physical values. Vector version
   */
  void revert(DataValues& v){
    preciceTrace("revert()");
    assertion(v.size()==(int)_weights.size());

    // scale residual
    for(int j=0; j<v.size(); j++){
      v[j] *= _invWeights[j];
    }
  }

  /**
     * @brief To transform physical values to balanced values. Matrix version
     */
    void apply(Eigen::MatrixXd& M){
      preciceTrace("apply()");
      assertion(M.rows()==(int)_weights.size(), M.rows(), (int)_weights.size());

      // scale matrix M
      for(int i=0; i<M.cols(); i++){
        for(int j=0; j<M.rows(); j++){
          M(j,i) *= _weights[j];
        }
      }
    }

    /**
     * @brief To transform physical values to balanced values. Vector version
     */
    void apply(Eigen::VectorXd& v){
      preciceTrace("apply()");

      assertion(v.size()==(int)_weights.size());

      // scale residual
      for(int j=0; j<v.size(); j++){
        v[j] *= _weights[j];
      }
    }

    /**
     * @brief To transform balanced values back to physical values. Matrix version
     */
    void revert(Eigen::MatrixXd& M){
      preciceTrace("revert()");

      assertion(M.rows()==(int)_weights.size());

      // scale matrix M
      for(int i=0; i<M.cols(); i++){
        for(int j=0; j<M.rows(); j++){
          M(j,i) *= _invWeights[j];
        }
      }
    }

    /**
     * @brief To transform balanced values back to physical values. Vector version
     */
    void revert(Eigen::VectorXd& v){
      preciceTrace("revert()");

      assertion(v.size()==(int)_weights.size());

      // scale residual
      for(int j=0; j<v.size(); j++){
        v[j] *= _invWeights[j];
      }
    }

  /**
   * @brief Update the scaling after every FSI iteration and require a new QR decomposition (if necessary)
   *
   * @param timestepComplete [IN] True if this FSI iteration also completed a timestep
   */
  void update(bool timestepComplete, const Eigen::VectorXd& oldValues, const Eigen::VectorXd& res){
    preciceTrace2(__func__, _nbNonConstTimesteps, _freezed);

    // if number of allowed non-const time steps is exceeded, do not update weights
    if(_freezed)
     return;

    // increment number of time steps that has been scaled with changing preconditioning weights
    if(timestepComplete){
     _nbNonConstTimesteps++;
     if(_nbNonConstTimesteps >= _maxNonConstTimesteps && _maxNonConstTimesteps > 0)
       _freezed = true;
    }

    // type specific update functionality
    _update_(timestepComplete, oldValues, res);
  }

  //@brief: returns true if a QR decomposition from scratch is necessary
  bool requireNewQR(){
    preciceTrace1("requireNewQR()", _requireNewQR);
    return _requireNewQR;
  }

  //@brief to tell the preconditioner that QR-decomposition has been recomputed
  void newQRfulfilled(){
    _requireNewQR = false;
  }

  std::vector<double>& getWeights()
  {
    return _weights;
  }

  bool isConst()
  {
    return _freezed;
  }

protected:

  //@brief weights used to scale the matrix V and the residual
  std::vector<double> _weights;

  //@brief inverse weights (for efficiency reasons)
  std::vector<double> _invWeights;

  //@brief dimension (scalar or vectorial) of each sub-vector
  std::vector<int> _dimensions;

  //@brief size of a scalar sub-vector (aka number of vertices)
  int _sizeOfSubVector;

  /** @brief maximum number of non-const time steps, i.e., after this number of time steps,
   *  the preconditioner is freezed with the current weights and becomes a constant preconditioner
   */
  int _maxNonConstTimesteps;

  /// @brief counts the number of completed time steps with a non-const weighting
  int _nbNonConstTimesteps;

  // true if a QR decomposition from scratch is necessary
  bool _requireNewQR;


  /// @brief true if _nbNonConstTimesteps >= _maxNonConstTimesteps, i.e., preconditioner is not updated any more.
  bool _freezed;


  /**
   * @brief Update the scaling after every FSI iteration and require a new QR decomposition (if necessary)
   *
   * @param timestepComplete [IN] True if this FSI iteration also completed a timestep
   */
  virtual void _update_(bool timestepComplete, const Eigen::VectorXd& oldValues, const Eigen::VectorXd& res) =0;


private:

  static logging::Logger _log;

};


}}} // namespace precice, cplscheme, impl
