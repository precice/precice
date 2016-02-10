/*
 * SVDFactorization.cpp
 *
 *  Created on: Feb 5, 2016
 *      Author: scheufks
 */

#include "SVDFactorization.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Globals.hpp"


namespace precice {
namespace cplscheme {
namespace impl {

tarch::logging::Log SVDFactorization::
      _log("precice::cplscheme::impl::SVDFactorization");


SVDFactorization::SVDFactorization(
    double eps,
    PtrPreconditioner preconditioner)
:
  _preconditioner(preconditioner),
  _parMatrixOps(nullptr),
  _psi(),
  _phi(),
  _sigma(),
  _rows(0),
  _cols(0),
  _globalRows(0),
  _truncationEps(eps),
  _preconditionerApplied(false),
  _initialized(false),
  _initialSVD(false),
  _infostream(),
  _fstream_set(false)
{}

void SVDFactorization::initialize(
    PtrParMatrixOps parOps,
    int globalRows)
{
  _parMatrixOps = parOps;
  _globalRows = globalRows;
  _initialized = true;
}


void SVDFactorization::applyPreconditioner()
{
  preciceTrace(__func__);

  if(_psi.size() > 0 && _phi.size() > 0){
    // apply preconditioner: \psi_i * P_i, corresponds to Wtil_i * P_i, local!
    _preconditioner->apply(_psi);
    // apply preconditioner: \phi_i^T * P_i^{-1}, corresponds to Z_i * P_i^{-1}, local!
    // here, \phi^T should be preconditioned from right with inv_weights, i.e., the columns
    // of \phi^T are scaled. This is identical to scaling the rows of \phi, i.e., applying
    // P^{-1} * \phi
    _preconditioner->revert(_phi);
  }
  _preconditionerApplied = true;
}

void SVDFactorization::revertPreconditioner()
{
  preciceTrace(__func__);

  if(_psi.size() > 0 && _phi.size() > 0){
    // revert preconditioner: \psi_i * P_i^{-1}, corresponds to Wtil_i * P_i^{-1}, local!
    _preconditioner->revert(_psi);
    // revert preconditioner: \phi_i^T * P_i, corresponds to Z_i * P_i, local!
    // here, \phi^T should be preconditioned from right with _weights, i.e., the columns
    // of \phi^T are scaled. This is identical to scaling the rows of \phi, i.e., applying
    // P * \phi
    _preconditioner->apply(_phi);
  }
  _preconditionerApplied = false;
}

void SVDFactorization::reset()
{
  _psi.resize(0,0);
  _phi.resize(0,0);
  _sigma.resize(0);
  _preconditionerApplied = false;
  _initialSVD = false;
}


SVDFactorization::Matrix& SVDFactorization::matrixPhi()
{
  assertion(_preconditionerApplied);
  return _phi;
}

SVDFactorization::Matrix& SVDFactorization::matrixPsi()
{
  assertion(_preconditionerApplied);
  return _psi;
}

SVDFactorization::Vector& SVDFactorization::singularValues()
{
  assertion(_preconditionerApplied);
  return _sigma;
}

void SVDFactorization::setPrecondApplied(bool b)
{
  _preconditionerApplied = b;
}

bool SVDFactorization::isPrecondApplied()
{
  return _preconditionerApplied;
}

bool SVDFactorization::isSVDinitialized(){
  return _initialSVD;
}

void SVDFactorization::setThreshold(double eps)
{
  _truncationEps = eps;
}

int SVDFactorization::cols()
{
  return _cols;
}

int SVDFactorization::rows()
{
  return _rows;
}

int SVDFactorization::rank()
{
  return _cols;
}


}}} // namespace precice, cplscheme, impl
