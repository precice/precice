/*
 * SVDFactorization.cpp
 *
 *  Created on: Feb 5, 2016
 *      Author: scheufks
 */

#include "SVDFactorization.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Globals.hpp"
#include "QRFactorization.hpp"


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

}

void SVDFactorization::revertPreconditioner()
{
  preciceTrace(__func__);

}

SVDFactorization::reset()
{
  _psi.resize(0,0);
  _phi.resize(0,0);
  _sigma.resize(0);
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
