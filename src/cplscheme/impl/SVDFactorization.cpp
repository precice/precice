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
  _psi(),
  _phi(),
  _sigma(),
  _rows(0),
  _cols(0),
  _truncationEps(eps),
  _preconditionerApplied(false),
  _infostream(),
  _fstream_set(false)
{}

void SVDFactorization::update(
    Matrix& A,
    Matrix& B)
{
  //TODO: implement
}


SVDFactorization::reset()
{
  _psi.resize(0,0);
  _phi.resize(0,0);
  _sigma.resize(0,0);
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

SVDFactorization::Matrix& SVDFactorization::matrixSigma()
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

int SVDFactorization::cols()
{
  return _cols;
}

int SVDFactorization::rows()
{
  return _rows;
}

}}} // namespace precice, cplscheme, impl
