#pragma once

#include "mapping/Mapping.hpp"
#include "boost/smart_ptr.hpp"
#include "tarch/la/DynamicMatrix.h"
#include "tarch/la/DynamicVector.h"
#include "tarch/la/LUDecomposition.h"
#include "tarch/la/TransposedMatrix.h"
#include "io/TXTWriter.hpp"
#include <limits>
#include <typeinfo>

#include <iostream>
using std::cout;
using std::endl;

#include "petnum.hpp"
#include "petscmat.h"
#include "petscksp.h"

namespace precice {
namespace mapping {

/**
 * @brief Mapping with radial basis functions.
 *
 * With help of the input data points and values an interpolant is constructed.
 * The interpolant is formed by a weighted sum of conditionally positive radial
 * basis functions and a (low order) polynomial, and evaluated at the output
 * data points.
 *
 * The radial basis function type has to be given as template parameter, and has
 * to be one of the defined types in this file.
 */
template<typename RADIAL_BASIS_FUNCTION_T>
class PetRadialBasisFctMapping : public Mapping
{
public:

  /**
   * @brief Constructor.
   *
   * @param constraint [IN] Specifies mapping to be consistent or conservative.
   * @param function [IN] Radial basis function used for mapping.
   */
  PetRadialBasisFctMapping (
    Constraint              constraint,
    RADIAL_BASIS_FUNCTION_T function );

  /**
   * @brief Destructor, empty.
   */
  virtual ~PetRadialBasisFctMapping() {}

  /**
   * @brief Computes the mapping coefficients from the in- and output mesh.
   */
  virtual void computeMapping();

  /**
   * @brief Returns true, if computeMapping() has been called.
   */
  virtual bool hasComputedMapping();

  /**
   * @brief Removes a computed mapping.
   */
  virtual void clear();

  /// @brief Maps input data to output data from input mesh to output mesh.
  virtual void map (
    int inputDataID,
    int outputDataID );

private:

  /// @brief Logging device.
  static tarch::logging::Log _log;

  bool _hasComputedMapping;

  /// @brief Radial basis function type used in interpolation.
  RADIAL_BASIS_FUNCTION_T _basisFunction;

  petsc::Matrix _matrixC;

  petsc::Matrix _matrixA;
};

/**
 * @brief Radial basis function with global support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 *
 * Evaluates to: radius^2 * log(radius).
 */
class PetThinPlateSplines
{
public:

  bool hasCompactSupport() const
  { return false; }

  double getSupportRadius() const
  { return std::numeric_limits<double>::max(); }

  double evaluate ( double radius ) const
  {
    double result = 0.0;
    if (tarch::la::greater(radius, 0.0)){
      result = std::log10(radius) * std::pow(radius, 2);
    }
    return result;
  }
};

/**
 * @brief Radial basis function with global support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 *
 * Evaluates to: sqrt(shape + radius^2).
 */
class PetMultiquadrics
{
public:

  PetMultiquadrics ( double c )
    : _cPow2(std::pow(c, 2)) {}

  bool hasCompactSupport() const
  { return false; }

  double getSupportRadius() const
  { return std::numeric_limits<double>::max(); }

  double evaluate ( double radius ) const
  {
    return std::sqrt(_cPow2 + std::pow(radius, 2));
  }

private:

  double _cPow2;
};

/**
 * @brief Radial basis function with global support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 * Takes a shape parameter (shape > 0.0) on construction.
 *
 * Evaluates to: 1 / (shape^2 + radius^2).
 */
class PetInverseMultiquadrics
{
public:

  PetInverseMultiquadrics ( double c )
    : _cPow2(std::pow(c, 2))
  {
    preciceCheck(tarch::la::greater(c, 0.0), "InverseMultiquadrics()",
                 "Shape parameter for radial-basis-function inverse multiquadric"
                 << " has to be larger than zero!");
  }

  bool hasCompactSupport() const
  { return false; }

  double getSupportRadius() const
  { return std::numeric_limits<double>::max(); }

  double evaluate ( double radius ) const
  {
    return 1.0 / std::sqrt(_cPow2 + std::pow(radius, 2));
  }

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  double _cPow2;
};

/**
 * @brief Radial basis function with global support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 *
 * Evaluates to: radius.
 */
class PetVolumeSplines
{
public:

  bool hasCompactSupport() const
  { return false; }

  double getSupportRadius() const
  { return std::numeric_limits<double>::max(); }

  double evaluate ( double radius ) const
  {
    return radius;
  }
};

/**
 * @brief Radial basis function with global support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 * Takes a shape parameter (shape > 0.0) on construction.
 *
 * Evaluates to: exp(-1 * shape * radius^2).
 */
class PetGaussian
{
public:

  PetGaussian ( double shape )
    : _shape(shape)
  {
    preciceCheck(tarch::la::greater(_shape, 0.0), "Gaussian()",
                 "Shape parameter for radial-basis-function gaussian"
                 << " has to be larger than zero!");
  }

  bool hasCompactSupport() const
  { return false; }

  double getSupportRadius() const
  { return std::numeric_limits<double>::max(); }

  double evaluate ( double radius ) const
  {
    return std::exp( - std::pow(_shape*radius,2.0) );
  }

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  double _shape;
};

/**
 * @brief Radial basis function with compact support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 * Takes the support radius (> 0.0) on construction.
 *
 *
 * Evaluates to: 1 - 30*rn^2 - 10*rn^3 + 45*rn^4 - 6*rn^5 - 60*log(rn^3),
 * where rn is the radius r normalized over the support radius sr: rn = r/sr.
 */
class PetCompactThinPlateSplinesC2
{
public:

  PetCompactThinPlateSplinesC2 ( double supportRadius )
    : _r(supportRadius)
  {
    preciceCheck(tarch::la::greater(_r, 0.0), "CompactThinPlateSplinesC2()",
                 "Support radius for radial-basis-function compact thin-plate-splines c2"
                 << " has to be larger than zero!");
  }

  bool hasCompactSupport() const
  { return true; }

  double getSupportRadius() const
  { return _r; }

  double evaluate ( double radius ) const
  {
    if (radius >= _r) return 0.0;
    double p = radius / _r;
    using std::pow;
    using std::log;
    return 1.0 - 30.0*pow(p,2.0) - 10.0*pow(p,3.0) + 45.0*pow(p,4.0)
      - 6.0*pow(p,5.0) - 60.0*log(pow(p,pow(p,3.0)));
  }

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  double _r;
};

/**
 * @brief Radial basis function with compact support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 * Takes the support radius (> 0.0) on construction.
 *
 *
 * Evaluates to: (1 - rn)^2,
 * where rn is the radius r normalized over the support radius sr: rn = r/sr.
 */
class PetCompactPolynomialC0
{
public:

  PetCompactPolynomialC0 ( double supportRadius )
    : _r(supportRadius)
  {
    preciceCheck(tarch::la::greater(_r, 0.0), "CompactPolynomialC0()",
                 "Support radius for radial-basis-function compact polynomial c0"
                 << " has to be larger than zero!");
  }

  bool hasCompactSupport() const
  { return true; }

  double getSupportRadius() const
  { return _r; }

  double evaluate ( double radius ) const
  {
    if (radius >= _r) return 0.0;
    return std::pow(1.0 - radius/_r, 2.0);
  }

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  double _r;
};

/**
 * @brief Radial basis function with compact support.
 *
 * To be used as template parameter for RadialBasisFctMapping.
 * Takes the support radius (> 0.0) on construction.
 *
 *
 * Evaluates to: (1 - rn)^8 * (32*rn^3 + 25*rn^2 + 8*rn + 1),
 * where rn is the radius r normalized over the support radius sr: rn = r/sr.
 */
class PetCompactPolynomialC6
{
public:

  PetCompactPolynomialC6 ( double supportRadius )
    : _r(supportRadius)
  {
    preciceCheck(tarch::la::greater(_r, 0.0), "CompactPolynomialC6()",
                 "Support radius for radial-basis-function compact polynomial c6"
                 << " has to be larger than zero!");
  }

  bool hasCompactSupport() const
  { return true; }

  double getSupportRadius() const
  { return _r; }

  double evaluate ( double radius ) const
  {
    if (radius >= _r) return 0.0;
    double p = radius / _r;
    using std::pow;
    return pow(1.0-p,8.0) * (32.0*pow(p,3.0) + 25.0*pow(p,2.0) + 8.0*p + 1.0);
  }

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  double _r;
};

// --------------------------------------------------- HEADER IMPLEMENTATIONS

template<typename RADIAL_BASIS_FUNCTION_T>
tarch::logging::Log PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::_log("precice::mapping::RadialBasisFctMapping");

template<typename RADIAL_BASIS_FUNCTION_T>
PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::PetRadialBasisFctMapping
(
  Constraint              constraint,
  RADIAL_BASIS_FUNCTION_T function )
  :
  Mapping ( constraint ),
  _hasComputedMapping ( false ),
  _basisFunction ( function ),
  _matrixC(PETSC_COMM_SELF, "C"),
  _matrixA(PETSC_COMM_SELF, "A")
{
  setInputRequirement(VERTEX);
  setOutputRequirement(VERTEX);
}

template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>::computeMapping()
{
  preciceTrace("computeMapping()");
  using namespace tarch::la;
  assertion2(input()->getDimensions() == output()->getDimensions(),
             input()->getDimensions(), output()->getDimensions());
  int dimensions = input()->getDimensions();
  mesh::PtrMesh inMesh;
  mesh::PtrMesh outMesh;
  if (getConstraint() == CONSERVATIVE){
    inMesh = output();
    outMesh = input();
  }
  else {
    inMesh = input();
    outMesh = output();
  }
  int inputSize = (int)inMesh->vertices().size();
  int outputSize = (int)outMesh->vertices().size();
  int polyparams = 1 + dimensions;
  PetscErrorCode ierr = 0;
  assertion1(inputSize >= 1 + polyparams, inputSize);
  int n = inputSize + polyparams; // Add linear polynom degrees

  _matrixC.reset(); 
  ierr = MatSetType(_matrixC.matrix, MATSBAIJ); CHKERRV(ierr); // create symmetric, block sparse matrix.
  ierr = MatSetSizes(_matrixC.matrix, PETSC_DECIDE, PETSC_DECIDE, n, n); CHKERRV(ierr);
  ierr = MatSetOption(_matrixC.matrix, MAT_SYMMETRY_ETERNAL, PETSC_TRUE); CHKERRV(ierr);
  ierr = MatSetUp(_matrixC.matrix); CHKERRV(ierr);

  _matrixA.reset();
  ierr = MatSetType(_matrixA.matrix, MATAIJ); CHKERRV(ierr); // create sparse matrix.
  ierr = MatSetSizes(_matrixA.matrix, PETSC_DECIDE, PETSC_DECIDE, outputSize, n); CHKERRV(ierr);
  ierr = MatSetUp(_matrixA.matrix); CHKERRV(ierr);

  // Fill upper right part (due to symmetry) of _matrixCLU with values
  int i = 0;
  utils::DynVector distance(dimensions);
  foreach (const mesh::Vertex& iVertex, inMesh->vertices()) {
    for (int j=iVertex.getID(); j < inputSize; j++) {
      distance = iVertex.getCoords() - inMesh->vertices()[j].getCoords();
      double coeff = _basisFunction.evaluate(norm2(distance));
      ierr = MatSetValue(_matrixC.matrix, i, j, coeff, INSERT_VALUES); CHKERRV(ierr); 
#     ifdef Asserts
      if (coeff == std::numeric_limits<double>::infinity()){
        preciceError("computeMapping()", "C matrix element has value inf. "
                     << "i = " << i << ", j = " << j
                     << ", coords i = " << iVertex.getCoords() << ", coords j = "
                     << inMesh->vertices()[j].getCoords() << ", dist = "
                     << distance << ", norm2 = " << norm2(distance) << ", rbf = "
                     << coeff
                     << ", rbf type = " << typeid(_basisFunction).name());
      }
#     endif
    }
    MatSetValue(_matrixC.matrix, i, inputSize, 1.0, INSERT_VALUES);
    for (int dim=0; dim < dimensions; dim++) {
      ierr = MatSetValue(_matrixC.matrix, i, inputSize+1+dim, iVertex.getCoords()[dim], INSERT_VALUES); CHKERRV(ierr);
    }
    i++;
  }

  // Petsc requires that all diagonal entries are set, even if set to zero.
  _matrixC.assemble(MAT_FLUSH_ASSEMBLY);
  petsc::Vector zeros(_matrixC);
  MatDiagonalSet(_matrixC.matrix, zeros.vector, ADD_VALUES);
  _matrixC.assemble(MAT_FINAL_ASSEMBLY);

  i = 0;
  foreach (const mesh::Vertex& iVertex, outMesh->vertices()){
    int j = 0;
    foreach (const mesh::Vertex& jVertex, inMesh->vertices()){
      distance = iVertex.getCoords() - jVertex.getCoords();
      double coeff = _basisFunction.evaluate(norm2(distance));
      ierr = MatSetValue(_matrixA.matrix, i, j, coeff, INSERT_VALUES); CHKERRV(ierr); 
#     ifdef Asserts
      if (coeff == std::numeric_limits<double>::infinity()){
        preciceError("computeMapping()", "A matrix element has value inf. "
                     << "i = " << i << ", j = " << j
                     << ", coords i = " << iVertex.getCoords() << ", coords j = "
                     << jVertex.getCoords() << ", dist = "
                     << distance << ", norm2 = " << norm2(distance) << ", rbf = "
                     << coeff
                     << ", rbf type = " << typeid(_basisFunction).name());
      }
#     endif
      j++;
    }
    ierr = MatSetValue(_matrixA.matrix, i, inputSize, 1.0, INSERT_VALUES); CHKERRV(ierr); 
    for (int dim=0; dim < dimensions; dim++){
      ierr = MatSetValue(_matrixA.matrix, i, inputSize+1+dim, iVertex.getCoords()[dim], INSERT_VALUES); CHKERRV(ierr); 
    }
    i++;
  }
  _matrixA.assemble();

  _hasComputedMapping = true;
}

template<typename RADIAL_BASIS_FUNCTION_T>
bool PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>:: hasComputedMapping()
{
  return _hasComputedMapping;
}

template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>:: clear()
{
  preciceTrace("clear()");
  _matrixC.reset();
  _matrixA.reset();
  _hasComputedMapping = false;
}

template<typename RADIAL_BASIS_FUNCTION_T>
void PetRadialBasisFctMapping<RADIAL_BASIS_FUNCTION_T>:: map
(
  int inputDataID,
  int outputDataID )
{
  preciceTrace2("map()", inputDataID, outputDataID);
  assertion(_hasComputedMapping);
  assertion2(input()->getDimensions() == output()->getDimensions(),
             input()->getDimensions(), output()->getDimensions());
  using namespace tarch::la;
  PetscErrorCode ierr = 0;
  utils::DynVector& inValues = input()->data(inputDataID)->values();
  utils::DynVector& outValues = output()->data(outputDataID)->values();

  int valueDim = input()->data(inputDataID)->getDimensions();
  assertion2(valueDim == output()->data(outputDataID)->getDimensions(),
             valueDim, output()->data(outputDataID)->getDimensions());
  int polyparams = 1 + input()->getDimensions();

  KSP solver;
  KSPCreate(PETSC_COMM_SELF, &solver);
  KSPSetOperators(solver, _matrixC.matrix, _matrixC.matrix);
  KSPSetFromOptions(solver);

  if (getConstraint() == CONSERVATIVE) {
    preciceDebug("Map conservative");
    static int mappingIndex = 0;
    petsc::Vector Au(_matrixC, "Au");
    petsc::Vector out(_matrixC, "out");
    petsc::Vector in(_matrixA, "in");

    for (int dim=0; dim < valueDim; dim++) {
      int size = in.getSize();
      for (int i=0; i < size; i++) { // Fill input data values
        in.setValue(i, inValues[i*valueDim + dim]);
      }
      in.assemble();

      ierr = MatMultTranspose(_matrixA.matrix, in.vector, Au.vector); CHKERRV(ierr);
      ierr = KSPSolve(solver, Au.vector, out.vector); CHKERRV(ierr);
      VecChop(out.vector, 1e-9);
      // Copy mapped data to output data values
      PetscScalar *outArray;
      ierr = VecGetArray(out.vector, &outArray);
      size = out.getSize();
      for (int i=0; i < size-polyparams; i++){
        outValues[i*valueDim + dim] = outArray[i];
      }
      VecRestoreArray(out.vector, &outArray);
    }
    mappingIndex++;
  }
  else { // Map consistent
    preciceDebug("Map consistent");
    petsc::Vector p(_matrixC, "p");
    petsc::Vector in(_matrixC, "in");
    petsc::Vector out(_matrixA, "out");
    PetscScalar *vecArray;

    // For every data dimension, perform mapping
    for (int dim=0; dim < valueDim; dim++){
      // Fill input from input data values (last polyparams entries remain zero)
      ierr = VecGetArray(in.vector, &vecArray);
      int size  = in.getSize();
      for (int i=0; i < size - polyparams; i++){
        vecArray[i] = inValues[i*valueDim + dim];
      }
      VecRestoreArray(in.vector, &vecArray);
      in.assemble();
      
      ierr = KSPSolve(solver, in.vector, p.vector); CHKERRV(ierr);
      ierr = MatMult(_matrixA.matrix, p.vector, out.vector); CHKERRV(ierr);

      VecChop(out.vector, 1e-9);
      // Copy mapped data to output data values
      ierr = VecGetArray(out.vector, &vecArray);
      size = out.getSize();
      for (int i=0; i < size; i++) {
        outValues[i*valueDim + dim] = vecArray[i];
      }
      VecRestoreArray(out.vector, &vecArray);
    }
  }
}

}} // namespace precice, mapping

