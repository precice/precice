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

  tarch::la::DynamicMatrix<double> _matrixCLU;

  tarch::la::DynamicVector<int> _pivotsCLU;

  tarch::la::DynamicMatrix<double> _matrixA;
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
  _matrixCLU (),
  _pivotsCLU (),
  _matrixA ()
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
  assertion1(inputSize >= 1 + polyparams, inputSize);
  int n = inputSize + polyparams; // Add linear polynom degrees
  _matrixCLU = DynamicMatrix<double>(n, n, 0.0);
  _pivotsCLU.clear();
  _pivotsCLU.append(n, 0);
  _matrixA = DynamicMatrix<double>(outputSize, n, 0.0);
  // Fill upper right part (due to symmetry) of _matrixCLU with values
  int i = 0;
  utils::DynVector distance(dimensions);
  foreach (const mesh::Vertex& iVertex, inMesh->vertices()){
    for (int j=iVertex.getID(); j < inputSize; j++){
      distance = iVertex.getCoords();
      distance -= inMesh->vertices()[j].getCoords();
      _matrixCLU(i,j) = _basisFunction.evaluate(norm2(distance));
#     ifdef Asserts
      if (_matrixCLU(i,j) == std::numeric_limits<double>::infinity()){
        preciceError("computeMapping()", "C matrix element has value inf. "
                     << "i = " << i << ", j = " << j
                     << ", coords i = " << iVertex.getCoords() << ", coords j = "
                     << inMesh->vertices()[j].getCoords() << ", dist = "
                     << distance << ", norm2 = " << norm2(distance) << ", rbf = "
                     << _basisFunction.evaluate(norm2(distance))
                     << ", rbf type = " << typeid(_basisFunction).name());
      }
#     endif
    }
    _matrixCLU(i,inputSize) = 1.0;
    for (int dim=0; dim < dimensions; dim++){
      _matrixCLU(i,inputSize+1+dim) = iVertex.getCoords()[dim];
    }
    i++;
  }
  // Copy values of upper right part of C to lower left part
  for (int i=0; i < n; i++){
    for (int j=i+1; j < n; j++){
      _matrixCLU(j,i) = _matrixCLU(i,j);
    }
  }

  // Fill _matrixA with values
  i = 0;
  foreach (const mesh::Vertex& iVertex, outMesh->vertices()){
    int j = 0;
    foreach (const mesh::Vertex& jVertex, inMesh->vertices()){
      distance = iVertex.getCoords();
      distance -= jVertex.getCoords();
      _matrixA(i,j) = _basisFunction.evaluate(norm2(distance));
#     ifdef Asserts
      if (_matrixA(i,j) == std::numeric_limits<double>::infinity()){
        preciceError("computeMapping()", "A matrix element has value inf. "
                     << "i = " << i << ", j = " << j
                     << ", coords i = " << iVertex.getCoords() << ", coords j = "
                     << jVertex.getCoords() << ", dist = "
                     << distance << ", norm2 = " << norm2(distance) << ", rbf = "
                     << _basisFunction.evaluate(norm2(distance))
                     << ", rbf type = " << typeid(_basisFunction).name());
      }
#     endif
      j++;
    }
    _matrixA(i,inputSize) = 1.0;
    for (int dim=0; dim < dimensions; dim++){
      _matrixA(i,inputSize+1+dim) = iVertex.getCoords()[dim];
    }
    i++;
  }

# ifdef PRECICE_STATISTICS
  static int computeIndex = 0;
  std::ostringstream streamC;
  if (getConstraint() == CONSERVATIVE){
    streamC << "conservative-matrixC-" << computeIndex << ".mat";
  }
  else {
    streamC << "consistent-matrixC-" << computeIndex << ".mat";
  }
  io::TXTWriter::write(_matrixCLU, streamC.str());
  std::ostringstream streamA;
  if (getConstraint() == CONSERVATIVE){
    streamA << "conservative-matrixA-" << computeIndex << ".mat";
  }
  else {
    streamA << "consistent-matrixA-" << computeIndex << ".mat";
  }
  io::TXTWriter::write(_matrixA, streamA.str());
  computeIndex++;
# endif // PRECICE_STATISTICS

//  preciceDebug ( "Matrix C = " << _matrixCLU );
  lu(_matrixCLU, _pivotsCLU);  // Compute LU decomposition

  int rankDeficiency = 0;
  for (int i=0; i < n; i++){
    if (equals(_matrixCLU(i,i), 0.0)){
      rankDeficiency++;
    }
  }
  if (rankDeficiency > 0){
    preciceWarning("computeMapping()", "Interpolation matrix C has rank "
                   << "deficiency of " << rankDeficiency);
  }
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
  _matrixCLU = tarch::la::DynamicMatrix<double>();
  _pivotsCLU.clear();
  _matrixA = tarch::la::DynamicMatrix<double>();
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

  utils::DynVector& inValues = input()->data(inputDataID)->values();
  utils::DynVector& outValues = output()->data(outputDataID)->values();
  int valueDim = input()->data(inputDataID)->getDimensions();
  assertion2(valueDim == output()->data(outputDataID)->getDimensions(),
             valueDim, output()->data(outputDataID)->getDimensions());
  int polyparams = 1 + input()->getDimensions();

  if (getConstraint() == CONSERVATIVE){
    preciceDebug("Map conservative");
    static int mappingIndex = 0;
    DynamicVector<double> Au(_matrixCLU.rows(), 0.0);
    DynamicVector<double> y(_matrixCLU.rows(), 0.0);
    DynamicVector<double> in(_matrixA.rows(), 0.0);
    DynamicVector<double> out(_matrixCLU.rows(), 0.0);

    preciceDebug("C rows=" << _matrixCLU.rows() << " cols=" << _matrixCLU.cols());
    preciceDebug("A rows=" << _matrixA.rows() << " cols=" << _matrixA.cols());
    preciceDebug("in size=" << in.size() << ", out size=" << out.size());

    for (int dim=0; dim < valueDim; dim++){
      for (int i=0; i < in.size(); i++){ // Fill input data values
        int index = i*valueDim + dim;
        in[i] = inValues[index];
      }
#     ifdef PRECICE_STATISTICS
      std::ostringstream stream;
      stream << "invec-dim" << dim << "-" << mappingIndex << ".mat";
      io::TXTWriter::write(in, stream.str());
#     endif

      multiply(transpose(_matrixA), in, Au); // Multiply by transposed of A
      // Account for pivoting in LU decomposition of C
      assertion2(Au.size() == _pivotsCLU.size(), in.size(), _pivotsCLU.size());
      for ( int i=0; i < Au.size(); i++ ){
        double temp = Au[i];
        Au[i] = Au[_pivotsCLU[i]];
        Au[_pivotsCLU[i]] = temp;
      }
      forwardSubstitution(_matrixCLU, Au, y);
      backSubstitution(_matrixCLU, y, out);
      // Copy mapped data to output data values
#     ifdef PRECICE_STATISTICS
      std::ostringstream stream2;
      stream2 << "outvec-dim" << dim << "-" << mappingIndex << ".mat";
      io::TXTWriter::write(out, stream2.str());
#     endif
      for (int i=0; i < out.size()-polyparams; i++){
        outValues[i*valueDim + dim] = out[i];
      }
    }
    mappingIndex++;
  }
  else { // Map consistent
    preciceDebug("Map consistent");
    DynamicVector<double> p(_matrixCLU.rows(), 0.0);
    DynamicVector<double> y(_matrixCLU.rows(), 0.0);
    DynamicVector<double> in(_matrixCLU.rows(), 0.0);
    DynamicVector<double> out(_matrixA.rows(), 0.0);
    // For every data dimension, perform mapping
    for (int dim=0; dim < valueDim; dim++){
      // Fill input from input data values (last polyparams entries remain zero)
      for (int i=0; i < in.size() - polyparams; i++){
        int index = i*valueDim + dim;
        in[i] = inValues[index];
      }
      // Account for pivoting in LU decomposition of C
      assertion2(in.size() == _pivotsCLU.size(), in.size(), _pivotsCLU.size());
      for (int i=0; i < in.size(); i++){
        double temp = in[i];
        in[i] = in[_pivotsCLU[i]];
        in[_pivotsCLU[i]] = temp;
      }
      forwardSubstitution(_matrixCLU, in, y);
      backSubstitution(_matrixCLU, y, p);
      multiply(_matrixA, p, out );
      // Copy mapped data to ouptut data values
      for (int i=0; i < out.size(); i++){
        outValues[i*valueDim + dim] = out[i];
      }
    }
  }
}

}} // namespace precice, mapping

