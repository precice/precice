#pragma once

#include <Eigen/Cholesky>
#include <Eigen/QR>
#include <Eigen/SVD>
#include <boost/range/adaptor/indexed.hpp>
#include <type_traits>
#include "impl/BasisFunctions.hpp"
#include "impl/RBFParameterTuner.hpp"
#include "mapping/MathHelper.hpp"
#include "mapping/RadialBasisFctSolver.hpp"
#include "mapping/config/MappingConfiguration.hpp"
#include "mapping/config/MappingConfigurationTypes.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "mesh/Mesh.hpp"
#include "profiling/Event.hpp"
#include "mapping/RadialBasisFctSolver.hpp"

namespace precice::mapping {

/**
 * This class assembles and solves an RBF system, given an input mesh and an output mesh with relevant vertex IDs.
 * The class uses a dense matrix decomposition in order to decompose the resulting system(s) and a backward substitution
 * in order to solve the system at runtime. The functionality uses Eigen and supports only serial execution. In case
 * the polynomial="separate" option is used, the polynomial system is solved using a QR decomposition.
 */
template <typename RADIAL_BASIS_FUNCTION_T>
class AutoTunedRBFSolver : public RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T> {
    mutable precice::logging::Logger _log{"mapping::AutoTunedRBFSolver"};
public:
    using BASIS_FUNCTION_T  = RADIAL_BASIS_FUNCTION_T;
    
    template <typename IndexContainer>
    Eigen::VectorXd solveConsistent(Eigen::VectorXd &inputData, const IndexContainer &inputIDs, Polynomial polynomial);

    template <typename IndexContainer>
    AutoTunedRBFSolver(RADIAL_BASIS_FUNCTION_T basisFunction, mesh::PtrMesh inputMesh, const IndexContainer &inputIDs,
                       mesh::PtrMesh outputMesh, const IndexContainer &outputIDs, std::vector<bool> deadAxis, Polynomial polynomial, MappingConfiguration::AutotuningParams rbfTunerConfig);

    AutoTunedRBFSolver() = default;
};


template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexContainer>
AutoTunedRBFSolver<RADIAL_BASIS_FUNCTION_T>::AutoTunedRBFSolver(RADIAL_BASIS_FUNCTION_T basisFunction, mesh::PtrMesh inputMesh, const IndexContainer &inputIDs,
                                                                    mesh::PtrMesh outputMesh, const IndexContainer &outputIDs, std::vector<bool> deadAxis, 
                                                                    Polynomial polynomial, MappingConfiguration::AutotuningParams rbfTunerConfig)
    : RadialBasisFctSolver<RADIAL_BASIS_FUNCTION_T>(basisFunction, inputMesh, inputIDs, outputMesh, outputIDs, deadAxis, polynomial, rbfTunerConfig) {}


// @todo: change the signature to Eigen::MatrixXd and process all components at once, the solve function of eigen can handle that
template <typename RADIAL_BASIS_FUNCTION_T>
template <typename IndexContainer>
Eigen::VectorXd AutoTunedRBFSolver<RADIAL_BASIS_FUNCTION_T>::solveConsistent(Eigen::VectorXd &inputData, const IndexContainer &inputIDs, Polynomial polynomial)
{
  PRECICE_ASSERT((this->_matrixQ.size() > 0 && polynomial == Polynomial::SEPARATE) || this->_matrixQ.size() == 0);
  Eigen::VectorXd polynomialContribution;
  // Solve polynomial QR and subtract it from the input data
  if (polynomial == Polynomial::SEPARATE) {
    polynomialContribution = this->_qrMatrixQ.solve(inputData);
    inputData -= (this->_matrixQ * polynomialContribution);
  }

  if (this->_rbfTunerConfig.autotuneShape) {
    if constexpr (RADIAL_BASIS_FUNCTION_T::hasCompactSupport()) {

      const bool isProbablyZeroData = this->_tuner->getLastOptimizationError() < 1e-15;
      // @todo: Data dimension not considered
      const bool isOptimizationInterval  = this->_rbfTunerConfig.iterationInterval == this->_iterSinceOptimization;
      const bool isOptimizedAndShouldOnce = this->_rbfTunerConfig.optimizeOnce() && this->_tuner->getLastOptimizationError() != std::numeric_limits<double>::max();

      // @todo: Add error criterion.
      if (isProbablyZeroData || isOptimizationInterval || !isOptimizedAndShouldOnce) {
        Sample sample        = this->_tuner->optimize(*this, this->_inputMesh.lock(), inputIDs, inputData);
        if (!this->_tuner->lastSampleWasOptimum()) {
          this->rebuildKernelDecomposition(this->_inputMesh.lock(), inputIDs, this->_tuner->getLastOptimizedRadius());
        }
        PRECICE_INFO("Using: radius={}, LOOCV={}", this->_tuner->getLastOptimizedRadius(), sample.error);
        this->_iterSinceOptimization = 1;
      } else {
        PRECICE_INFO("Using radius={} again", this->_tuner->getLastOptimizedRadius());
        this->_iterSinceOptimization++;
      }
    } else {
      // Should not happen, since an error will already be thrown in the configuration.
      PRECICE_ERROR("Parameter optimization is not available for RBFs without a \"support-radius\" configuration.");
    }
  }

  // Integrated polynomial (and separated)
  PRECICE_ASSERT(inputData.size() == this->_matrixA.cols());
  Eigen::VectorXd p = this->_decMatrixC->solve(inputData);

  // if (polynomial != Polynomial::ON && computeCrossValidation) {
  //   precice::profiling::Event e("map.rbf.evaluateLOOCV");
  //   PRECICE_INFO("Cross validation error (LOOCV): {}", evaluateRippaLOOCVerror(p));
  // }
  PRECICE_ASSERT(p.size() == this->_matrixA.cols());

  Eigen::VectorXd out;
  if (this->_rbfTunerConfig.autotuneShape) {
    if constexpr (RADIAL_BASIS_FUNCTION_T::hasCompactSupport()) {
      const Eigen::Index outSize    = this->_matrixA.rows();
      const Eigen::Index inSize     = inputData.size();
      const Eigen::Index polyParams = this->_matrixA.cols() - inSize;

      // not yet necessary: integrated polynomial is only supported for SPD functions
      out = applyKernelToMatrix<RADIAL_BASIS_FUNCTION_T>(this->_matrixA.block(0, 0, outSize, inSize), this->_tuner->getLastOptimizedRadius()) * p.segment(0, inSize);
      out += this->_matrixA.block(0, inSize, outSize, polyParams) * p.segment(inSize, polyParams);
    } else {
      // Should not happen, since an error will already be thrown in the configuration.
      PRECICE_ERROR("Parameter optimization is not available for RBFs without a \"support-radius\" configuration.");
    }
  } else {
    out = this->_matrixA * p;
  }

  // Add the polynomial part again for separated polynomial
  if (polynomial == Polynomial::SEPARATE) {
    out += (this->_matrixV * polynomialContribution);
  }
  return out;
}

};