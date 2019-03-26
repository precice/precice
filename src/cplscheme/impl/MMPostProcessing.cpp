#include "MMPostProcessing.hpp"
#include <Eigen/Dense>
#include "QRFactorization.hpp"
#include "com/Communication.hpp"
#include "cplscheme/CouplingData.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "utils/EigenHelperFunctions.hpp"
#include "utils/Helpers.hpp"
#include "utils/MasterSlave.hpp"
//#include "utils/NumericalCompare.hpp"

namespace precice
{
namespace cplscheme
{
namespace impl
{

/* ----------------------------------------------------------------------------
 *     Constructor
 * ----------------------------------------------------------------------------
 */
MMPostProcessing::MMPostProcessing(
    impl::PtrPostProcessing coarseModelOptimization,
    int                     maxIterationsUsed,
    int                     timestepsReused,
    int                     filter,
    double                  singularityLimit,
    bool                    estimateJacobian,
    std::vector<int>        fineDataIDs,
    std::vector<int>        coarseDataIDs,
    //std::map<int, double> scalings,
    PtrPreconditioner preconditioner)
    : PostProcessing(),
      _coarseModelOptimization(coarseModelOptimization),
      _preconditioner(preconditioner),
      _maxIterationsUsed(maxIterationsUsed),
      _timestepsReused(timestepsReused),
      _singularityLimit(singularityLimit),
      _designSpecification(),
      _coarseModel_designSpecification(),
      _fineDataIDs(fineDataIDs),
      _coarseDataIDs(coarseDataIDs),
      _estimateJacobian(estimateJacobian),
      _maxIterCoarseModelOpt(maxIterationsUsed),
      _filter(filter)
{
  CHECK(_maxIterationsUsed > 0,
        "Maximal iterations used for MM post-processing has to be larger than zero!");
  CHECK(_maxIterCoarseModelOpt > 0,
        "Maximal iterations used for coarse model optimization for MM post-processing has to "
            << "be larger than zero!");
  CHECK(_timestepsReused >= 0,
        "Number of old timesteps to be reused for MM post-processing has to be >= 0!");
}

/** ---------------------------------------------------------------------------------------------
 *         initialize()
 *
 * @brief: Initializes all the needed variables and data
 *  ---------------------------------------------------------------------------------------------
 */
void MMPostProcessing::initialize(DataMap &cplData)
{
  TRACE(cplData.size());
  size_t              entries       = 0;
  size_t              coarseEntries = 0;
  std::vector<size_t> subVectorSizes; //needed for preconditioner

  assertion(_fineDataIDs.size() == _coarseDataIDs.size(), _fineDataIDs.size(), _coarseDataIDs.size());
  assertion(_dataIDs.empty(), _dataIDs.size());

  _dataIDs.insert(_dataIDs.end(), _fineDataIDs.begin(), _fineDataIDs.end());
  _dataIDs.insert(_dataIDs.end(), _coarseDataIDs.begin(), _coarseDataIDs.end());

  for (auto &elem : _fineDataIDs) {
    CHECK(utils::contained(elem, cplData),
          "Data with ID " << elem << " is not contained in data given at initialization!");
    entries += cplData[elem]->values->size();
    subVectorSizes.push_back(cplData[elem]->values->size());
  }

  for (auto &elem : _coarseDataIDs) {
    CHECK(utils::contained(elem, cplData),
          "Data with ID " << elem << " is not contained in data given at initialization!");
    coarseEntries += cplData[elem]->values->size();
  }

  /**
   * initialize the coarse model optimization method
   */
  // view on coarse coupling data only (otherwise problem with secondary data: mixed up with fine data)
  DataMap coarseCplData;
  for (int id : _coarseDataIDs) {
    DataMap::value_type pair = std::make_pair(id, cplData[id]);
    coarseCplData.insert(pair);
  }
  _coarseModelOptimization->initialize(coarseCplData);

  // the coarse model also uses the fine mesh (only evaluation in solver is on coarse model)
  assertion(entries == coarseEntries, entries, coarseEntries);

  _matrixCols.push_front(0);
  _firstIteration = true;
  _firstTimeStep  = true;

  assertion(_coarseOldResiduals.size() == 0);
  assertion(_fineOldResiduals.size() == 0);
  _coarseOldResiduals = Eigen::VectorXd::Zero(entries);
  _fineOldResiduals   = Eigen::VectorXd::Zero(entries);
  _fineResiduals      = Eigen::VectorXd::Zero(entries);
  _coarseResiduals    = Eigen::VectorXd::Zero(entries);
  _outputFineModel    = Eigen::VectorXd::Zero(entries);
  _outputCoarseModel  = Eigen::VectorXd::Zero(entries);
  _input_Xstar        = Eigen::VectorXd::Zero(entries);

  // if design specifiaction not initialized yet
  if (not(_designSpecification.size() > 0)) {
    _designSpecification = Eigen::VectorXd::Zero(entries);
  }
  _coarseModel_designSpecification = Eigen::VectorXd::Zero(entries);

  /**
   *  make dimensions public to all procs,
   *  last entry _dimOffsets[MasterSlave::_size] holds the global dimension, global,n
   */
  if (utils::MasterSlave::_masterMode || utils::MasterSlave::_slaveMode) {
    assertion(utils::MasterSlave::_communication.get() != NULL);
    assertion(utils::MasterSlave::_communication->isConnected());

    /** provide vertex offset information for all processors
     *  mesh->getVertexOffsets() provides an array that stores the number of mesh vertices on each processor
     *  This information needs to be gathered for all meshes. To get the number of respective unknowns of a specific processor
     *  we need to multiply the number of vertices with the dimensionality of the vector-valued data for each coupling data.
     */
    _dimOffsets.resize(utils::MasterSlave::_size + 1);
    _dimOffsets[0] = 0;
    for (size_t i = 0; i < _dimOffsets.size() - 1; i++) {
      int accumulatedNumberOfUnknowns = 0;
      for (auto &elem : _fineDataIDs) {
        auto &offsets = cplData[elem]->mesh->getVertexOffsets();
        accumulatedNumberOfUnknowns += offsets[i] * cplData[elem]->dimension;
      }
      _dimOffsets[i + 1] = accumulatedNumberOfUnknowns;
    }

    // test that the computed number of unknown per proc equals the number of entries actually present on that proc
    size_t unknowns = _dimOffsets[utils::MasterSlave::_rank + 1] - _dimOffsets[utils::MasterSlave::_rank];
    assertion(entries == unknowns, entries, unknowns);
  }

  if (_estimateJacobian) {
    _MMMappingMatrix = Eigen::MatrixXd::Zero(getLSSystemRows(), entries);
    // do not initialize Tkprev (MMMappingMatrix_prev), as we need a different constructing rule
    // in the first step (see computeCoarseModelDesignSpecifiaction).

    //_MMMappingMatrix_prev = Eigen::MatrixXd::Identity(getLSSystemRows(), entries);
  }

  /**
   // Fetch secondary data IDs, to be relaxed with same coefficients from IQN-ILS
   for (DataMap::value_type& pair : cplData){
   if (not utils::contained(pair.first, _dataIDs)) {
   _secondaryDataIDs.push_back(pair.first);
   int secondaryEntries = pair.second->values->size();
   //      _secondaryOldXTildes[pair.first] Eigen::VectorXd::Zero(secondaryEntries);
   _secondaryResiduals[pair.first] = Eigen::VectorXd::Zero(secondaryEntries);
   }
   }
   */

  // Append old value columns, if not done outside of post-processing already
  for (DataMap::value_type &pair : cplData) {
    int cols = pair.second->oldValues.cols();
    if (cols < 1) { // Add only, if not already done
      //assertion(pair.second->values->size() > 0, pair.first);
      utils::append(pair.second->oldValues, (Eigen::VectorXd) Eigen::VectorXd::Zero(pair.second->values->size()));
    }
  }

  _preconditioner->initialize(subVectorSizes);
}

/** ---------------------------------------------------------------------------------------------
 *         registerSolutionCoarseModelOptimization()
 *
 * @brief: registers the current solution from the coarse model optimization problem as input
 *         for the fine model evaluation step.
 *  ---------------------------------------------------------------------------------------------
 */
void MMPostProcessing::registerSolutionCoarseModelOptimization(
    DataMap &cplData)
{
  TRACE();
  // extract new solution x_star from coarse model optimization problem from coarse cplData
  int off = 0;
  for (int id : _coarseDataIDs) {
    int   size       = cplData[id]->values->size();
    auto &valuesPart = *(cplData[id]->values);
    for (int i = 0; i < size; i++) {
      // the coarse model optimization reverts its own scaling, hence valuesPart is not scaled, can be copied.
      _input_Xstar[i + off] = valuesPart[i];
    }
    off += size;
  }

  // register new solution x_star from coarse model optimization problem as input
  // to fine model evaluation.
  off = 0;
  for (int id : _fineDataIDs) {
    int   size       = cplData[id]->values->size();
    auto &valuesPart = *(cplData[id]->values);
    for (int i = 0; i < size; i++) {
      // write new coarse model solution back as input data for the fine model evaluation
      // _input_xStar needs to be updated in each iteration
      // we want valuesPart not to be scaled, hence pure copying
      valuesPart[i] = _input_Xstar[i + off];
    }
    off += size;
  }
  // in the context of mannifold mapping post processing we want _input_Xstar to be scaled according to the scaling factors
  //  scale(_input_Xstar, cplData);
}

/** ---------------------------------------------------------------------------------------------
 *         setDesignSpecification()
 *
 * @brief: sets a design specification for the fine model optimization problem
 *         i.e., x_star = argmin_x || f(x) - q ||
 *  ---------------------------------------------------------------------------------------------
 */
void MMPostProcessing::setDesignSpecification(
    Eigen::VectorXd &q)
{
  TRACE();
  assertion(q.size() == _fineResiduals.size(), q.size(), _fineResiduals.size());
  _designSpecification = (q.norm() <= 1.0e-15) ? Eigen::VectorXd::Zero(_fineResiduals.size()) : q;

  // only in the first step, the coarse model design specification equals the design specification
  // for the overall objective function (initial coarse solution)
  if (_firstTimeStep)
    _coarseModel_designSpecification = _designSpecification;
}

/** ---------------------------------------------------------------------------------------------
 *         getDesignSpecification()
 *
 * @brief: Returns the design specification corresponding to the given coupling data that is updated in every
 *         manifold mapping cycle. This information is needed for convergence measurements in the
 *         coupling scheme.
 *  ---------------------------------------------------------------------------------------------
 */
std::map<int, Eigen::VectorXd> MMPostProcessing::getDesignSpecification(
    DataMap &cplData)
{
  std::map<int, Eigen::VectorXd> designSpecifications;
  int                            off = 0;
  for (int id : _fineDataIDs) {
    int             size = cplData[id]->values->size();
    Eigen::VectorXd q    = Eigen::VectorXd::Zero(size);
    for (int i = 0; i < size; i++) {
      q(i) = _designSpecification(i + off);
    }
    off += size;
    std::map<int, Eigen::VectorXd>::value_type pair = std::make_pair(id, q);
    designSpecifications.insert(pair);
  }
  off = 0;
  for (int id : _coarseDataIDs) {
    int             size = cplData[id]->values->size();
    Eigen::VectorXd q    = Eigen::VectorXd::Zero(size);
    for (int i = 0; i < size; i++) {
      q(i) = _coarseModel_designSpecification(i + off);
    }
    off += size;
    std::map<int, Eigen::VectorXd>::value_type pair = std::make_pair(id, q);
    designSpecifications.insert(pair);
  }
  return designSpecifications;
}

/** ---------------------------------------------------------------------------------------------
 *         updateDifferenceMatrices()
 *
 * @brief: computes the current coarse and fine model residual, computes the differences and
 *         updates the difference matrices F and C. Also stores the residuals
 *  ---------------------------------------------------------------------------------------------
 */
void MMPostProcessing::updateDifferenceMatrices(
    DataMap &cplData)
{
  TRACE();

  /**
   * Compute current residual: vertex-data - oldData
   */
  _fineResiduals   = _outputFineModel - _input_Xstar;
  _coarseResiduals = _outputCoarseModel - _input_Xstar;

  /**
   * Update matrices C, F with newest information
   */
  if (not _firstIteration) {
    DEBUG("   Update Difference Matrices C and F with coarse and fine model responses");
    assertion(_matrixF.cols() == _matrixC.cols(), _matrixF.cols(), _matrixC.cols());
    assertion(getLSSystemCols() <= _maxIterationsUsed, getLSSystemCols(), _maxIterationsUsed);

    if (2 * getLSSystemCols() >= getLSSystemRows())
      WARN(
          "The number of columns in the least squares system exceeded half the number of unknowns at the interface. The system will probably become bad or ill-conditioned and the quasi-Newton post processing may not converge. Maybe the number of allowed columns (maxIterationsUsed) should be limited.");

    Eigen::VectorXd colF = _fineResiduals - _fineOldResiduals;
    Eigen::VectorXd colC = _coarseResiduals - _coarseOldResiduals;

    bool columnLimitReached = getLSSystemCols() == _maxIterationsUsed;
    bool overdetermined     = getLSSystemCols() <= getLSSystemRows();
    if (not columnLimitReached && overdetermined) {

      utils::appendFront(_matrixF, colF);
      utils::appendFront(_matrixC, colC);

      _matrixCols.front()++;
    } else {
      utils::shiftSetFirst(_matrixF, colF);
      utils::shiftSetFirst(_matrixC, colC);

      _matrixCols.front()++;
      _matrixCols.back()--;
      if (_matrixCols.back() == 0) {
        _matrixCols.pop_back();
      }
    }
  }

  /**
   *  Store residuals
   */
  _fineOldResiduals   = _fineResiduals;
  _coarseOldResiduals = _coarseResiduals;
}

/** ---------------------------------------------------------------------------------------------
 *         performPostProcessing()
 *
 * @brief: performs one iteration of the manifold mapping post processing. It steers the execution
 *         of fine and coarse model evaluations and also calls the coarse model optimization.
 *  ---------------------------------------------------------------------------------------------
 */
void MMPostProcessing::performPostProcessing(
    DataMap &cplData)
{
  TRACE(_dataIDs.size(), cplData.size());

  //assertion(_fineDataIDs.size() == _scalings.size(), _fineDataIDs.size(), _scalings.size());
  assertion(_fineOldResiduals.size() == _fineResiduals.size(), _fineOldResiduals.size(), _fineResiduals.size());
  assertion(_coarseResiduals.size() == _fineResiduals.size(), _coarseResiduals.size(), _fineResiduals.size());
  assertion(_coarseOldResiduals.size() == _fineResiduals.size(), _coarseOldResiduals.size(), _fineResiduals.size());
  assertion(_outputFineModel.size() == _fineResiduals.size(), _outputFineModel.size(), _fineResiduals.size());
  assertion(_input_Xstar.size() == _fineResiduals.size(), _input_Xstar.size(), _fineResiduals.size());

  /**
   * Manifold Mapping cycle, computing the new design specification for the coarse model
   * using input and output datafrom the coarse and the fine model of the previous iterates (time steps)
   * Also updating mapping matrix (if jacobian_estimation = enabled)
   */
  if (not(*_isCoarseModelOptimizationActive)) {

    /**
     * assume the coarse model and the fine model has been evaluated for the new coarse model
     * solution _input_Xstar, obtained in the coarse model optimization step.
     */

    // view on coarse coupling data only
    DataMap coarseCplData;
    for (int id : _coarseDataIDs) {
      DataMap::value_type pair = std::make_pair(id, cplData[id]);
      coarseCplData.insert(pair);
    }
    // assume that we always start with a coarse model optimization step at the very beginning.
    // every time we get here, the coarse model optimization has just converged.
    _coarseModelOptimization->iterationsConverged(coarseCplData);
    _iterCoarseModelOpt = 0;

    // update the difference matrices with the newest residual deltas
    concatenateCouplingData(cplData);
    updateDifferenceMatrices(cplData);

    /**
     *  === update and apply preconditioner ===
     *
     * IQN-ILS would also work without W and xUpdate scaling, IQN-IMVJ unfortunately not
     * Note: here, the _residuals are H(x)- x - q, i.e., residual of the fixed-point iteration
     *       minus the design specification of the optimization problem (!= null if MM is used)
     */
    Eigen::VectorXd objective = _fineResiduals - _designSpecification;
    _preconditioner->update(false, _outputFineModel, objective);
    /// @todo evaluate whether the pure residual should be used for updating the preconditioner or residual - design specification
    if (getLSSystemCols() > 0) {
      _preconditioner->apply(_matrixF);
      _preconditioner->apply(_matrixC);
    }
    _preconditioner->apply(_fineResiduals);
    _preconditioner->apply(_coarseResiduals);
    if (isSet(_designSpecification))
      _preconditioner->apply(_designSpecification);

    if (_estimateJacobian && _MMMappingMatrix_prev.rows() > 0) {
      _preconditioner->apply(_MMMappingMatrix_prev, false);
      _preconditioner->revert(_MMMappingMatrix_prev, true);
    }

    /** compute the new design specification for the coarse model optimization
     *  updates: _coarseModel_designSpecification
     *           i.e., qk = c(xk) - Tk * ( f(xk) - q )
     *  updates: _MMMappingMatrix, i.e.,
     *           Tk = Tkprev + (C - Tkprev * F) * pseudoInv_F,        iff Tkprev exists
     *           Tk = C * pseudoInv_F + (I - Uc*Uc^T)*(I - Uf*Uf^T),  else (in the first step or after rescaling)
     */
    computeCoarseModelDesignSpecifiaction();

    assertion(isSet(_coarseModel_designSpecification)); // the coarse model design specification is computed within the MM cycle and should therefore be set and valid

    // undo preconditioning
    if (_estimateJacobian && _MMMappingMatrix_prev.rows() > 0) {
      _preconditioner->revert(_MMMappingMatrix_prev, false);
      _preconditioner->apply(_MMMappingMatrix_prev, true);
    }

    if (getLSSystemCols() > 0) {
      _preconditioner->revert(_matrixF);
      _preconditioner->revert(_matrixC);
    }
    _preconditioner->revert(_fineResiduals);
    _preconditioner->revert(_coarseResiduals);
    _preconditioner->revert(_designSpecification);
    // The coarse model design specification is computed with scaled data and needs to be re-scaled to normal.
    // It is to be scaled again in the coarse model optimization scheme.
    _preconditioner->revert(_coarseModel_designSpecification);
    //unscale(_coarseModel_designSpecification, cplData);

    /**
     * now, the difference matrices for the MM mapping as well as the design specification for the coarse
     * model optimization problem are updated (also Jacobian of MM mapping matrix if required).
     * next step: coarse model optimization, set the steering variable accordingly
     */
    (*_isCoarseModelOptimizationActive) = true;

    /** Undo of cplData scaling is not necessary, as we only read information from the cpl data.
     * The write back step is done in registerSolutionCoarseModelOptimization
     */

    // one MM iteration completed
    its++;
    _firstIteration = false;
  }

  /**
    * coarse model optimization cycle for the problem x_star = argmin_x|| c(x) - q_k ||
    */
  if (*_isCoarseModelOptimizationActive) {
    // view on coarse coupling data only
    DataMap coarseCplData;
    for (int id : _coarseDataIDs) {
      DataMap::value_type pair = std::make_pair(id, cplData[id]);
      coarseCplData.insert(pair);
    }

    // no preconditioning, i.e. scaling of input/output data. This is done in the coarse optimization routine, only for the coarse cplData.
    // the _coarseModel_designSpecification is scaled back at this point

    /** perform the coarse model optimization, determine x_star
     *          x_k+1 = argmin_x || c(x) - q_k ||
     *        ------------------------------------
     *  In the first time step, the coarse model design specification is equal to the
     *  design specification of the overall objective function. Here, a initial coarse
     *  model solution is obtained for a initial guess.
     */
    _coarseModelOptimization->optimize(coarseCplData, _coarseModel_designSpecification);

    /**
     * after the coarse model optimization has converged successfully,
     * the fine and the coarse model has to be evaluated for the new solution x_star.
     * Hence, x_star needs to be copied to the fine model input values.
     */
    registerSolutionCoarseModelOptimization(cplData);

    _iterCoarseModelOpt++;
    // if coarse model optimization exceeds max iteration count, print warning and break coarse model optimization iteration
    if (_iterCoarseModelOpt >= _maxIterCoarseModelOpt) {
      //(*_isCoarseModelOptimizationActive)  = false;
      _notConvergedWithinMaxIter = true;
      WARN("The coarse model optimization in coupling iteration " << its
                                                                  << " exceeds maximal number of optimization cycles (" << _maxIterCoarseModelOpt << " without convergence!");
    }
  }

  if (_notConvergedWithinMaxIter) {
    if (std::isnan(utils::MasterSlave::l2norm(_input_Xstar))) {
      ERROR("The coupling iteration in time step " << tSteps << " failed to converge and NaN values occured throughout the coupling process. "
                                                   << "This is most likely due to the fact that the coarse model failed to converge within "
                                                   << "the given maximum number of allowed iterations: " << _maxIterCoarseModelOpt);
    }
  }

  DEBUG("  * Manifold Mapping Iterations: " << its);
  DEBUG("  * Coarse Model Optimization Iterations: " << _iterCoarseModelOpt);
}

/** ---------------------------------------------------------------------------------------------
 *         computeCoarseModelDesignSpecifiaction()
 *
 * @brief: computes/updates the design specification for the coarse model optimization problem
 *     	   i. e., q_k = c(x_k) - T_k * (f(x_k) - q), q = 0 is the fine model design specification
 *  ---------------------------------------------------------------------------------------------
 */
void MMPostProcessing::computeCoarseModelDesignSpecifiaction()
{
  TRACE();

  /** update design specification
   *  alpha = (f(x) - q),
   *  q_k   = c(x)
   */
  Eigen::VectorXd alpha            = _fineResiduals - _designSpecification;
  _coarseModel_designSpecification = _coarseResiduals;

  // if residual differences are available for fine and coarse model
  // (either from previous iterations or from previous time steps or both)
  if (getLSSystemCols() > 0) {
    // compute SVDs of _matrixF and _matriC
    Eigen::VectorXd S_F, S_C;
    Eigen::MatrixXd V_F, U_F, V_C, U_C, Sigma_F, pseudoSigma_F;

    // Remove dependent columns of _matrixC and _matrixF
    int nbRemoveCols = 1;
    while (nbRemoveCols > 0) {
      nbRemoveCols = 0;
      if (getLSSystemCols() == 0)
        break;

      // Calculate singular value decomposition with Eigen
      Eigen::JacobiSVD<Eigen::MatrixXd> svd(_matrixF, Eigen::ComputeThinU | Eigen::ComputeThinV);
      Eigen::VectorXd                   singularValues = svd.singularValues();

      for (int i = 0; i < singularValues.rows(); i++) {
        if (std::abs(singularValues(i)) <= _singularityLimit) {
          std::cout << "singular value: " << singularValues(i) << '\n';

          // Remove the column from _matrixC and _matrixF
          removeMatrixColumn(i - nbRemoveCols);
          nbRemoveCols++;
        }
      }
      if (nbRemoveCols)
        DEBUG("Manifold mapping: remove " << nbRemoveCols << " columns from the Jacobian matrices");
    }

    assert(_matrixF.cols() == _matrixC.cols());

    if (getLSSystemCols() > 0) {
      // Calculate singular value decomposition with Eigen
      Eigen::JacobiSVD<Eigen::MatrixXd> svd_C(_matrixC, Eigen::ComputeThinU | Eigen::ComputeThinV);
      Eigen::JacobiSVD<Eigen::MatrixXd> svd_F(_matrixF, Eigen::ComputeThinU | Eigen::ComputeThinV);

      Eigen::MatrixXd pseudoSigma_F = svd_F.singularValues().asDiagonal();

      for (int i = 0; i < pseudoSigma_F.cols(); i++)
        pseudoSigma_F(i, i) = 1.0 / pseudoSigma_F(i, i);

      Eigen::MatrixXd pseudoMatrixF = svd_F.matrixV() * pseudoSigma_F * svd_F.matrixU().transpose();

      U_F = svd_F.matrixU();
      U_C = svd_C.matrixU();

      // if Jacobian matrix of MM mapping matrix is not set up explicitly, perform
      // matrix-vector update
      if (not _estimateJacobian) {
        Eigen::VectorXd beta = U_F * (U_F.transpose() * alpha);

        _coarseModel_designSpecification -= alpha;
        _coarseModel_designSpecification -= _matrixC * (pseudoMatrixF * alpha);
        _coarseModel_designSpecification += U_C * (U_C.transpose() * (alpha - beta));
        _coarseModel_designSpecification += beta;
      }

      // Jacobian matrix of MM mapping matrix is estimated and set up explicitly
      // multi-vector method for update of Jacobian with implicit incorporation of
      // information from previous time steps.
      if (_estimateJacobian) {
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(_matrixF.rows(), _matrixF.rows());

        // if previous Jacobian exists, i.e., no re-scaling and not first estimation
        if (_MMMappingMatrix_prev.rows() == getLSSystemRows()) {
          _MMMappingMatrix = _MMMappingMatrix_prev + (_matrixC - _MMMappingMatrix_prev * _matrixF) * pseudoMatrixF;

          // if no previous Jacobian exists, set up Jacobian with IQN-ILS update rule + stabilization term
        } else {
          _MMMappingMatrix = _matrixC * pseudoMatrixF + (I - U_C * U_C.transpose()) * (I - U_F * U_F.transpose());
        }

        // compute new design specification for coarse model optimization: qk = c(x) - Tk( f(x) - q )
        _coarseModel_designSpecification -= _MMMappingMatrix * alpha;
      }
    }
  }

  // if no residual differences for the fine and coarse model are given so far
  if ((_firstIteration && _firstTimeStep) || getLSSystemCols() <= 0) {
    assertion(getLSSystemCols() <= 0, getLSSystemCols());
    if (_estimateJacobian && (_MMMappingMatrix_prev.rows() == getLSSystemRows())) {
      _coarseModel_designSpecification -= _MMMappingMatrix_prev * alpha;
    } else {
      _coarseModel_designSpecification -= alpha;
    }
  }
}

void MMPostProcessing::concatenateCouplingData(
    DataMap &cplData)
{
  TRACE();

  int offset = 0;
  int k      = 0;
  assertion(_fineDataIDs.size() == _coarseDataIDs.size(), _fineDataIDs.size(), _coarseDataIDs.size());
  for (int id : _fineDataIDs) {
    int         size            = cplData[id]->values->size();
    auto &      values          = *cplData[id]->values;
    auto &      coarseValues    = *cplData[_coarseDataIDs.at(k)]->values;
    const auto &coarseOldValues = cplData[_coarseDataIDs.at(k)]->oldValues.col(0);
    assertion(values.size() == coarseValues.size(), values.size(), coarseValues.size());
    assertion(values.size() == coarseOldValues.size(), values.size(), coarseOldValues.size());
    for (int i = 0; i < size; i++) {
      _outputFineModel[i + offset] = values[i];
      // ignore input from fine model as it must be exactly the
      // same as the input for the coarse model, if the fine model is evaluated
      _outputCoarseModel[i + offset] = coarseValues[i];
      _input_Xstar[i + offset]       = coarseOldValues[i];
    }
    offset += size;
    k++;
  }
}

/** -----------------------------------------------------------------------------------
 *            isSet()
 *   @brief indicates whether the design specification has been set and is valid
 *  -----------------------------------------------------------------------------------
 */
bool MMPostProcessing::isSet(Eigen::VectorXd &designSpec)
{
  // design specification is considered to be set and active if
  // 1. its size is larger then zero (i. e., it must be equal to the number of unknowns)
  // 2. its l2-norm is larger then 1.0e-15
  bool set((designSpec.size() > 0)); // && (designSpec.norm() > 1.0e-15));
  if (set)
    assertion(designSpec.size() == _fineResiduals.size(), designSpec.size(), _fineResiduals.size());
  return set;
}

/** ---------------------------------------------------------------------------------------------
 *         iterationsConverged()
 *
 * @brief: Is called when the convergence criterion for the coupling is fullfilled and finalizes
 *         the manifold mapping post processing. Stores new differences in F and C, clears or
 *         updates F and C according to the number of reused time steps
 *  ---------------------------------------------------------------------------------------------
 */
void MMPostProcessing::iterationsConverged(
    DataMap &cplData)
{
  TRACE();

  its = 0;
  tSteps++;
  deletedColumns = 0;

  /** in each advance() cycle in the coupling scheme, either MMPostProcessing::performPostProcessing() or
   *  MMPostProcessing::iterationsConverged() is called. Both, if and only if the coarse model has converged.
   *  Hence, we need to call iterationsConverged() for the coarse model.
   */
  _coarseModelOptimization->iterationsConverged(cplData);
  _iterCoarseModelOpt = 0;

  // the most recent differences for the F, C matrices have not been added so far
  // this has to be done in iterations converged, as PP won't be called any more if
  // convergence was achieved
  concatenateCouplingData(cplData);
  updateDifferenceMatrices(cplData);

  /**
   * Difference matrices and Jacobian updated, MM cycle completed, start with coarse model
   * optimization in next cycle.
   * next step: coarse model optimization, set the steering variable accordingly
   */
  (*_isCoarseModelOptimizationActive) = true;

  // reset the coarse model design specification
  _coarseModel_designSpecification = _designSpecification;

  // update the preconditioner
  Eigen::VectorXd objective = _fineResiduals - _designSpecification;
  _preconditioner->update(false, _outputFineModel, objective);

  // if the multi-vector generalized broyden like update for the manifold matrix estimation process is used
  // store the estimated matrix from the last time step.
  if (_estimateJacobian && _MMMappingMatrix.cols() > 0) {
    _MMMappingMatrix_prev = _MMMappingMatrix;
    _preconditioner->revert(_MMMappingMatrix_prev, false);
    _preconditioner->apply(_MMMappingMatrix_prev, true);
  }

  _firstTimeStep = false;
  if (_matrixCols.front() == 0) { // Did only one iteration
    _matrixCols.pop_front();
  }

#ifndef NDEBUG
  std::ostringstream stream;
  stream << "Matrix column counters: ";
  for (int cols : _matrixCols) {
    stream << cols << ", ";
  }
  DEBUG(stream.str());
#endif // Debug

  if (_timestepsReused == 0) {
    _matrixF.resize(0, 0);
    _matrixC.resize(0, 0);
    _matrixCols.clear();
  } else if ((int) _matrixCols.size() > _timestepsReused) {
    int toRemove = _matrixCols.back();
    assertion(toRemove > 0, toRemove);
    DEBUG("Removing " << toRemove << " cols from mannifold mapping least-squares system with " << getLSSystemCols() << " cols");
    assertion(_matrixF.cols() == _matrixC.cols(), _matrixF.cols(), _matrixC.cols());
    assertion(getLSSystemCols() > toRemove, getLSSystemCols(), toRemove);

    // remove columns
    for (int i = 0; i < toRemove; i++) {
      utils::removeColumnFromMatrix(_matrixF, _matrixF.cols() - 1);
      utils::removeColumnFromMatrix(_matrixC, _matrixC.cols() - 1);
    }
    _matrixCols.pop_back();
  }

  _matrixCols.push_front(0);
  _firstIteration = true;
}

/** ---------------------------------------------------------------------------------------------
 *         removeMatrixColumn()
 *
 * @brief: removes a column from the least squares system, i. e., from the matrices F and C
 *  ---------------------------------------------------------------------------------------------
 */
void MMPostProcessing::removeMatrixColumn(
    int columnIndex)
{
  TRACE(columnIndex, _matrixF.cols());

  // debugging information, can be removed
  deletedColumns++;

  assertion(_matrixF.cols() > 1);
  utils::removeColumnFromMatrix(_matrixF, columnIndex);
  utils::removeColumnFromMatrix(_matrixC, columnIndex);

  // Reduce column count
  std::deque<int>::iterator iter = _matrixCols.begin();
  int                       cols = 0;
  while (iter != _matrixCols.end()) {
    cols += *iter;
    if (cols > columnIndex) {
      assertion(*iter > 0);
      *iter -= 1;
      if (*iter == 0) {
        _matrixCols.erase(iter);
      }
      break;
    }
    iter++;
  }
}

void MMPostProcessing::exportState(
    io::TXTWriter &writer)
{
}

void MMPostProcessing::importState(
    io::TXTReader &reader)
{
}

int MMPostProcessing::getDeletedColumns()
{
  return deletedColumns;
}

int MMPostProcessing::getLSSystemCols()
{
  int cols = 0;
  for (int col : _matrixCols) {
    cols += col;
  }
  //if(_hasNodesOnInterface){
  //	assertion(cols == _matrixF.cols(), cols, _matrixF.cols(), _matrixCols);
  //	assertion(cols == _matrixC.cols(), cols, _matrixC.cols());
  //}

  return cols;
}

int MMPostProcessing::getLSSystemRows()
{
  if (utils::MasterSlave::_masterMode || utils::MasterSlave::_slaveMode) {
    return _dimOffsets.back();
  }
  return _fineResiduals.size();
  //return _matrixF.rows();
}
}
}
} // namespace precice, cplscheme, impl
