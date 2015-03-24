// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_ITESTS_COUPLINGINTERFACETEST_HPP_
#define PRECICE_ITESTS_COUPLINGINTERFACETEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"
#include "precice/SolverInterface.hpp"
#include <string>
#include <vector>

namespace precice {
namespace tests {

/**
 * @brief Provides tests for class SolverInterface.
 */
class SolverInterfaceTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor, takes relative/absolute path to this directory.
   */
  SolverInterfaceTest();

  /**
   * @brief Destructor.
   */
  virtual ~SolverInterfaceTest() {}

  /**
   * @brief Retrieves path to test directory.
   */
  virtual void setUp();

  /**
   * @brief Runs all tests.
   */
  virtual void run();

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  // @brief Path to this directory.
  std::string _pathToTests;

# ifndef PRECICE_NO_MPI

  /**
   * @brief As SolverInterface::configure(), but without changing logging config.
   */
  void configureSolverInterface (
    const std::string& configFilename,
    SolverInterface&   interface );

  /**
   * @brief Test reading of a full features coupling configuration file.
   */
  void testConfiguration();

  /**
   * @brief Test to run simple "do nothing" coupling between two solvers.
   */
  void testExplicit();

  /**
   * @brief Test to run a simple "do nothing" coupling with subcycling solvers.
   */
  void testExplicitWithSubcycling();

  /**
   * @brief One solver uses incremental position set, read/write methods.
   */
  void testExplicitWithDataExchange();

  /**
   * @brief The second solver initializes the data of the first.
   *
   * A mapping is employed for the second solver, i.e., at the end of
   * initializeData(), the mapping needs to be invoked.
   */
  void testExplicitWithDataInitialization();

  /**
   * @brief One solver uses block set/get/read/write methods.
   */
  void testExplicitWithBlockDataExchange();

  /**
   * @brief Runs a coupled simulation where one solver supplies a geometry.
   *
   * SolverOne only reads the displacements of the geometry and checks whether
   * they are equals to the coordinates of SolverTwo. SolverTwo creates and
   * displaces the coordinates.
   */
  void testExplicitWithSolverGeometry();

  /**
   * @brief Runs a coupled simulation where one solver displaces a geometry.
   */
  void testExplicitWithDisplacingGeometry();

  /**
   * @brief Runs a coupled simulation with one solver doing pre-coupled timesteps.
   */
  //   void testCoupledSimulationWithPreCoupledTimesteps ();

  /**
   * @brief Runs a coupled sim. with data scaling applied.
   *
   * SolverOne writes vector data on a cube geometry. The data values are defined
   * and stay constant over the coupling cycles. SolverTwo has a scaling of the
   * values activated and reads the scaled values.
   */
  void testExplicitWithDataScaling();

  /**
   * @brief Explicit coupled checkpointing with stationary data mapping.
   */
  void testExplicitWithCheckpointingStatMapping();

  /**
   * @brief Test simple coupled simulation with coupling iterations.
   */
  void testImplicit();

  /**
   * @brief Tests an implicit coupled simulation where one of the solvers
   *        performs a subcycling.
   */
  void testImplicitWithSubcycling();

  /**
   * @brief Implicit coupled checkpointing with stationary data mapping.
   */
  void testImplicitWithCheckpointingMappingStat();

  /**
   * @brief Initializes, runs and finalizes a coupled simulation.
   *
   * Uses the timestep prescribed by the coupling interface.
   *
   * @param solverName [IN] Name of the solver accessing the coupling interface.
   * @param configurationFileName [IN] Name of file holding xml configuration.
   */
  void runSolver (
    const std::string& solverName,
    const std::string& configurationFileName,
    int&               timestepsComputed,
    double&            timeComputed );



  /**
   * @brief Tests various distributed communication schemes.
   */
  void testDistributedCommunications();


  /**
   * @brief Tests stationary mapping with solver provided meshes.
   */
  void testStationaryMappingWithSolverMesh();

  /**
   * @brief Buggy simulation setup of FSI coupling between Flite and Calculix.
   *
   * Bug: after first call of advance by Flite the mapped forces are value NaN.
   *
   * Some information about the coupling:
   * - explicit coupling scheme
   * - Flite (incompressible Navier-Stokes) starts simulation
   * - Mapping is done on Flite side with RBF
   */
  void testBug();

  /**
   * @brief Three solvers are coupled in a fork S2 <-> S1 <-> S3.
   *
   * Both couplings are explicit, solver 1 provides the mesh to the other two
   * solvers.
   */
  void testThreeSolvers();

  /**
   * @brief Four solvers are multi-coupled.
   *
   */
  void testMultiCoupling();

  /**
   * @brief All NASTIN meshes are restarted. One provided as well as one "from"
   * are tested.
   *
   */
  void testNASTINMeshRestart();

  void runThreeSolvers (
      const std::string&      configFilename,
      const std::vector<int>& expectedCallsOfAdvance );

# endif // defined( not PRECICE_NO_MPI )
};

}} // namespace precice, tests

#endif /* PRECICE_ITESTS_COUPLINGINTERFACETEST_HPP_ */
