#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>
#include <precice/impl/SolverInterfaceImpl.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Parallel)
BOOST_AUTO_TEST_CASE(TestFinalize)
{
  PRECICE_TEST("SolverOne"_on(2_ranks), "SolverTwo"_on(2_ranks));

  if (context.isNamed("SolverOne")) {
    precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
    auto                     meshName = "MeshOne";
    double                   xCoord   = 0.0 + context.rank;
    interface.setMeshVertex(meshName, Eigen::Vector3d(xCoord, 0.0, 0.0).data());
    interface.initialize();
    BOOST_TEST(precice::testing::WhiteboxAccessor::impl(interface).mesh("MeshOne").vertices().size() == 1);
    BOOST_TEST(precice::testing::WhiteboxAccessor::impl(interface).mesh("MeshTwo").vertices().size() == 1);
    interface.finalize();
  } else {
    BOOST_TEST(context.isNamed("SolverTwo"));
    precice::SolverInterface interface(context.name, context.config(), context.rank, context.size);
    auto                     meshName = "MeshTwo";
    double                   xCoord   = 0.0 + context.rank;
    interface.setMeshVertex(meshName, Eigen::Vector3d(xCoord, 0.0, 0.0).data());
    interface.initialize();
    BOOST_TEST(precice::testing::WhiteboxAccessor::impl(interface).mesh("MeshTwo").vertices().size() == 1);
    interface.finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // Integration
BOOST_AUTO_TEST_SUITE_END() // Parallel

#endif // PRECICE_NO_MPI
