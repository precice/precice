#ifndef PRECICE_NO_GINKGO

#include "mapping/GinkgoRadialBasisFctSolver.hpp"
#include "mapping/impl/BasisFunctions.hpp"
#include "mapping/tests/RadialBasisFctHelper.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Utils.hpp"
#include "mesh/Vertex.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

using namespace precice;
using namespace precice::mesh;
using namespace precice::mapping;
using namespace precice::testing;
using precice::testing::TestContext;

BOOST_AUTO_TEST_SUITE(MappingTests)
BOOST_AUTO_TEST_SUITE(GinkgoRadialBasisFunctionSolver)

#undef doLocalCode
#define doLocalCode(Type, function, polynomial, EXECUTOR, SOLVER)                                                                                                                                                    \
  {                                                                                                                                                                                                                  \
    MappingConfiguration::GinkgoParameter gpm;                                                                                                                                                                       \
    gpm.executor          = EXECUTOR;                                                                                                                                                                                \
    gpm.deviceId          = 0;                                                                                                                                                                                \
    gpm.solver            = SOLVER;                                                                                                                                                                                  \
    gpm.maxIterations     = 100;                                                                                                                                                                                     \
    gpm.usePreconditioner = false;                                                                                                                                                                                   \
    RadialBasisFctMapping<GinkgoRadialBasisFctSolver<Type>, MappingConfiguration::GinkgoParameter> consistentMap2D(Mapping::CONSISTENT, 2, function, {{false, false, false}}, polynomial, gpm);                      \
    perform2DTestConsistentMapping(consistentMap2D);                                                                                                                                                                 \
    RadialBasisFctMapping<GinkgoRadialBasisFctSolver<Type>, MappingConfiguration::GinkgoParameter> consistentMap2DVector(Mapping::CONSISTENT, 2, function, {{false, false, false}}, polynomial, gpm);                \
    perform2DTestConsistentMappingVector(consistentMap2DVector);                                                                                                                                                     \
    RadialBasisFctMapping<GinkgoRadialBasisFctSolver<Type>, MappingConfiguration::GinkgoParameter> consistentMap3D(Mapping::CONSISTENT, 3, function, {{false, false, false}}, polynomial, gpm);                      \
    perform3DTestConsistentMapping(consistentMap3D);                                                                                                                                                                 \
    RadialBasisFctMapping<GinkgoRadialBasisFctSolver<Type>, MappingConfiguration::GinkgoParameter> scaledConsistentMap2D(Mapping::SCALED_CONSISTENT_SURFACE, 2, function, {{false, false, false}}, polynomial, gpm); \
    perform2DTestScaledConsistentMapping(scaledConsistentMap2D);                                                                                                                                                     \
    RadialBasisFctMapping<GinkgoRadialBasisFctSolver<Type>, MappingConfiguration::GinkgoParameter> scaledConsistentMap3D(Mapping::SCALED_CONSISTENT_SURFACE, 3, function, {{false, false, false}}, polynomial, gpm); \
    perform3DTestScaledConsistentMapping(scaledConsistentMap3D);                                                                                                                                                     \
    RadialBasisFctMapping<GinkgoRadialBasisFctSolver<Type>, MappingConfiguration::GinkgoParameter> conservativeMap2D(Mapping::CONSERVATIVE, 2, function, {{false, false, false}}, polynomial, gpm);                  \
    perform2DTestConservativeMapping(conservativeMap2D);                                                                                                                                                             \
    RadialBasisFctMapping<GinkgoRadialBasisFctSolver<Type>, MappingConfiguration::GinkgoParameter> conservativeMap2DVector(Mapping::CONSERVATIVE, 2, function, {{false, false, false}}, polynomial, gpm);            \
    perform2DTestConservativeMappingVector(conservativeMap2DVector);                                                                                                                                                 \
    RadialBasisFctMapping<GinkgoRadialBasisFctSolver<Type>, MappingConfiguration::GinkgoParameter> conservativeMap3D(Mapping::CONSERVATIVE, 3, function, {{false, false, false}}, polynomial, gpm);                  \
    perform3DTestConservativeMapping(conservativeMap3D);                                                                                                                                                             \
  }

#define TEST_FOR_ALL_RBFS(EXECUTOR, SOLVER)                                              \
  BOOST_AUTO_TEST_CASE(MapThinPlateSplines)                                              \
  {                                                                                      \
    PRECICE_TEST(1_rank);                                                                \
    ThinPlateSplines fct;                                                                \
    doLocalCode(ThinPlateSplines, fct, Polynomial::SEPARATE, EXECUTOR, SOLVER);          \
  }                                                                                      \
  BOOST_AUTO_TEST_CASE(MapMultiquadrics)                                                 \
  {                                                                                      \
    PRECICE_TEST(1_rank);                                                                \
    Multiquadrics fct(1e-3);                                                             \
    doLocalCode(Multiquadrics, fct, Polynomial::SEPARATE, EXECUTOR, SOLVER);             \
  }                                                                                      \
  BOOST_AUTO_TEST_CASE(MapInverseMultiquadrics)                                          \
  {                                                                                      \
    PRECICE_TEST(1_rank);                                                                \
    InverseMultiquadrics fct(1e-3);                                                      \
    doLocalCode(InverseMultiquadrics, fct, Polynomial::SEPARATE, EXECUTOR, SOLVER);      \
  }                                                                                      \
  BOOST_AUTO_TEST_CASE(MapVolumeSplines)                                                 \
  {                                                                                      \
    PRECICE_TEST(1_rank);                                                                \
    VolumeSplines fct;                                                                   \
    doLocalCode(VolumeSplines, fct, Polynomial::SEPARATE, EXECUTOR, SOLVER);             \
  }                                                                                      \
  BOOST_AUTO_TEST_CASE(MapGaussian)                                                      \
  {                                                                                      \
    PRECICE_TEST(1_rank);                                                                \
    Gaussian fct(1.0);                                                                   \
    doLocalCode(Gaussian, fct, Polynomial::SEPARATE, EXECUTOR, SOLVER);                  \
  }                                                                                      \
  BOOST_AUTO_TEST_CASE(MapCompactThinPlateSplinesC2)                                     \
  {                                                                                      \
    PRECICE_TEST(1_rank);                                                                \
    double                    supportRadius = 1.2;                                       \
    CompactThinPlateSplinesC2 fct(supportRadius);                                        \
    doLocalCode(CompactThinPlateSplinesC2, fct, Polynomial::SEPARATE, EXECUTOR, SOLVER); \
  }                                                                                      \
  BOOST_AUTO_TEST_CASE(MapCompactPolynomialC0)                                           \
  {                                                                                      \
    PRECICE_TEST(1_rank);                                                                \
    double              supportRadius = 1.2;                                             \
    CompactPolynomialC0 fct(supportRadius);                                              \
    doLocalCode(CompactPolynomialC0, fct, Polynomial::SEPARATE, EXECUTOR, SOLVER);       \
  }                                                                                      \
  BOOST_AUTO_TEST_CASE(MapCompactPolynomialC2)                                           \
  {                                                                                      \
    PRECICE_TEST(1_rank);                                                                \
    double              supportRadius = 1.2;                                             \
    CompactPolynomialC2 fct(supportRadius);                                              \
    doLocalCode(CompactPolynomialC2, fct, Polynomial::SEPARATE, EXECUTOR, SOLVER);       \
  }                                                                                      \
  BOOST_AUTO_TEST_CASE(MapCompactPolynomialC4)                                           \
  {                                                                                      \
    PRECICE_TEST(1_rank);                                                                \
    double              supportRadius = 1.2;                                             \
    CompactPolynomialC4 fct(supportRadius);                                              \
    doLocalCode(CompactPolynomialC4, fct, Polynomial::SEPARATE, EXECUTOR, SOLVER);       \
  }                                                                                      \
  BOOST_AUTO_TEST_CASE(MapCompactPolynomialC6)                                           \
  {                                                                                      \
    PRECICE_TEST(1_rank);                                                                \
    double              supportRadius = 1.2;                                             \
    CompactPolynomialC6 fct(supportRadius);                                              \
    doLocalCode(CompactPolynomialC6, fct, Polynomial::SEPARATE, EXECUTOR, SOLVER);       \
  }

BOOST_AUTO_TEST_SUITE(Reference)
TEST_FOR_ALL_RBFS("reference-executor", "gmres-solver");
BOOST_AUTO_TEST_SUITE_END()

#ifdef PRECICE_WITH_OMP
BOOST_AUTO_TEST_SUITE(OpenMP)
TEST_FOR_ALL_RBFS("omp-executor", "gmres-solver");
BOOST_AUTO_TEST_SUITE_END()
#endif

#ifdef PRECICE_WITH_CUDA
BOOST_AUTO_TEST_SUITE(Cuda)
TEST_FOR_ALL_RBFS("cuda-executor", "gmres-solver");
BOOST_AUTO_TEST_SUITE_END()
#endif

#ifdef PRECICE_WITH_CUDA
BOOST_AUTO_TEST_SUITE(cuSolver)
TEST_FOR_ALL_RBFS("cuda-executor", "qr-solver");
BOOST_AUTO_TEST_SUITE_END()
#endif

#ifdef PRECICE_WITH_HIP
BOOST_AUTO_TEST_SUITE(Hip)
TEST_FOR_ALL_RBFS("hip-executor", "gmres-solver");
BOOST_AUTO_TEST_SUITE_END()
#endif

#ifdef PRECICE_WITH_HIP
BOOST_AUTO_TEST_SUITE(hipSolver)
TEST_FOR_ALL_RBFS("hip-executor", "qr-solver");
BOOST_AUTO_TEST_SUITE_END()
#endif

#undef TEST_FOR_ALL_RBFS
#undef doLocalCode

BOOST_AUTO_TEST_SUITE_END() // RadialBasisFunctionMapping
BOOST_AUTO_TEST_SUITE_END()

#endif
