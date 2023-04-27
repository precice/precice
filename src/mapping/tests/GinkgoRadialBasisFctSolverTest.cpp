#ifndef PRECICE_NO_GINKGO

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
#define doLocalCode(Type, function, polynomial, EXECUTOR, SOLVER)                                                                                           \
  {                                                                                                                                                         \
    bool                                  useEigen = false;                                                                                                 \
    MappingConfiguration::GinkgoParameter gpm;                                                                                                              \
    gpm.executor      = EXECUTOR;                                                                                                                           \
    gpm.solver        = SOLVER;                                                                                                                             \
    gpm.maxIterations = 1e2;                                                                                                                                \
    RadialBasisFctMapping<Type> consistentMap2D(Mapping::CONSISTENT, 2, function, {{false, false, false}}, polynomial, useEigen, gpm);                      \
    perform2DTestConsistentMapping(consistentMap2D);                                                                                                        \
    RadialBasisFctMapping<Type> consistentMap2DVector(Mapping::CONSISTENT, 2, function, {{false, false, false}}, polynomial, useEigen, gpm);                \
    perform2DTestConsistentMappingVector(consistentMap2DVector);                                                                                            \
    RadialBasisFctMapping<Type> consistentMap3D(Mapping::CONSISTENT, 3, function, {{false, false, false}}, polynomial, useEigen, gpm);                      \
    perform3DTestConsistentMapping(consistentMap3D);                                                                                                        \
    RadialBasisFctMapping<Type> scaledConsistentMap2D(Mapping::SCALED_CONSISTENT_SURFACE, 2, function, {{false, false, false}}, polynomial, useEigen, gpm); \
    perform2DTestScaledConsistentMapping(scaledConsistentMap2D);                                                                                            \
    RadialBasisFctMapping<Type> scaledConsistentMap3D(Mapping::SCALED_CONSISTENT_SURFACE, 3, function, {{false, false, false}}, polynomial, useEigen, gpm); \
    perform3DTestScaledConsistentMapping(scaledConsistentMap3D);                                                                                            \
    RadialBasisFctMapping<Type> conservativeMap2D(Mapping::CONSERVATIVE, 2, function, {{false, false, false}}, polynomial, useEigen, gpm);                  \
    perform2DTestConservativeMapping(conservativeMap2D);                                                                                                    \
    RadialBasisFctMapping<Type> conservativeMap2DVector(Mapping::CONSERVATIVE, 2, function, {{false, false, false}}, polynomial, useEigen, gpm);            \
    perform2DTestConservativeMappingVector(conservativeMap2DVector);                                                                                        \
    RadialBasisFctMapping<Type> conservativeMap3D(Mapping::CONSERVATIVE, 3, function, {{false, false, false}}, polynomial, useEigen, gpm);                  \
    perform3DTestConservativeMapping(conservativeMap3D);                                                                                                    \
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

TEST_FOR_ALL_RBFS("reference-executor", "cg-solver");

BOOST_AUTO_TEST_SUITE_END()

#undef TEST_FOR_ALL_RBFS
#undef doLocalCode

BOOST_AUTO_TEST_SUITE_END() // RadialBasisFunctionMapping
BOOST_AUTO_TEST_SUITE_END()

#endif