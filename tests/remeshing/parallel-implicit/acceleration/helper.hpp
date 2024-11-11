#pragma once

#include <string>

#include <precice/Participant.hpp>
#include "testing/QuickTest.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

namespace precice::tests::remesh::parallelImplicit {

using testing::QuickTest;
using testing::operator""_mesh;
using testing::operator""_read;
using testing::operator""_write;

namespace acceleration {

// The data written in the acceleration tests are always 40, 30, 15
// first, second, third, and forth are read data from the first, second, third and forth iteration
// The end of a timewindow and the end of the simulation are never accelerated and thus 15

inline void runResetA(testing::TestContext &context, double first, double second, double third, double forth)
{
  constexpr double y = 0.0;

  BOOST_REQUIRE(context.size == 1);
  BOOST_REQUIRE(context.rank == 0);
  Participant p{context.name, context.config(), 0, 1};

  // A - Adaptive Geometry
  if (context.isNamed("A")) {
    QuickTest(p, "MA"_mesh, "DB"_read, "DA"_write)
        .setVertices({1.0, y, 2.0, y})
        .initialize()
        .expectWriteCheckpoint()
        .expectAll(00.00)
        .writeAll(40.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(first)
        .writeAll(30.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(second)
        .writeAll(15.00)
        .advance()

        .expectWriteCheckpoint()
        .expectAll(15)
        .writeAll(40.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(third)
        .writeAll(30.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(forth)
        .writeAll(15.00)
        .advance()

        .expectWriteCheckpoint()
        .expectAll(15)
        .resetMesh()
        .setVertices({2.0, y})
        .writeAll(40.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(first)
        .writeAll(30.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(second)
        .writeAll(15.00)
        .advance()

        .expectWriteCheckpoint()
        .expectAll(15)
        .writeAll(40.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(third)
        .writeAll(30.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(forth)
        .writeAll(15.00)
        .advance()

        .expectCouplingCompleted()
        .expectAll(15)
        .finalize();
  }
  // B - Changing Geometry
  if (context.isNamed("B")) {
    QuickTest(p, "MB"_mesh, "DA"_read, "DB"_write)
        .setVertices({1.0, y, 2.0, y})
        .initialize()
        .expectWriteCheckpoint()
        .expectAll(00.00)
        .writeAll(40.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(first)
        .writeAll(30.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(second)
        .writeAll(15.00)
        .advance()

        .expectWriteCheckpoint()
        .expectAll(15.00)
        .writeAll(40.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(third)
        .writeAll(30.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(forth)
        .writeAll(15.00)
        .advance()

        .expectWriteCheckpoint()
        .expectAll(15.00)
        .writeAll(40.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(first)
        .writeAll(30.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(second)
        .writeAll(15.00)
        .advance()

        .expectWriteCheckpoint()
        .expectAll(15.00)
        .writeAll(40.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(third)
        .writeAll(30.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(forth)
        .writeAll(15.00)
        .advance()

        .expectCouplingCompleted()
        .expectAll(15)
        .finalize();
  }
}

inline void runResetAConstant(testing::TestContext &context)
{
  runResetA(context, 20, 25, 27.5, 28.75);
}

inline void runResetAAitken(testing::TestContext &context)
{
  runResetA(context, 20, 26.666666666, 27.5, 28.888888888);
}

inline void runResetAIQNILS(testing::TestContext &context)
{
  runResetA(context, 4, 11.4285714285714, 23.280254777, 26.324041811);
}

inline void runResetAIQNIMVJ(testing::TestContext &context)
{
  runResetA(context, 4, 11.428571428, 6.560509554, -120.204081632);
}

inline void runResetBoth(testing::TestContext &context, double first, double second, double third, double forth)
{
  constexpr double y = 0.0;

  BOOST_REQUIRE(context.size == 1);
  BOOST_REQUIRE(context.rank == 0);
  Participant p{context.name, context.config(), 0, 1};

  // A - Adaptive Geometry
  if (context.isNamed("A")) {
    QuickTest(p, "MA"_mesh, "DB"_read, "DA"_write)
        .setVertices({1.0, y, 2.0, y})
        .initialize()
        .expectWriteCheckpoint()
        .expectAll(00.00)
        .writeAll(40.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(first)
        .writeAll(30.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(second)
        .writeAll(15.00)
        .advance()

        .expectWriteCheckpoint()
        .expectAll(15.00) // end of lw
        .writeAll(40.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(third)
        .writeAll(30.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(forth)
        .writeAll(15.00)
        .advance()

        .expectWriteCheckpoint()
        .expectAll(15.00) // end of lw
        .resetMesh()
        .setVertices({2.0, y})
        .writeAll(40.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(first)
        .writeAll(30.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(second)
        .writeAll(15.00)
        .advance()

        .expectWriteCheckpoint()
        .expectAll(15.00) // end sample
        .writeAll(40.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(third)
        .writeAll(30.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(forth)
        .writeAll(15.00)
        .advance()

        .expectCouplingCompleted()
        .expectAll(15.00) // end sample
        .finalize();
  }
  // B - Changing Geometry
  if (context.isNamed("B")) {
    QuickTest(p, "MB"_mesh, "DA"_read, "DB"_write)
        .setVertices({1.0, y, 2.0, y})
        .initialize()
        .expectWriteCheckpoint()
        .expectAll(00.00)
        .writeAll(40.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(first)
        .writeAll(30.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(second)
        .writeAll(15.00)
        .advance()

        .expectWriteCheckpoint()
        .expectAll(15.00) // end of tw
        .writeAll(40.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(third)
        .writeAll(30.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(forth)
        .writeAll(15.00)
        .advance()

        .expectWriteCheckpoint()
        .expectAll(15.00) // end of tw
        .resetMesh()
        .setVertices({2.0, y, 3.0, y})
        .writeAll(40.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(first)
        .writeAll(30.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(second)
        .writeAll(15.00)
        .advance()

        .expectWriteCheckpoint()
        .expectAll(15.00) // end of tw
        .writeAll(40.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(third)
        .writeAll(30.00)
        .advance()

        .expectReadCheckpoint()
        .expectAll(forth)
        .writeAll(15.00)
        .advance()

        .expectCouplingCompleted()
        .expectAll(15.00) // end sample
        .finalize();
  }
}

inline void runResetBothConstant(testing::TestContext &context)
{
  runResetBoth(context, 20, 25, 27.5, 28.75);
}

inline void runResetBothAitken(testing::TestContext &context)
{
  runResetBoth(context, 20, 26.666666666, 27.5, 28.888888888);
}

inline void runResetBothIQNILS(testing::TestContext &context)
{
  runResetBoth(context, 4, 11.4285714285714, 23.280254777, 26.324041811);
}

inline void runResetBothIQNIMVJ(testing::TestContext &context)
{
  runResetBoth(context, 4, 11.4285714285714, 6.560509554, -120.204081632);
}
} // namespace acceleration
} // namespace precice::tests::remesh::parallelImplicit
