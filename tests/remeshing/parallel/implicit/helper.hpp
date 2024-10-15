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

/* data format 12.34
 * 1 = rank
 * 2 = vertex x coordinate
 * 3 = time window starting with 0
 * 4 = iteration starting with 0
 */

namespace noop {

inline void runResetA(testing::TestContext &context)
{
  constexpr double y = 0.0;

  Participant p{context.name, context.config(), context.rank, context.size};

  // A - Static Geometry
  if (context.isNamed("A")) {
    if (context.isPrimary()) {
      QuickTest(p, "MA"_mesh, "DB"_read, "DA"_write)
          .setVertices({1.0, y, 2.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({01.00, 02.00})
          .advance()

          .expectReadCheckpoint()
          .expect({01.00, 02.00})
          .write({01.01, 02.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({01.01, 02.01})
          .resetMesh()
          .setVertices({1.0, y, 2.0, y})
          .write({01.10, 02.10})
          .advance()

          .expectReadCheckpoint()
          .expect({01.10, 02.10})
          .write({01.11, 02.11})
          .advance()

          .expectCouplingCompleted()
          .expect({01.11, 02.11})
          .finalize();
    } else {
      QuickTest(p, "MA"_mesh, "DB"_read, "DA"_write)
          .setVertices({3.0, y, 4.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({13.00, 14.00})
          .advance()

          .expectReadCheckpoint()
          .expect({13.00, 14.00})
          .write({13.01, 14.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({13.01, 14.01})
          .resetMesh()
          .setVertices({3.0, y, 4.0, y})
          .write({13.10, 14.10})
          .advance()

          .expectReadCheckpoint()
          .expect({13.10, 14.10})
          .write({13.11, 14.11})
          .advance()

          .expectCouplingCompleted()
          .expect({13.11, 14.11})
          .finalize();
    }
  }
  // B - Adaptive Geometry
  if (context.isNamed("B")) {
    if (context.isPrimary()) {
      QuickTest(p, "MB"_mesh, "DA"_read, "DB"_write)
          .setVertices({1.0, y, 2.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({01.00, 02.00})
          .advance()

          .expectReadCheckpoint()
          .expect({01.00, 02.00})
          .write({01.01, 02.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({01.01, 02.01})
          .write({01.10, 02.10})
          .advance()

          .expectReadCheckpoint()
          .expect({01.10, 02.10})
          .write({01.11, 02.11})
          .advance()

          .expectCouplingCompleted()
          .expect({01.11, 02.11})
          .finalize();
    } else {
      QuickTest(p, "MB"_mesh, "DA"_read, "DB"_write)
          .setVertices({3.0, y, 4.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({13.00, 14.00})
          .advance()

          .expectReadCheckpoint()
          .expect({13.00, 14.00})
          .write({13.01, 14.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({13.01, 14.01})
          .write({13.10, 14.10})
          .advance()

          .expectReadCheckpoint()
          .expect({13.10, 14.10})
          .write({13.11, 14.11})
          .advance()

          .expectCouplingCompleted()
          .expect({13.11, 14.11})
          .finalize();
    }
  }
}

inline void runResetBoth(testing::TestContext &context)
{
  constexpr double y = 0.0;

  Participant p{context.name, context.config(), context.rank, context.size};

  auto qt = context.isNamed("A") ? QuickTest(p, "MA"_mesh, "DB"_read, "DA"_write) : QuickTest(p, "MB"_mesh, "DA"_read, "DB"_write);

  if (context.isPrimary()) {
    qt.setVertices({1.0, y, 2.0, y})
        .initialize()
        .expectWriteCheckpoint()
        .expect({00.00, 00.00})
        .write({01.00, 02.00})
        .advance()

        .expectReadCheckpoint()
        .expect({01.00, 02.00})
        .write({01.01, 02.01})
        .advance()

        .expectWriteCheckpoint()
        .expect({01.01, 02.01})
        .resetMesh()
        .setVertices({1.0, y, 2.0, y})
        .write({01.10, 02.10})
        .advance()

        .expectReadCheckpoint()
        .expect({01.10, 02.10})
        .write({01.11, 02.11})
        .advance()

        .expectCouplingCompleted()
        .expect({01.11, 02.11})
        .finalize();
  } else {
    qt.setVertices({3.0, y, 4.0, y})
        .initialize()
        .expectWriteCheckpoint()
        .expect({00.00, 00.00})
        .write({13.00, 14.00})
        .advance()

        .expectReadCheckpoint()
        .expect({13.00, 14.00})
        .write({13.01, 14.01})
        .advance()

        .expectWriteCheckpoint()
        .expect({13.01, 14.01})
        .resetMesh()
        .setVertices({3.0, y, 4.0, y})
        .write({13.10, 14.10})
        .advance()

        .expectReadCheckpoint()
        .expect({13.10, 14.10})
        .write({13.11, 14.11})
        .advance()

        .expectCouplingCompleted()
        .expect({13.11, 14.11})
        .finalize();
  }
}
} // namespace noop

namespace changemapping {

inline void runResetA(testing::TestContext &context)
{
  constexpr double y = 0.0;

  Participant p{context.name, context.config(), context.rank, context.size};

  // A - Static Geometry
  if (context.isNamed("A")) {
    if (context.isPrimary()) {
      QuickTest(p, "MA"_mesh, "DB"_read, "DA"_write)
          .setVertices({1.0, y, 2.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({01.00, 02.00})
          .advance()

          .expectReadCheckpoint()
          .expect({01.00, 02.00})
          .write({01.01, 02.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({01.01, 02.01})
          .resetMesh()
          .setVertices({2.0, y})
          .write({02.10})
          .advance()

          .expectReadCheckpoint()
          .expect({02.10})
          .write({02.11})
          .advance()

          .expectCouplingCompleted()
          .expect({02.11})
          .finalize();
    } else {
      QuickTest(p, "MA"_mesh, "DB"_read, "DA"_write)
          .setVertices({3.0, y, 4.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({13.00, 14.00})
          .advance()

          .expectReadCheckpoint()
          .expect({13.00, 14.00})
          .write({13.01, 14.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({13.01, 14.01})
          .resetMesh()
          .setVertices({3.0, y})
          .write({13.10})
          .advance()

          .expectReadCheckpoint()
          .expect({13.10})
          .write({13.11})
          .advance()

          .expectCouplingCompleted()
          .expect({13.11})
          .finalize();
    }
  }
  // B - Adaptive Geometry
  if (context.isNamed("B")) {
    if (context.isPrimary()) {
      QuickTest(p, "MB"_mesh, "DA"_read, "DB"_write)
          .setVertices({1.0, y, 2.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({01.00, 02.00})
          .advance()

          .expectReadCheckpoint()
          .expect({01.00, 02.00})
          .write({01.01, 02.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({01.01, 02.01})
          .write({01.10, 02.10})
          .advance()

          .expectReadCheckpoint()
          .expect({02.10, 02.10})
          .write({01.11, 02.11})
          .advance()

          .expectCouplingCompleted()
          .expect({02.11, 02.11})
          .finalize();
    } else {
      QuickTest(p, "MB"_mesh, "DA"_read, "DB"_write)
          .setVertices({3.0, y, 4.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({13.00, 14.00})
          .advance()

          .expectReadCheckpoint()
          .expect({13.00, 14.00})
          .write({13.01, 14.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({13.01, 14.01})
          .write({13.10, 14.10})
          .advance()

          .expectReadCheckpoint()
          .expect({13.10, 13.10})
          .write({13.11, 14.11})
          .advance()

          .expectCouplingCompleted()
          .expect({13.11, 13.11})
          .finalize();
    }
  }
}

inline void runResetBoth(testing::TestContext &context)
{
  constexpr double y = 0.0;

  Participant p{context.name, context.config(), context.rank, context.size};

  if (context.isNamed("A")) {
    QuickTest qt(p, "MA"_mesh, "DB"_read, "DA"_write);
    if (context.isPrimary()) {
      qt.setVertices({1.0, y, 2.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({01.00, 02.00})
          .advance()

          .expectReadCheckpoint()
          .expect({01.00, 02.00})
          .write({01.01, 02.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({01.01, 02.01})
          .resetMesh()
          .setVertices({1.0, y})
          .write({01.10})
          .advance()

          .expectReadCheckpoint()
          .expect({02.10})
          .write({01.11})
          .advance()

          .expectCouplingCompleted()
          .expect({02.11})
          .finalize();
    } else {
      qt.setVertices({4.0, y, 5.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({14.00, 15.00})
          .advance()

          .expectReadCheckpoint()
          .expect({14.00, 15.00})
          .write({14.01, 15.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({14.01, 15.01})
          .resetMesh()
          .setVertices({4.0, y})
          .write({14.10})
          .advance()

          .expectReadCheckpoint()
          .expect({15.10})
          .write({14.11})
          .advance()

          .expectCouplingCompleted()
          .expect({15.11})
          .finalize();
    }
  } else {
    QuickTest qt(p, "MB"_mesh, "DA"_read, "DB"_write);
    if (context.isPrimary()) {
      qt.setVertices({1.0, y, 2.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({01.00, 02.00})
          .advance()

          .expectReadCheckpoint()
          .expect({01.00, 02.00})
          .write({01.01, 02.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({01.01, 02.01})
          .resetMesh()
          .setVertices({2.0, y})
          .write({02.10})
          .advance()

          .expectReadCheckpoint()
          .expect({01.10})
          .write({02.11})
          .advance()

          .expectCouplingCompleted()
          .expect({01.11})
          .finalize();
    } else {
      qt.setVertices({4.0, y, 5.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({14.00, 15.00})
          .advance()

          .expectReadCheckpoint()
          .expect({14.00, 15.00})
          .write({14.01, 15.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({14.01, 15.01})
          .resetMesh()
          .setVertices({5.0, y})
          .write({15.10})
          .advance()

          .expectReadCheckpoint()
          .expect({14.10})
          .write({15.11})
          .advance()

          .expectCouplingCompleted()
          .expect({14.11})
          .finalize();
    }
  }
}

} // namespace changemapping

namespace changepartition {

/// Changes partitioning from 12|34 to 1|234 and 123|4
inline void runOverlapBoth(testing::TestContext &context)
{
  constexpr double y = 0.0;

  Participant p{context.name, context.config(), context.rank, context.size};

  // A - Static Geometry
  if (context.isNamed("A")) {
    if (context.isPrimary()) {
      QuickTest(p, "MA"_mesh, "DB"_read, "DA"_write)
          .setVertices({1.0, y, 2.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({01.00, 02.00})
          .advance()

          .expectReadCheckpoint()
          .expect({01.00, 02.00})
          .write({01.01, 02.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({01.01, 02.01})
          .resetMesh()
          .setVertices({1.0, y, 2.0, y, 3.0, y})
          .write({01.10, 02.10, 03.10})
          .advance()

          .expectReadCheckpoint()
          .expect({01.10, 12.10, 13.10})
          .write({01.11, 02.11, 03.11})
          .advance()

          .expectCouplingCompleted()
          .expect({01.11, 12.11, 13.11})
          .finalize();
    } else {
      QuickTest(p, "MA"_mesh, "DB"_read, "DA"_write)
          .setVertices({3.0, y, 4.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({13.00, 14.00})
          .advance()

          .expectReadCheckpoint()
          .expect({13.00, 14.00})
          .write({13.01, 14.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({13.01, 14.01})
          .resetMesh()
          .setVertices({4.0, y})
          .write({14.10})
          .advance()

          .expectReadCheckpoint()
          .expect({14.10})
          .write({14.11})
          .advance()

          .expectCouplingCompleted()
          .expect({14.11})
          .finalize();
    }
  }
  // B - Adaptive Geometry
  if (context.isNamed("B")) {
    if (context.isPrimary()) {
      QuickTest(p, "MB"_mesh, "DA"_read, "DB"_write)
          .setVertices({1.0, y, 2.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({01.00, 02.00})
          .advance()

          .expectReadCheckpoint()
          .expect({01.00, 02.00})
          .write({01.01, 02.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({01.01, 02.01})
          .resetMesh()
          .setVertices({1.0, y})
          .write({01.10})
          .advance()

          .expectReadCheckpoint()
          .expect({01.10})
          .write({01.11})
          .advance()

          .expectCouplingCompleted()
          .expect({01.11})
          .finalize();
    } else {
      QuickTest(p, "MB"_mesh, "DA"_read, "DB"_write)
          .setVertices({3.0, y, 4.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({13.00, 14.00})
          .advance()

          .expectReadCheckpoint()
          .expect({13.00, 14.00})
          .write({13.01, 14.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({13.01, 14.01})
          .resetMesh()
          .setVertices({2.0, y, 3.0, y, 4.0, y})
          .write({12.10, 13.10, 14.10})
          .advance()

          .expectReadCheckpoint()
          .expect({02.10, 03.10, 14.10})
          .write({12.11, 13.11, 14.11})
          .advance()

          .expectCouplingCompleted()
          .expect({02.11, 03.11, 14.11})
          .finalize();
    }
  }
}

/// Changes partitioning of A from 12|34 to 34|12
inline void runSwapA(testing::TestContext &context)
{
  constexpr double y = 0.0;

  Participant p{context.name, context.config(), context.rank, context.size};

  // A - Static Geometry
  if (context.isNamed("A")) {
    if (context.isPrimary()) {
      QuickTest(p, "MA"_mesh, "DB"_read, "DA"_write)
          .setVertices({1.0, y, 2.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({01.00, 02.00})
          .advance()

          .expectReadCheckpoint()
          .expect({01.00, 02.00})
          .write({01.01, 02.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({01.01, 02.01})
          .resetMesh()
          .setVertices({3.0, y, 4.0, y})
          .write({03.10, 04.10})
          .advance()

          .expectReadCheckpoint()
          .expect({13.10, 14.10})
          .write({03.11, 04.11})
          .advance()

          .expectCouplingCompleted()
          .expect({13.11, 14.11})
          .finalize();
    } else {
      QuickTest(p, "MA"_mesh, "DB"_read, "DA"_write)
          .setVertices({3.0, y, 4.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({13.00, 14.00})
          .advance()

          .expectReadCheckpoint()
          .expect({13.00, 14.00})
          .write({13.01, 14.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({13.01, 14.01})
          .resetMesh()
          .setVertices({1.0, y, 2.0, y})
          .write({11.10, 12.10})
          .advance()

          .expectReadCheckpoint()
          .expect({01.10, 02.10})
          .write({11.11, 12.11})
          .advance()

          .expectCouplingCompleted()
          .expect({01.11, 02.11})
          .finalize();
    }
  }
  // B - Adaptive Geometry
  if (context.isNamed("B")) {
    if (context.isPrimary()) {
      QuickTest(p, "MB"_mesh, "DA"_read, "DB"_write)
          .setVertices({1.0, y, 2.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({01.00, 02.00})
          .advance()

          .expectReadCheckpoint()
          .expect({01.00, 02.00})
          .write({01.01, 02.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({01.01, 02.01})
          .write({01.10, 02.10})
          .advance()

          .expectReadCheckpoint()
          .expect({11.10, 12.10})
          .write({01.11, 02.11})
          .advance()

          .expectCouplingCompleted()
          .expect({11.11, 12.11})
          .finalize();
    } else {
      QuickTest(p, "MB"_mesh, "DA"_read, "DB"_write)
          .setVertices({3.0, y, 4.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({13.00, 14.00})
          .advance()

          .expectReadCheckpoint()
          .expect({13.00, 14.00})
          .write({13.01, 14.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({13.01, 14.01})
          .write({13.10, 14.10})
          .advance()

          .expectReadCheckpoint()
          .expect({03.10, 04.10})
          .write({13.11, 14.11})
          .advance()

          .expectCouplingCompleted()
          .expect({03.11, 04.11})
          .finalize();
    }
  }
}

/// Changes partitioning of A from 12|34 to 13|24
inline void runScatterA(testing::TestContext &context)
{
  constexpr double y = 0.0;

  Participant p{context.name, context.config(), context.rank, context.size};

  // A - Static Geometry
  if (context.isNamed("A")) {
    if (context.isPrimary()) {
      QuickTest(p, "MA"_mesh, "DB"_read, "DA"_write)
          .setVertices({1.0, y, 2.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({01.00, 02.00})
          .advance()

          .expectReadCheckpoint()
          .expect({01.00, 02.00})
          .write({01.01, 02.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({01.01, 02.01})
          .resetMesh()
          .setVertices({1.0, y, 3.0, y})
          .write({01.10, 03.10})
          .advance()

          .expectReadCheckpoint()
          .expect({01.10, 13.10})
          .write({01.11, 03.11})
          .advance()

          .expectCouplingCompleted()
          .expect({01.11, 13.11})
          .finalize();
    } else {
      QuickTest(p, "MA"_mesh, "DB"_read, "DA"_write)
          .setVertices({3.0, y, 4.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({13.00, 14.00})
          .advance()

          .expectReadCheckpoint()
          .expect({13.00, 14.00})
          .write({13.01, 14.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({13.01, 14.01})
          .resetMesh()
          .setVertices({2.0, y, 4.0, y})
          .write({12.10, 14.10})
          .advance()

          .expectReadCheckpoint()
          .expect({02.10, 14.10})
          .write({12.11, 14.11})
          .advance()

          .expectCouplingCompleted()
          .expect({02.11, 14.11})
          .finalize();
    }
  }
  // B - Adaptive Geometry
  if (context.isNamed("B")) {
    if (context.isPrimary()) {
      QuickTest(p, "MB"_mesh, "DA"_read, "DB"_write)
          .setVertices({1.0, y, 2.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({01.00, 02.00})
          .advance()

          .expectReadCheckpoint()
          .expect({01.00, 02.00})
          .write({01.01, 02.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({01.01, 02.01})
          .write({01.10, 02.10})
          .advance()

          .expectReadCheckpoint()
          .expect({01.10, 12.10})
          .write({01.11, 02.11})
          .advance()

          .expectCouplingCompleted()
          .expect({01.11, 12.11})
          .finalize();
    } else {
      QuickTest(p, "MB"_mesh, "DA"_read, "DB"_write)
          .setVertices({3.0, y, 4.0, y})
          .initialize()
          .expectWriteCheckpoint()
          .expect({00.00, 00.00})
          .write({13.00, 14.00})
          .advance()

          .expectReadCheckpoint()
          .expect({13.00, 14.00})
          .write({13.01, 14.01})
          .advance()

          .expectWriteCheckpoint()
          .expect({13.01, 14.01})
          .write({13.10, 14.10})
          .advance()

          .expectReadCheckpoint()
          .expect({03.10, 14.10})
          .write({13.11, 14.11})
          .advance()

          .expectCouplingCompleted()
          .expect({03.11, 14.11})
          .finalize();
    }
  }
}

} // namespace changepartition

namespace convergence {

inline void runResetA(testing::TestContext &context)
{
  constexpr double y = 0.0;

  BOOST_REQUIRE(context.size == 1);
  BOOST_REQUIRE(context.rank == 0);
  Participant p{context.name, context.config(), 0, 1};

  // The data format uses the following format:
  //  0 xpos . (00 | 10 | 11 )

  // A - Adaptive Geometry
  if (context.isNamed("A")) {
    QuickTest(p, "MA"_mesh, "DB"_read, "DA"_write)
        .setVertices({1.0, y, 2.0, y})
        .initialize()
        .expectWriteCheckpoint()
        .expect({00.00, 00.00})
        .write({01.00, 02.00})
        .advance()

        .expectReadCheckpoint()
        .expect({01.00, 02.00})
        .write({01.10, 02.10})
        .advance()

        .expectReadCheckpoint()
        .expect({01.10, 02.10})
        .write({01.11, 02.11})
        .advance()

        .expectWriteCheckpoint()
        .expect({01.11, 02.11})
        .resetMesh()
        .setVertices({2.0, y})
        .write({02.00})
        .advance()

        .expectReadCheckpoint()
        .expect({02.00})
        .write({02.10})
        .advance()

        .expectReadCheckpoint()
        .expect({02.10})
        .write({02.11})
        .advance()

        .expectCouplingCompleted()
        .expect({02.11})
        .finalize();
  }
  // B - Changing Geometry
  if (context.isNamed("B")) {
    QuickTest(p, "MB"_mesh, "DA"_read, "DB"_write)
        .setVertices({1.0, y, 2.0, y})
        .initialize()
        .expectWriteCheckpoint()
        .expect({00.00, 00.00})
        .write({01.00, 02.00})
        .advance()

        .expectReadCheckpoint()
        .expect({01.00, 02.00})
        .write({01.10, 02.10})
        .advance()

        .expectReadCheckpoint()
        .expect({01.10, 02.10})
        .write({01.11, 02.11})
        .advance()

        .expectWriteCheckpoint()
        .expect({01.11, 02.11})
        .write({01.00, 02.00})
        .advance()

        .expectReadCheckpoint()
        .expect({02.00, 02.00})
        .write({01.10, 02.10})
        .advance()

        .expectReadCheckpoint()
        .expect({02.10, 02.10})
        .write({01.11, 02.11})
        .advance()

        .expectCouplingCompleted()
        .expect({02.11, 02.11})
        .finalize();
  }
}

inline void runResetBoth(testing::TestContext &context)
{
  constexpr double y = 0.0;

  BOOST_REQUIRE(context.size == 1);
  BOOST_REQUIRE(context.rank == 0);
  Participant p{context.name, context.config(), 0, 1};

  // The data format uses the following format:
  //  0 xpos . (00 | 10 | 11 )

  // A - Adaptive Geometry
  if (context.isNamed("A")) {
    QuickTest(p, "MA"_mesh, "DB"_read, "DA"_write)
        .setVertices({1.0, y, 2.0, y})
        .initialize()
        .expectWriteCheckpoint()
        .expect({00.00, 00.00})
        .write({01.00, 02.00})
        .advance()

        .expectReadCheckpoint()
        .expect({01.00, 02.00})
        .write({01.10, 02.10})
        .advance()

        .expectReadCheckpoint()
        .expect({01.10, 02.10})
        .write({01.11, 02.11})
        .advance()

        .expectWriteCheckpoint()
        .expect({01.11, 02.11})
        .resetMesh()
        .setVertices({2.0, y})
        .write({02.00})
        .advance()

        .expectReadCheckpoint()
        .expect({02.00})
        .write({02.10})
        .advance()

        .expectReadCheckpoint()
        .expect({02.10})
        .write({02.11})
        .advance()

        .expectCouplingCompleted()
        .expect({02.11})
        .finalize();
  }
  // B - Changing Geometry
  if (context.isNamed("B")) {
    QuickTest(p, "MB"_mesh, "DA"_read, "DB"_write)
        .setVertices({1.0, y, 2.0, y})
        .initialize()
        .expectWriteCheckpoint()
        .expect({00.00, 00.00})
        .write({01.00, 02.00})
        .advance()

        .expectReadCheckpoint()
        .expect({01.00, 02.00})
        .write({01.10, 02.10})
        .advance()

        .expectReadCheckpoint()
        .expect({01.10, 02.10})
        .write({01.11, 02.11})
        .advance()

        .expectWriteCheckpoint()
        .expect({01.11, 02.11})
        .resetMesh()
        .setVertices({1., y, 2., y, 3., y})
        .write({01.00, 02.00, 03.00})
        .advance()

        .expectReadCheckpoint()
        .expect({02.00, 02.00, 02.00})
        .write({01.10, 02.10, 03.10})
        .advance()

        .expectReadCheckpoint()
        .expect({02.10, 02.10, 02.10})
        .write({01.11, 02.11, 03.11})
        .advance()

        .expectCouplingCompleted()
        .expect({02.11, 02.11, 02.11})
        .finalize();
  }
}

} // namespace convergence

} // namespace precice::tests::remesh::parallelImplicit
