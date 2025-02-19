#pragma once

#include <precice/Participant.hpp>
#include <string>
#include <vector>
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

namespace precice::testing {

struct QuickTest {

  struct Mesh {
    std::string name;
  };

  struct ReadData {
    std::string name;
  };

  struct WriteData {
    std::string name;
  };

  QuickTest(Participant &p, Mesh m, ReadData r)
      : interface(p), meshName(m.name), readDataName(r.name)
  {
  }

  QuickTest(Participant &p, Mesh m, WriteData w)
      : interface(p), meshName(m.name), writeDataName(w.name)
  {
  }

  QuickTest(Participant &p, Mesh m, ReadData r, WriteData w)
      : interface(p), meshName(m.name), readDataName(r.name), writeDataName(w.name)
  {
  }

  QuickTest &setVertices(const std::vector<double> &pos)
  {
    auto n = pos.size() / interface.getMeshDimensions(meshName);
    vertexIDs.resize(n, -1);
    interface.setMeshVertices(meshName, pos, vertexIDs);
    return *this;
  }

  QuickTest &initialize()
  {
    interface.initialize();
    return *this;
  }

  QuickTest &resetMesh()
  {
    BOOST_TEST_MESSAGE("Remeshing");
    interface.resetMesh(meshName);
    return *this;
  }

  QuickTest &readCheckpoint()
  {
    interface.requiresReadingCheckpoint();
    return *this;
  }

  QuickTest &writeCheckpoint()
  {
    interface.requiresWritingCheckpoint();
    return *this;
  }

  void finalize()
  {
    interface.finalize();
  }

  QuickTest &advance()
  {
    interface.advance(interface.getMaxTimeStepSize());
    return *this;
  }

  QuickTest &advance(double dt)
  {
    interface.advance(dt);
    return *this;
  }

  QuickTest &write(const std::vector<double> &data)
  {
    interface.writeData(meshName, writeDataName, vertexIDs, data);
    return *this;
  }

  QuickTest &writeAll(double d)
  {
    std::vector<double> data(vertexIDs.size() * interface.getDataDimensions(meshName, writeDataName), d);
    interface.writeData(meshName, writeDataName, vertexIDs, data);
    return *this;
  }

  std::vector<double> read()
  {
    auto                n = vertexIDs.size() * interface.getDataDimensions(meshName, readDataName);
    std::vector<double> result(n, -0.0);
    interface.readData(meshName, readDataName, vertexIDs, interface.getMaxTimeStepSize(), result);
    return result;
  }

  QuickTest &expect(const std::vector<double> &expected)
  {
    auto data = read();
    BOOST_TEST(data == expected, boost::test_tools::per_element());
    return *this;
  }

  QuickTest &expectAll(double e)
  {
    auto                data = read();
    std::vector<double> expected(vertexIDs.size() * interface.getDataDimensions(meshName, readDataName), e);
    BOOST_TEST(data == expected, boost::test_tools::per_element());
    return *this;
  }

  QuickTest &expectReadCheckpoint()
  {
    BOOST_TEST(interface.requiresReadingCheckpoint());
    BOOST_TEST(!interface.requiresWritingCheckpoint());
    return *this;
  }

  QuickTest &expectWriteCheckpoint()
  {
    BOOST_TEST(!interface.requiresReadingCheckpoint());
    BOOST_TEST(interface.requiresWritingCheckpoint());
    return *this;
  }

  QuickTest &expectCouplingCompleted()
  {
    BOOST_TEST(!interface.isCouplingOngoing());
    return *this;
  }

  Participant &         interface;
  std::string           meshName;
  std::string           readDataName  = "unused";
  std::string           writeDataName = "unused";
  std::vector<VertexID> vertexIDs;
};

inline QuickTest::Mesh operator""_mesh(const char *name, std::size_t)
{
  return {name};
}

inline QuickTest::ReadData operator""_read(const char *name, std::size_t)
{
  return {name};
}

inline QuickTest::WriteData operator""_write(const char *name, std::size_t)
{
  return {name};
}

} // namespace precice::testing
