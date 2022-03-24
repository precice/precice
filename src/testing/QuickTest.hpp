#pragma once

#include <boost/test/unit_test.hpp>
#include <cmath>
#include <precice/SolverInterface.hpp>
#include <string>
#include <vector>

namespace precice {
namespace testing {

struct QuickTest {

  /// Represents named data and its type
  struct Data {
    std::string name;
    bool        vectorial;
  };

  /// Creates quicktest based on one mesh
  QuickTest(SolverInterface &si, const std::string &mesh)
      : interface(&si), dims(si.getDimensions()), meshID(si.getMeshID(mesh))
  {
  }

  /** Sets mesh vertices based on the given coordinates while saving ids for later access.
   *
   * Will automatically deduce the amount of vertices to set.
   *
   * @see precice::SolverInterface::setMeshVertices()
   */
  QuickTest &setVertices(const std::vector<double> &pos)
  {
    int n = static_cast<int>(pos.size()) / dims;
    BOOST_REQUIRE(n > 0);
    vertexIDs.resize(n, -1);
    interface->setMeshVertices(meshID, n, pos.data(), vertexIDs.data());
    BOOST_REQUIRE(std::count(vertexIDs.begin(), vertexIDs.end(), -1) == 0);
    return *this;
  }

  /** Wrapper around initialize() ignoring the return dt.
   *
   * @see precice::SolverInterface::initialize()
   */
  QuickTest &initialize()
  {
    interface->initialize();
    return *this;
  }

  /** Wrapper around finalize.
   *
   * @see precice::SolverInterface::finalize()
   */
  void finalize()
  {
    interface->finalize();
  }

  /** Wrapper around advance(dt) ignoring the return dt.
   *
   * @see precice::SolverInterface::advance()
   */
  QuickTest &advance(double dt = 1.0)
  {
    interface->advance(dt);
    return *this;
  }

  /// Initializes data
  QuickTest &initializeData()
  {
    auto action = precice::constants::actionWriteInitialData();
    if (interface->isActionRequired(action)) {
      interface->markActionFulfilled(action);
    }
    interface->initializeData();
    return *this;
  }

  /// Writes the data to vertices defined in \ref setVertices()
  QuickTest &write(Data d, const std::vector<double> &data)
  {
    auto dataID = interface->getDataID(d.name, meshID);
    if (d.vectorial) {
      auto n = data.size() / dims;
      interface->writeBlockVectorData(dataID, n, vertexIDs.data(), data.data());
    } else {
      interface->writeBlockScalarData(dataID, data.size(), vertexIDs.data(), data.data());
    }
    return *this;
  }

  /// Reads the data from vertices defined in \ref setVertices() and returns a vector
  QuickTest &read(Data d)
  {
    auto dataID = interface->getDataID(d.name, meshID);
    auto n      = vertexIDs.size();

    std::vector<double> result;
    if (d.vectorial) {
      result.resize(n * dims, -0.0);
      interface->readBlockVectorData(dataID, n, vertexIDs.data(), result.data());
    } else {
      result.resize(n, -0.0);
      interface->readBlockScalarData(dataID, n, vertexIDs.data(), result.data());
    }
    reads.push_back(std::move(result));
    return *this;
  }

  std::vector<double> &last()
  {
    return reads.back();
  }

  SolverInterface *                interface;
  int                              dims;
  int                              meshID;
  std::vector<int>                 vertexIDs;
  std::vector<std::vector<double>> reads;
};

inline QuickTest::Data operator""_scalar(const char *name, std::size_t)
{
  return {name, false};
}

inline QuickTest::Data operator""_vector(const char *name, std::size_t)
{
  return {name, true};
}

} // namespace testing
} // namespace precice
