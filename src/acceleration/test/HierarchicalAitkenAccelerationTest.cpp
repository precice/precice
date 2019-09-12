#include <Eigen/Core>
#include "cplscheme/CouplingData.hpp"
#include "acceleration/HierarchicalAitkenAcceleration.hpp"
#include "acceleration/Acceleration.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "testing/Testing.hpp"
#include "utils/EigenHelperFunctions.hpp"

BOOST_AUTO_TEST_SUITE(AccelerationTests)

using namespace precice;
using namespace cplscheme;
using namespace acceleration;

BOOST_AUTO_TEST_CASE(HierarchicalAitkenAccelerationTest)
{
  Acceleration::DataMap       dataMap;
  int                           dataID = 0;
  std::vector<int>              dataIDs;
  dataIDs.push_back(dataID);

  Eigen::VectorXd highF(5), midF(5), lowF(5), values(5), temp(5);
  highF << 0.0, -1.0, 0.0, 1.0, 0.0;
  midF << 1.0, 1.0, 0.0, -1.0, -1.0;
  lowF << 1.0, 1.0, 1.0, 1.0, 1.0;
  values = highF;

  temp = midF;
  temp *= 2.0;
  values += temp;
  temp = lowF;
  temp *= 4.0;
  values += temp;
  bool            initializeValues = false;
  mesh::PtrMesh   dummyMesh(new mesh::Mesh("DummyMesh", 3, false));
  PtrCouplingData ptrCplData = PtrCouplingData(new CouplingData(
      &values, dummyMesh, initializeValues, 1));
  temp                       = values;
  temp *= 2.0;
  utils::appendFront(ptrCplData->oldValues, temp);
  dataMap.insert(std::make_pair(dataID, ptrCplData));

  double                                 initRelaxation = 1.0;
  HierarchicalAitkenAcceleration hierarchAitken(initRelaxation, dataIDs);
  hierarchAitken.initialize(dataMap);
  hierarchAitken.performAcceleration(dataMap);

  dataMap[dataID]->oldValues.col(0) = *dataMap[dataID]->values;
  temp                              = highF;
  temp *= 0.5;
  *dataMap[dataID]->values -= temp;
  temp = midF;
  temp *= 0.1;
  *dataMap[dataID]->values -= temp;
  temp = lowF;
  temp *= 1.5;
  *dataMap[dataID]->values -= temp;
  hierarchAitken.performAcceleration(dataMap);
  //TODO BOOST_TEST missing. -> https://github.com/precice/precice/issues/88
}

BOOST_AUTO_TEST_SUITE_END()
