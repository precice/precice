#include "../impl/PostProcessing.hpp"
#include "../impl/HierarchicalAitkenPostProcessing.hpp"
#include "../CouplingData.hpp"
#include "utils/Globals.hpp"
#include "utils/Parallel.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Mesh.hpp"
#include <Eigen/Core>
#include "utils/EigenHelperFunctions.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(CplSchemeTests)

using namespace precice;
using namespace cplscheme;

BOOST_AUTO_TEST_CASE(HierarchicalAitkenPostProcessingTest) {
  impl::PostProcessing::DataMap dataMap;
  int dataID = 0;
  std::vector<int> dataIDs;
  dataIDs.push_back(dataID);

  Eigen::VectorXd highF(5), midF(5), lowF(5), values(5), temp(5);
  highF <<  0.0, -1.0, 0.0,  1.0,  0.0;
  midF << 1.0,  1.0, 0.0, -1.0, -1.0;
  lowF << 1.0,  1.0, 1.0,  1.0,  1.0;
  values = highF;

  temp = midF; temp *= 2.0; values += temp;
  temp = lowF; temp *= 4.0; values += temp;
  bool initializeValues = false;
  mesh::PtrMesh dummyMesh ( new mesh::Mesh("dummyMesh", 3, false) );
  PtrCouplingData ptrCplData = PtrCouplingData(new CouplingData(
      &values,dummyMesh,initializeValues, 1));
  temp = values;
  temp *= 2.0;
  utils:: appendFront(ptrCplData->oldValues, temp);
  dataMap.insert ( std::make_pair(dataID, ptrCplData));

  double initRelaxation = 1.0;
  impl::HierarchicalAitkenPostProcessing hierarchAitken ( initRelaxation, dataIDs);
  hierarchAitken.initialize ( dataMap );
  hierarchAitken.performPostProcessing ( dataMap );

  dataMap[dataID]->oldValues.col(0) = *dataMap[dataID]->values;
  temp = highF; temp *= 0.5; *dataMap[dataID]->values -= temp;
  temp = midF; temp *= 0.1; *dataMap[dataID]->values -= temp;
  temp = lowF; temp *= 1.5; *dataMap[dataID]->values -= temp;
  hierarchAitken.performPostProcessing ( dataMap );
  //TODO BOOST_TEST missing. Nothing is tested in this test.
}

BOOST_AUTO_TEST_SUITE_END()
