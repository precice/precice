// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "HierarchicalAitkenPostProcessingTest.hpp"
#include "../impl/PostProcessing.hpp"
#include "../impl/HierarchicalAitkenPostProcessing.hpp"
#include "../CouplingData.hpp"
#include "utils/Globals.hpp"
#include "utils/Parallel.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::cplscheme::tests::HierarchicalAitkenPostProcessingTest)

namespace precice {
namespace cplscheme {
namespace tests {

tarch::logging::Log HierarchicalAitkenPostProcessingTest::
  _log ( "precice::cplscheme::tests::HierarchicalAitkenPostProcessingTest" );

HierarchicalAitkenPostProcessingTest:: HierarchicalAitkenPostProcessingTest ()
:
  TestCase ( "precice::cplscheme::tests::HierarchicalAitkenPostProcessingTest" )
{}

void HierarchicalAitkenPostProcessingTest:: run ()
{
  PRECICE_MASTER_ONLY {
    preciceTrace ( "run()" );
    impl::PostProcessing::DataMap dataMap;
    int dataID = 0;
    std::vector<int> dataIDs;
    dataIDs.push_back(dataID);

    utils::DynVector highF ( 5 );
    utils::DynVector midF ( 5 );
    utils::DynVector lowF ( 5 );
    utils::DynVector values ( 5 );
    utils::DynVector temp ( 5 );
    assignList(highF) = 0.0, -1.0, 0.0,  1.0,  0.0;
    assignList(midF)  = 1.0,  1.0, 0.0, -1.0, -1.0;
    assignList(lowF)  = 1.0,  1.0, 1.0,  1.0,  1.0;
    values = highF;
    temp = midF; temp *= 2.0; values += temp;
    temp = lowF; temp *= 4.0; values += temp;
    bool initializeValues = false;
    PtrCouplingData ptrCplData = PtrCouplingData(new CouplingData(
                  &values,initializeValues));
    temp = values;
    temp *= 2.0;
    ptrCplData->oldValues.appendFront ( temp );
    dataMap.insert ( std::make_pair(dataID, ptrCplData));

    double initRelaxation = 1.0;
    impl::HierarchicalAitkenPostProcessing hierarchAitken ( initRelaxation, dataIDs);
    hierarchAitken.initialize ( dataMap );
    hierarchAitken.performPostProcessing ( dataMap );

    dataMap[dataID]->oldValues.column(0) = *dataMap[dataID]->values;
    temp = highF; temp *= 0.5; *dataMap[dataID]->values -= temp;
    temp = midF; temp *= 0.1; *dataMap[dataID]->values -= temp;
    temp = lowF; temp *= 1.5; *dataMap[dataID]->values -= temp;
    hierarchAitken.performPostProcessing ( dataMap );
  }
}

}}} // namespace precice, cplscheme, tests
