// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "PropertyContainerTest.hpp"
#include "mesh/PropertyContainer.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::mesh::tests::PropertyContainerTest)

namespace precice {
namespace mesh {
namespace tests {

tarch::logging::Log PropertyContainerTest::
  _log ( "precice::mesh::PropertyContainerTest" );

PropertyContainerTest:: PropertyContainerTest()
:
  TestCase ( "container::PropertyContainerTest" )
{}

void PropertyContainerTest:: run ()
{
  PRECICE_MASTER_ONLY {
    testMethod ( testSinglePropertyContainer );
    testMethod ( testHierarchicalPropertyContainers );
    testMethod ( testMultipleParents );
  }
}

void PropertyContainerTest:: testSinglePropertyContainer ()
{
  preciceTrace  ( "testSinglePropertyContainer()" );

  PropertyContainer propertyContainer;
  int propertyIndex = 0;
  int integerValue = 1;
  double doubleValue = 2.0;
  utils::Vector3D vectorValue (0.0, 1.0, 2.0);

  validate ( not propertyContainer.hasProperty(propertyIndex) );

  propertyContainer.setProperty (propertyIndex, integerValue);
  validate ( propertyContainer.hasProperty(propertyIndex) );
  integerValue = 0;
  validateEquals ( propertyContainer.getProperty<int>(propertyIndex), 1 );

  propertyContainer.setProperty ( propertyIndex, doubleValue );
  validateNumericalEquals ( propertyContainer.getProperty<double>(propertyIndex), 2.0 );

  propertyContainer.setProperty ( propertyIndex, vectorValue );
  for ( int dim=0; dim < 3; dim++ ) {
    validateNumericalEquals (
      propertyContainer.getProperty<utils::Vector3D>(propertyIndex)(dim),
      static_cast<double>(dim) );
  }

  validate ( propertyContainer.deleteProperty (propertyIndex) );
  validate ( not propertyContainer.hasProperty(propertyIndex) );
}

void PropertyContainerTest:: testHierarchicalPropertyContainers ()
{
  preciceTrace ( "testHierarchicalPropertyContainers" );

  PropertyContainer parent, child;
  int index = 0;

  child.addParent ( parent );
  validate (! child.hasProperty(index));
  parent.setProperty ( index, 1 );
  std::vector<int> properties;
  child.getProperties ( index, properties );
  validateEquals ( properties.size(), 1 );
  validateEquals ( properties[0], 1 );
  child.setProperty ( index, 2 );
  validateEquals ( child.getProperty<int>(index), 2 );
  child.getProperties ( index, properties );
  validateEquals ( properties.size(), 2 );
  validateEquals ( properties[1], 2 );
}

void PropertyContainerTest:: testMultipleParents ()
{
  preciceTrace ( "testMultipleParents()" );

  PropertyContainer parent0, parent1, child;

  child.addParent ( parent0 );
  child.addParent ( parent1 );
  int id = 0;
  validate ( ! child.hasProperty(id) );
  parent0.setProperty ( id, 0 );
  std::vector<int> properties;
  child.getProperties ( id, properties );
  validateEquals ( properties.size(), 1 );
  validateEquals ( properties[0], 0 );

  parent1.setProperty ( id, 1 );
  properties.clear ();
  child.getProperties ( id, properties );
  validateEquals ( properties.size(), 2 );
  validateEquals ( properties[0], 0 );
  validateEquals ( properties[1], 1 );
}

}}} // namespace precice, mesh, tests
