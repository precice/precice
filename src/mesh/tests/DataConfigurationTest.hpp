// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_MESH_TESTS_DATACONFIGURATIONTEST_HPP_
#define PRECICE_MESH_TESTS_DATACONFIGURATIONTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"
#include <string>

namespace precice {
namespace mesh {
namespace tests {

/**
 * @brief Provides tests for class precice::mesh::config::DataConfiguration.
 */
class DataConfigurationTest : public tarch::tests::TestCase
{
public:

   DataConfigurationTest ();

   virtual ~DataConfigurationTest () {};

   virtual void setUp () {}

   virtual void run ();

private:

   static logging::Logger _log;
};

}}} // namespace precice, mesh, tests

#endif /* PRECICE_MESH_TESTS_DATACONFIGURATIONTEST_HPP_ */
