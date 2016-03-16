// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IO_TESTS_ImportVRMLTest_HPP_
#define PRECICE_IO_TESTS_ImportVRMLTest_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"
#include <string>

namespace precice {
namespace io {
namespace tests {

/**
 *
 */
class ImportVRMLTest : public tarch::tests::TestCase
{
public:

  /**
   * Constructor.
   *
   * @param pathToSrc [IN] Relative path to src directory of preCICE.
   */
  ImportVRMLTest ();

  /**
   * Destructor.
   */
  virtual ~ImportVRMLTest() {};

  /**
   * @brief Retrieves path to sources.
   */
  virtual void setUp ();

  /**
   * Runs all testcases.
   */
  virtual void run ();

private:

  // @brief Logging device.
  static logging::Logger _log;

  // @brief Relative path to src directory of preCICE.
  std::string _pathToTests;

  /**
   * @brief 2D test case for importing a VRML file.
   *
   * The import of data values is tested.
   */
  void testImportSquare ();

  /**
   * @brief 3D test case for importing a VRML file.
   */
  void testImportCube ();

  /**
   * @brief 3D test case for importing a VRML file.
   */
  void testImportSphere ();

  /**
   * @brief 3D test case for importing a VRML file.
   */
  void testImportApe ();

  void testImportBunny ();

  void testImportDragon ();

  void testImportReactorPipe ();
};

}}} // namespace precice, io, tests

#endif //PRECICE_IO_TESTS_ImportVRMLTest_HPP_
