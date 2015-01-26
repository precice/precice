// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IO_TESTS_TXTTABLEWRITER_HPP_
#define PRECICE_IO_TESTS_TXTTABLEWRITER_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace io {
namespace tests {

class TXTTableWriterTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  TXTTableWriterTest ();

  /**
   * @brief Destructor, empty.
   */
  virtual ~TXTTableWriterTest () {}

  /**
   * @brief Empty.
   */
  virtual void setUp () {}

  /**
   * @brief Calls all tests.
   */
  virtual void run ();

private:

  static tarch::logging::Log _log;

  /**
   * @brief Tests writing table data to a file.
   */
  void test ();
};

}}} // namespace precice, io, tests

#endif // PRECICE_IO_TESTS_TXTTABLEWRITER_HPP_
