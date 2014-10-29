// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_UTILS_MASTER_SLAVE_HPP_
#define PRECICE_UTILS_MASTER_SLAVE_HPP_

#include "tarch/logging/Log.h"
#include "com/SharedPointer.hpp"
#include "utils/Globals.hpp"
#include "utils/Dimensions.hpp"

namespace precice {
namespace utils {

/**
 * @brief Utility class for managing Master-Slave operations.
 */
class MasterSlave
{
public:

  static int _rank;
  static int _size;

  static bool _masterMode;
  static bool _slaveMode;

  static com::PtrCommunication _communication;

  static void configure(int rank, int size);

  static double l2norm(const DynVector& vec);

  static double dot(const DynVector& vec1, const DynVector& vec2);

private:

  static tarch::logging::Log _log;

};


}} // namespace precice, utils

#endif /* PRECICE_UTILS_MASTER_SLAVE_HPP_ */
