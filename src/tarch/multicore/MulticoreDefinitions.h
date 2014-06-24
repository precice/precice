// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_MULTICORE_DEFINITIONS_H_
#define _TARCH_MULTICORE_DEFINITIONS_H_

#if defined(SharedOMP) || defined(SharedTBB) || defined(SharedCobra)
  #define SharedMemoryParallelisation
#endif

#endif
