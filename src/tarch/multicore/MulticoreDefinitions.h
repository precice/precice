#ifndef _TARCH_MULTICORE_DEFINITIONS_H_
#define _TARCH_MULTICORE_DEFINITIONS_H_

#if defined(SharedOMP) || defined(SharedTBB) || defined(SharedCobra)
  #define SharedMemoryParallelisation
#endif

#endif
