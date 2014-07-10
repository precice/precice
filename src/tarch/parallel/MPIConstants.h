// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano

#include "tarch/la/Vector.h"
#include <sstream>

#ifndef _TARCH_PARALLEL_MPI_CONSTANTS_H_
#define _TARCH_PARALLEL_MPI_CONSTANTS_H_


#include "tarch/compiler/CompilerSpecificSettings.h"

#ifdef __INTEL_COMPILER
#define CompilerICC
#endif

#if !defined(CompilerDefinesMPIMaxNameString)
#define MPI_MAX_NAME_STRING            80
#endif


#define MPI_MAX_NAME_STRING_ADDED_ONE  (MPI_MAX_NAME_STRING+1)

// Note: implementation was moved to header due to strange linker error on supermuc

namespace tarch {
  namespace parallel {
    class StringTools {
    public:
      static std::string convert(const tarch::la::Vector<MPI_MAX_NAME_STRING_ADDED_ONE,short>& value ) {
        std::ostringstream result;

        std::string::size_type i=0;
        while (i<MPI_MAX_NAME_STRING && value( static_cast<int>(i) )!=0 ) {
          result << static_cast<char>(value( static_cast<int>(i) ));
          i++;
        }

        return result.str();
      }


      static tarch::la::Vector<MPI_MAX_NAME_STRING_ADDED_ONE,short> convert( const std::string& value ) {
        tarch::la::Vector<MPI_MAX_NAME_STRING_ADDED_ONE,short> result;
        assertion( value.length() <= MPI_MAX_NAME_STRING );
        std::string::size_type i=0;
        while (i<value.length()) {
          result( static_cast<int>(i) ) = value.at( static_cast<int>(i) );
          i++;
        }
        result( static_cast<int>(i) ) = 0;
        return result;
      }
    };
}}


#endif
