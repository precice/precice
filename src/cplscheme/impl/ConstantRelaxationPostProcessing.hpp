// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_CONSTANTRELAXATIONPOSTPROCESSING_HPP_
#define PRECICE_CPLSCHEME_CONSTANTRELAXATIONPOSTPROCESSING_HPP_

#include "PostProcessing.hpp"
#include "tarch/logging/Log.h"
#include <map>

namespace precice {
namespace cplscheme {
namespace impl {

class ConstantRelaxationPostProcessing : public PostProcessing
{
public:

   ConstantRelaxationPostProcessing (
      double relaxation,
      int    dataID );

   virtual ~ConstantRelaxationPostProcessing() {}

   virtual int getDataID () const
   {
      return _dataID;
   }

   virtual void initialize ( DataMap & cplData );

   virtual void performPostProcessing ( DataMap & cplData );

   virtual void iterationsConverged ( DataMap & cplData )
   {}

private:

   static tarch::logging::Log _log;

   double _relaxation;

   int _dataID;
};

}}} // namespace precice, cplscheme, impl

#endif /* PRECICE_CPLSCHEME_CONSTANTRELAXATIONPOSTPROCESSING_HPP_ */
