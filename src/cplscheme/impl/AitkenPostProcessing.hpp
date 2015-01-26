// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_AIKTENPOSTPROCESSING_HPP_
#define PRECICE_CPLSCHEME_AIKTENPOSTPROCESSING_HPP_

#include "PostProcessing.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/logging/Log.h"
#include <map>

namespace precice {
namespace cplscheme {
namespace impl {

class AitkenPostProcessing : public PostProcessing
{
public:

   AitkenPostProcessing (
      double initialRelaxationFactor,
      std::vector<int>    dataIDs );

   virtual ~AitkenPostProcessing() {}

   virtual std::vector<int> getDataIDs() const
   {
      return _dataIDs;
   }

   virtual void initialize ( DataMap& cpldata );

   virtual void performPostProcessing ( DataMap& cpldata );

   virtual void iterationsConverged ( DataMap& cpldata );

private:

   static tarch::logging::Log _log;

   double _initialRelaxation;

   std::vector<int> _dataIDs;

   double _aitkenFactor;

   int _iterationCounter;

   utils::DynVector _residuals;
};

}}} // namespace precice, cplscheme, impl

#endif /* PRECICE_CPLSCHEME_AIKTENPOSTPROCESSING_HPP_ */
