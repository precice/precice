// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_IMPL__SHAREDPOINTER_HPP_
#define PRECICE_CPLSCHEME_IMPL__SHAREDPOINTER_HPP_

#include "boost/smart_ptr.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

class ConvergenceMeasure;
typedef boost::shared_ptr<ConvergenceMeasure> PtrConvergenceMeasure;

class PostProcessing;
typedef boost::shared_ptr<PostProcessing> PtrPostProcessing;


}}} // namespace precice, cplscheme, impl

#endif /* PRECICE_CPLSCHEME_IMPL__SHAREDPOINTER_HPP_ */
