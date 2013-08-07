// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_SHAREDPOINTER_HPP_
#define PRECICE_CPLSCHEME_SHAREDPOINTER_HPP_

#include "boost/smart_ptr.hpp"


namespace precice {
namespace cplscheme {


class CouplingScheme;
typedef boost::shared_ptr<CouplingScheme> PtrCouplingScheme;

class CouplingSchemeConfiguration;
typedef boost::shared_ptr<CouplingSchemeConfiguration> PtrCouplingSchemeConfiguration;

class PostProcessingConfiguration;
typedef boost::shared_ptr<PostProcessingConfiguration> PtrPostProcessingConfiguration;

struct CouplingData;
typedef boost::shared_ptr<CouplingData> PtrCouplingData;

}} // namespace precice, cplscheme

#endif /* PRECICE_CPLSCHEME_SHAREDPOINTER_HPP_ */
