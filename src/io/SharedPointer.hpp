// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IO_SHAREDPOINTER_HPP_
#define PRECICE_IO_SHAREDPOINTER_HPP_

#include "boost/smart_ptr.hpp"

namespace precice {
namespace io {

class Export;
typedef boost::shared_ptr<Export> PtrExport;

class ExportConfiguration;
typedef boost::shared_ptr<ExportConfiguration> PtrExportConfiguration;


}} // namespace precice, io

#endif /* PRECICE_IO_SHAREDPOINTER_HPP_ */
