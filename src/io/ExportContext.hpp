// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IO_EXPORTCONTEXT_HPP_
#define PRECICE_IO_EXPORTCONTEXT_HPP_

#include "io/SharedPointer.hpp"
#include <string>

namespace precice {
namespace io {

struct ExportContext
{
  // @brief Exporters performing the actual export.
  io::PtrExport exporter;

  // @brief Path to export location.
  std::string location;

  // @brief Exporting timestep interval (equals -1 when not set).
  int timestepInterval;

  // @brief Flag for synchronuous triggering of solver plots.
  bool triggerSolverPlot;

  // @brief Flag to decide whether neighbors of mapping shall be plotted.
  //bool plotNeighbors;

  // @brief Flag to decide whether normals of geometry elements shall be plotted.
  //bool plotNormals;

  bool exportSpacetree;

  // @brief If true, export is done in every iteration (also implicit).
  bool everyIteration;

  /**
   * @brief Constructor.
   */
  ExportContext()
  : exporter(),
    location(),
    timestepInterval(-1),
    triggerSolverPlot(false),
    //plotNeighbors ( false ),
    //plotNormals(true),
    exportSpacetree(false),
    everyIteration(false)
  {}
};

}} // namespace precice, impl

#endif /* PRECICE_IO_EXPORTCONTEXT_HPP_ */
