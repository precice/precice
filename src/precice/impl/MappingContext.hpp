// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IMPL_MAPPINGCONTEXT_HPP_
#define PRECICE_IMPL_MAPPINGCONTEXT_HPP_

#include "mesh/SharedPointer.hpp"
#include "mapping/SharedPointer.hpp"
#include "mapping/config/MappingConfiguration.hpp"

namespace precice {
namespace impl {

/**
 * @brief Holds a data mapping and related information.
 */
struct MappingContext
{
  // @brief Data mapping.
  mapping::PtrMapping mapping;

  // @brief id of mesh from which is mapped
  int fromMeshID;

  // @brief id of mesh to which is mapped
  int toMeshID;

  // @brief Time of execution of mapping.
  mapping::MappingConfiguration::Timing timing;

  // @brief True, if computation and mapping is done repeatedly for single values.
  //bool isIncremental;

  // @brief True, if data has been mapped already.
  bool hasMappedData;

  /**
   * @brief Constructor.
   */
  MappingContext()
  : mapping(),
    fromMeshID(-1),
    toMeshID(-1),
    timing(mapping::MappingConfiguration::INITIAL),
    //isIncremental(false),
    hasMappedData(false)
  {}
};

}} // namespace precice, impl

#endif /* PRECICE_IMPL_MAPPINGCONTEXT_HPP_ */
