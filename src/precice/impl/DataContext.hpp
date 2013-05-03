// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IMPL_DATACONTEXT_HPP_
#define PRECICE_IMPL_DATACONTEXT_HPP_

#include "MappingContext.hpp"
#include "mesh/SharedPointer.hpp"
#include "mapping/SharedPointer.hpp"

namespace precice {
namespace impl {

/**
 * @brief Stores one Data object with related context.
 */
struct DataContext
{
  bool used;

  mesh::PtrData data;

  mesh::PtrData localData;

  mesh::PtrMesh mesh;

  MappingContext mappingContext;

  DataContext():
    used(false),
    data(),
    localData(),
    mesh(),
    mappingContext()
  {}
};

}} // namespace precice, impl

#endif /* PRECICE_IMPL_DATACONTEXT_HPP_ */
