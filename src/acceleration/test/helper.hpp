#pragma once

#include <memory>

#include "cplscheme/CouplingData.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice::testing {
inline cplscheme::PtrCouplingData makeCouplingData(mesh::PtrData data, mesh::PtrMesh mesh, bool exchangeSubsteps)
{
  bool requiresInitialization = false;
  return std::make_shared<cplscheme::CouplingData>(std::move(data), std::move(mesh), requiresInitialization, exchangeSubsteps, cplscheme::CouplingData::Direction::Send);
}
} // namespace precice::testing
