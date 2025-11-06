#pragma once

#include <mapping/MathHelper.hpp>
#include <mapping/config/MappingConfigurationTypes.hpp>
#include <mesh/Mesh.hpp>
#include "mapping/impl/RBFParameterTuner.hpp"

namespace precice::mapping {

class BisectionRBFTuner {

  mutable logging::Logger _log{"mapping::BisectionRBFTuner"};


};


} // namespace precice::mapping
