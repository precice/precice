#pragma once

#include <memory>

namespace precice::cplscheme::impl {

class ConvergenceMeasure;
class ParallelMatrixOperations;

using PtrConvergenceMeasure = std::shared_ptr<ConvergenceMeasure>;
} // namespace precice::cplscheme::impl
