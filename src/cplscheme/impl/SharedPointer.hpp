#pragma once

#include <memory>

namespace precice {
namespace cplscheme {
namespace impl {

class ConvergenceMeasure;
class ParallelMatrixOperations;

using PtrConvergenceMeasure = std::shared_ptr<ConvergenceMeasure>;
} // namespace impl
} // namespace cplscheme
} // namespace precice
