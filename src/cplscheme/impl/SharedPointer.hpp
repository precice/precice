#pragma once

#include <memory>

namespace precice
{
namespace cplscheme
{
namespace impl
{

class ConvergenceMeasure;
class PostProcessing;
class Preconditioner;
class ParallelMatrixOperations;

using PtrConvergenceMeasure = std::shared_ptr<ConvergenceMeasure>;
using PtrPostProcessing     = std::shared_ptr<PostProcessing>;
using PtrPreconditioner     = std::shared_ptr<Preconditioner>;
using PtrParMatrixOps       = std::shared_ptr<ParallelMatrixOperations>;
}
}
} // namespace precice, cplscheme, impl
