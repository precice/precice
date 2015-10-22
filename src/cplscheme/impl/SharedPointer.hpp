#pragma once

#include <memory>

namespace precice {
namespace cplscheme {
namespace impl {

class ConvergenceMeasure;
class PostProcessing;
class Preconditioner;

using PtrConvergenceMeasure = std::shared_ptr<ConvergenceMeasure>;
using PtrPostProcessing     = std::shared_ptr<PostProcessing>;
using PtrPreconditioner     = std::shared_ptr<Preconditioner>;

}}} // namespace precice, cplscheme, impl

