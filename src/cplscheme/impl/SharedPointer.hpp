#pragma once

#include <memory>

namespace precice {
namespace cplscheme {
namespace impl {

class ConvergenceMeasure;
class PostProcessing;

using PtrConvergenceMeasure = std::shared_ptr<ConvergenceMeasure>;
using PtrPostProcessing     = std::shared_ptr<PostProcessing>;

}}} // namespace precice, cplscheme, impl

