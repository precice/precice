#pragma once

#include <memory>

namespace precice {
namespace acceleration {
namespace impl {
class ParallelMatrixOperations;
class Preconditioner;

using PtrParMatrixOps   = std::shared_ptr<ParallelMatrixOperations>;
using PtrPreconditioner = std::shared_ptr<Preconditioner>;
} // namespace impl
} // namespace acceleration
} // namespace precice
