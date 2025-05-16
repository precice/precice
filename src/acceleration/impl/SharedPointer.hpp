#pragma once

#include <memory>

namespace precice::acceleration::impl {
class ParallelMatrixOperations;
class Preconditioner;

using PtrParMatrixOps   = std::shared_ptr<ParallelMatrixOperations>;
using PtrPreconditioner = std::shared_ptr<Preconditioner>;
} // namespace precice::acceleration::impl
