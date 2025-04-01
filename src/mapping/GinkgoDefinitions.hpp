#pragma once

#include <ginkgo/ginkgo.hpp>

// Every class uses Ginkgo's default_precision = double
// Ginkgo Data Structures
using GinkgoVector = gko::matrix::Dense<>;
using GinkgoMatrix = gko::matrix::Dense<>;
using GinkgoScalar = gko::matrix::Dense<>;
// Ginkgo Solver
using cg         = gko::solver::Cg<>;
using gmres      = gko::solver::Gmres<>;
using triangular = gko::solver::UpperTrs<>;
// Ginkgo Preconditioner
using jacobi   = gko::preconditioner::Jacobi<>;
using cholesky = gko::preconditioner::Ic<>;
