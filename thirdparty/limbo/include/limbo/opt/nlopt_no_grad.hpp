//| Copyright Inria May 2015
//| This project has received funding from the European Research Council (ERC) under
//| the European Union's Horizon 2020 research and innovation programme (grant
//| agreement No 637972) - see http://www.resibots.eu
//|
//| Contributor(s):
//|   - Jean-Baptiste Mouret (jean-baptiste.mouret@inria.fr)
//|   - Antoine Cully (antoinecully@gmail.com)
//|   - Konstantinos Chatzilygeroudis (konstantinos.chatzilygeroudis@inria.fr)
//|   - Federico Allocati (fede.allocati@gmail.com)
//|   - Vaios Papaspyros (b.papaspyros@gmail.com)
//|   - Roberto Rama (bertoski@gmail.com)
//|
//| This software is a computer library whose purpose is to optimize continuous,
//| black-box functions. It mainly implements Gaussian processes and Bayesian
//| optimization.
//| Main repository: http://github.com/resibots/limbo
//| Documentation: http://www.resibots.eu/limbo
//|
//| This software is governed by the CeCILL-C license under French law and
//| abiding by the rules of distribution of free software.  You can  use,
//| modify and/ or redistribute the software under the terms of the CeCILL-C
//| license as circulated by CEA, CNRS and INRIA at the following URL
//| "http://www.cecill.info".
//|
//| As a counterpart to the access to the source code and  rights to copy,
//| modify and redistribute granted by the license, users are provided only
//| with a limited warranty  and the software's author,  the holder of the
//| economic rights,  and the successive licensors  have only  limited
//| liability.
//|
//| In this respect, the user's attention is drawn to the risks associated
//| with loading,  using,  modifying and/or developing or reproducing the
//| software by the user in light of its specific status of free software,
//| that may mean  that it is complicated to manipulate,  and  that  also
//| therefore means  that it is reserved for developers  and  experienced
//| professionals having in-depth computer knowledge. Users are therefore
//| encouraged to load and test the software's suitability as regards their
//| requirements in conditions enabling the security of their systems and/or
//| data to be ensured and,  more generally, to use and operate it in the
//| same conditions as regards security.
//|
//| The fact that you are presently reading this means that you have had
//| knowledge of the CeCILL-C license and that you accept its terms.
//|
#ifndef LIMBO_OPT_NLOPT_NO_GRAD_HPP
#define LIMBO_OPT_NLOPT_NO_GRAD_HPP

#ifndef USE_NLOPT
#warning No NLOpt
#else
#include <limbo/opt/nlopt_base.hpp>

namespace limbo {
    namespace defaults {
        struct opt_nloptnograd {
            /// @ingroup opt_defaults
            /// number of calls to the optimized function
            BO_PARAM(int, iterations, 500);
            /// @ingroup opt_defaults
            /// tolerance for convergence: stop when an optimization step (or an
            /// estimate of the optimum) changes the objective function value by
            /// less than tol multiplied by the absolute value of the function
            /// value.
            /// IGNORED if negative
            BO_PARAM(double, fun_tolerance, -1);
            /// @ingroup opt_defaults
            /// tolerance for convergence: stop when an optimization step (or an
            /// estimate of the optimum) changes all the parameter values by
            /// less than tol multiplied by the absolute value of the parameter
            /// value.
            /// IGNORED if negative
            BO_PARAM(double, xrel_tolerance, -1);
        };
    } // namespace defaults
    namespace opt {
        /**
          @ingroup opt
        Binding to gradient-free NLOpt algorithms.
         See: http://ab-initio.mit.edu/wiki/index.php/NLopt_Algorithms

         Algorithms:
         - GN_DIRECT
         - GN_DIRECT_L, [default]
         - GN_DIRECT_L_RAND
         - GN_DIRECT_NOSCAL
         - GN_DIRECT_L_NOSCAL
         - GN_DIRECT_L_RAND_NOSCAL
         - GN_ORIG_DIRECT
         - GN_ORIG_DIRECT_L
         - GN_CRS2_LM
         - GN_MLSL
         - GN_MLSL_LDS
         - GN_ISRES
         - LN_COBYLA
         - LN_AUGLAG_EQ
         - LN_BOBYQA
         - LN_NEWUOA
         - LN_NEWUOA_BOUND
         - LN_PRAXIS
         - LN_NELDERMEAD
         - LN_SBPLX
         - LN_AUGLAG

         Parameters:
         - int iterations
         - double fun_tolerance
         - double xrel_tolerance
        */
        template <typename Params, nlopt::algorithm Algorithm = nlopt::GN_DIRECT_L_RAND>
        struct NLOptNoGrad : public NLOptBase<Params, Algorithm> {
        public:
            void initialize(int dim) override
            {
                // Assert that the algorithm is non-gradient
                // TO-DO: Add support for MLSL (Multi-Level Single-Linkage)
                // TO-DO: Add better support for ISRES (Improved Stochastic Ranking Evolution Strategy)
                // clang-format off
                static_assert(Algorithm == nlopt::LN_COBYLA || Algorithm == nlopt::LN_BOBYQA ||
                    Algorithm == nlopt::LN_NEWUOA || Algorithm == nlopt::LN_NEWUOA_BOUND ||
                    Algorithm == nlopt::LN_PRAXIS || Algorithm == nlopt::LN_NELDERMEAD ||
                    Algorithm == nlopt::LN_SBPLX || Algorithm == nlopt::GN_DIRECT ||
                    Algorithm == nlopt::GN_DIRECT_L || Algorithm == nlopt::GN_DIRECT_L_RAND ||
                    Algorithm == nlopt::GN_DIRECT_NOSCAL || Algorithm == nlopt::GN_DIRECT_L_NOSCAL ||
                    Algorithm == nlopt::GN_DIRECT_L_RAND_NOSCAL || Algorithm == nlopt::GN_ORIG_DIRECT ||
                    Algorithm == nlopt::GN_ORIG_DIRECT_L || Algorithm == nlopt::GN_CRS2_LM ||
                    Algorithm == nlopt::LN_AUGLAG || Algorithm == nlopt::LN_AUGLAG_EQ ||
                    Algorithm == nlopt::GN_ISRES || Algorithm == nlopt::GN_ESCH, "NLOptNoGrad accepts gradient free nlopt algorithms only");
                // clang-format on

                NLOptBase<Params, Algorithm>::initialize(dim);

                this->_opt.set_maxeval(Params::opt_nloptnograd::iterations());
                this->_opt.set_ftol_rel(Params::opt_nloptnograd::fun_tolerance());
                this->_opt.set_xtol_rel(Params::opt_nloptnograd::xrel_tolerance());
            }
        };
    } // namespace opt
} // namespace limbo

#endif
#endif
