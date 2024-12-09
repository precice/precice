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
#ifndef LIMBO_OPT_GRADIENT_ASCENT_HPP
#define LIMBO_OPT_GRADIENT_ASCENT_HPP

#include <algorithm>

#include <Eigen/Core>

#include <limbo/opt/optimizer.hpp>
#include <limbo/tools/macros.hpp>
#include <limbo/tools/math.hpp>

namespace limbo {
    namespace defaults {
        struct opt_gradient_ascent {
            /// @ingroup opt_defaults
            /// number of max iterations
            BO_PARAM(int, iterations, 300);

            /// @ingroup opt_defaults
            /// alpha - learning rate
            BO_PARAM(double, alpha, 0.001);

            /// @ingroup opt_defaults
            /// gamma - for momentum
            BO_PARAM(double, gamma, 0.0);

            /// @ingroup opt_defaults
            /// nesterov momentum; turn on/off
            BO_PARAM(bool, nesterov, false);

            /// @ingroup opt_defaults
            /// norm epsilon for stopping
            BO_PARAM(double, eps_stop, 0.0);
        };
    } // namespace defaults
    namespace opt {
        /// @ingroup opt
        /// Gradient Ascent with or without momentum (Nesterov or simple)
        /// Equations from: http://ruder.io/optimizing-gradient-descent/index.html#gradientdescentoptimizationalgorithms
        /// (I changed a bit the notation; η to α)
        ///
        /// Parameters:
        /// - int iterations
        /// - double alpha
        /// - double gamma
        /// - bool nesterov
        /// - double eps_stop
        template <typename Params>
        struct GradientAscent {
            template <typename F>
            Eigen::VectorXd operator()(const F& f, const Eigen::VectorXd& init, bool bounded) const
            {
                assert(Params::opt_gradient_ascent::gamma() >= 0. && Params::opt_gradient_ascent::gamma() < 1.);
                assert(Params::opt_gradient_ascent::alpha() >= 0.);

                size_t param_dim = init.size();
                double gamma = Params::opt_gradient_ascent::gamma();
                double alpha = Params::opt_gradient_ascent::alpha();
                double stop = Params::opt_gradient_ascent::eps_stop();
                bool is_nesterov = Params::opt_gradient_ascent::nesterov();

                Eigen::VectorXd v = Eigen::VectorXd::Zero(param_dim);

                Eigen::VectorXd params = init;

                if (bounded) {
                    for (int j = 0; j < params.size(); j++) {
                        if (params(j) < 0)
                            params(j) = 0;
                        if (params(j) > 1)
                            params(j) = 1;
                    }
                }

                for (int i = 0; i < Params::opt_gradient_ascent::iterations(); ++i) {
                    Eigen::VectorXd prev_params = params;
                    Eigen::VectorXd query_params = params;
                    // if Nesterov momentum, change query parameters
                    if (is_nesterov) {
                        query_params.array() += gamma * v.array();

                        // make sure that the parameters are still in bounds, if needed
                        if (bounded) {
                            for (int j = 0; j < query_params.size(); j++) {
                                if (query_params(j) < 0)
                                    query_params(j) = 0;
                                if (query_params(j) > 1)
                                    query_params(j) = 1;
                            }
                        }
                    }
                    auto perf = opt::eval_grad(f, query_params);

                    Eigen::VectorXd grad = opt::grad(perf);
                    v = gamma * v.array() + alpha * grad.array();

                    params.array() += v.array();

                    if (bounded) {
                        for (int j = 0; j < params.size(); j++) {
                            if (params(j) < 0)
                                params(j) = 0;
                            if (params(j) > 1)
                                params(j) = 1;
                        }
                    }

                    if ((prev_params - params).norm() < stop)
                        break;
                }

                return params;
            }
        };
    } // namespace opt
} // namespace limbo

#endif
