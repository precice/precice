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
#ifndef LIMBO_OPT_ADAM_HPP
#define LIMBO_OPT_ADAM_HPP

#include <algorithm>

#include <Eigen/Core>

#include <limbo/opt/optimizer.hpp>
#include <limbo/tools/macros.hpp>
#include <limbo/tools/math.hpp>

namespace limbo {
    namespace defaults {
        struct opt_adam {
            /// @ingroup opt_defaults
            /// number of max iterations
            BO_PARAM(int, iterations, 300);

            /// @ingroup opt_defaults
            /// alpha - learning rate
            BO_PARAM(double, alpha, 0.001);

            /// @ingroup opt_defaults
            /// β1
            BO_PARAM(double, b1, 0.9);

            /// @ingroup opt_defaults
            /// β2
            BO_PARAM(double, b2, 0.999);

            /// @ingroup opt_defaults
            /// norm epsilon for stopping
            BO_PARAM(double, eps_stop, 0.0);
        };
    } // namespace defaults
    namespace opt {
        /// @ingroup opt
        /// Adam optimizer
        /// Equations from: http://ruder.io/optimizing-gradient-descent/index.html#gradientdescentoptimizationalgorithms
        /// (I changed a bit the notation; η to α)
        ///
        /// Parameters:
        /// - int iterations
        /// - double alpha
        /// - double b1
        /// - double b2
        /// - double eps_stop
        template <typename Params>
        struct Adam {
            template <typename F>
            Eigen::VectorXd operator()(const F& f, const Eigen::VectorXd& init, bool bounded) const
            {
                assert(Params::opt_adam::b1() >= 0. && Params::opt_adam::b1() < 1.);
                assert(Params::opt_adam::b2() >= 0. && Params::opt_adam::b2() < 1.);
                assert(Params::opt_adam::alpha() >= 0.);

                size_t param_dim = init.size();
                double b1 = Params::opt_adam::b1();
                double b2 = Params::opt_adam::b2();
                double b1_t = b1;
                double b2_t = b2;
                double alpha = Params::opt_adam::alpha();
                double stop = Params::opt_adam::eps_stop();
                double epsilon = 1e-8;

                Eigen::VectorXd m = Eigen::VectorXd::Zero(param_dim);
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

                for (int i = 0; i < Params::opt_adam::iterations(); ++i) {
                    Eigen::VectorXd prev_params = params;
                    auto perf = opt::eval_grad(f, params);

                    Eigen::VectorXd grad = opt::grad(perf);
                    m = b1 * m.array() + (1. - b1) * grad.array();
                    v = b2 * v.array() + (1. - b2) * grad.array().square();

                    Eigen::VectorXd m_hat = m.array() / (1. - b1_t);
                    Eigen::VectorXd v_hat = v.array() / (1. - b2_t);

                    params.array() += alpha * m_hat.array() / (v_hat.array().sqrt() + epsilon);

                    b1_t *= b1;
                    b2_t *= b2;

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
