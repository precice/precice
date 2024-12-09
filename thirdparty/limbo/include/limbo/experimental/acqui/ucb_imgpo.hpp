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
#ifndef LIMBO_EXPERIMENTAL_ACQUI_UCB_HPP
#define LIMBO_EXPERIMENTAL_ACQUI_UCB_HPP

#include <Eigen/Core>

#include <limbo/tools/macros.hpp>

namespace limbo {
    namespace defaults {
        struct acqui_ucb_imgpo {
            BO_PARAM(double, nu, 0.05);
        };
    }
    namespace acqui {
        namespace experimental {
            template <typename Params, typename Model>
            class UCB_IMGPO {
            public:
                UCB_IMGPO(const Model& model, size_t M = 1) : _model(model), _M(M) {}

                size_t dim_in() const { return _model.dim_in(); }

                size_t dim_out() const { return _model.dim_out(); }

                template <typename AggregatorFunction>
                double operator()(const Eigen::VectorXd& v, const AggregatorFunction& afun) const
                {
                    Eigen::VectorXd mu;
                    double sigma;
                    std::tie(mu, sigma) = _model.query(v);
                    // UCB - nu = 0.05
                    // sqrt(2*log(pi^2*M^2/(12*nu)))
                    double gp_varsigma = std::sqrt(2.0 * std::log(std::pow(M_PI, 2.0) * std::pow(_M, 2.0) / (12.0 * Params::acqui_ucb_imgpo::nu())));
                    return (afun(mu) + (gp_varsigma + 0.2) * std::sqrt(sigma));
                }

            protected:
                const Model& _model;
                size_t _M;
            };
        }
    }
}

#endif
