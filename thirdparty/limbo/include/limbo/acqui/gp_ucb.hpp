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
#ifndef LIMBO_ACQUI_GP_UCB_HPP
#define LIMBO_ACQUI_GP_UCB_HPP

#include <Eigen/Core>

#include <limbo/opt/optimizer.hpp>
#include <limbo/tools/macros.hpp>

namespace limbo {
    namespace defaults {
        struct acqui_gpucb {
            /// @ingroup acqui_defaults
            BO_PARAM(double, delta, 0.1);
        };
    }
    namespace acqui {
        /** @ingroup acqui
        \rst
        GP-UCB (Upper Confidence Bound). See :cite:`brochu2010tutorial`, p. 15. See also: http://arxiv.org/abs/0912.3995

        .. math::
          UCB(x) = \mu(x) + \kappa \sigma(x).

        with:

        .. math::
          \kappa = \sqrt{2 \log{(\frac{n^{D/2+2}\pi^2}{3 \delta})}}

        where :math:`n` is the number of past evaluations of the objective function and :math:`D` the dimensionality of the parameters (dim_in).

        Parameters:
          - `double delta` (a small number in [0,1], e.g. 0.1)
        \endrst
        */
        template <typename Params, typename Model>
        class GP_UCB {
        public:
            GP_UCB(const Model& model, int iteration) : _model(model)
            {
                double nt = std::pow(iteration, dim_in() / 2.0 + 2.0);
                static constexpr double delta3 = Params::acqui_gpucb::delta() * 3;
                static constexpr double pi2 = M_PI * M_PI;
                _beta = std::sqrt(2.0 * std::log(nt * pi2 / delta3));
            }

            size_t dim_in() const { return _model.dim_in(); }

            size_t dim_out() const { return _model.dim_out(); }

            template <typename AggregatorFunction>
            opt::eval_t operator()(const Eigen::VectorXd& v, const AggregatorFunction& afun, bool gradient) const
            {
                assert(!gradient);
                Eigen::VectorXd mu;
                double sigma;
                std::tie(mu, sigma) = _model.query(v);
                return opt::no_grad(afun(mu) + _beta * std::sqrt(sigma));
            }

        protected:
            const Model& _model;
            double _beta;
        };
    }
}
#endif
