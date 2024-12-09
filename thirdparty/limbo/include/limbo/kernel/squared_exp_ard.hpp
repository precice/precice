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
#ifndef LIMBO_KERNEL_SQUARED_EXP_ARD_HPP
#define LIMBO_KERNEL_SQUARED_EXP_ARD_HPP

#include <limbo/kernel/kernel.hpp>

namespace limbo {
    namespace defaults {
        struct kernel_squared_exp_ard {
            /// @ingroup kernel_defaults
            BO_PARAM(int, k, 0); //equivalent to the standard exp ARD
            /// @ingroup kernel_defaults
            BO_PARAM(double, sigma_sq, 1);
        };
    } // namespace defaults

    namespace kernel {
        /**
        @ingroup kernel
        \rst

        Squared exponential covariance function with automatic relevance detection (to be used with a likelihood optimizer)
        Computes the squared exponential covariance like this:

        .. math::
            k_{SE}(v1, v2) = \sigma^2 \exp \Big(-\frac{1}{2}(v1-v2)^TM(v1-v2)\Big),

        with :math:`M = \Lambda\Lambda^T + diag(l_1^{-2}, \dots, l_n^{-2})` being the characteristic length scales and :math:`\alpha` describing the variability of the latent function. The parameters :math:`l_1^2, \dots, l_n^2, \Lambda,\sigma^2` are expected in this order in the parameter array. :math:`\Lambda` is a :math:`D\times k` matrix with :math:`k<D`.

        Parameters:
           - ``double sigma_sq`` (initial signal variance)
           - ``int k`` (number of columns of :math:`\Lambda` matrix)

        Reference: :cite:`Rasmussen2006`, p. 106 & :cite:`brochu2010tutorial`, p. 10
        \endrst
        */
        template <typename Params>
        struct SquaredExpARD : public BaseKernel<Params, SquaredExpARD<Params>> {
            SquaredExpARD(int dim = 1) : _ell(dim), _A(dim, Params::kernel_squared_exp_ard::k()), _input_dim(dim)
            {
                Eigen::VectorXd p = Eigen::VectorXd::Zero(_ell.size() + _ell.size() * Params::kernel_squared_exp_ard::k() + 1);
                p(p.size() - 1) = std::log(std::sqrt(Params::kernel_squared_exp_ard::sigma_sq()));
                this->set_params(p);
            }

            size_t params_size() const { return _ell.size() + _ell.size() * Params::kernel_squared_exp_ard::k() + 1; }

            // Return the hyper parameters in log-space
            Eigen::VectorXd params() const { return _h_params; }

            // We expect the input parameters to be in log-space
            void set_params(const Eigen::VectorXd& p)
            {
                _h_params = p;
                for (size_t i = 0; i < _input_dim; ++i)
                    _ell(i) = std::exp(p(i));
                for (size_t j = 0; j < (unsigned int)Params::kernel_squared_exp_ard::k(); ++j)
                    for (size_t i = 0; i < _input_dim; ++i)
                        _A(i, j) = p((j + 1) * _input_dim + i);
                _sf2 = std::exp(2.0 * p(params_size() - 1));
            }

            Eigen::VectorXd gradient(const Eigen::VectorXd& x1, const Eigen::VectorXd& x2) const
            {
                if (Params::kernel_squared_exp_ard::k() > 0) {
                    Eigen::VectorXd grad = Eigen::VectorXd::Zero(this->params_size());
                    Eigen::MatrixXd K = (_A * _A.transpose());
                    K.diagonal() += (Eigen::MatrixXd)(_ell.array().inverse().square());
                    double z = ((x1 - x2).transpose() * K * (x1 - x2));
                    double k = _sf2 * std::exp(-0.5 * z);

                    grad.head(_input_dim) = (x1 - x2).cwiseQuotient(_ell).array().square() * k;

                    for (size_t j = 0; j < (unsigned int)Params::kernel_squared_exp_ard::k(); ++j) {
                        Eigen::VectorXd G = -((x1 - x2).transpose() * _A.col(j))(0) * (x1 - x2) * k;
                        grad.segment((j + 1) * _input_dim, _input_dim) = G;
                    }

                    grad(grad.size() - 1) = 2 * k;

                    return grad;
                }
                else {
                    Eigen::VectorXd grad(this->params_size());
                    Eigen::VectorXd z = (x1 - x2).cwiseQuotient(_ell).array().square();
                    double k = _sf2 * std::exp(-0.5 * z.sum());
                    grad.head(_input_dim) = z * k;

                    grad(grad.size() - 1) = 2 * k;
                    return grad;
                }
            }

            double kernel(const Eigen::VectorXd& x1, const Eigen::VectorXd& x2) const
            {
                assert(x1.size() == _ell.size());
                double z;
                if (Params::kernel_squared_exp_ard::k() > 0) {
                    Eigen::MatrixXd K = (_A * _A.transpose());
                    K.diagonal() += (Eigen::MatrixXd)(_ell.array().inverse().square());
                    z = ((x1 - x2).transpose() * K * (x1 - x2));
                }
                else {
                    z = (x1 - x2).cwiseQuotient(_ell).squaredNorm();
                }
                return _sf2 * std::exp(-0.5 * z);
            }

            const Eigen::VectorXd& ell() const { return _ell; }

        protected:
            double _sf2;
            Eigen::VectorXd _ell;
            Eigen::MatrixXd _A;
            size_t _input_dim;
            Eigen::VectorXd _h_params;
        };
    } // namespace kernel
} // namespace limbo

#endif
