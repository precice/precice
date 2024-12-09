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
#ifndef LIMBO_MODEL_GP_HPP
#define LIMBO_MODEL_GP_HPP

#include <cassert>
#include <iostream>
#include <limits>
#include <vector>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/LU>

// Quick hack for definition of 'I' in <complex.h>
#undef I

#include <limbo/kernel/matern_five_halves.hpp>
#include <limbo/kernel/squared_exp_ard.hpp>
#include <limbo/mean/constant.hpp>
#include <limbo/mean/data.hpp>
#include <limbo/model/gp/kernel_lf_opt.hpp>
#include <limbo/model/gp/no_lf_opt.hpp>
#include <limbo/tools.hpp>

namespace limbo {
    namespace model {
        /// @ingroup model
        /// A classic Gaussian process.
        /// It is parametrized by:
        /// - a kernel function
        /// - a mean function
        /// - [optional] an optimizer for the hyper-parameters
        template <typename Params, typename KernelFunction = kernel::MaternFiveHalves<Params>, typename MeanFunction = mean::Data<Params>, typename HyperParamsOptimizer = gp::NoLFOpt<Params>>
        class GP {
        public:
            /// useful because the model might be created before knowing anything about the process
            GP() : _dim_in(-1), _dim_out(-1), _inv_kernel_updated(false) {}

            /// useful because the model might be created before having samples
            GP(int dim_in, int dim_out)
                : _dim_in(dim_in), _dim_out(dim_out), _kernel_function(dim_in), _mean_function(dim_out), _inv_kernel_updated(false) {}

            /// Compute the GP from samples and observations. This call needs to be explicit!
            void compute(const std::vector<Eigen::VectorXd>& samples,
                const std::vector<Eigen::VectorXd>& observations, bool compute_kernel = true)
            {
                assert(samples.size() != 0);
                assert(observations.size() != 0);
                assert(samples.size() == observations.size());

                if (_dim_in != samples[0].size()) {
                    _dim_in = samples[0].size();
                    _kernel_function = KernelFunction(_dim_in); // the cost of building a functor should be relatively low
                }

                if (_dim_out != observations[0].size()) {
                    _dim_out = observations[0].size();
                    _mean_function = MeanFunction(_dim_out); // the cost of building a functor should be relatively low
                }

                _samples = samples;

                _observations.resize(observations.size(), _dim_out);
                for (int i = 0; i < _observations.rows(); ++i)
                    _observations.row(i) = observations[i];

                _mean_observation = _observations.colwise().mean();

                this->_compute_obs_mean();
                if (compute_kernel)
                    this->_compute_full_kernel();
            }

            /// Do not forget to call this if you use hyper-parameters optimization!!
            void optimize_hyperparams()
            {
                _hp_optimize(*this);
            }

            /// add sample and update the GP. This code uses an incremental implementation of the Cholesky
            /// decomposition. It is therefore much faster than a call to compute()
            void add_sample(const Eigen::VectorXd& sample, const Eigen::VectorXd& observation)
            {
                if (_samples.empty()) {
                    if (_dim_in != sample.size()) {
                        _dim_in = sample.size();
                        _kernel_function = KernelFunction(_dim_in); // the cost of building a functor should be relatively low
                    }
                    if (_dim_out != observation.size()) {
                        _dim_out = observation.size();
                        _mean_function = MeanFunction(_dim_out); // the cost of building a functor should be relatively low
                    }
                }
                else {
                    assert(sample.size() == _dim_in);
                    assert(observation.size() == _dim_out);
                }

                _samples.push_back(sample);

                _observations.conservativeResize(_observations.rows() + 1, _dim_out);
                _observations.bottomRows<1>() = observation.transpose();

                _mean_observation = _observations.colwise().mean();

                this->_compute_obs_mean();
                this->_compute_incremental_kernel();
            }

            /**
             \\rst
             return :math:`\mu`, :math:`\sigma^2` (un-normalized). If there is no sample, return the value according to the mean function. Using this method instead of separate calls to mu() and sigma() is more efficient because some computations are shared between mu() and sigma().
             \\endrst
            */
            std::tuple<Eigen::VectorXd, double> query(const Eigen::VectorXd& v) const
            {
                if (_samples.size() == 0)
                    return std::make_tuple(_mean_function(v, *this),
                        _kernel_function(v, v) + _kernel_function.noise());

                Eigen::VectorXd k = _compute_k(v);
                return std::make_tuple(_mu(v, k), _sigma(v, k) + _kernel_function.noise());
            }

            /**
             \\rst
             return :math:`\mu` (un-normalized). If there is no sample, return the value according to the mean function.
             \\endrst
            */
            Eigen::VectorXd mu(const Eigen::VectorXd& v) const
            {
                if (_samples.size() == 0)
                    return _mean_function(v, *this);
                return _mu(v, _compute_k(v));
            }

            /**
             \\rst
             return :math:`\sigma^2` (un-normalized). If there is no sample, return the max :math:`\sigma^2`.
             \\endrst
            */
            double sigma(const Eigen::VectorXd& v) const
            {
                if (_samples.size() == 0)
                    return _kernel_function(v, v) + _kernel_function.noise();
                return _sigma(v, _compute_k(v)) + _kernel_function.noise();
            }

            /// return the number of dimensions of the input
            int dim_in() const
            {
                assert(_dim_in != -1); // need to compute first!
                return _dim_in;
            }

            /// return the number of dimensions of the output
            int dim_out() const
            {
                assert(_dim_out != -1); // need to compute first!
                return _dim_out;
            }

            const KernelFunction& kernel_function() const { return _kernel_function; }

            KernelFunction& kernel_function() { return _kernel_function; }

            const MeanFunction& mean_function() const { return _mean_function; }

            MeanFunction& mean_function() { return _mean_function; }

            /// return the maximum observation (only call this if the output of the GP is of dimension 1)
            Eigen::VectorXd max_observation() const
            {
                if (_observations.cols() > 1)
                    std::cout << "WARNING max_observation with multi dimensional "
                                 "observations doesn't make sense"
                              << std::endl;
                return tools::make_vector(_observations.maxCoeff());
            }

            /// return the mean observation (only call this if the output of the GP is of dimension 1)
            Eigen::VectorXd mean_observation() const
            {
                assert(_dim_out > 0);
                return _samples.size() > 0 ? _mean_observation
                                           : Eigen::VectorXd::Zero(_dim_out);
            }

            const Eigen::MatrixXd& mean_vector() const { return _mean_vector; }

            const Eigen::MatrixXd& obs_mean() const { return _obs_mean; }

            /// return the number of samples used to compute the GP
            int nb_samples() const { return _samples.size(); }

            ///  recomputes the GP
            void recompute(bool update_obs_mean = true, bool update_full_kernel = true)
            {
                assert(!_samples.empty());

                if (update_obs_mean)
                    this->_compute_obs_mean();

                if (update_full_kernel)
                    this->_compute_full_kernel();
                else
                    this->_compute_alpha();
            }

            void compute_inv_kernel()
            {
                size_t n = _obs_mean.rows();
                // K^{-1} using Cholesky decomposition
                _inv_kernel = Eigen::MatrixXd::Identity(n, n);

                _matrixL.template triangularView<Eigen::Lower>().solveInPlace(_inv_kernel);
                _matrixL.template triangularView<Eigen::Lower>().transpose().solveInPlace(_inv_kernel);

                _inv_kernel_updated = true;
            }

            /// compute and return the log likelihood
            double compute_log_lik()
            {
                size_t n = _obs_mean.rows();

                // --- cholesky ---
                // see:
                // http://xcorr.net/2008/06/11/log-determinant-of-positive-definite-matrices-in-matlab/
                long double logdet = 2 * _matrixL.diagonal().array().log().sum();

                double a = (_obs_mean.transpose() * _alpha)
                               .trace(); // generalization for multi dimensional observation

                _log_lik = -0.5 * a - 0.5 * logdet - 0.5 * n * std::log(2 * M_PI);

                return _log_lik;
            }

            /// compute and return the gradient of the log likelihood wrt to the kernel parameters
            Eigen::VectorXd compute_kernel_grad_log_lik()
            {
                size_t n = _obs_mean.rows();

                // compute K^{-1} only if needed
                if (!_inv_kernel_updated) {
                    compute_inv_kernel();
                }
                Eigen::MatrixXd w = _inv_kernel;

                // alpha * alpha.transpose() - K^{-1}
                w = _alpha * _alpha.transpose() - w;

                // only compute half of the matrix (symmetrical matrix)
                Eigen::VectorXd grad = Eigen::VectorXd::Zero(_kernel_function.h_params_size());
                for (size_t i = 0; i < n; ++i) {
                    for (size_t j = 0; j <= i; ++j) {
                        Eigen::VectorXd g = _kernel_function.grad(_samples[i], _samples[j], i, j);
                        if (i == j)
                            grad += w(i, j) * g * 0.5;
                        else
                            grad += w(i, j) * g;
                    }
                }

                return grad;
            }

            /// compute and return the gradient of the log likelihood wrt to the mean parameters
            Eigen::VectorXd compute_mean_grad_log_lik()
            {
                size_t n = _obs_mean.rows();

                // compute K^{-1} only if needed
                if (!_inv_kernel_updated) {
                    compute_inv_kernel();
                }

                Eigen::VectorXd grad = Eigen::VectorXd::Zero(_mean_function.h_params_size());
                for (int i_obs = 0; i_obs < _dim_out; ++i_obs)
                    for (size_t n_obs = 0; n_obs < n; n_obs++) {
                        grad += _obs_mean.col(i_obs).transpose() * _inv_kernel.col(n_obs) * _mean_function.grad(_samples[n_obs], *this).row(i_obs);
                    }

                return grad;
            }

            /// return the likelihood (do not compute it -- return last computed)
            double get_log_lik() const { return _log_lik; }

            /// set the log likelihood (e.g. computed from outside)
            void set_log_lik(double log_lik) { _log_lik = log_lik; }

            /// compute and return the log probability of LOO CV
            double compute_log_loo_cv()
            {
                // compute K^{-1} only if needed
                if (!_inv_kernel_updated) {
                    compute_inv_kernel();
                }

                Eigen::VectorXd inv_diag = _inv_kernel.diagonal().array().inverse();

                _log_loo_cv = (((-0.5 * (_alpha.array().square().array().colwise() * inv_diag.array())).array().colwise() - 0.5 * inv_diag.array().log().array()) - 0.5 * std::log(2 * M_PI)).colwise().sum().sum();

                return _log_loo_cv;
            }

            /// compute and return the gradient of the log probability of LOO CV wrt to the kernel parameters
            Eigen::VectorXd compute_kernel_grad_log_loo_cv()
            {
                size_t n = _obs_mean.rows();
                size_t n_params = _kernel_function.h_params_size();

                // compute K^{-1} only if needed
                if (!_inv_kernel_updated) {
                    compute_inv_kernel();
                }

                Eigen::VectorXd grad = Eigen::VectorXd::Zero(n_params);
                Eigen::MatrixXd grads = Eigen::MatrixXd::Zero(n_params, _dim_out);

                // only compute half of the matrix (symmetrical matrix)
                // TO-DO: Make it better
                std::vector<std::vector<Eigen::VectorXd>> full_dk;
                for (size_t i = 0; i < n; i++) {
                    full_dk.push_back(std::vector<Eigen::VectorXd>());
                    for (size_t j = 0; j <= i; j++)
                        full_dk[i].push_back(_kernel_function.grad(_samples[i], _samples[j], i, j));
                    for (size_t j = i + 1; j < n; j++)
                        full_dk[i].push_back(Eigen::VectorXd::Zero(n_params));
                }
                for (size_t i = 0; i < n; i++)
                    for (size_t j = 0; j < i; ++j)
                        full_dk[j][i] = full_dk[i][j];

                Eigen::VectorXd inv_diag = _inv_kernel.diagonal().array().inverse();

                for (int j = 0; j < grad.size(); j++) {
                    Eigen::MatrixXd dKdTheta_j = Eigen::MatrixXd::Zero(n, n);
                    for (size_t i = 0; i < n; i++) {
                        for (size_t k = 0; k < n; k++)
                            dKdTheta_j(i, k) = full_dk[i][k](j);
                    }
                    Eigen::MatrixXd Zeta_j = _inv_kernel * dKdTheta_j;
                    Eigen::MatrixXd Zeta_j_alpha = Zeta_j * _alpha;
                    Eigen::MatrixXd Zeta_j_K = Zeta_j * _inv_kernel;

                    grads.row(j) = ((_alpha.array() * Zeta_j_alpha.array() - 0.5 * ((1. + _alpha.array().square().array().colwise() * inv_diag.array()).array().colwise() * Zeta_j_K.diagonal().array())).array().colwise() * inv_diag.array()).colwise().sum();

                    // for (size_t i = 0; i < n; i++)
                    //     grads.row(j).array() += (_alpha.row(i).array() * Zeta_j_alpha.row(i).array() - 0.5 * (1. + _alpha.row(i).array().square() / _inv_kernel.diagonal()(i)) * Zeta_j_K.diagonal()(i)) / _inv_kernel.diagonal()(i);
                }

                grad = grads.rowwise().sum();

                return grad;
            }

            /// return the LOO-CV log probability (do not compute it -- return last computed)
            double get_log_loo_cv() const { return _log_loo_cv; }

            /// set the LOO-CV log probability (e.g. computed from outside)
            void set_log_loo_cv(double log_loo_cv) { _log_loo_cv = log_loo_cv; }

            /// LLT matrix (from Cholesky decomposition)
            const Eigen::MatrixXd& matrixL() const { return _matrixL; }

            const Eigen::MatrixXd& alpha() const { return _alpha; }

            /// return the list of samples
            const std::vector<Eigen::VectorXd>& samples() const { return _samples; }

            /// return the list of observations
            std::vector<Eigen::VectorXd> observations() const
            {
                std::vector<Eigen::VectorXd> observations;
                for (int i = 0; i < _observations.rows(); i++) {
                    observations.push_back(_observations.row(i));
                }

                return observations;
            }

            /// return the observations (in matrix form)
            /// (NxD), where N is the number of points and D is the dimension output
            const Eigen::MatrixXd& observations_matrix() const
            {
                return _observations;
            }

            bool inv_kernel_computed() { return _inv_kernel_updated; }

            /// save the parameters and the data for the GP to the archive (text or binary)
            template <typename A>
            void save(const std::string& directory) const
            {
                A archive(directory);
                save(archive);
            }

            /// save the parameters and the data for the GP to the archive (text or binary)
            template <typename A>
            void save(const A& archive) const
            {
                if (_kernel_function.h_params_size() > 0) {
                    archive.save(_kernel_function.h_params(), "kernel_params");
                }
                if (_mean_function.h_params_size() > 0) {
                    archive.save(_mean_function.h_params(), "mean_params");
                }
                archive.save(_samples, "samples");
                archive.save(_observations, "observations");
                archive.save(_matrixL, "matrixL");
                archive.save(_alpha, "alpha");
            }

            /// load the parameters and the data for the GP from the archive (text or binary)
            /// if recompute is true, we do not read the kernel matrix
            /// but we recompute it given the data and the hyperparameters
            template <typename A>
            void load(const std::string& directory, bool recompute = true)
            {
                A archive(directory);
                load(archive, recompute);
            }

            /// load the parameters and the data for the GP from the archive (text or binary)
            /// if recompute is true, we do not read the kernel matrix
            /// but we recompute it given the data and the hyperparameters
            template <typename A>
            void load(const A& archive, bool recompute = true)
            {
                _samples.clear();
                archive.load(_samples, "samples");

                archive.load(_observations, "observations");

                _dim_in = _samples[0].size();
                _kernel_function = KernelFunction(_dim_in);

                if (_kernel_function.h_params_size() > 0) {
                    Eigen::VectorXd h_params;
                    archive.load(h_params, "kernel_params");
                    assert(h_params.size() == (int)_kernel_function.h_params_size());
                    _kernel_function.set_h_params(h_params);
                }

                _dim_out = _observations.cols();
                _mean_function = MeanFunction(_dim_out);

                if (_mean_function.h_params_size() > 0) {
                    Eigen::VectorXd h_params;
                    archive.load(h_params, "mean_params");
                    assert(h_params.size() == (int)_mean_function.h_params_size());
                    _mean_function.set_h_params(h_params);
                }

                _mean_observation = _observations.colwise().mean();

                if (recompute)
                    this->recompute(true, true);
                else {
                    archive.load(_matrixL, "matrixL");
                    archive.load(_alpha, "alpha");
                }
            }

        protected:
            int _dim_in;
            int _dim_out;

            KernelFunction _kernel_function;
            MeanFunction _mean_function;

            std::vector<Eigen::VectorXd> _samples;
            Eigen::MatrixXd _observations;
            Eigen::MatrixXd _mean_vector;
            Eigen::MatrixXd _obs_mean;

            Eigen::MatrixXd _alpha;
            Eigen::VectorXd _mean_observation;

            Eigen::MatrixXd _kernel, _inv_kernel;

            Eigen::MatrixXd _matrixL;

            double _log_lik, _log_loo_cv;
            bool _inv_kernel_updated;

            HyperParamsOptimizer _hp_optimize;

            void _compute_obs_mean()
            {
                assert(!_samples.empty());
                _mean_vector.resize(_samples.size(), _dim_out);
                for (int i = 0; i < _mean_vector.rows(); i++) {
                    assert(_samples[i].cols() == 1);
                    assert(_samples[i].rows() != 0);
                    assert(_samples[i].rows() == _dim_in);
                    _mean_vector.row(i) = _mean_function(_samples[i], *this);
                }
                _obs_mean = _observations - _mean_vector;
            }

            void _compute_full_kernel()
            {
                size_t n = _samples.size();
                _kernel.resize(n, n);

                // O(n^2) [should be negligible]
                for (size_t i = 0; i < n; i++)
                    for (size_t j = 0; j <= i; ++j)
                        _kernel(i, j) = _kernel_function(_samples[i], _samples[j], i, j);

                for (size_t i = 0; i < n; i++)
                    for (size_t j = 0; j < i; ++j)
                        _kernel(j, i) = _kernel(i, j);

                // O(n^3)
                _matrixL = Eigen::LLT<Eigen::MatrixXd>(_kernel).matrixL();

                this->_compute_alpha();

                // notify change of kernel
                _inv_kernel_updated = false;
            }

            void _compute_incremental_kernel()
            {
                // Incremental LLT
                // This part of the code is inspired from the Bayesopt Library (cholesky_add_row function).
                // However, the mathematical foundations can be easily retrieved by detailing the equations of the
                // extended L matrix that produces the desired kernel.

                size_t n = _samples.size();
                _kernel.conservativeResize(n, n);

                for (size_t i = 0; i < n; ++i) {
                    _kernel(i, n - 1) = _kernel_function(_samples[i], _samples[n - 1], i, n - 1);
                    _kernel(n - 1, i) = _kernel(i, n - 1);
                }

                _matrixL.conservativeResizeLike(Eigen::MatrixXd::Zero(n, n));

                double L_j;
                for (size_t j = 0; j < n - 1; ++j) {
                    L_j = _kernel(n - 1, j) - (_matrixL.block(j, 0, 1, j) * _matrixL.block(n - 1, 0, 1, j).transpose())(0, 0);
                    _matrixL(n - 1, j) = (L_j) / _matrixL(j, j);
                }

                L_j = _kernel(n - 1, n - 1) - (_matrixL.block(n - 1, 0, 1, n - 1) * _matrixL.block(n - 1, 0, 1, n - 1).transpose())(0, 0);
                _matrixL(n - 1, n - 1) = sqrt(L_j);

                this->_compute_alpha();

                // notify change of kernel
                _inv_kernel_updated = false;
            }

            void _compute_alpha()
            {
                // alpha = K^{-1} * this->_obs_mean;
                Eigen::TriangularView<Eigen::MatrixXd, Eigen::Lower> triang = _matrixL.template triangularView<Eigen::Lower>();
                _alpha = triang.solve(_obs_mean);
                triang.adjoint().solveInPlace(_alpha);
            }

            Eigen::VectorXd _mu(const Eigen::VectorXd& v, const Eigen::VectorXd& k) const
            {
                return (k.transpose() * _alpha) + _mean_function(v, *this).transpose();
            }

            double _sigma(const Eigen::VectorXd& v, const Eigen::VectorXd& k) const
            {
                Eigen::VectorXd z = _matrixL.triangularView<Eigen::Lower>().solve(k);
                double res = _kernel_function(v, v) - z.dot(z);

                return (res <= std::numeric_limits<double>::epsilon()) ? 0 : res;
            }

            Eigen::VectorXd _compute_k(const Eigen::VectorXd& v) const
            {
                Eigen::VectorXd k(_samples.size());
                for (int i = 0; i < k.size(); i++)
                    k[i] = _kernel_function(_samples[i], v);
                return k;
            }
        };
        /// GPBasic is a GP with a "mean data" mean function, Exponential kernel,
        /// and NO hyper-parameter optimization
        template <typename Params>
        using GPBasic = GP<Params, kernel::MaternFiveHalves<Params>, mean::Data<Params>, gp::NoLFOpt<Params>>;

        /// GPOpt is a GP with a "mean data" mean function, Exponential kernel with Automatic Relevance
        /// Determination (ARD), and hyper-parameter optimization based on Rprop
        template <typename Params>
        using GPOpt = GP<Params, kernel::SquaredExpARD<Params>, mean::Data<Params>, gp::KernelLFOpt<Params>>;
    } // namespace model
} // namespace limbo

#endif
