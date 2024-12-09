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

#ifndef LIMBO_MODEL_SPGP_HPP
#define LIMBO_MODEL_SPGP_HPP

#include <cassert>
#include <iostream>
#include <limits>
#include <vector>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/LU>

#include <limbo/opt.hpp>
#include <limbo/tools.hpp>

namespace limbo {
    namespace defaults {
        struct model_spgp {
            BO_PARAM(double, jitter, 0.000001);
            BO_PARAM(double, samples_percent, 10);
            BO_PARAM(int, min_m, 1);

            /// kernel parameters
            BO_PARAM(double, sig, 0.01);
            BO_PARAM(double, pred_kernel_sigma_sq, 0.5);
            BO_PARAM(double, pred_kernel_l, 0.5);
        };
    }

    namespace model {
        /// @ingroup model
        /// Sparse Gaussian process.
        /// Reference and credits for this code:
        /// Sparse Gaussian Processes using Pseudo-inputs (NIPS, 2005)
        /// Edward Snelson and Zoubin Ghahramani
        /// It is parametrized by:
        /// - a mean function
        /// - [optionnal] an optimizer for the hyper-parameters
        template <typename Params, typename KernelFunction, typename MeanFunction, class HyperParamsOptimizer = opt::NLOptGrad<Params, nlopt::LD_LBFGS>>
        class SPGP {
        public:
            struct HyperParams {
                Eigen::MatrixXd b;
                double c;
                double sig;
                Eigen::MatrixXd xb;

                HyperParams(const Eigen::VectorXd& w, const size_t& m, const size_t& dim_in)
                {
                    b = w.segment(m * dim_in, dim_in).transpose().array().exp();
                    c = std::exp(w[(m + 1) * dim_in]);
                    sig = std::exp(w[(m + 1) * dim_in + 1]);
                    xb = w.segment(0, m * dim_in);
                    xb.resize(m, dim_in);
                }
            };

            /// useful because the model might be created before knowing anything about the process
            SPGP() : _dim_in(-1), _dim_out(-1) {}

            /// useful because the model might be created before having samples
            SPGP(int dim_in, int dim_out) : _dim_in(dim_in), _dim_out(dim_out)
            {
                _mean_function = MeanFunction(_dim_out);
                _kernel_function = KernelFunction(_dim_in);
                _b = Eigen::VectorXd::Ones(_dim_in) * Params::model_spgp::pred_kernel_l();
                _c = Params::model_spgp::pred_kernel_sigma_sq();
                _sig = Params::model_spgp::sig();
            }

            /// useful to construct without optimizing
            SPGP(const std::vector<Eigen::VectorXd>& samples,
                const std::vector<Eigen::VectorXd>& observations)
            {
                _init(samples, observations);
            }

            /// execute the hyperparameters optimization
            void optimize_hyperparams()
            {
                _optimize_init = true;
                _optimize_hyperparams();
            }

            /// Compute the SPGP from samples and observations. This call needs to be explicit!
            void compute(const std::vector<Eigen::VectorXd>& samples,
                const std::vector<Eigen::VectorXd>& observations)
            {
                compute(_to_matrix(samples), _to_matrix(observations));
            }
            void compute(const std::vector<Eigen::VectorXd>& samples,
                const Eigen::MatrixXd& observations)
            {
                compute(_to_matrix(samples), observations);
            }
            void compute(const Eigen::MatrixXd& samples,
                const Eigen::MatrixXd& observations)
            {
                assert(samples.rows() != 0);
                assert(observations.rows() != 0);
                assert(samples.rows() == observations.rows());

                _optimize_init = true;
                _init(samples, observations);
                _compute();
            }

            /// add sample and recompute the SPGP
            void add_sample(const Eigen::VectorXd& sample, const Eigen::VectorXd& observation)
            {
                if (_samples.rows() == 0) {
                    _dim_in = sample.size();
                    _kernel_function = KernelFunction(_dim_in); // the cost of building a functor should be relatively low

                    _dim_out = observation.size();
                    _mean_function = MeanFunction(_dim_out); // the cost of building a functor should be relatively low

                    // default values for the kernel
                    _b = Eigen::VectorXd::Ones(_dim_in) * Params::model_spgp::pred_kernel_l();
                    _c = Params::model_spgp::pred_kernel_sigma_sq();
                    _sig = Params::model_spgp::sig();
                }
                else {
                    assert(sample.size() == _dim_in);
                    assert(observation.size() == _dim_out);
                }

                // add the new sample
                _samples.conservativeResize(_samples.rows() + 1, _dim_in);
                _samples.row(_samples.rows() - 1) = sample.transpose();
                _update_m();

                // add the new observation and update
                _observations.conservativeResize(_observations.rows() + 1, _dim_out);
                _observations.row(_observations.rows() - 1) = observation.transpose();
                _compute_observations_zm(); // TODO: if the mean function doesn't use _obs_mean we can do a partial calculation

                _optimize_init = true; // maybe we can add one random pseudo-input before optimizing to avoid re-initialization
                _compute();
            }

            /**
             \\rst
             return :math:`\mu`, :math:`\sigma^2` (unormalized). If there is no sample, return the value according to the mean function. Using this method instead of separate calls to mu() and sigma() is more efficient because some computations are shared between mu() and sigma().
             \\endrst
	  		*/
            std::tuple<Eigen::VectorXd, double> query(const Eigen::VectorXd& v) const
            {
                std::pair<Eigen::MatrixXd, Eigen::MatrixXd>&& result = predict(v.transpose());
                return std::make_tuple(result.first, result.second(0, 0));
            }

            /**
             \\rst
             return :math:`\mu`, :math:`\sigma^2`. Predict a bunch of points.
             \\endrst
            */
            std::pair<Eigen::MatrixXd, Eigen::MatrixXd> predict(const Eigen::MatrixXd& xt) const
            {
                return _predict(xt, true, true);
            }

            /**
             \\rst
             return :math:`\mu` (unormalized). If there is no sample, return the value according to the mean function.
             \\endrst
	  		*/
            Eigen::MatrixXd mu(const Eigen::MatrixXd& v) const
            {
                return _predict(v, true, false).first;
            }
            std::vector<Eigen::VectorXd> mu_mult(const Eigen::MatrixXd& v) const
            {
                return _to_vector(_predict(v, true, false).first);
            }

            /**
             \\rst
             return :math:`\sigma^2` (unormalized). If there is no sample, return the max :math:`\sigma^2`.
             \\endrst
	  		*/
            double sigma(const Eigen::VectorXd& v) const
            {
                return (_predict(v, false, true).second)(0, 0);
            }
            // TODO: Deberia machear la interfaz de la GP
            Eigen::VectorXd sigma_mult(const Eigen::MatrixXd& v) const
            {
                return (_predict(v, false, true).second);
            }

            /// return the number of dimensions of the input
            int dim_in() const
            {
                assert(_dim_in != -1); // need to compute first !
                return _dim_in;
            }

            /// return the number of dimensions of the output
            int dim_out() const
            {
                assert(_dim_out != -1); // need to compute first !
                return _dim_out;
            }

            // NOTE: This is worth keeping?
            // const KernelFunction& kernel_function() const { return _kernel_function; }
            // KernelFunction& kernel_function() { return _kernel_function; }

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

            /// return the mean observation
            Eigen::VectorXd mean_observation() const
            {
                return _samples.size() > 0 ? _obs_mean
                                           : Eigen::VectorXd::Zero(_dim_out);
            }

            /// return the number of samples used to compute the SPGP
            int nb_samples() const { return _samples.rows(); }

            /// return the number of pseudo samples used to compute the SPGP
            int nb_pseudo_samples() const { return _pseudo_samples.rows(); }

            ///  recomputes the SPGP
            void recompute(bool update_obs_mean = true)
            {
                _optimize_init = true;
                _compute();
            }

            /// return the list of samples that have been tested so far
            std::vector<Eigen::VectorXd> samples() const { return _to_vector(_samples); }

            /// return the list of pseudo-samples beign used
            std::vector<Eigen::VectorXd> pseudo_samples() const { return _to_vector(_pseudo_samples); }

            // test the likelihood function given a dataset
            void test_likelihood(const Eigen::MatrixXd& data,
                const Eigen::MatrixXd& samples,
                const Eigen::MatrixXd& observations)
            {
                size_t m = samples.rows() / 10;
                _init(samples, observations);
                int size_w = (m + 1) * samples.cols() + 2;
                for (int i = 0; i < data.rows(); i++) {
                    Eigen::VectorXd w = data.row(i).segment(0, size_w);
                    double fw = data.row(i).segment(size_w, 1)(0);
                    Eigen::VectorXd dfw = data.row(i).segment(size_w + 1, size_w);

                    opt::eval_t result = _likelihood(w, true);
                    double rfw = std::get<0>(result);
                    Eigen::VectorXd rdfw = std::get<1>(result).get();

                    std::cout << "[" << i + 1 << "] "
                              << "fw: " << std::abs(rfw + fw)
                              << " dfw: " << (rdfw - dfw).norm() << std::endl;
                }
            }

        protected:
            /// set on initialization
            size_t _m;
            size_t _dim_in = -1;
            size_t _dim_out = -1;
            double _del;
            Eigen::MatrixXd _samples;
            Eigen::MatrixXd _observations;
            Eigen::MatrixXd _observations_zm;
            Eigen::VectorXd _obs_mean;
            MeanFunction _mean_function;
            KernelFunction _kernel_function;

            /// set after the hyperparameters calculation
            Eigen::MatrixXd _pseudo_samples;
            Eigen::VectorXd _b; // lengthscales
            double _c; // noise
            double _sig; // signal variance
            // Eigen::VectorXd _noises;

            /// calculated in _compute
            Eigen::MatrixXd _Lm;
            Eigen::MatrixXd _bet;
            Eigen::MatrixXd _matrixL;

            /// auxiliary variables
            bool _optimize_init = true;
            bool _optimized = false;
            Eigen::VectorXd _w_init;
            HyperParamsOptimizer _hp_optimize;

            void _init(const std::vector<Eigen::VectorXd>& samples, const std::vector<Eigen::VectorXd>& observations)
            {
                _init(_to_matrix(samples), _to_matrix(observations));
            }
            void _init(const Eigen::MatrixXd& samples, const Eigen::MatrixXd& observations)
            {
                _samples = samples;
                _observations = observations;

                _dim_in = _samples.cols();
                _dim_out = _observations.cols();
                _mean_function = MeanFunction(_dim_out);
                _kernel_function = KernelFunction(_dim_in);

                _compute_observations_zm();
                _update_m();

                _del = Params::model_spgp::jitter();

                _optimize_init = true;
                srand(time(NULL));
            }

            void _compute_observations_zm()
            {
                _obs_mean = _observations.colwise().mean();
                _observations_zm.resize(_observations.rows(), _observations.cols());
                for (int i = 0; i < _observations.rows(); i++) {
                    _observations_zm.row(i) = _observations.row(i) - (_mean_function(_samples.row(i), *this)).transpose();
                }
            }

            void _update_m()
            {
                // NOTE: We could add the option to use a function to change the m dinamycally based on the samples.
                _m = Params::model_spgp::samples_percent() * _samples.rows() / 100;
                if (_m < Params::model_spgp::min_m())
                    _m = Params::model_spgp::min_m();
            }

            void _compute(bool optimize = true)
            {
                if (optimize)
                    _optimize_hyperparams();

                Eigen::MatrixXd kernel_m = _compute_kernel_matrix(_pseudo_samples, _pseudo_samples) + Eigen::MatrixXd::Identity(_pseudo_samples.rows(), _pseudo_samples.rows()) * Params::model_spgp::jitter();
                _matrixL = kernel_m.llt().matrixL();
                Eigen::MatrixXd kernel_mn = _compute_kernel_matrix(_pseudo_samples, _samples);

                Eigen::MatrixXd V = (_matrixL.template triangularView<Eigen::Lower>()).solve(kernel_mn);
                Eigen::MatrixXd ep = Eigen::MatrixXd::Ones(_observations_zm.rows(), 1).array() + (_k_diag(_observations_zm).array() - V.array().pow(2).colwise().sum().transpose()) / _sig;

                Eigen::MatrixXd ep_sqrt = ep.array().sqrt();
                V = V.cwiseQuotient(ep_sqrt.transpose().replicate(_m, 1));
                Eigen::MatrixXd y = _observations_zm.cwiseQuotient(ep_sqrt);

                _Lm = (_sig * Eigen::MatrixXd::Identity(_m, _m) + V * V.transpose()).llt().matrixL();
                _bet = (_Lm.template triangularView<Eigen::Lower>()).solve(V * y);
            }

            void _optimize_hyperparams()
            {
                if (_optimize_init) {

                    // Initialize parameter vector
                    _update_m();
                    _w_init = Eigen::VectorXd((_m + 1) * _dim_in + 2);

                    // Initialize pseudo-inputs to a random subset of training inputs
                    Eigen::VectorXd positions = Eigen::VectorXd::LinSpaced(_samples.rows(), 0, _samples.rows() - 1);
                    std::random_shuffle(positions.data(), positions.data() + positions.size());
                    for (size_t i = 0; i < _m; ++i)
                        _w_init.segment(i * _dim_in, _dim_in) = _samples.row(positions[i]);

                    // Initialize hyperparameters sensibly in log space
                    _w_init.segment(_m * _dim_in, _dim_in) = -2 * ((_samples.colwise().maxCoeff() - _samples.colwise().minCoeff()).array() / 2).log(); // -2*log((max(x)-min(x))'/2)
                    _w_init[(_m + 1) * _dim_in] = std::log(_observations_zm.array().pow(2).mean());
                    _w_init[(_m + 1) * _dim_in + 1] = std::log((_observations_zm.array().pow(2) / 4).mean());

                    _optimize_init = false;
                }

                auto optimize = [&](const Eigen::VectorXd& x, bool g) {
                    return this->_likelihood(x, g);
                };
                Eigen::VectorXd result = _hp_optimize(optimize, _w_init, false);

                HyperParams hp(result, _m, _dim_in);
                _pseudo_samples = hp.xb;
                _b = hp.b.transpose();
                _c = hp.c;
                _sig = hp.sig;
                // _noises = hp.sig * Eigen::VectorXd::Ones(_m);

                _optimized = true;
            }

            opt::eval_t _likelihood(const Eigen::VectorXd& input, bool eval_grad = false) const // opt wrapper
            {
                // Unzip parameters
                HyperParams hp(input, _m, _dim_in);
                return _likelihood_wp(hp.xb, hp.b, hp.c, hp.sig, _samples, _observations_zm, eval_grad, true);
            }

            opt::eval_t _likelihood_wp(Eigen::MatrixXd& xb, const Eigen::MatrixXd& b, const double& c, const double& sig,
                const Eigen::MatrixXd& samples, const Eigen::MatrixXd& observations_zm, bool eval_grad, bool inverse) const
            {
                // Prepare xb and x
                int n = samples.rows();
                Eigen::MatrixXd idm = Eigen::MatrixXd::Identity(_m, _m);
                Eigen::MatrixXd b_sqrt = b.array().sqrt();
                xb = xb.cwiseProduct(b_sqrt.replicate(_m, 1));
                Eigen::MatrixXd x = samples.cwiseProduct(b_sqrt.replicate(n, 1));

                // Construct Q
                Eigen::MatrixXd Q = xb * xb.transpose();
                Q = Eigen::MatrixXd(Q.diagonal().replicate(1, _m) + Q.diagonal().transpose().replicate(_m, 1) - 2 * Q);
                Q = (Q.array() * -0.5).exp() * c;
                Q += _del * idm;

                // NOTE: Need to be changed for the generic kernel construction
                Eigen::MatrixXd k_m = x.cwiseProduct(x).rowwise().sum().transpose().replicate(_m, 1);
                Eigen::MatrixXd k_n = xb.cwiseProduct(xb).rowwise().sum().replicate(1, n);
                Eigen::MatrixXd K = -2 * xb * x.transpose() + k_m + k_n;
                K = (K.array() * -0.5).exp() * c;

                // Resolve
                Eigen::MatrixXd L = Q.llt().matrixL();
                Eigen::TriangularView<Eigen::MatrixXd, Eigen::Lower> L_t = L.template triangularView<Eigen::Lower>();
                Eigen::MatrixXd V = L_t.solve(K);

                Eigen::MatrixXd ep = 1 + (c - V.array().pow(2).colwise().sum().transpose()) / sig; // ep = 1 + (c-sum(V.^2)')/sig;
                K = K.array() / ep.array().sqrt().transpose().replicate(_m, 1);
                V = V.array() / ep.array().sqrt().transpose().replicate(_m, 1);
                Eigen::MatrixXd y = observations_zm.array() / ep.array().sqrt();

                Eigen::MatrixXd Lm = (sig * idm + V * V.transpose()).llt().matrixL();
                Eigen::TriangularView<Eigen::MatrixXd, Eigen::Lower> Lm_t = Lm.template triangularView<Eigen::Lower>();
                Eigen::MatrixXd invLmV = Lm_t.solve(V);
                Eigen::MatrixXd bet = invLmV * y;

                // Calculate likelihood
                Eigen::MatrixXd fw = Lm.diagonal().array().log().colwise().sum() + (n - _m) / 2 * std::log(sig) + (y.transpose() * y - bet.transpose() * bet).array() / (2 * sig) + ep.array().log().sum() / 2 + 0.5 * n * std::log(2 * M_PI);

                if (!eval_grad) {
                    if (inverse)
                        return opt::no_grad(-fw(0));
                    else
                        return opt::no_grad(fw(0));
                }

                // ** Derivates calculation - precomputations
                Eigen::MatrixXd Lt = L * Lm;
                Eigen::MatrixXd B1 = Lt.transpose().lu().solve(invLmV); // Lt is psd but Lt' it's not. We must use lu here. NOTE: Maybe it can be changed to solve the system with psd.
                Eigen::MatrixXd b1 = Lt.transpose().lu().solve(bet);

                Eigen::MatrixXd invLV = L.transpose().lu().solve(V);
                Eigen::MatrixXd invL = L.inverse();
                Eigen::MatrixXd invQ = invL.transpose() * invL; // delete invL
                Eigen::MatrixXd invLt = Lt.inverse();
                Eigen::MatrixXd invA = invLt.transpose() * invLt; // delete invLt

                //NOTE: Throws exception after solving a transpose if multiplied directly with V, that's why its beign casted.
                Eigen::MatrixXd mu = (((Eigen::MatrixXd)Lm.transpose().lu().solve(bet).transpose()) * V).transpose(); // NOTE: Check if the cast is correct here
                Eigen::MatrixXd sumVsq = V.array().pow(2).colwise().sum().transpose();

                //bigsum = y.*(bet'*invLmV)'/sig - sum(invLmV.*invLmV)'/2 - (y.^2+mu.^2)/2/sig + 0.5
                Eigen::MatrixXd bigsum = (((Eigen::MatrixXd)y.cwiseProduct((bet.transpose() * invLmV).transpose()) / sig) - ((Eigen::MatrixXd)invLmV.array().pow(2).colwise().sum().transpose() / 2) - ((Eigen::MatrixXd)(y.array().pow(2) + mu.array().pow(2)) / (2 * sig))).array() + 0.5;

                Eigen::MatrixXd TT = invLV * (invLV.transpose().cwiseProduct(bigsum.replicate(1, _m)));

                // Pseudo inputs and lengthscales
                Eigen::MatrixXd dfxb(_m, _dim_in);
                Eigen::MatrixXd dfb(_dim_in, 1);
                for (size_t i = 0; i < _dim_in; i++) {
                    // dnnQ = dist(xb(:,i),xb(:,i)).*Q;
                    Eigen::MatrixXd r = xb.block(0, i, xb.rows(), 1);

                    Eigen::MatrixXd dnnQ = _dist(xb.block(0, i, xb.rows(), 1), xb.block(0, i, xb.rows(), 1)).cwiseProduct(Q);
                    // dNnK = dist(-xb(:,i),-x(:,i)).*K;
                    Eigen::MatrixXd dNnK = _dist(xb.block(0, i, xb.rows(), 1) * -1, x.block(0, i, x.rows(), 1) * -1).cwiseProduct(K);

                    // epdot = (dNnK.*invLV)*(-2/sig); epPmod = -sum(epdot)';
                    Eigen::MatrixXd epdot = dNnK.cwiseProduct(invLV) * (-2 / sig);
                    Eigen::MatrixXd epPmod = -1 * epdot.colwise().sum().transpose();

                    // dfxb(:,i) = - b1.*(dNnK*(y-mu)/sig + dnnQ*b1) + sum((invQ - invA*sig).*dnnQ,2) + epdot*bigsum - 2/sig*sum(dnnQ.*TT,2);
                    dfxb.block(0, i, dfxb.rows(), 1) = b1.cwiseProduct((dNnK * (y - mu)) / sig + dnnQ * b1) * -1 + (invQ - invA * sig).cwiseProduct(dnnQ).rowwise().sum() + epdot * bigsum - (dnnQ.cwiseProduct(TT).rowwise().sum()) * (2 / sig);

                    // dfb(i,1) = (((y-mu)'.*(b1'*dNnK))/sig + (epPmod.*bigsum)')*x(:,i);
                    dfb.block(i, 0, 1, 1) = ((y - mu).transpose().cwiseProduct(b1.transpose() * dNnK) / sig + epPmod.cwiseProduct(bigsum).transpose()) * x.block(0, i, x.rows(), 1);

                    // Overwrite dNnK and write derivates
                    dNnK = Eigen::MatrixXd(dNnK.cwiseProduct(B1));
                    dfxb.block(0, i, dfxb.rows(), 1) += dNnK.rowwise().sum(); // dfxb(:,i) = dfxb(:,i) + sum(dNnK,2);
                    dfb.block(i, 0, 1, 1) -= dNnK.colwise().sum() * x.block(0, i, x.rows(), 1); // dfb(i,1) = dfb(i,1) - sum(dNnK,1)*x(:,i);

                    dfxb.block(0, i, dfxb.rows(), 1) *= std::sqrt(b(i)); // dfxb(:,i) = dfxb(:,i)*sqrt(b(i));

                    dfb.block(i, 0, 1, 1) /= std::sqrt(b(i)); // dfb(i,1) = dfb(i,1)/sqrt(b(i));
                    dfb.block(i, 0, 1, 1) += (dfxb.block(0, i, dfxb.rows(), 1).transpose() * // dfb(i,1) = dfb(i,1) + (dfxb(:,i)'*xb(:,i))/b(i);
                                                 xb.block(0, i, xb.rows(), 1))
                        / b(i);
                    dfb.block(i, 0, 1, 1) *= std::sqrt(b(i)) / 2; // dfb(i,1) = dfb(i,1)*sqrt(b(i))/2;
                }

                // Size
                // epc = (c./ep - sumVsq - del*sum((invLV).^2)')/sig
                Eigen::MatrixXd epc = (ep.cwiseInverse() * c - sumVsq - _del * (((Eigen::MatrixXd)invLV.array().pow(2)).colwise().sum().transpose())) / sig;

                // dfc = (m + del*trace(invQ-sig*invA) - sig*sum(sum(invA.*Q')))/2 -
                //          mu'*(y-mu)/sig + (b1'*(Q-del*eye(m))*b1)/2 + epc'*bigsum;
                Eigen::MatrixXd dfc = (_m + _del * (invQ - sig * invA).trace() - sig * (invA.cwiseProduct(Q.transpose()).sum())) / 2 - ((mu.transpose() * (y - mu)) / sig).array() + ((b1.transpose() * (Q - _del * idm) * b1) / 2).array() + (epc.transpose() * bigsum).array();

                // Noise - dfsig = sum(bigsum./ep);
                Eigen::MatrixXd dfsig = bigsum.cwiseQuotient(ep).colwise().sum();

                // Assemble the result
                Eigen::VectorXd dfw(dfxb.size() + dfb.size() + dfc.size() + dfsig.size());
                dfw.segment(0, dfxb.size()) = Eigen::Map<Eigen::MatrixXd>(dfxb.data(), 1, dfxb.size()).transpose();
                dfw.segment(dfxb.size(), dfb.size()) = dfb;
                dfw.segment(dfxb.size() + dfb.size(), dfc.size()) = dfc;
                dfw.segment(dfxb.size() + dfb.size() + dfc.size(), dfsig.size()) = dfsig;

                if (inverse) {
                    Eigen::VectorXd rdfw = dfw.array() * -1;
                    return {-fw(0), rdfw}; // limbo maximizes instead of minimizing
                }
                else {
                    return {fw(0), dfw};
                }
            }

            std::pair<Eigen::MatrixXd, Eigen::MatrixXd> _predict(const Eigen::MatrixXd& xt, bool calc_mu = true, bool calc_s2 = true) const
            {
                Eigen::MatrixXd mu(xt.rows(), _dim_out);
                Eigen::MatrixXd s2(xt.rows(), 1);

                if (_samples.rows() == 0) {
                    for (size_t i = 0; i < (size_t)xt.rows(); ++i) {
                        if (calc_mu)
                            mu.row(i) = _mean_function(xt.row(i), *this);
                        if (calc_s2)
                            s2(i, 0) = _kernel_function(xt.row(i), xt.row(i));
                    }
                    return {mu, s2.array()};
                }

                Eigen::MatrixXd K = _compute_kernel_matrix(_pseudo_samples, xt);
                Eigen::MatrixXd lst = (_matrixL.template triangularView<Eigen::Lower>()).solve(K);
                Eigen::MatrixXd lmst = (_Lm.template triangularView<Eigen::Lower>()).solve(lst);

                if (calc_mu) {
                    for (size_t i = 0; i < (size_t)xt.rows(); ++i)
                        mu.row(i) = _mean_function(xt.row(i), *this);
                    mu += (_bet.transpose() * lmst).transpose();
                }

                if (calc_s2)
                    s2 = _k_diag(xt).array() - lst.array().pow(2).colwise().sum().transpose() + _sig * lmst.array().pow(2).colwise().sum().transpose() + (_optimized ? _sig : 0);
                return {mu, s2.array()};
            }

            Eigen::MatrixXd _compute_kernel_matrix(const Eigen::MatrixXd& p_1, const Eigen::MatrixXd& p_2) const
            {
                Eigen::MatrixXd K;
                int n_1 = p_1.rows();
                int n_2 = p_2.rows();

                Eigen::MatrixXd bx_1 = _b.transpose().array().sqrt().replicate(n_1, 1);
                Eigen::MatrixXd bx_2 = _b.transpose().array().sqrt().replicate(n_2, 1);

                Eigen::MatrixXd x_1 = p_1.cwiseProduct(bx_1);
                Eigen::MatrixXd x_2 = p_2.cwiseProduct(bx_2);

                K = (-2 * x_1 * x_2.transpose()).array() + x_2.array().pow(2).rowwise().sum().transpose().replicate(n_1, 1) + x_1.array().pow(2).rowwise().sum().replicate(1, n_2);
                K = _c * (K * -0.5).array().exp();

                return K;
            }

            Eigen::MatrixXd _k_diag(const Eigen::MatrixXd& x) const
            {
                // NOTE: This calculation it's only for the ExpArd kernel
                return Eigen::MatrixXd::Constant(x.rows(), 1, _c);
            }

            // Compute pairwise distance matrix from two column vectors x0 and x1
            Eigen::MatrixXd _dist(const Eigen::MatrixXd& x_0, const Eigen::MatrixXd& x_1) const
            {
                return x_0.replicate(1, x_1.size()) - x_1.transpose().replicate(x_0.size(), 1);
            }

            // Construct a matrixXd from a std::vector of eigen::vectorXds
            // NOTE: Maybe this could be done better?
            Eigen::MatrixXd _to_matrix(const std::vector<Eigen::VectorXd>& xs) const
            {
                Eigen::MatrixXd result(xs.size(), xs[0].size());
                for (size_t i = 0; i < (size_t)result.rows(); ++i) {
                    result.row(i) = xs[i];
                }
                return result;
            }
            Eigen::MatrixXd _to_matrix(std::vector<Eigen::VectorXd>& xs) const { return _to_matrix(xs); }

            // Construct a std::vector of eigen::vectorXds from a matrixXd
            // NOTE: Maybe this could be done better?
            std::vector<Eigen::VectorXd> _to_vector(const Eigen::MatrixXd& m) const
            {
                std::vector<Eigen::VectorXd> result(m.rows());
                for (size_t i = 0; i < result.size(); ++i) {
                    result[i] = m.row(i);
                }
                return result;
            }
            std::vector<Eigen::VectorXd> _to_vector(Eigen::MatrixXd& m) const { return _to_vector(m); }

            /*
            void _compute_kernel_m()
            {
                _kernel_m.resize(m, m);

                // O(m^2) [should be negligible]
                for (int i = 0; i < m; i++)
                    for (int j = 0; j <= i; ++j)
                        _kernel_m(i, j) = _kernel_function(_pseudo_samples.row(i), _pseudo_samples.row(j)) + ((i == j) ? _noises[i] : 0); // noise only on the diagonal

                for (int i = 0; i < m; i++)
                    for (int j = 0; j < i; ++j)
                        _kernel_m(j, i) = _kernel_m(i, j);

                // O(m^3)
                _matrixL = Eigen::LLT<Eigen::MatrixXd>(_kernel_m).matrixL();
            }

            void _compute_kernel_mn()
            {
                int n = _samples.size();
                _kernel_mn.resize(m, n);

                // O(m*n)
                for (int i = 0; i < m; i++)
                    for (int j = 0; j < n; ++j)
                        _kernel_mn(i, j) = _kernel_function(_pseudo_samples.row(i), _samples.row(j));
            }

            Eigen::MatrixXd _compute_kernel(Eigen::MatrixXd& x_1, Eigen::MatrixXd& x_2)
            {
                int m = x_1.size();
                int n = x_2.size();
                Eigen::MatrixXd kernel(m, n);

                for (int i = 0; i < m; i++)
                    for (int j = 0; j < n; ++j)
                        kernel(i, j) = _kernel_function(x_1.row(i), x_2.row(j));

                return kernel;
            }

            double _kernel_function(const Eigen::VectorXd& v1, const Eigen::VectorXd& v2) const
            {
                // SquaredExpARD kernel
                const Eigen::VectorXd _ell = _b;
                double _sf2 = _c;
                double z = (v1 - v2).cwiseQuotient(_ell).squaredNorm();
                return _sf2 * std::exp(-0.5 * z);
            }
            */
        };
    }
}

#endif
