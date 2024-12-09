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
#ifndef LIMBO_BAYES_OPT_IMGPO_HPP
#define LIMBO_BAYES_OPT_IMGPO_HPP

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iterator>

#include <boost/parameter/aux_/void.hpp>

#include <Eigen/Core>

#include <limbo/bayes_opt/bo_base.hpp>
#include <limbo/tools/macros.hpp>
#include <limbo/tools/random_generator.hpp>

namespace limbo {
    namespace defaults {
        struct bayes_opt_imgpo {
            BO_PARAM(bool, hp_opt, false);
        };
    }

    namespace bayes_opt {
        namespace experimental {
            struct TreeNode {
                std::vector<Eigen::VectorXd> x_max, x_min, x, f;
                std::vector<bool> leaf, samp;
            };
            // clang-format off
        template <class Params,
          class A1 = boost::parameter::void_,
          class A2 = boost::parameter::void_,
          class A3 = boost::parameter::void_,
          class A4 = boost::parameter::void_,
          class A5 = boost::parameter::void_>
            // clang-format on
            // Bayesian Optimization with Exponential Convergence (NIPS 2015 paper)
            // Kenji Kawaguchi, Leslie Pack Kaelbling, Tomás Lozano-Pérez
            // http://papers.nips.cc/paper/5715-bayesian-optimization-with-exponential-convergence
            class IMGPO : public limbo::bayes_opt::BoBase<Params, A1, A2, A3, A4, A5> {
            public:
                using base_t = BoBase<Params, A1, A2, A3, A4, A5>;
                using model_t = typename base_t::model_t;
                using acquisition_function_t = typename base_t::acquisition_function_t;

                template <typename StateFunction, typename AggregatorFunction = FirstElem>
                void optimize(const StateFunction& sfun, const AggregatorFunction& afun = AggregatorFunction(), bool reset = true)
                {
                    this->_init(sfun, afun, reset);
                    size_t h_upper = 1000;
                    // Init tree
                    if (this->_total_iterations == 0)
                        _init_tree(h_upper);

                    // Init model
                    _model = model_t(StateFunction::dim_in(), StateFunction::dim_out());

                    // Init root
                    _tree[0].x_max.push_back(Eigen::VectorXd::Ones(StateFunction::dim_in()));
                    _tree[0].x_min.push_back(Eigen::VectorXd::Zero(StateFunction::dim_in()));
                    _tree[0].x.push_back(Eigen::VectorXd::Ones(StateFunction::dim_in()) * 0.5);
                    _tree[0].f.push_back(sfun(_tree[0].x[0]));
                    _tree[0].leaf.push_back(true);
                    _tree[0].samp.push_back(true);

                    this->add_new_sample(_tree[0].x[0], _tree[0].f[0]);
                    double LB = afun(_tree[0].f[0]);
                    double inf = std::numeric_limits<double>::infinity();

                    size_t depth_T = 0, M = 1;
                    double rho_avg = 0, rho_bar = 0;
                    size_t t = 0, XI_max = 4, split_n = 1;
                    int xi_max = 0;
                    double XI = 1;
                    double LB_old = LB;
                    // N = N+1
                    this->_current_iteration++;
                    this->_total_iterations++;

                    while (this->_samples.size() == 0 || !this->_stop(*this, afun)) {
                        std::vector<int> i_max(depth_T + 1, -1);
                        std::vector<double> b_max(depth_T + 1, -inf);
                        size_t h_max = depth_T + 1;
                        double b_hi_max = -inf;
                        t = t + 1;

                        // Steps (i)-(ii)
                        for (size_t h = 0; h <= depth_T; h++) {
                            if (h >= h_max)
                                break;

                            while (true) {
                                for (size_t i = 0; i < _tree[h].x.size(); i++) {
                                    if (_tree[h].leaf[i] == true) {
                                        double b_hi = afun(_tree[h].f[i]);
                                        if (b_hi > b_hi_max) {
                                            b_hi_max = b_hi;
                                            i_max[h] = i;
                                            b_max[h] = b_hi;
                                        }
                                    }
                                }
                                if (i_max[h] == -1)
                                    break;

                                if (_tree[h].samp[i_max[h]] == true) {
                                    break;
                                }
                                else {
                                    Eigen::VectorXd xxx = _tree[h].x[i_max[h]];

                                    auto tmp_sample = sfun(xxx);
                                    _tree[h].samp[i_max[h]] = true;
                                    this->add_new_sample(xxx, tmp_sample);

                                    // N = N+1
                                    this->_current_iteration++;
                                    this->_total_iterations++;
                                }
                            }
                        }

                        // Steps (iii)
                        for (size_t h = 0; h <= depth_T; h++) {
                            if (h >= h_max)
                                break;
                            if (i_max[h] != -1) {
                                int xi = -1;
                                for (size_t h2 = h + 1; h2 <= std::min(depth_T, h + std::min((size_t)std::ceil(XI), XI_max)); h2++) {
                                    if (i_max[h2] != -1) {
                                        xi = h2 - h;
                                        break;
                                    }
                                }

                                double z_max = -inf;
                                size_t M2 = M;
                                if (xi != -1) {
                                    std::vector<TreeNode> tmp_tree = std::vector<TreeNode>(h + xi + 1);
                                    tmp_tree[h].x_max.push_back(_tree[h].x_max[i_max[h]]);
                                    tmp_tree[h].x_min.push_back(_tree[h].x_min[i_max[h]]);
                                    tmp_tree[h].x.push_back(_tree[h].x[i_max[h]]);

                                    M2 = M;
                                    for (size_t h2 = h; h2 <= h + xi - 1; h2++) {
                                        for (size_t ii = 0; ii < std::pow(3, h2 - h); ii++) {
                                            if (ii >= tmp_tree[h].x.size())
                                                break;
                                            Eigen::VectorXd xx = tmp_tree[h].x[ii];
                                            Eigen::VectorXd to_split = tmp_tree[h2].x_max[ii].array() - tmp_tree[h2].x_min[ii].array();
                                            size_t tmp, splitd;
                                            to_split.maxCoeff(&splitd, &tmp);
                                            auto x_g = xx, x_d = xx;
                                            x_g(splitd) = (5 * tmp_tree[h2].x_min[ii](splitd) + tmp_tree[h2].x_max[ii](splitd)) / 6.0;
                                            x_d(splitd) = (tmp_tree[h2].x_min[ii](splitd) + 5 * tmp_tree[h2].x_max[ii](splitd)) / 6.0;

                                            // TO-DO: Properly handle bl_samples etc
                                            _model.compute(this->_samples, this->_observations);
                                            acquisition_function_t acqui_g(_model, M2);
                                            z_max = std::max(z_max, acqui_g(x_g, afun));
                                            M2++;

                                            _model.compute(this->_samples, this->_observations);
                                            acquisition_function_t acqui_d(_model, M2);
                                            z_max = std::max(z_max, acqui_d(x_d, afun));
                                            M2++;

                                            if (z_max >= b_max[h + xi])
                                                break;

                                            tmp_tree[h2 + 1].x.push_back(x_g);
                                            Eigen::VectorXd newmin = tmp_tree[h2].x_min[ii];
                                            tmp_tree[h2 + 1].x_min.push_back(newmin);
                                            Eigen::VectorXd newmax = tmp_tree[h2].x_max[ii];
                                            newmax(splitd) = (2 * tmp_tree[h2].x_min[ii](splitd) + tmp_tree[h2].x_max[ii](splitd)) / 3.0;
                                            tmp_tree[h2 + 1].x_max.push_back(newmax);

                                            tmp_tree[h2 + 1].x.push_back(x_d);
                                            Eigen::VectorXd newmax2 = tmp_tree[h2].x_max[ii];
                                            tmp_tree[h2 + 1].x_max.push_back(newmax2);
                                            Eigen::VectorXd newmin2 = tmp_tree[h2].x_min[ii];
                                            newmin2(splitd) = (tmp_tree[h2].x_min[ii](splitd) + 2 * tmp_tree[h2].x_max[ii](splitd)) / 3.0;
                                            tmp_tree[h2 + 1].x_min.push_back(newmin2);

                                            tmp_tree[h2 + 1].x.push_back(xx);
                                            Eigen::VectorXd newmin3 = tmp_tree[h2].x_min[ii];
                                            newmin3(splitd) = (2 * tmp_tree[h2].x_min[ii](splitd) + tmp_tree[h2].x_max[ii](splitd)) / 3.0;
                                            tmp_tree[h2 + 1].x_min.push_back(newmin3);
                                            Eigen::VectorXd newmax3 = tmp_tree[h2].x_max[ii];
                                            newmax3(splitd) = (tmp_tree[h2].x_min[ii](splitd) + 2 * tmp_tree[h2].x_max[ii](splitd)) / 3.0;
                                            tmp_tree[h2 + 1].x_max.push_back(newmax3);
                                        }
                                        if (z_max >= b_max[h + xi])
                                            break;
                                    }
                                }

                                if (xi != -1 && z_max < b_max[h + xi]) {
                                    M = M2;
                                    i_max[h] = -1;
                                    xi_max = std::max(xi, xi_max);
                                }
                            }
                        }

                        // Steps (iv)-(v)
                        double b_hi_max_2 = -inf, rho_t = 0.0;
                        for (size_t h = 0; h <= depth_T; h++) {
                            if (h >= h_max)
                                break;
                            if (i_max[h] != -1 && b_max[h] > b_hi_max_2) {
                                rho_t += 1.0;
                                depth_T = std::max(depth_T, h + 1);
                                split_n++;
                                _tree[h].leaf[i_max[h]] = 0;

                                Eigen::VectorXd xx = _tree[h].x[i_max[h]];
                                Eigen::VectorXd to_split = _tree[h].x_max[i_max[h]].array() - _tree[h].x_min[i_max[h]].array();
                                size_t tmp, splitd;
                                to_split.maxCoeff(&splitd, &tmp);
                                auto x_g = xx, x_d = xx;
                                x_g(splitd) = (5 * _tree[h].x_min[i_max[h]](splitd) + _tree[h].x_max[i_max[h]](splitd)) / 6.0;
                                x_d(splitd) = (_tree[h].x_min[i_max[h]](splitd) + 5 * _tree[h].x_max[i_max[h]](splitd)) / 6.0;

                                // left node
                                _tree[h + 1].x.push_back(x_g);
                                // TO-DO: Properly handle bl_samples etc
                                _model.compute(this->_samples, this->_observations);
                                acquisition_function_t acqui_g(_model, M);
                                double UCB = acqui_g(x_g, afun);
                                Eigen::VectorXd fsample_g;
                                if ((UCB - LB) < 1e-6) {
                                    Eigen::VectorXd mu;
                                    double sigma;
                                    std::tie(mu, sigma) = _model.query(x_g);
                                    // TO-DO: we should fix this somehow for general acquis
                                    double var = (std::sqrt(2.0 * std::log(std::pow(M_PI, 2.0) * std::pow(M, 2.0) / (12.0 * 0.05))) + 0.2) * std::sqrt(sigma);
                                    M++;
                                    fsample_g = mu.array() + var;
                                    _tree[h + 1].samp.push_back(false);
                                }
                                else {
                                    fsample_g = sfun(x_g);
                                    _tree[h + 1].samp.push_back(true);

                                    this->add_new_sample(x_g, fsample_g);

                                    b_hi_max_2 = std::max(b_hi_max_2, fsample_g(0));

                                    // N = N+1
                                    this->_current_iteration++;
                                    this->_total_iterations++;
                                }

                                _tree[h + 1].f.push_back(fsample_g);

                                Eigen::VectorXd newmin = _tree[h].x_min[i_max[h]];
                                _tree[h + 1].x_min.push_back(newmin);
                                Eigen::VectorXd newmax = _tree[h].x_max[i_max[h]];
                                newmax(splitd) = (2 * _tree[h].x_min[i_max[h]](splitd) + _tree[h].x_max[i_max[h]](splitd)) / 3.0;
                                _tree[h + 1].x_max.push_back(newmax);
                                _tree[h + 1].leaf.push_back(true);

                                // right node
                                _tree[h + 1].x.push_back(x_d);
                                // TO-DO: Properly handle bl_samples etc
                                _model.compute(this->_samples, this->_observations);
                                acquisition_function_t acqui_d(_model, M);
                                double UCB2 = acqui_d(x_d, afun);
                                Eigen::VectorXd fsample_d;
                                if ((UCB2 - LB) < 1e-6) {
                                    Eigen::VectorXd mu;
                                    double sigma;
                                    std::tie(mu, sigma) = _model.query(x_d);
                                    // TO-DO: we should fix this somehow for general acquis
                                    double var = (std::sqrt(2.0 * std::log(std::pow(M_PI, 2.0) * std::pow(M, 2.0) / (12.0 * 0.05))) + 0.2) * std::sqrt(sigma);
                                    M++;
                                    fsample_d = mu.array() + var;
                                    _tree[h + 1].samp.push_back(false);
                                }
                                else {
                                    fsample_d = sfun(x_d);
                                    _tree[h + 1].samp.push_back(true);

                                    this->add_new_sample(x_d, fsample_d);

                                    b_hi_max_2 = std::max(b_hi_max_2, fsample_d(0));

                                    // N = N+1
                                    this->_current_iteration++;
                                    this->_total_iterations++;
                                }

                                _tree[h + 1].f.push_back(fsample_d);

                                Eigen::VectorXd newmax2 = _tree[h].x_max[i_max[h]];
                                _tree[h + 1].x_max.push_back(newmax2);
                                Eigen::VectorXd newmin2 = _tree[h].x_min[i_max[h]];
                                newmin2(splitd) = (_tree[h].x_min[i_max[h]](splitd) + 2 * _tree[h].x_max[i_max[h]](splitd)) / 3.0;
                                _tree[h + 1].x_min.push_back(newmin2);
                                _tree[h + 1].leaf.push_back(true);

                                // central node
                                _tree[h + 1].x.push_back(xx);
                                _tree[h + 1].f.push_back(_tree[h].f[i_max[h]]);
                                _tree[h + 1].samp.push_back(true);
                                Eigen::VectorXd newmin3 = _tree[h].x_min[i_max[h]];
                                Eigen::VectorXd newmax3 = _tree[h].x_max[i_max[h]];
                                newmin3(splitd) = (2 * _tree[h].x_min[i_max[h]](splitd) + _tree[h].x_max[i_max[h]](splitd)) / 3.0;
                                newmax3(splitd) = (_tree[h].x_min[i_max[h]](splitd) + 2 * _tree[h].x_max[i_max[h]](splitd)) / 3.0;
                                _tree[h + 1].x_min.push_back(newmin3);
                                _tree[h + 1].x_max.push_back(newmax3);
                                _tree[h + 1].leaf.push_back(true);

                                // update LB
                                LB = afun(this->best_observation());

                                // std::cout << t << " N= " << this->_current_iteration << " n= " << split_n << " fmax_hat= " << LB << " rho= " << rho_bar << " xi= " << xi_max << " h= " << depth_T << std::endl;
                            }
                        }

                        // Finalizing iteration
                        rho_avg = (rho_avg * (t - 1) + rho_t) / t;
                        rho_bar = std::max(rho_bar, rho_avg);
                        // update XI
                        if (std::abs(LB_old - LB) < 1e-6)
                            XI = std::max(XI - 0.5, 1.0);
                        else
                            XI = XI + 4.0;
                        LB_old = LB;
                    }

                    if (Params::bayes_opt_imgpo::hp_opt())
                        _model.optimize_hyperparams();
                }

                template <typename AggregatorFunction = FirstElem>
                const Eigen::VectorXd& best_observation(const AggregatorFunction& afun = AggregatorFunction()) const
                {
                    auto rewards = std::vector<double>(this->_observations.size());
                    std::transform(this->_observations.begin(), this->_observations.end(), rewards.begin(), afun);
                    auto max_e = std::max_element(rewards.begin(), rewards.end());
                    return this->_observations[std::distance(rewards.begin(), max_e)];
                }

                template <typename AggregatorFunction = FirstElem>
                const Eigen::VectorXd& best_sample(const AggregatorFunction& afun = AggregatorFunction()) const
                {
                    auto rewards = std::vector<double>(this->_observations.size());
                    std::transform(this->_observations.begin(), this->_observations.end(), rewards.begin(), afun);
                    auto max_e = std::max_element(rewards.begin(), rewards.end());
                    return this->_samples[std::distance(rewards.begin(), max_e)];
                }

                const model_t& model() const { return _model; }

            protected:
                void _init_tree(size_t h_max = 1000)
                {
                    _tree.clear();
                    _tree = std::vector<TreeNode>(h_max);
                }

                model_t _model;
                std::vector<TreeNode> _tree;
            };
        }
    }
}

#endif
