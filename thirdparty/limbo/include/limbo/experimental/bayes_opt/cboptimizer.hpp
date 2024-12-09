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
#ifndef LIMBO_BAYES_OPT_CBOPTIMIZER_HPP
#define LIMBO_BAYES_OPT_CBOPTIMIZER_HPP

#include <algorithm>
#include <iostream>
#include <iterator>

#include <boost/parameter/aux_/void.hpp>

#include <Eigen/Core>

#include <limbo/bayes_opt/bo_base.hpp>
#include <limbo/tools/macros.hpp>
#include <limbo/tools/random_generator.hpp>
#ifdef USE_NLOPT
#include <limbo/opt/nlopt_no_grad.hpp>
#elif defined USE_LIBCMAES
#include <limbo/opt/cmaes.hpp>
#else
#include <limbo/opt/grid_search.hpp>
#endif

namespace limbo {
    namespace defaults {
        struct bayes_opt_cboptimizer {
            BO_PARAM(int, hp_period, -1);
            BO_PARAM(bool, bounded, true);
        };
    }

    namespace experimental {
        BOOST_PARAMETER_TEMPLATE_KEYWORD(constraint_modelfun)

        namespace bayes_opt {

            using cboptimizer_signature = boost::parameter::parameters<boost::parameter::optional<limbo::tag::acquiopt>,
                boost::parameter::optional<limbo::tag::statsfun>,
                boost::parameter::optional<limbo::tag::initfun>,
                boost::parameter::optional<limbo::tag::acquifun>,
                boost::parameter::optional<limbo::tag::stopcrit>,
                boost::parameter::optional<limbo::tag::modelfun>,
                boost::parameter::optional<limbo::experimental::tag::constraint_modelfun>>;

            // clang-format off
        /**
        The classic Bayesian optimization algorithm.

        \\rst
        References: :cite:`brochu2010tutorial,Mockus2013`
        \\endrst

        This class takes the same template parameters as BoBase. It adds:
        \\rst
        +---------------------+------------+----------+---------------+
        |type                 |typedef     | argument | default       |
        +=====================+============+==========+===============+
        |acqui. optimizer     |acqui_opt_t | acquiopt | see below     |
        +---------------------+------------+----------+---------------+
        \\endrst

        The default value of acqui_opt_t is:
        - ``opt::Cmaes<Params>`` if libcmaes was found in `waf configure`
        - ``opt::NLOptNoGrad<Params, nlopt::GN_DIRECT_L_RAND>`` if NLOpt was found but libcmaes was not found
        - ``opt::GridSearch<Params>`` otherwise (please do not use this: the algorithm will not work at all!)
        */
        template <class Params,
          class A1 = boost::parameter::void_,
          class A2 = boost::parameter::void_,
          class A3 = boost::parameter::void_,
          class A4 = boost::parameter::void_,
          class A5 = boost::parameter::void_,
          class A6 = boost::parameter::void_,
          class A7 = boost::parameter::void_>
            // clang-format on
            class CBOptimizer : public limbo::bayes_opt::BoBase<Params, A1, A2, A3, A4, A5, A6> {
            public:
                // defaults
                struct defaults {
#ifdef USE_NLOPT
                    using acquiopt_t = limbo::opt::NLOptNoGrad<Params, nlopt::GN_DIRECT_L_RAND>;
#elif defined(USE_LIBCMAES)
                    using acquiopt_t = limbo::opt::Cmaes<Params>;
#else
#warning NO NLOpt, and NO Libcmaes: the acquisition function will be optimized by a grid search algorithm (which is usually bad). Please install at least NLOpt or libcmaes to use limbo!.
                    using acquiopt_t = limbo::opt::GridSearch<Params>;
#endif
                    using kf_t = limbo::kernel::Exp<Params>;
                    using mean_t = limbo::mean::Constant<Params>;
                    using constraint_model_t = limbo::model::GP<Params, kf_t, mean_t>;
                };
                /// link to the corresponding BoBase (useful for typedefs)
                using base_t = limbo::bayes_opt::BoBase<Params, A1, A2, A3, A4, A5, A6>;
                using model_t = typename base_t::model_t;

                using acquisition_function_t = typename base_t::acquisition_function_t;
                // extract the types
                using args = typename cboptimizer_signature::bind<A1, A2, A3, A4, A5, A6, A7>::type;
                using acqui_optimizer_t = typename boost::parameter::binding<args, limbo::tag::acquiopt, typename defaults::acquiopt_t>::type;

                using constraint_model_t = typename boost::parameter::binding<args, limbo::experimental::tag::constraint_modelfun, typename defaults::constraint_model_t>::type;

                /// The main function (run the Bayesian optimization algorithm)
                template <typename StateFunction, typename AggregatorFunction = FirstElem>
                void optimize(const StateFunction& sfun, const AggregatorFunction& afun = AggregatorFunction(), bool reset = true)
                {
                    _nb_constraints = StateFunction::nb_constraints();
                    _dim_out = StateFunction::dim_out();

                    this->_init(sfun, afun, reset);

                    if (!this->_observations.empty()) {
                        _split_observations();
                        _model.compute(this->_samples, _obs[0]);
                        if (_nb_constraints > 0)
                            _constraint_model.compute(this->_samples, _obs[1]);
                    }
                    else {
                        _model = model_t(StateFunction::dim_in(), StateFunction::dim_out());
                        if (_nb_constraints > 0)
                            _constraint_model = constraint_model_t(StateFunction::dim_in(), _nb_constraints);
                    }

                    acqui_optimizer_t acqui_optimizer;

                    while (!this->_stop(*this, afun)) {
                        acquisition_function_t acqui(_model, _constraint_model, this->_current_iteration);

                        auto acqui_optimization =
                            [&](const Eigen::VectorXd& x, bool g) { return acqui(x, afun, g); };
                        Eigen::VectorXd starting_point = tools::random_vector(StateFunction::dim_in(), Params::bayes_opt_cboptimizer::bounded());
                        Eigen::VectorXd new_sample = acqui_optimizer(acqui_optimization, starting_point, Params::bayes_opt_cboptimizer::bounded());
                        this->eval_and_add(sfun, new_sample);

                        this->_update_stats(*this, afun);

                        _model.add_sample(this->_samples.back(), _obs[0].back());
                        if (_nb_constraints > 0)
                            _constraint_model.add_sample(this->_samples.back(), _obs[1].back());

                        if (Params::bayes_opt_cboptimizer::hp_period() > 0
                            && (this->_current_iteration + 1) % Params::bayes_opt_cboptimizer::hp_period() == 0) {
                            _model.optimize_hyperparams();
                            if (_nb_constraints > 0)
                                _constraint_model.optimize_hyperparams();
                        }

                        this->_current_iteration++;
                        this->_total_iterations++;
                    }
                }

                /// return the best observation so far (i.e. max(f(x)))
                template <typename AggregatorFunction = FirstElem>
                const Eigen::VectorXd& best_observation(const AggregatorFunction& afun = AggregatorFunction()) const
                {
                    _dim_out = _model.dim_out();
                    _split_observations();

                    std::vector<Eigen::VectorXd> obs = _feasible_observations();
                    if (obs.size() > 0)
                        _obs[0] = obs;

                    auto rewards = std::vector<double>(_obs[0].size());
                    std::transform(_obs[0].begin(), _obs[0].end(), rewards.begin(), afun);
                    auto max_e = std::max_element(rewards.begin(), rewards.end());
                    return _obs[0][std::distance(rewards.begin(), max_e)];
                }

                /// return the best sample so far (i.e. the argmax(f(x)))
                template <typename AggregatorFunction = FirstElem>
                const Eigen::VectorXd& best_sample(const AggregatorFunction& afun = AggregatorFunction()) const
                {
                    _dim_out = _model.dim_out();
                    _split_observations();

                    std::vector<Eigen::VectorXd> obs = _feasible_observations();
                    if (obs.size() > 0)
                        _obs[0] = obs;

                    auto rewards = std::vector<double>(_obs[0].size());
                    std::transform(_obs[0].begin(), _obs[0].end(), rewards.begin(), afun);
                    auto max_e = std::max_element(rewards.begin(), rewards.end());
                    return this->_samples[std::distance(rewards.begin(), max_e)];
                }

                const model_t& model() const { return _model; } // returns model for the objective function
            protected:
                model_t _model;
                constraint_model_t _constraint_model;

                size_t _nb_constraints;
                mutable size_t _dim_out;
                mutable std::vector<std::vector<Eigen::VectorXd>> _obs;

                std::vector<Eigen::VectorXd> _feasible_observations() const
                {
                    std::vector<Eigen::VectorXd> feasible_obs;
                    for (size_t i = 0; i < _obs[0].size(); ++i) {
                        if (_obs[1][i].prod() > 0)
                            feasible_obs.push_back(_obs[0][i]);
                    }

                    return feasible_obs;
                }

                void _split_observations() const
                {
                    _obs.clear();
                    _obs.resize(2);

                    for (size_t i = 0; i < this->_observations.size(); ++i) {
                        assert(size_t(this->_observations[i].size()) == _dim_out + _nb_constraints);

                        Eigen::VectorXd vec_obj(_dim_out);
                        for (size_t j = 0; j < _dim_out; ++j)
                            vec_obj(j) = this->_observations[i](j);
                        _obs[0].push_back(vec_obj);

                        if (_nb_constraints > 0) {
                            Eigen::VectorXd vec_con(_nb_constraints);
                            for (int j = _dim_out, ind = 0; j < this->_observations[i].size(); ++j, ++ind)
                                vec_con(ind) = this->_observations[i](j);
                            _obs[1].push_back(vec_con);
                        }
                    }
                }
            };
        }
    }
}

#endif
