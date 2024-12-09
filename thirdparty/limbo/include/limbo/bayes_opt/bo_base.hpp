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
#ifndef LIMBO_BAYES_OPT_BO_BASE_HPP
#define LIMBO_BAYES_OPT_BO_BASE_HPP

#include <exception>
#include <iostream>
#include <limits>
#include <vector>

// Quick hack for definition of 'I' in <complex.h>
#undef I
#include <boost/fusion/include/accumulate.hpp>
#include <boost/fusion/include/for_each.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/parameter.hpp>
#define BOOST_NO_SCOPED_ENUMS
#include <boost/filesystem.hpp>

#include <Eigen/Core>

// we need everything to have the defaults
#include <limbo/acqui/ucb.hpp>
#include <limbo/init/random_sampling.hpp>
#include <limbo/kernel/exp.hpp>
#include <limbo/mean/data.hpp>
#include <limbo/model/gp.hpp>
#include <limbo/stat/aggregated_observations.hpp>
#include <limbo/stat/console_summary.hpp>
#include <limbo/stat/samples.hpp>
#include <limbo/stop/chain_criteria.hpp>
#include <limbo/stop/max_iterations.hpp>
#include <limbo/tools/macros.hpp>
#include <limbo/tools/math.hpp>
#include <limbo/tools/sys.hpp>

namespace limbo {
    namespace defaults {
        struct bayes_opt_bobase {
            BO_PARAM(bool, stats_enabled, true);
            BO_PARAM(bool, bounded, true);
        };
    }
    template <typename BO, typename AggregatorFunction>
    struct RefreshStat_f {
        RefreshStat_f(BO& bo, const AggregatorFunction& afun)
            : _bo(bo), _afun(afun) {}

        BO& _bo;
        const AggregatorFunction& _afun;

        template <typename T>
        void operator()(T& x) const { x(_bo, _afun); }
    };

    struct FirstElem {
        using result_type = double;
        double operator()(const Eigen::VectorXd& x) const
        {
            return x(0);
        }
    };
    class EvaluationError : public std::exception {
    };

    // we use optimal named template parameters
    // see:
    // http://www.boost.org/doc/libs/1_55_0/libs/parameter/doc/html/index.html#parameter-enabled-class-templates

    BOOST_PARAMETER_TEMPLATE_KEYWORD(initfun)
    BOOST_PARAMETER_TEMPLATE_KEYWORD(acquifun)
    BOOST_PARAMETER_TEMPLATE_KEYWORD(modelfun)
    BOOST_PARAMETER_TEMPLATE_KEYWORD(statsfun)
    BOOST_PARAMETER_TEMPLATE_KEYWORD(stopcrit)

    namespace bayes_opt {

        using bobase_signature = boost::parameter::parameters<boost::parameter::optional<tag::statsfun>,
            boost::parameter::optional<tag::initfun>,
            boost::parameter::optional<tag::acquifun>,
            boost::parameter::optional<tag::stopcrit>,
            boost::parameter::optional<tag::modelfun>>;

        // clang-format off
        template <class Params,
          class A1 = boost::parameter::void_,
          class A2 = boost::parameter::void_,
          class A3 = boost::parameter::void_,
          class A4 = boost::parameter::void_,
          class A5 = boost::parameter::void_,
          class A6 = boost::parameter::void_>
        // clang-format on
        /**
        \rst

        Base class for Bayesian optimizers

        Parameters:
          - ``bool Params::bayes_opt_bobase::stats_enabled``: activate / deactivate the statistics

        This class is templated by several types with default values (thanks to boost::parameters).

        +----------------+---------+---------+---------------+
        |type            |typedef  | argument| default       |
        +================+=========+=========+===============+
        |init. func.     |init_t   | initfun | RandomSampling|
        +----------------+---------+---------+---------------+
        |model           |model_t  | modelfun| GP<...>       |
        +----------------+---------+---------+---------------+
        |acquisition fun.|aqui_t   | acquifun| GP_UCB        |
        +----------------+---------+---------+---------------+
        |statistics      | stat_t  | statfun | see below     |
        +----------------+---------+---------+---------------+
        |stopping crit.  | stop_t  | stopcrit| MaxIterations |
        +----------------+---------+---------+---------------+

        \endrst

        For GP, the default value is: ``model::GP<Params, kf_t, mean_t, opt_t>>``,
          - with ``kf_t = kernel::SquaredExpARD<Params>``
          - with ``mean_t = mean::Data<Params>``
          - with ``opt_t = model::gp::KernelLFOpt<Params>``

         (meaning: kernel with automatic relevance determination and mean equals to the mean of the input data, that is, center the data automatically)

        For Statistics, the default value is: ``boost::fusion::vector<stat::Samples<Params>, stat::AggregatedObservations<Params>, stat::ConsoleSummary<Params>>``

        Example of customization:
          - ``using Kernel_t = kernel::MaternFiveHalves<Params>;``
          - ``using Mean_t = mean::Data<Params>;``
          - ``using GP_t = model::GP<Params, Kernel_t, Mean_t>;``
          - ``using Acqui_t = acqui::UCB<Params, GP_t>;``
          - ``bayes_opt::BOptimizer<Params, modelfun<GP_t>, acquifun<Acqui_t>> opt;``

        */
        class BoBase {
        public:
            using params_t = Params;
            // defaults
            struct defaults {
                using init_t = init::RandomSampling<Params>; // 1
                using model_t = model::GP<Params>; // 2
                // WARNING: you have to specify the acquisition  function
                // if you use a custom model
                using acqui_t = acqui::UCB<Params, model_t>; // 3
                using stat_t = boost::fusion::vector<stat::Samples<Params>, stat::AggregatedObservations<Params>, stat::ConsoleSummary<Params>>; // 4
                using stop_t = boost::fusion::vector<stop::MaxIterations<Params>>; // 5
            };

            // extract the types
            using args = typename bobase_signature::bind<A1, A2, A3, A4, A5, A6>::type;
            using init_function_t = typename boost::parameter::binding<args, tag::initfun, typename defaults::init_t>::type;
            using acquisition_function_t = typename boost::parameter::binding<args, tag::acquifun, typename defaults::acqui_t>::type;
            using model_t = typename boost::parameter::binding<args, tag::modelfun, typename defaults::model_t>::type;
            using Stat = typename boost::parameter::binding<args, tag::statsfun, typename defaults::stat_t>::type;
            using StoppingCriteria = typename boost::parameter::binding<args, tag::stopcrit, typename defaults::stop_t>::type;

            using stopping_criteria_t = typename boost::mpl::if_<boost::fusion::traits::is_sequence<StoppingCriteria>, StoppingCriteria, boost::fusion::vector<StoppingCriteria>>::type;
            using stat_t = typename boost::mpl::if_<boost::fusion::traits::is_sequence<Stat>, Stat, boost::fusion::vector<Stat>>::type;

            /// default constructor
            BoBase() : _total_iterations(0) { _make_res_dir(); }

            /// copy is disabled (dangerous and useless)
            BoBase(const BoBase& other) = delete;
            /// copy is disabled (dangerous and useless)
            BoBase& operator=(const BoBase& other) = delete;

            /// return true if the statitics are enabled (they can be disabled to avoid dumping data, e.g. for unit tests)
            bool stats_enabled() const { return Params::bayes_opt_bobase::stats_enabled(); }

            /// return the name of the directory in which results (statistics) are written
            const std::string& res_dir() const { return _res_dir; }

            /// return the vector of points of observations (observations can be multi-dimensional, hence the VectorXd) -- f(x)
            const std::vector<Eigen::VectorXd>& observations() const { return _observations; }

            /// return the list of the points that have been evaluated so far (x)
            const std::vector<Eigen::VectorXd>& samples() const { return _samples; }

            /// return the current iteration number
            int current_iteration() const { return _current_iteration; }

            int total_iterations() const { return _total_iterations; }

            /// Add a new sample / observation pair
            /// - does not update the model!
            /// - we don't add NaN and inf observations
            void add_new_sample(const Eigen::VectorXd& s, const Eigen::VectorXd& v)
            {
                if (tools::is_nan_or_inf(v))
                    throw EvaluationError();
                _samples.push_back(s);
                _observations.push_back(v);
            }

            /// Evaluate a sample and add the result to the 'database' (sample / observations vectors) -- it does not update the model
            template <typename StateFunction>
            void eval_and_add(const StateFunction& seval, const Eigen::VectorXd& sample)
            {
                this->add_new_sample(sample, seval(sample));
            }

        protected:
            template <typename StateFunction, typename AggregatorFunction>
            void _init(const StateFunction& seval, const AggregatorFunction& afun, bool reset = true)
            {
                this->_current_iteration = 0;
                if (reset) {
                    this->_total_iterations = 0;
                    this->_samples.clear();
                    this->_observations.clear();
                }

                if (this->_total_iterations == 0)
                    init_function_t()(seval, afun, *this);
            }

            template <typename BO, typename AggregatorFunction>
            bool _stop(const BO& bo, const AggregatorFunction& afun) const
            {
                stop::ChainCriteria<BO, AggregatorFunction> chain(bo, afun);
                return boost::fusion::accumulate(_stopping_criteria, false, chain);
            }

            template <typename BO, typename AggregatorFunction>
            void _update_stats(BO& bo, const AggregatorFunction& afun)
            { // not const, because some stat class
                // modify the optimizer....
                boost::fusion::for_each(_stat, RefreshStat_f<BO, AggregatorFunction>(bo, afun));
            }

            void _make_res_dir()
            {
                if (!Params::bayes_opt_bobase::stats_enabled())
                    return;
                _res_dir = tools::hostname() + "_" + tools::date() + "_" + tools::getpid();
                boost::filesystem::path my_path(_res_dir);
                boost::filesystem::create_directory(my_path);
            }

            std::string _res_dir;
            int _current_iteration;
            int _total_iterations;
            stopping_criteria_t _stopping_criteria;
            stat_t _stat;

            std::vector<Eigen::VectorXd> _observations;
            std::vector<Eigen::VectorXd> _samples;
        };
    }
}

#endif
