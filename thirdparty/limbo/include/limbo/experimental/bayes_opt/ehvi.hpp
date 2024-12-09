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
#ifndef LIMBO_EXPERIMENTAL_BAYES_OPT_EHVI_HPP
#define LIMBO_EXPERIMENTAL_BAYES_OPT_EHVI_HPP

#include <algorithm>

#include <ehvi/ehvi_calculations.h>
#include <ehvi/ehvi_sliceupdate.h>

#include <limbo/tools/macros.hpp>

#include <limbo/experimental/acqui/ehvi.hpp>
#include <limbo/experimental/bayes_opt/bo_multi.hpp>

namespace limbo {
    namespace defaults {
        struct bayes_opt_ehvi {
            BO_PARAM(double, x_ref, -11);
            BO_PARAM(double, y_ref, -11);
        };
    }

    namespace experimental {
        namespace bayes_opt {

            BOOST_PARAMETER_TEMPLATE_KEYWORD(acquiopt)

            using ehvi_signature = boost::parameter::parameters<boost::parameter::optional<tag::acquiopt>>;

            template <class Params,
                class A1 = boost::parameter::void_,
                class A2 = boost::parameter::void_,
                class A3 = boost::parameter::void_,
                class A4 = boost::parameter::void_,
                class A5 = boost::parameter::void_,
                class A6 = boost::parameter::void_>

            class Ehvi : public BoMulti<Params, A1, A2, A3, A4, A5, A6> {
            public:
                struct defaults {
#ifdef USE_NLOPT
                    using acquiopt_t = opt::NLOptNoGrad<Params, nlopt::GN_DIRECT_L_RAND>;
#elif defined(USE_LIBCMAES)
                    using acquiopt_t = opt::Cmaes<Params>;
#else
#warning NO NLOpt, and NO Libcmaes: the acquisition function will be optimized by a grid search algorithm (which is usually bad). Please install at least NLOpt or libcmaes to use limbo!.
                    using acquiopt_t = opt::GridSearch<Params>;
#endif
                };

                using args = typename ehvi_signature::bind<A1, A2, A3, A4, A5, A6>::type;
                using acqui_optimizer_t = typename boost::parameter::binding<args, tag::acquiopt, typename defaults::acquiopt_t>::type;

                using pareto_point_t = std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>;
                using base_t = limbo::experimental::bayes_opt::BoMulti<Params, A1, A2, A3, A4, A5, A6>;
                using model_t = typename base_t::model_t;

                template <typename EvalFunction>
                void optimize(const EvalFunction& feval, bool reset = true)
                {
                    this->_init(feval, FirstElem(), reset);

                    acqui_optimizer_t inner_opt;

                    while (this->_samples.size() == 0 || !this->_stop(*this, FirstElem())) {
                        this->template update_pareto_model<EvalFunction::dim_in()>();
                        this->update_pareto_data();

                        // copy in the ehvi structure to compute expected improvement
                        std::deque<individual*> pop;
                        for (auto x : this->pareto_data()) {
                            individual* ind = new individual;
                            ind->f[0] = std::get<1>(x)(0);
                            ind->f[1] = std::get<1>(x)(1);
                            ind->f[2] = 0;
                            pop.push_back(ind);
                        }

                        auto acqui = acqui::Ehvi<Params, model_t>(
                            this->_models, pop,
                            Eigen::Vector3d(Params::bayes_opt_ehvi::x_ref(), Params::bayes_opt_ehvi::y_ref(), 0));

                        // maximize with inner opt
                        using pair_t = std::pair<Eigen::VectorXd, double>;
                        pair_t init(Eigen::VectorXd::Zero(1), -std::numeric_limits<float>::max());

                        auto body = [&](int i) -> pair_t {
                            auto x = this->pareto_data()[i];

                            auto acqui_optimization =
                                [&](const Eigen::VectorXd& x, bool g) { return opt::no_grad(acqui(x)); };

                            Eigen::VectorXd s = inner_opt(acqui_optimization, std::get<0>(x), true);
                            double hv = acqui(s);

                            return std::make_pair(s, hv);
                        };

                        auto comp = [](const pair_t& v1, const pair_t& v2) {
                            return v1.second > v2.second;
                        };

                        auto m = tools::par::max(init, this->pareto_data().size(), body, comp);

                        // take the best
                        Eigen::VectorXd new_sample = m.first;

                        // delete pop
                        for (auto x : pop)
                            delete x;

                        // add sample
                        this->add_new_sample(new_sample, feval(new_sample));
                        this->_update_stats(*this, FirstElem());
                        this->_current_iteration++;
                        this->_total_iterations++;
                    }
                }
            };
        }
    }
}

#endif
