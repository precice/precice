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
#ifndef LIMBO_BAYES_OPT_NSBO_HPP
#define LIMBO_BAYES_OPT_NSBO_HPP

#include <algorithm>

#include <limbo/experimental/bayes_opt/bo_multi.hpp>

namespace limbo {
    namespace experimental {
        namespace bayes_opt {

            template <class Params,
                // clang-format off
              class A2 = boost::parameter::void_,
              class A3 = boost::parameter::void_,
              class A4 = boost::parameter::void_,
              class A5 = boost::parameter::void_,
              class A6 = boost::parameter::void_>
            // clang-format on
            class Nsbo : public BoMulti<Params, A2, A3, A4, A5, A6> {
            public:
                using pareto_point_t = std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>;

                template <typename EvalFunction>
                void optimize(const EvalFunction& feval, bool reset = true)
                {
                    this->_init(feval, FirstElem(), reset);

                    while (this->_samples.size() == 0 || !this->_stop(*this, FirstElem())) {
                        std::cout << "updating pareto model...";
                        std::cout.flush();
                        this->template update_pareto_model<EvalFunction::dim_in()>();
                        std::cout << "ok" << std::endl;
                        auto pareto = this->pareto_model();

                        // Pareto front of the variances
                        auto p_variance = pareto::pareto_set<2>(pareto);
                        auto best = p_variance[rand() % p_variance.size()];
                        Eigen::VectorXd best_v = std::get<0>(best);

                        this->add_new_sample(best_v, feval(best_v));
                        std::cout << this->_current_iteration << " | " << best_v.transpose() << "-> "
                                  << this->_observations.back().transpose()
                                  << " (expected:" << this->_models[0].mu(best_v) << " "
                                  << this->_models[1].mu(best_v) << ")"
                                  << " sigma:" << this->_models[0].sigma(best_v) << " "
                                  << this->_models[1].sigma(best_v) << std::endl;
                        this->_update_stats(*this, FirstElem());
                        this->_current_iteration++;
                        this->_total_iterations++;
                    }
                }

            protected:
                // former way to deal with the template of RefreshStat
                // void _update_stats()
                // {
                //     boost::fusion::for_each(this->_stat, RefreshStat_f<Nsbo>(*this));
                // }
            };
        }
    }
}

#endif
