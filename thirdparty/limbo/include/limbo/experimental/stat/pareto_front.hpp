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
#ifndef LIMBO_STAT_PARETO_FRONT_HPP
#define LIMBO_STAT_PARETO_FRONT_HPP

#include <limbo/experimental/tools/pareto.hpp>
#include <limbo/stat/stat_base.hpp>

namespace limbo {
    namespace experimental {
        namespace stat {
            template <typename Params>
            struct ParetoFront : public limbo::stat::StatBase<Params> {
                // point, obj, sigma
                using pareto_point_t = std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>;
                using pareto_t = std::vector<pareto_point_t>;

                template <typename BO, typename AggregatorFunction>
                void operator()(const BO& bo, const AggregatorFunction&)
                {
                    if (!bo.stats_enabled() || bo.observations().empty())
                        return;
                    std::string fname = bo.res_dir() + "/" + "pareto_front_" + std::to_string(bo.current_iteration()) + ".dat";
                    std::ofstream ofs(fname.c_str());
                    auto pareto = _pareto_data(bo);
                    for (auto x : pareto) {
                        ofs << std::get<0>(x).transpose() << " "
                            << std::get<1>(x).transpose() << std::endl;
                    }
                }

            protected:
                template <typename BO>
                pareto_t _pareto_data(const BO& bo)
                {
                    std::vector<Eigen::VectorXd> v(bo.samples().size());
                    size_t dim = bo.observations().size();
                    std::fill(v.begin(), v.end(), Eigen::VectorXd::Zero(dim));
                    return pareto::pareto_set<1>(
                        _pack_data(bo.samples(), bo.observations(), v));
                }
                pareto_t _pack_data(const std::vector<Eigen::VectorXd>& points,
                    const std::vector<Eigen::VectorXd>& objs,
                    const std::vector<Eigen::VectorXd>& sigma) const
                {
                    assert(points.size() == objs.size());
                    assert(sigma.size() == objs.size());
                    pareto_t p(points.size());
                    tools::par::loop(0, p.size(), [&](size_t k) {
                        // clang-format off
                    p[k] = std::make_tuple(points[k], objs[k], sigma[k]);
                        // clang-format on
                    });
                    return p;
                }
            };
        }
    }
}

#endif
