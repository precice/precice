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
#ifndef LIMBO_STAT_GP_HPP
#define LIMBO_STAT_GP_HPP

#include <cmath>

#include <limbo/stat/stat_base.hpp>

namespace limbo {
    namespace stat {
        /// @ingroup stat
        /// filename: `gp_<iteration>.dat`
        template <typename Params>
        struct GP : public StatBase<Params> {
            template <typename BO, typename AggregatorFunction>
            void operator()(const BO& bo, const AggregatorFunction& afun)
            {
                if (!bo.stats_enabled())
                    return;
                std::string fname = bo.res_dir() + "/" + "gp_" + std::to_string(bo.total_iterations()) + ".dat";
                std::ofstream ofs(fname.c_str());
                int gp_in = bo.model().dim_in();
                int gp_out = bo.model().dim_out();
                ofs << "#Point[" << gp_in << "d] mu[" << gp_out << "d] sigma[1d] acquisition[1d]" << std::endl;
                _explore(0, ofs, bo, afun, Eigen::VectorXd::Constant(bo.model().dim_in(), 0));
            }

        protected:
            // recursively explore all the dimensions
            template <typename BO, typename AggregatorFunction>
            void _explore(int dim_in, std::ofstream& ofs,
                const BO& bo, const AggregatorFunction& afun,
                const Eigen::VectorXd& current) const
            {
                for (double x = 0; x <= 1.0f; x += 1.0f / (double)Params::stat_gp::bins()) {
                    Eigen::VectorXd point = current;
                    point[dim_in] = x;
                    if (dim_in == current.size() - 1) {
                        auto q = bo.model().query(point);
                        double acqui = opt::fun(typename BO::acquisition_function_t(bo.model(), bo.current_iteration())(point, afun, false));
                        ofs << point.transpose() << " "
                            << std::get<0>(q).transpose() << " "
                            << std::get<1>(q) << " "
                            << acqui << std::endl;
                    }
                    else {
                        _explore(dim_in + 1, ofs, bo, afun, point);
                    }
                }
            }
        };
    }
}

#endif
