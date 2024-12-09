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
#ifndef LIMBO_OPT_GRID_SEARCH_HPP
#define LIMBO_OPT_GRID_SEARCH_HPP

#include <limits>

#include <Eigen/Core>

#include <limbo/opt/optimizer.hpp>
#include <limbo/tools/macros.hpp>

namespace limbo {
    namespace defaults {
        struct opt_gridsearch {
            /// @ingroup opt_defaults
            /// number of bins for each dimension
            BO_PARAM(int, bins, 5);
        };
    }
    namespace opt {
        /// @ingroup opt
        /// Grid search
        ///
        /// Parameters:
        /// - int bins
        template <typename Params>
        struct GridSearch {
        public:
            template <typename F>
            Eigen::VectorXd operator()(const F& f, const Eigen::VectorXd& init, bool bounded) const
            {
                // Grid search does not support unbounded search
                assert(bounded);
                size_t dim = init.size();
                return _inner_search(f, 0, Eigen::VectorXd::Constant(dim, 0.5));
            }

        protected:
            template <typename F>
            Eigen::VectorXd _inner_search(const F& f, size_t depth, const Eigen::VectorXd& current) const
            {
                size_t dim = current.size();
                double step_size = 1.0 / (double)Params::opt_gridsearch::bins();
                double upper_lim = 1.0 + step_size;
                double best_fit = -std::numeric_limits<double>::max();
                Eigen::VectorXd current_result(dim);
                for (double x = 0; x < upper_lim; x += step_size) {
                    Eigen::VectorXd new_point = current;
                    new_point[depth] = x;
                    double val;
                    if (depth == dim - 1) {
                        val = eval(f, new_point);
                        if (val > best_fit) {
                            best_fit = val;
                            current_result = new_point;
                        }
                    }
                    else {
                        Eigen::VectorXd temp_result = _inner_search(f, depth + 1, new_point);
                        val = eval(f, temp_result);
                        if (val > best_fit) {
                            best_fit = val;
                            current_result = temp_result;
                        }
                    }
                }
                return current_result;
            }
        };
    }
}

#endif
