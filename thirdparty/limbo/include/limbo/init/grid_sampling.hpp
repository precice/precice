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
#ifndef LIMBO_INIT_GRID_SAMPLING_HPP
#define LIMBO_INIT_GRID_SAMPLING_HPP

#include <Eigen/Core>

#include <limbo/tools/macros.hpp>

namespace limbo {
    namespace defaults {
        struct init_gridsampling {
            ///@ingroup init_defaults
            BO_PARAM(int, bins, 5);
        };
    }
    namespace init {
        /** @ingroup init
          \rst
          Grid sampling.

          Parameter:
            - ``int bins`` (number of bins)
          \endrst
        */
        template <typename Params>
        struct GridSampling {
            template <typename StateFunction, typename AggregatorFunction, typename Opt>
            void operator()(const StateFunction& seval, const AggregatorFunction&, Opt& opt) const
            {
                _explore(0, seval, Eigen::VectorXd::Constant(StateFunction::dim_in(), 0), opt);
            }

        private:
            // recursively explore all the dimensions
            template <typename StateFunction, typename Opt>
            void _explore(int dim_in, const StateFunction& seval, const Eigen::VectorXd& current,
                Opt& opt) const
            {
                for (double x = 0; x <= 1.0f; x += 1.0f / (double)Params::init_gridsampling::bins()) {
                    Eigen::VectorXd point = current;
                    point[dim_in] = x;
                    if (dim_in == current.size() - 1) {
                        opt.eval_and_add(seval, point);
                    }
                    else {
                        _explore(dim_in + 1, seval, point, opt);
                    }
                }
            }
        };
    }
}
#endif
