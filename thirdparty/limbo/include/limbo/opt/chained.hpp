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
#ifndef LIMBO_OPT_CHAINED_HPP
#define LIMBO_OPT_CHAINED_HPP

#include <algorithm>

#include <Eigen/Core>

#include <limbo/opt/optimizer.hpp>

namespace limbo {
    namespace opt {

        // Needed for the variadic data structure
        template <typename Params, typename... Optimizers>
        struct Chained {
        };

        // Base case: just 1 optimizer to call
        template <typename Params, typename Optimizer>
        struct Chained<Params, Optimizer> {
            template <typename F>
            Eigen::VectorXd operator()(const F& f, const Eigen::VectorXd& init, bool bounded) const
            {
                return Optimizer()(f, init, bounded);
            };
        };

        // Recursive case: call current optimizer, and pass result as init value for the next one
        template <typename Params, typename Optimizer, typename... Optimizers>
        struct Chained<Params, Optimizer, Optimizers...> : Chained<Params, Optimizers...> {
            template <typename F>
            Eigen::VectorXd operator()(const F& f, const Eigen::VectorXd& init, bool bounded) const
            {
                return Chained<Params, Optimizers...>::operator()(f, Optimizer()(f, init, bounded), bounded);
            };
        };
    }
}

#endif
