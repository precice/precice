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
#ifndef LIMBO_OPT_OPTIMIZER_HPP
#define LIMBO_OPT_OPTIMIZER_HPP

#include <tuple>

#include <Eigen/Core>

#include <boost/optional.hpp>

namespace limbo {
    ///\defgroup opt_tools
    namespace opt {

        ///@ingroup opt_tools
        /// return type of the function to optimize
        using eval_t = std::pair<double, boost::optional<Eigen::VectorXd>>;

        ///@ingroup opt_tools
        ///return with opt::no_grad(your_val) if no gradient is available (to be used in functions to be optimized)
        inline eval_t no_grad(double x) { return eval_t{x, boost::optional<Eigen::VectorXd>{}}; }

        ///@ingroup opt_tools
        /// get the gradient from a function evaluation (eval_t)
        inline const Eigen::VectorXd& grad(const eval_t& fg)
        {
            assert(std::get<1>(fg).is_initialized());
            return std::get<1>(fg).get();
        }

        ///@ingroup opt_tools
        /// get the value from a function evaluation (eval_t)
        inline double fun(const eval_t& fg)
        {
            return std::get<0>(fg);
        }

        ///@ingroup opt_tools
        /// Evaluate f without gradient (to be called from the optimization algorithms that do not use the gradient)
        template <typename F>
        inline double eval(const F& f, const Eigen::VectorXd& x)
        {
            return std::get<0>(f(x, false));
        }

        ///@ingroup opt_tools
        /// Evaluate f with gradient (to be called from the optimization algorithms that use the gradient)
        template <typename F>
        inline eval_t eval_grad(const F& f, const Eigen::VectorXd& x)
        {
            return f(x, true);
        }
    }
}

#endif
