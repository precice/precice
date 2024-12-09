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
#ifndef LIMBO_OPT_NLOPT_BASE_HPP
#define LIMBO_OPT_NLOPT_BASE_HPP

#ifndef USE_NLOPT
#warning No NLOpt
#else
#include <Eigen/Core>

#include <vector>

#include <nlopt.hpp>

#include <limbo/opt/optimizer.hpp>
#include <limbo/tools/macros.hpp>

namespace limbo {
    namespace opt {
        /**
          @ingroup opt
        Base class for NLOpt wrappers
        */
        template <typename Params, nlopt::algorithm Algorithm>
        struct NLOptBase {
        public:
            virtual void initialize(int dim)
            {
                _opt = nlopt::opt(Algorithm, dim);
                _initialized = true;
            }

            template <typename F>
            Eigen::VectorXd operator()(const F& f, const Eigen::VectorXd& init, bool bounded)
            {
                int dim = init.size();
                if (!_initialized)
                    initialize(dim);

                _opt.set_max_objective(nlopt_func<F>, (void*)&f);

                std::vector<double> x(dim);
                Eigen::VectorXd::Map(&x[0], dim) = init;

                if (bounded) {
                    _opt.set_lower_bounds(std::vector<double>(dim, 0));
                    _opt.set_upper_bounds(std::vector<double>(dim, 1));
                }

                double max;

                try {
                    _opt.optimize(x, max);
                }
                catch (nlopt::roundoff_limited& e) {
                    // In theory it's ok to ignore this error
                    std::cerr << "[NLOptNoGrad]: " << e.what() << std::endl;
                }
                catch (std::invalid_argument& e) {
                    // In theory it's ok to ignore this error
                    std::cerr << "[NLOptNoGrad]: " << e.what() << std::endl;
                }
                catch (std::runtime_error& e) {
                    // In theory it's ok to ignore this error
                    std::cerr << "[NLOptGrad]: " << e.what() << std::endl;
                }

                return Eigen::VectorXd::Map(x.data(), x.size());
            }

            // Inequality constraints of the form f(x) <= 0
            template <typename F>
            void add_inequality_constraint(const F& f, double tolerance = 1e-8)
            {
                if (_initialized)
                    _opt.add_inequality_constraint(nlopt_func<F>, (void*)&f, tolerance);
                else
                    std::cerr << "[NLOptNoGrad]: Trying to set an inequality constraint without having initialized the optimizer! Nothing will be done!" << std::endl;
            }

            // Equality constraints of the form f(x) = 0
            template <typename F>
            void add_equality_constraint(const F& f, double tolerance = 1e-8)
            {
                if (_initialized)
                    _opt.add_equality_constraint(nlopt_func<F>, (void*)&f, tolerance);
                else
                    std::cerr << "[NLOptNoGrad]: Trying to set an equality constraint without having initialized the optimizer! Nothing will be done!" << std::endl;
            }

        protected:
            nlopt::opt _opt;
            bool _initialized = false;

            template <typename F>
            static double nlopt_func(const std::vector<double>& x, std::vector<double>& grad, void* my_func_data)
            {
                F* f = (F*)(my_func_data);
                Eigen::VectorXd params = Eigen::VectorXd::Map(x.data(), x.size());
                double v;
                if (!grad.empty()) {
                    auto r = eval_grad(*f, params);
                    v = opt::fun(r);
                    Eigen::VectorXd g = opt::grad(r);
                    Eigen::VectorXd::Map(&grad[0], g.size()) = g;
                }
                else {
                    v = eval(*f, params);
                }
                return v;
            }
        };
    } // namespace opt
} // namespace limbo

#endif
#endif
