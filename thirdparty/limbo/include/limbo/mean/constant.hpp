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
#ifndef LIMBO_MEAN_CONSTANT_HPP
#define LIMBO_MEAN_CONSTANT_HPP

#include <limbo/mean/mean.hpp>

namespace limbo {
    namespace defaults {
        struct mean_constant {
            ///@ingroup mean_defaults
            BO_PARAM(double, constant, 1);
        };
    } // namespace defaults

    namespace mean {
        /** @ingroup mean
          A constant mean (the traditionnal choice for Bayesian optimization)

          Parameter:
            - ``double constant`` (the value of the constant)
        */
        template <typename Params>
        struct Constant : public BaseMean<Params> {
            Constant(size_t dim_out = 1) : _dim_out(dim_out), _constant(Params::mean_constant::constant()) {}

            template <typename GP>
            Eigen::VectorXd operator()(const Eigen::VectorXd& v, const GP&) const
            {
                return Eigen::VectorXd::Constant(_dim_out, _constant);
            }

            template <typename GP>
            Eigen::MatrixXd grad(const Eigen::VectorXd& x, const GP& gp) const
            {
                return Eigen::MatrixXd::Ones(_dim_out, 1);
            }

            size_t h_params_size() const { return 1; }

            Eigen::VectorXd h_params() const
            {
                Eigen::VectorXd p(1);
                p << _constant;
                return p;
            }

            void set_h_params(const Eigen::VectorXd& p)
            {
                _constant = p(0);
            }

        protected:
            size_t _dim_out;
            double _constant;
        };
    } // namespace mean
} // namespace limbo

#endif
