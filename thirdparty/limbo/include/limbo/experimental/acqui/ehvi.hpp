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
#ifndef LIMBO_EXPERIMENTAL_ACQUI_EHVI_HPP
#define LIMBO_EXPERIMENTAL_ACQUI_EHVI_HPP

#include <vector>

#include <Eigen/Core>

#include <ehvi/ehvi_calculations.h>

namespace limbo {
    namespace acqui {
        // only work in 2D for now
        template <typename Params, typename Model>
        class Ehvi {
        public:
            Ehvi(const std::vector<Model>& models, const std::deque<individual*>& pop,
                const Eigen::VectorXd& ref_point)
                : _models(models), _pop(pop), _ref_point(ref_point)
            {
                assert(_models.size() == 2);
            }

            size_t dim() const { return _models[0].dim(); }

            template <typename AggregatorFunction = FirstElem>
            double operator()(const Eigen::VectorXd& v, const AggregatorFunction& afun = AggregatorFunction()) const
            {
                assert(_models.size() == 2);
                double r[3] = {_ref_point(0), _ref_point(1), _ref_point(2)};
                double mu[3] = {afun(_models[0].mu(v)), afun(_models[1].mu(v)), 0};
                double s[3] = {_models[0].sigma(v), _models[1].sigma(v), 0};
                double ehvi = ehvi2d(_pop, r, mu, s);
                return ehvi;
            }

        protected:
            const std::vector<Model>& _models;
            const std::deque<individual*>& _pop;
            Eigen::VectorXd _ref_point;
        };
    }
}

#endif
