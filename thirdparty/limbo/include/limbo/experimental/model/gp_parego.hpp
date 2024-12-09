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
#ifndef LIMBO_MODEL_GP_PAREGO_HPP
#define LIMBO_MODEL_GP_PAREGO_HPP

#include <cassert>
#include <iostream>
#include <limits>
#include <vector>

#include <Eigen/Cholesky>
#include <Eigen/Core>
#include <Eigen/LU>

#include <limbo/model/gp/no_lf_opt.hpp>

namespace limbo {
    namespace experimental {
        namespace defaults {
            struct model_gp_parego {
                BO_PARAM(double, rho, 0.05);
            };
        }
        namespace model {

            /// this is the model used in Parego
            /// reference: Knowles, J. (2006). ParEGO: A hybrid algorithm
            /// with on-line landscape approximation for expensive multiobjective
            /// optimization problems.
            /// IEEE Transactions On Evolutionary Computation, 10(1), 50-66.
            /// Main idea:
            /// - this models aggregates all the objective values with the Tchebycheff distance
            /// - objectives are weighted using a random vector
            /// - a single model is built
            template <typename Params, typename Model>
            class GPParego : public Model {
            public:
                GPParego() {}
                GPParego(int dim_in, int dim_out) : Model(dim_in, 1), _nb_objs(dim_out) {}
                void compute(const std::vector<Eigen::VectorXd>& samples,
                    const std::vector<Eigen::VectorXd>& observations)
                {
                    _raw_observations = observations;
                    _nb_objs = observations[0].size();
                    auto new_observations = _scalarize_obs(observations);
                    Model::compute(samples, new_observations);
                }
                /// add sample will NOT be incremental (we call compute each time)
                void add_sample(const Eigen::VectorXd& sample, const Eigen::VectorXd& observation)
                {
                    this->_samples.push_back(sample);
                    _raw_observations.push_back(observation);

                    this->compute(this->_samples, _raw_observations);
                }

            protected:
                size_t _nb_objs;
                std::vector<Eigen::VectorXd> _raw_observations;
                std::vector<Eigen::VectorXd> _scalarize_obs(const std::vector<Eigen::VectorXd>& observations)
                {
                    Eigen::VectorXd lambda = tools::random_vector(_nb_objs);
                    double sum = lambda.sum();
                    lambda = lambda / sum;
                    // scalarize (Tchebycheff)
                    std::vector<Eigen::VectorXd> scalarized;
                    for (auto x : observations) {
                        double y = (lambda.array() * x.array()).maxCoeff();
                        double s = (lambda.array() * x.array()).sum();
                        auto v = tools::make_vector(y + Params::model_gp_parego::rho() * s);
                        scalarized.push_back(v);
                    }
                    return scalarized;
                }
            };
        }
    }
}

#endif
