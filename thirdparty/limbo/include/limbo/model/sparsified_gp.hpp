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
#ifndef LIMBO_MODEL_SPARSE_GP_HPP
#define LIMBO_MODEL_SPARSE_GP_HPP

#include <algorithm>
#include <mutex>

#include <limbo/model/gp.hpp>
#include <limbo/tools/parallel.hpp>

namespace limbo {
    namespace defaults {
        struct model_sparse_gp {
            BO_PARAM(int, max_points, 200);
        };
    } // namespace defaults

    namespace model {
        /// @ingroup model
        /// A sparse Gaussian process model.
        /// It is parametrized by:
        /// - a kernel function
        /// - a mean function
        /// - [optional] an optimizer for the hyper-parameters
        /// A sparsification based on the density of points is performed
        /// until a desired number of points is reached
        template <typename Params, typename KernelFunction = kernel::MaternFiveHalves<Params>, typename MeanFunction = mean::Data<Params>, typename HyperParamsOptimizer = gp::NoLFOpt<Params>>
        class SparsifiedGP : public GP<Params, KernelFunction, MeanFunction, HyperParamsOptimizer> {
        public:
            using base_gp_t = GP<Params, KernelFunction, MeanFunction, HyperParamsOptimizer>;

            /// useful because the model might be created before knowing anything about the process
            SparsifiedGP() : base_gp_t() {}

            /// useful because the model might be created before having samples
            SparsifiedGP(int dim_in, int dim_out)
                : base_gp_t(dim_in, dim_out) {}

            /// Compute the GP from samples and observations. This call needs to be explicit!
            void compute(const std::vector<Eigen::VectorXd>& samples,
                const std::vector<Eigen::VectorXd>& observations, bool compute_kernel = true)
            {
                /// if the number of samples is less or equal than the desired
                /// compute the normal GP
                if (samples.size() <= Params::model_sparse_gp::max_points())
                    base_gp_t::compute(samples, observations, compute_kernel);
                /// otherwise, sparsify the samples
                else {
                    std::vector<Eigen::VectorXd> samp, obs;

                    std::tie(samp, obs) = _sparsify(samples, observations);

                    /// now compute the normal GP with less points
                    base_gp_t::compute(samp, obs, compute_kernel);
                }
            }

            /// add sample and update the GP. If the number of samples is bigger than
            /// the desired maximum points, we re-sparsify and re-compute the GP
            void add_sample(const Eigen::VectorXd& sample, const Eigen::VectorXd& observation)
            {
                base_gp_t::add_sample(sample, observation);
                /// if we surpassed the maximum points, re-sparsify
                /// and recompute
                if (this->_samples.size() > Params::model_sparse_gp::max_points()) {
                    /// get observations in appropriate format
                    std::vector<Eigen::VectorXd> observations;
                    for (size_t i = 0; i < this->_samples.size(); i++) {
                        observations.push_back(this->_observations.row(i));
                    }

                    compute(this->_samples, observations, true);
                }
            }

        protected:
            /// get the densest point in a list of samples
            /// D is the dimensionality of the samples
            /// N is the number of samples
            /// distances is an NxN matrix where element (i,j) contains
            /// the (pre)computed distance between the ith and the jth samples
            int _get_most_dense_point(int D, int N, const Eigen::MatrixXd& distances) const
            {
                std::mutex update_mutex;
                double min_dist = std::numeric_limits<double>::max();
                int denser = -1;

                tools::par::loop(0u, N, [&](size_t i) {
                    double dist = 0.;
                    std::vector<double> neighbors(N);
                    Eigen::VectorXd::Map(neighbors.data(), neighbors.size()) = distances.row(i);
                    /// remove self distance
                    neighbors.erase(neighbors.begin() + i);

                    std::partial_sort(neighbors.begin(), neighbors.begin() + D, neighbors.end());

                    for (int j = 0; j < D; j++) {
                        dist += neighbors[j];
                    }

                    update_mutex.lock();
                    if (dist < min_dist) {
                        min_dist = dist;
                        denser = i;
                    }
                    update_mutex.unlock();
                });

                return denser;
            }

            /// get sparsified data
            std::pair<std::vector<Eigen::VectorXd>, std::vector<Eigen::VectorXd>> _sparsify(const std::vector<Eigen::VectorXd>& samples, const std::vector<Eigen::VectorXd>& observations) const
            {
                /// compute big distance matrix to avoid re-computations
                Eigen::MatrixXd distances(samples.size(), samples.size());
                tools::par::loop(0u, samples.size(), [&](size_t i) {
                    for (size_t j = 0; j < samples.size(); j++) {
                        if (i == j)
                            continue;
                        distances(i, j) = (samples[i] - samples[j]).norm();
                    }
                });

                std::vector<Eigen::VectorXd> samp = samples, obs = observations;
                while (samp.size() > Params::model_sparse_gp::max_points()) {
                    int k = _get_most_dense_point(samp[0].size(), samp.size(), distances);
                    /// sanity check
                    if (k < 0)
                        break;
                    samp.erase(samp.begin() + k);
                    obs.erase(obs.begin() + k);

                    _remove_column(distances, k);
                    _remove_row(distances, k);
                }

                return std::make_pair(samp, obs);
            }

            /// remove row from an Eigen matrix
            void _remove_row(Eigen::MatrixXd& matrix, unsigned int rowToRemove) const
            {
                unsigned int numRows = matrix.rows() - 1;
                unsigned int numCols = matrix.cols();

                if (rowToRemove < numRows)
                    matrix.block(rowToRemove, 0, numRows - rowToRemove, numCols) = matrix.block(rowToRemove + 1, 0, numRows - rowToRemove, numCols);

                matrix.conservativeResize(numRows, numCols);
            }

            /// remove column from an Eigen matrix
            void _remove_column(Eigen::MatrixXd& matrix, unsigned int colToRemove) const
            {
                unsigned int numRows = matrix.rows();
                unsigned int numCols = matrix.cols() - 1;

                if (colToRemove < numCols)
                    matrix.block(0, colToRemove, numRows, numCols - colToRemove) = matrix.block(0, colToRemove + 1, numRows, numCols - colToRemove);

                matrix.conservativeResize(numRows, numCols);
            }
        };
    } // namespace model
} // namespace limbo

#endif
