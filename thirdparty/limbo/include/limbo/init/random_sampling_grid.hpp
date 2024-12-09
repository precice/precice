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
#ifndef LIMBO_INIT_RANDOM_SAMPLING_GRID_HPP
#define LIMBO_INIT_RANDOM_SAMPLING_GRID_HPP

#include <Eigen/Core>

#include <limbo/tools/macros.hpp>
#include <limbo/tools/random_generator.hpp>

namespace limbo {
    namespace defaults {
        struct init_randomsamplinggrid {
            ///@ingroup init_defaults
            BO_PARAM(int, samples, 10);
            ///@ingroup init_defaults
            BO_PARAM(int, bins, 5);
        };
    }
    namespace init {
        /** @ingroup init
          \rst
          Grid-based random sampling: in effect, define a grid and takes random samples from this grid. The number of bins in the grid is ``bins``^number_of_dimensions

          For instance, if bins = 5 and there are 3 inputs dimensions, then the grid is 5*5*5=125, and random points will all be points of this grid.

          Parameters:
            - ``int samples`` (total number of samples)
            - ``int bins`` (number of bins for each dimensions)
          \endrst
        */
        template <typename Params>
        struct RandomSamplingGrid {
            template <typename StateFunction, typename AggregatorFunction, typename Opt>
            void operator()(const StateFunction& seval, const AggregatorFunction&, Opt& opt) const
            {
                // Only works with bounded BO
                assert(Params::bayes_opt_bobase::bounded());

                tools::rgen_int_t rgen(0, Params::init_randomsamplinggrid::bins());
                for (int i = 0; i < Params::init_randomsamplinggrid::samples(); i++) {
                    Eigen::VectorXd new_sample(StateFunction::dim_in());
                    for (size_t i = 0; i < StateFunction::dim_in(); i++)
                        new_sample[i] = rgen.rand() / double(Params::init_randomsamplinggrid::bins());
                    opt.eval_and_add(seval, new_sample);
                }
            }
        };
    }
}

#endif
