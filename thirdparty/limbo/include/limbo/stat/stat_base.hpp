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
#ifndef LIMBO_STAT_STAT_BASE_HPP
#define LIMBO_STAT_STAT_BASE_HPP

#include <fstream>
#include <string>

#include <memory>

namespace limbo {
    namespace stat {
        /**
          Base class for statistics

          The only method provided is protected :

          \rst
          .. code-block:: cpp

            template <typename BO>
            void _create_log_file(const BO& bo, const std::string& name)


          This method allocates an attribute `_log_file` (type: `std::shared_ptr<std::ofstream>`) if it has not been created yet, and does nothing otherwise. This method is designed so that you can safely call it in operator() while being 'guaranteed' that the file exists. Using this method is not mandatory for a statistics class.
          \endrst
        */
        template <typename Params>
        struct StatBase {
            StatBase() {}

            /// main method (to be written in derived classes)
            template <typename BO>
            void operator()(const BO& bo)
            {
                assert(false);
            }

        protected:
            std::shared_ptr<std::ofstream> _log_file;

            template <typename BO>
            void _create_log_file(const BO& bo, const std::string& name)
            {
                if (!_log_file && bo.stats_enabled()) {
                    std::string log = bo.res_dir() + "/" + name;
                    _log_file = std::make_shared<std::ofstream>(log.c_str());
                    assert(_log_file->good());
                }
            }
        };
    }
}

#endif
