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
#ifndef LIMBO_BAYES_OPT_PAREGO_HPP
#define LIMBO_BAYES_OPT_PAREGO_HPP

#include <algorithm>

#include <limbo/bayes_opt/bo_base.hpp>
#include <limbo/bayes_opt/boptimizer.hpp>
#include <limbo/experimental/model/gp_parego.hpp>
#include <limbo/tools/macros.hpp>

namespace limbo {
    namespace experimental {
        namespace bayes_opt {

            BOOST_PARAMETER_TEMPLATE_KEYWORD(parego_modelfun)

            using parego_signature = boost::parameter::parameters<boost::parameter::optional<tag::parego_modelfun>>;

            // clang-format off
            template <class Params,
                      class A1 = boost::parameter::void_,
                      class A2 = boost::parameter::void_,
                      class A3 = boost::parameter::void_,
                      class A4 = boost::parameter::void_,
                      class A5 = boost::parameter::void_>
            // we find the model a wrap it into a GPParego
            // YOU NEED TO PASS A parego_modelfun and not a modelfun !
            class Parego : public limbo::bayes_opt::BOptimizer<
              Params,
              modelfun<
                  typename model::GPParego<Params,
                    typename boost::parameter::binding<
                         typename parego_signature::bind<A1, A2, A3, A4, A5>::type,
                         tag::parego_modelfun,
                         typename limbo::bayes_opt::BoBase<Params, A1, A2, A3, A4, A5>::defaults::model_t
                    > ::type // end binding
                  > // end GPParego
              > // end model_fun
            , A1, A2, A3, A4, A5> {//pass the remaining arguments
              // nothing here !
            };
            // clang-format on
        }
    }
}

#endif
