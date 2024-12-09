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
#ifndef LIMBO_BAYES_OPT_BO_MULTI_HPP
#define LIMBO_BAYES_OPT_BO_MULTI_HPP

#include <Eigen/Core>

#include <limbo/bayes_opt/bo_base.hpp>
#include <limbo/experimental/tools/pareto.hpp>

#ifndef USE_SFERES
#warning No sferes
#else
#ifndef USE_TBB
#define NO_PARALLEL
#endif
// Quick hack for definition of 'I' in <complex.h>
#undef I
#include <sferes/gen/evo_float.hpp>
#include <sferes/phen/parameters.hpp>
#ifdef USE_TBB
#include <sferes/eval/parallel.hpp>
#endif
#include <sferes/ea/nsga2.hpp>
#include <sferes/modif/dummy.hpp>
#endif

namespace limbo {
    namespace experimental {
        namespace bayes_opt {
            namespace multi {
#ifdef USE_SFERES
                struct SferesParams {
                    struct evo_float {
                        using mutation_t = sferes::gen::evo_float::mutation_t;
                        using cross_over_t = sferes::gen::evo_float::cross_over_t;
                        SFERES_CONST float cross_rate = 0.5f;
                        SFERES_CONST float mutation_rate = 0.1f;
                        SFERES_CONST float eta_m = 15.0f;
                        SFERES_CONST float eta_c = 10.0f;
                        SFERES_CONST mutation_t mutation_type = sferes::gen::evo_float::polynomial;
                        SFERES_CONST cross_over_t cross_over_type = sferes::gen::evo_float::sbx;
                    };

                    struct pop {
                        SFERES_CONST unsigned size = 100;
                        SFERES_CONST unsigned nb_gen = 1000;
                        SFERES_CONST int dump_period = -1;
                        SFERES_CONST int initial_aleat = 1;
                    };

                    struct parameters {
                        SFERES_CONST float min = 0.0f;
                        SFERES_CONST float max = 1.0f;
                    };
                };

                SFERES_FITNESS(SferesFitBase, sferes::fit::Fitness){
                    template <typename Indiv>
                    void eval(const Indiv& indiv){}};

                template <typename M>
                class SferesFit : public SferesFitBase<> {
                public:
                    SferesFit(const std::vector<M>& models) : _models(models) {}
                    SferesFit() {}

                    const std::vector<float>& objs() const { return _objs; }

                    float obj(size_t i) const { return _objs[i]; }

                    template <typename Indiv>
                    void eval(const Indiv& indiv)
                    {
                        this->_objs.resize(_models.size());
                        Eigen::VectorXd v(indiv.size());
                        for (size_t j = 0; j < indiv.size(); ++j)
                            v[j] = indiv.data(j);
                        // we protect against overestimation because this has some spurious effect
                        for (size_t i = 0; i < _models.size(); ++i)
                            this->_objs[i] = std::min(_models[i].mu(v)(0), _models[i].max_observation()(0));
                    }

                protected:
                    std::vector<M> _models;
                };
#endif
            }

            // to be removed once moved out of experimental?
            BOOST_PARAMETER_TEMPLATE_KEYWORD(initfun)
            BOOST_PARAMETER_TEMPLATE_KEYWORD(acquifun)
            BOOST_PARAMETER_TEMPLATE_KEYWORD(modelfun)
            BOOST_PARAMETER_TEMPLATE_KEYWORD(statsfun)
            BOOST_PARAMETER_TEMPLATE_KEYWORD(stopcrit)

            using bo_multi_signature = boost::parameter::parameters<boost::parameter::optional<tag::statsfun>,
                boost::parameter::optional<tag::initfun>,
                boost::parameter::optional<tag::acquifun>,
                boost::parameter::optional<tag::stopcrit>,
                boost::parameter::optional<tag::modelfun>>;

            template <class Params,
                class A1 = boost::parameter::void_,
                class A2 = boost::parameter::void_,
                class A3 = boost::parameter::void_,
                class A4 = boost::parameter::void_,
                class A5 = boost::parameter::void_,
                class A6 = boost::parameter::void_>
            class BoMulti : public limbo::bayes_opt::BoBase<Params, A1, A2, A3, A4, A5, A6> {
            public:
                using args = typename bo_multi_signature::bind<A1, A2, A3, A4, A5, A6>::type;

                using base_t = limbo::bayes_opt::BoBase<Params, A1, A2, A3, A4, A5, A6>;
                using model_t = typename base_t::model_t;
                using acquisition_function_t = typename base_t::acquisition_function_t;
                // point, obj, sigma
                using pareto_point_t = std::tuple<Eigen::VectorXd, Eigen::VectorXd, Eigen::VectorXd>;
                using pareto_t = std::vector<pareto_point_t>;

                size_t nb_objs() const { return this->_observations[0].size(); }

                const pareto_t& pareto_model() const { return _pareto_model; }

                const pareto_t& pareto_data() const { return _pareto_data; }

                const std::vector<model_t>& models() const { return _models; }

                // will be called at the end of the algo
                void update_pareto_data()
                {
                    std::vector<Eigen::VectorXd> v(this->_samples.size());
                    size_t dim = this->_observations[0].size();
                    std::fill(v.begin(), v.end(), Eigen::VectorXd::Zero(dim));
                    _pareto_data = pareto::pareto_set<1>(
                        _pack_data(this->_samples, this->_observations, v));
                }

                // will be called at the end of the algo
                template <int D>
                void update_pareto_model()
                {
                    this->_update_models();
#ifdef USE_SFERES
                    using gen_t = sferes::gen::EvoFloat<D, multi::SferesParams>;
                    using phen_t = sferes::phen::Parameters<gen_t, multi::SferesFit<model_t>, multi::SferesParams>;
                    using eval_t = sferes::eval::Parallel<multi::SferesParams>;
                    using stat_t = boost::fusion::vector<>;
                    using modifier_t = sferes::modif::Dummy<>;
                    using nsga2_t = sferes::ea::Nsga2<phen_t, eval_t, stat_t, modifier_t, multi::SferesParams>;

                    // commented to remove a dependency to a particular version of sferes
                    nsga2_t ea;
                    ea.set_fit_proto(multi::SferesFit<model_t>(_models));
                    ea.run();
                    auto pareto_front = ea.pareto_front();
                    tools::par::sort(pareto_front.begin(), pareto_front.end(), sferes::fit::compare_objs_lex());
                    _pareto_model.resize(pareto_front.size());
                    Eigen::VectorXd point(D), objs(nb_objs()), sigma(nb_objs());
                    for (size_t p = 0; p < pareto_front.size(); ++p) {
                        for (size_t i = 0; i < pareto_front[p]->size(); ++i)
                            point(i) = pareto_front[p]->data(i);
                        for (size_t i = 0; i < nb_objs(); ++i) {
                            objs(i) = pareto_front[p]->fit().obj(i);
                            sigma(i) = _models[i].sigma(point);
                        }
                        _pareto_model[p] = std::make_tuple(point, objs, sigma);
                    }
#endif
                }

            protected:
                std::vector<model_t> _models;
                pareto_t _pareto_model;
                pareto_t _pareto_data;

                pareto_t _pack_data(const std::vector<Eigen::VectorXd>& points,
                    const std::vector<Eigen::VectorXd>& objs,
                    const std::vector<Eigen::VectorXd>& sigma) const
                {
                    assert(points.size() == objs.size());
                    assert(sigma.size() == objs.size());
                    pareto_t p(points.size());
                    tools::par::loop(0, p.size(), [&](size_t k) {
                        p[k] = std::make_tuple(points[k], objs[k], sigma[k]);
                    });
                    return p;
                }

                void _update_models()
                {
                    size_t dim = this->_samples[0].size();
                    std::vector<std::vector<Eigen::VectorXd>> uni_obs(nb_objs());
                    for (size_t i = 0; i < this->_observations.size(); ++i)
                        for (int j = 0; j < this->_observations[i].size(); ++j)
                            uni_obs[j].push_back(Eigen::VectorXd::Constant(1, this->_observations[i][j]));
                    std::vector<model_t> models(nb_objs(), model_t(dim, 1));
                    _models = models;
                    for (size_t i = 0; i < uni_obs.size(); ++i) {
                        _models[i].compute(this->_samples, uni_obs[i]);
                    }
                }
            };
        }
    }
}

#endif
