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
#ifndef LIMBO_TOOLS_PARALLEL_HPP
#define LIMBO_TOOLS_PARALLEL_HPP

#include <algorithm>
#include <vector>

#ifdef USE_TBB
// Quick hack for definition of 'I' in <complex.h>
#undef I
#include <tbb/blocked_range.h>
#include <tbb/concurrent_vector.h>
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>
#include <tbb/parallel_reduce.h>
#include <tbb/parallel_sort.h>

#ifndef USE_TBB_ONEAPI
#include <tbb/task_scheduler_init.h>
#else
#include <oneapi/tbb/global_control.h>
using namespace oneapi;
#endif

#endif

///@defgroup par_tools

namespace limbo {
    namespace tools {
        namespace par {
#ifdef USE_TBB
#ifdef __GXX_EXPERIMENTAL_CXX0X__
            template <typename X> // old fashion way to create template alias (for GCC
            // 4.6...)
            struct vector {
                typedef tbb::concurrent_vector<X> type;
            };
#else
            template <typename X>
            using vector = tbb::concurrent_vector<X>; // Template alias (for GCC 4.7 and later)
#endif
            /// @ingroup par_tools
            /// convert a std::vector to something else (e.g. a std::list)
            template <typename V>
            inline std::vector<typename V::value_type> convert_vector(const V& v)
            {
                std::vector<typename V::value_type> v2(v.size());
                std::copy(v.begin(), v.end(), v2.begin());
                return v2;
            }
#else
#ifdef __GXX_EXPERIMENTAL_CXX0X__
            template <typename X> // old fashion way to create template alias (for GCC
            // 4.6...)
            struct vector {
                typedef std::vector<X> type;
            };
#else
            template <typename X>
            using vector = std::vector<X>; // Template alias (for GCC 4.7 and later)
#endif

            template <typename V>
            inline V convert_vector(const V& v)
            {
                return v;
            }

#endif

#ifdef USE_TBB
            inline void init(int threads = -1)
            {
#ifndef USE_TBB_ONEAPI
            static tbb::task_scheduler_init init(threads);
#else
            if (threads < 0)
                threads = tbb::info::default_concurrency();
            static tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, threads);

#endif
            }
#else
            /// @ingroup par_tools
            /// init TBB (if activated) for multi-core computing
            inline void init(int threads = -1)
            {
            }
#endif

            ///@ingroup par_tools
            /// parallel for
            template <typename F>
            inline void loop(size_t begin, size_t end, const F& f)
            {
#ifdef USE_TBB
                tbb::parallel_for(size_t(begin), end, size_t(1), [&](size_t i) {
                    // clang-format off
                f(i);
                    // clang-format on
                });
#else
                for (size_t i = begin; i < end; ++i)
                    f(i);
#endif
            }

            /// @ingroup par_tools
            /// parallel for_each
            template <typename Iterator, typename F>
            inline void for_each(Iterator begin, Iterator end, const F& f)
            {
#ifdef USE_TBB
                tbb::parallel_for_each(begin, end, f);
#else
                for (Iterator i = begin; i != end; ++i)
                    f(*i);
#endif
            }

            /// @ingroup par_tools
            /// parallel max
            template <typename T, typename F, typename C>
            inline T max(const T& init, int num_steps, const F& f, const C& comp)
            {
#ifdef USE_TBB
                auto body = [&](const tbb::blocked_range<size_t>& r, T current_max) -> T {
                    // clang-format off
            for (size_t i = r.begin(); i != r.end(); ++i)
            {
                T v = f(i);
                if (comp(v, current_max))
                  current_max = v;
            }
            return current_max;
                    // clang-format on
                };
                auto joint = [&](const T& p1, const T& p2) -> T {
                    // clang-format off
            if (comp(p1, p2))
                return p1;
            return p2;
                    // clang-format on
                };
                return tbb::parallel_reduce(tbb::blocked_range<size_t>(0, num_steps), init,
                    body, joint);
#else
                T current_max = init;
                for (int i = 0; i < num_steps; ++i) {
                    T v = f(i);
                    if (comp(v, current_max))
                        current_max = v;
                }
                return current_max;
#endif
            }
            /// @ingroup par_tools
            /// parallel sort
            template <typename T1, typename T2, typename T3>
            inline void sort(T1 i1, T2 i2, T3 comp)
            {
#ifdef USE_TBB
                tbb::parallel_sort(i1, i2, comp);
#else
                std::sort(i1, i2, comp);
#endif
            }

            /// @ingroup par_tools
            /// replicate a function nb times
            template <typename F>
            inline void replicate(size_t nb, const F& f)
            {
#ifdef USE_TBB
                tbb::parallel_for(size_t(0), nb, size_t(1), [&](size_t i) {
                    // clang-format off
                f();
                    // clang-format on
                });
#else
                for (size_t i = 0; i < nb; ++i)
                    f();
#endif
            }
        }
    }
}

#endif
