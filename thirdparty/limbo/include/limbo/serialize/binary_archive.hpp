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
#ifndef LIMBO_SERIALIZE_BINARY_ARCHIVE_HPP
#define LIMBO_SERIALIZE_BINARY_ARCHIVE_HPP

#include <cassert>
#include <iostream>
#include <sstream>
#include <string>

// Quick hack for definition of 'I' in <complex.h>
#undef I
#include <boost/filesystem.hpp>

#include <Eigen/Core>

namespace limbo {
    namespace serialize {

        class BinaryArchive {
        public:
            BinaryArchive(const std::string& dir_name) : _dir_name(dir_name) {}

            /// write an Eigen::Matrix*
            void save(const Eigen::MatrixXd& v, const std::string& object_name) const
            {
                _create_directory();

                std::ofstream out(fname(object_name).c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
                _write_binary(out, v);
                out.close();
            }

            /// write a vector of Eigen::Vector*
            template <typename T>
            void save(const std::vector<T>& v, const std::string& object_name) const
            {
                _create_directory();

                std::stringstream s;

                int size = v.size();
                s.write((char*)(&size), sizeof(int));
                for (auto& x : v) {
                    _write_binary(s, x);
                }

                std::ofstream out(fname(object_name).c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
                out << s.rdbuf();
                out.close();
            }

            /// load an Eigen matrix (or vector)
            template <typename M>
            void load(M& m, const std::string& object_name) const
            {
                std::ifstream in(fname(object_name).c_str(), std::ios::in | std::ios::binary);
                _read_binary(in, m);
                in.close();
            }

            /// load a vector of Eigen::Vector*
            template <typename V>
            void load(std::vector<V>& m_list, const std::string& object_name) const
            {
                m_list.clear();

                std::ifstream in(fname(object_name).c_str(), std::ios::in | std::ios::binary);

                int size;
                in.read((char*)(&size), sizeof(int));

                for (int i = 0; i < size; i++) {
                    V v;
                    _read_binary(in, v);
                    m_list.push_back(v);
                }
                in.close();
                assert(!m_list.empty());
            }

            std::string fname(const std::string& object_name) const
            {
                return _dir_name + "/" + object_name + ".bin";
            }

            const std::string& directory() const
            {
                return _dir_name;
            }

        protected:
            std::string _dir_name;

            void _create_directory() const
            {
                boost::filesystem::path my_path(_dir_name);
                boost::filesystem::create_directories(my_path);
            }

            template <class Matrix, class Stream>
            void _write_binary(Stream& out, const Matrix& matrix) const
            {
                typename Matrix::Index rows = matrix.rows(), cols = matrix.cols();
                out.write((char*)(&rows), sizeof(typename Matrix::Index));
                out.write((char*)(&cols), sizeof(typename Matrix::Index));
                out.write((char*)matrix.data(), rows * cols * sizeof(typename Matrix::Scalar));
            }

            template <class Matrix, class Stream>
            void _read_binary(Stream& in, Matrix& matrix) const
            {
                typename Matrix::Index rows = 0, cols = 0;
                in.read((char*)(&rows), sizeof(typename Matrix::Index));
                in.read((char*)(&cols), sizeof(typename Matrix::Index));
                matrix.resize(rows, cols);
                in.read((char*)matrix.data(), rows * cols * sizeof(typename Matrix::Scalar));
            }
        };
    } // namespace serialize
} // namespace limbo

#endif