// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_UTILS_VALIDATOR_HPP_
#define PRECICE_UTILS_VALIDATOR_HPP_

#include <string>

namespace precice {
namespace utils {


template<typename VALUE_T>
class Validator
{
public:

   Validator() {}

   virtual ~Validator() {}

   virtual bool validateValue ( const VALUE_T& value ) =0;

   virtual Validator<VALUE_T>& clone() const =0;

   virtual std::string getErrorMessage() const =0;

   virtual std::string getDocumentation() const =0;

private:

   Validator ( const Validator<VALUE_T>& rhs );

   Validator<VALUE_T>& operator= ( const Validator<VALUE_T>& rhs );
};

}} // namespace precice, utils

#endif /* PRECICE_UTILS_VALIDATOR_HPP_ */
