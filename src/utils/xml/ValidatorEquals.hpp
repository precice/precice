// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_UTILS_VALIDATOREQUALS_HPP_
#define PRECICE_UTILS_VALIDATOREQUALS_HPP_

#include "Validator.hpp"
#include "tarch/logging/Log.h"
#include "utils/Helpers.hpp"

namespace precice {
namespace utils {

template<typename VALUE_T>
class ValidatorEquals;


template<typename VALUE_T>
class ValidatorEquals : public Validator<VALUE_T>
{
public:

   ValidatorEquals ( VALUE_T valueToEqual )
   :
     Validator<VALUE_T>(),
     _valueToEqual(valueToEqual)
   {}

   virtual ~ValidatorEquals() {}

   virtual bool validateValue ( const VALUE_T& value )
   {
     preciceTrace1("validateValue()", value);
     return value == _valueToEqual;
   }

   virtual Validator<VALUE_T>& clone() const
   {
      ValidatorEquals<VALUE_T>* validator =
          new ValidatorEquals<VALUE_T>(_valueToEqual);
      return *validator;
   }

   virtual std::string getErrorMessage() const
   {
      std::ostringstream stream;
      stream << _valueToEqual;
      return std::string("value must be \"" + stream.str() + "\"");
   }

   virtual std::string getDocumentation() const
   {
     std::ostringstream stream;
     stream << "'" << _valueToEqual  << "'";
     return stream.str();
   }

private:

   static tarch::logging::Log _log;

   ValidatorEquals ( const ValidatorEquals<VALUE_T>& rhs );

   ValidatorEquals<VALUE_T>& operator= ( const ValidatorEquals<VALUE_T>& rhs );

   VALUE_T _valueToEqual;
};

template<typename VALUE_T>
tarch::logging::Log precice::utils::ValidatorEquals<VALUE_T>::
   _log = tarch::logging::Log ("precice::utils::ValidatorEquals");

}} // namespace precice, utils

#endif /* PRECICE_UTILS_VALIDATOREQUALS_HPP_ */
