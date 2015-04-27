// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_UTILS_VALIDATOROR_HPP_
#define PRECICE_UTILS_VALIDATOROR_HPP_

#include "Validator.hpp"

namespace precice {
  namespace utils {

    template<typename VALUE_T>
    class ValidatorOr;

    template<typename VALUE_T>
    const Validator<VALUE_T>& operator||
    (
      const Validator<VALUE_T>& lhs,
      const Validator<VALUE_T>& rhs );
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace utils {

template<typename VALUE_T>
class ValidatorOr
:
   public Validator<VALUE_T>
{
public:

   ValidatorOr (
     const Validator<VALUE_T>& lhs,
     const Validator<VALUE_T>& rhs );

   virtual ~ValidatorOr()
   {
     delete _lhs;
     delete _rhs;
   }

   virtual bool validateValue ( const VALUE_T& value );

   virtual Validator<VALUE_T>& clone() const;

   virtual std::string getErrorMessage() const;

   virtual std::string getDocumentation() const;

private:

   ValidatorOr ( const ValidatorOr<VALUE_T> & rhs );

   ValidatorOr<VALUE_T> & operator= ( const ValidatorOr<VALUE_T> & rhs );

   Validator<VALUE_T> * _lhs;

   Validator<VALUE_T> * _rhs;
};

#include "ValidatorOr.cpph"

}} // namespace precice, utils

#endif /* PRECICE_UTILS_VALIDATOROR_HPP_ */
