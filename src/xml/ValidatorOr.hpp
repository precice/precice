#pragma once

#include "Validator.hpp"

namespace precice {
  namespace xml {

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
namespace xml {

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

}} // namespace precice, xml
