#pragma once

#include "Validator.hpp"
#include "logging/Logger.hpp"

namespace precice
{
namespace xml
{

template <typename VALUE_T>
class ValidatorEquals;

template <typename VALUE_T>
class ValidatorEquals : public Validator<VALUE_T>
{
public:
  ValidatorEquals(VALUE_T valueToEqual)
      : Validator<VALUE_T>(),
        _valueToEqual(valueToEqual)
  {
  }

  virtual ~ValidatorEquals() {}

  virtual bool validateValue(const VALUE_T &value)
  {
    TRACE(value);
    return value == _valueToEqual;
  }

  virtual Validator<VALUE_T> &clone() const
  {
    ValidatorEquals<VALUE_T> *validator =
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
    stream << "'" << _valueToEqual << "'";
    return stream.str();
  }

private:
  logging::Logger _log{"xml::ValidatorEquals"};

  ValidatorEquals(const ValidatorEquals<VALUE_T> &rhs);

  ValidatorEquals<VALUE_T> &operator=(const ValidatorEquals<VALUE_T> &rhs);

  VALUE_T _valueToEqual;
};

}
} // namespace precice, xml
