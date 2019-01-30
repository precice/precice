#pragma once

#include "Validator.hpp"
#include "logging/Logger.hpp"

namespace precice
{
namespace xml
{

template <typename VALUE_T>
class ValidatorEquals : public Validator<VALUE_T>
{
public:
  ValidatorEquals(VALUE_T valueToEqual)
      : Validator<VALUE_T>(),
        _valueToEqual(valueToEqual)
  {
  }

  ValidatorEquals(const ValidatorEquals &other) = delete;

  ValidatorEquals& operator=(const ValidatorEquals &other) = delete;

  ~ValidatorEquals() override {}

  bool validateValue(const VALUE_T &value) override
  {
    TRACE(value);
    return value == _valueToEqual;
  }

  std::unique_ptr<Validator<VALUE_T>> clone() const override
  {
    return std::unique_ptr<ValidatorEquals>(new ValidatorEquals(_valueToEqual));
  }

  std::string getErrorMessage() const override
  {
    std::ostringstream stream;
    stream << _valueToEqual;
    return std::string("value must be \"" + stream.str() + "\"");
  }

  std::string getDocumentation() const override
  {
    std::ostringstream stream;
    stream << "'" << _valueToEqual << "'";
    return stream.str();
  }

private:
  logging::Logger _log{"xml::ValidatorEquals"};

  VALUE_T _valueToEqual;
};

template <typename VALUE_T>
std::unique_ptr<Validator<VALUE_T>> makeValidatorEquals(const VALUE_T& value)
{
    using VAL = ValidatorEquals<VALUE_T>;
    return std::unique_ptr<VAL>(new VAL(value));
}

template <typename VALUE_T>
std::unique_ptr<Validator<VALUE_T>> makeValidatorEquals(VALUE_T&& value)
{
    using VAL = ValidatorEquals<VALUE_T>;
    return std::unique_ptr<VAL>(new VAL(std::forward<VALUE_T>(value)));
}

} // namespace xml
} // namespace precice
