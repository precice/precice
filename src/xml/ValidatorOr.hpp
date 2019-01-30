#pragma once

#include "Validator.hpp"

namespace precice
{
namespace xml
{

template <typename VALUE_T>
class ValidatorOr
    : public Validator<VALUE_T>
{
public:
  ValidatorOr(
      const std::unique_ptr<Validator<VALUE_T>> &lhs,
      const std::unique_ptr<Validator<VALUE_T>> &rhs)
:
    _lhs(lhs->clone()),
    _rhs(rhs->clone())
  {
  }

  ValidatorOr(
      std::unique_ptr<Validator<VALUE_T>> &&lhs,
      const std::unique_ptr<Validator<VALUE_T>> &rhs):
    _lhs(lhs),
    _rhs(rhs->clone())
  {
  }

  ValidatorOr(
      const std::unique_ptr<Validator<VALUE_T>> &lhs,
      std::unique_ptr<Validator<VALUE_T>> &&rhs):
    _lhs(lhs->clone()),
    _rhs(rhs)
  {
  }

  ValidatorOr(
      std::unique_ptr<Validator<VALUE_T>> &&lhs,
      std::unique_ptr<Validator<VALUE_T>> &&rhs):
    _lhs(lhs),
    _rhs(rhs)
  {
  }


  ~ValidatorOr() override {}

  ValidatorOr(const ValidatorOr &other) = delete;

  ValidatorOr &operator=(const ValidatorOr &other) = delete;

  bool validateValue(const VALUE_T &value) override
  {
    return _lhs->validateValue(value) || _rhs->validateValue(value);
  }

  std::unique_ptr<Validator<VALUE_T>> clone() const override
  {
    return std::unique_ptr<ValidatorOr>(new ValidatorOr(_lhs, _rhs));
  }

  std::string getErrorMessage() const override
  {
    return _lhs->getErrorMessage() + " or " + _rhs->getErrorMessage();
  }

  std::string getDocumentation() const override
  {
    return _lhs->getDocumentation() + " or " + _rhs->getDocumentation();
  }

private:
  std::unique_ptr<Validator<VALUE_T>> _lhs;
  std::unique_ptr<Validator<VALUE_T>> _rhs;
};

template <typename VALUE_T>
const std::unique_ptr<Validator<VALUE_T>> operator||(
    const std::unique_ptr<Validator<VALUE_T>> &lhs,
    const std::unique_ptr<Validator<VALUE_T>> &rhs)
{
    using VAL = ValidatorOr<VALUE_T>;
    return std::unique_ptr<VAL>(new VAL(lhs, rhs));
}

} // namespace xml
} // namespace precice
