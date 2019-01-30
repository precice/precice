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
      const Validator<VALUE_T> &lhs,
      const Validator<VALUE_T> &rhs)
  {
    _lhs = &lhs.clone();
    _rhs = &rhs.clone();
  }

  ~ValidatorOr() override {};

  ValidatorOr(const ValidatorOr &other) = delete;

  ValidatorOr &operator=(const ValidatorOr &other) = delete;

  bool validateValue(const VALUE_T &value) override
  {
    return _lhs->validateValue(value) || _rhs->validateValue(value);
  }

  ValidatorPtr &clone() const override
  {
    return std::unique_ptr<ValidatorOr>(*_lhs, *_rhs)
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
  ValidatorPtr _lhs;
  ValidatorPtr _rhs;
};

template <typename VALUE_T>
const ValidatorPtr<VALUE_T> operator||(
    const Validator<VALUE_T> &lhs,
    const Validator<VALUE_T> &rhs)
{
  return std::unique_ptr<ValidatorOr<VALUE_T>>(lhs, rhs);
}

} // namespace precice, xml
