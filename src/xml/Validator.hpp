#pragma once

#include <string>

namespace precice
{
namespace xml
{

template <typename VALUE_T>
class Validator
{
public:
  Validator() {}

  virtual ~Validator() {}

  virtual bool validateValue(const VALUE_T &value) = 0;

  virtual Validator<VALUE_T> &clone() const = 0;

  virtual std::string getErrorMessage() const = 0;

  virtual std::string getDocumentation() const = 0;

private:
  Validator(const Validator<VALUE_T> &rhs);

  Validator<VALUE_T> &operator=(const Validator<VALUE_T> &rhs);
};
}
} // namespace precice, xml
