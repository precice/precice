#pragma once

#include <string>
#include <memory>

namespace precice
{
namespace xml
{

template <typename VALUE_T>
class Validator
{
public:
  Validator() = default;

  Validator(const Validator &other) = delete;

  Validator& operator=(const Validator &other) = delete;

  virtual ~Validator() {}

  virtual bool validateValue(const VALUE_T &value) = 0;

  virtual std::unique_ptr<Validator> clone() const = 0;

  virtual std::string getErrorMessage() const = 0;

  virtual std::string getDocumentation() const = 0;
};
}
} // namespace precice, xml
