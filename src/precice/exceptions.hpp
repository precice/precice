#pragma once

#include <stdexcept>

namespace precice {

class Error : public std::runtime_error {
public:
  Error(const std::string &what_arg)
      : std::runtime_error(what_arg){};
};

class ActionError : public precice::Error {
public:
  ActionError(const std::string &what_arg)
      : precice::Error(what_arg){};
};

class AccelerationError : public precice::Error {
public:
  AccelerationError(const std::string &what_arg)
      : precice::Error(what_arg){};
};

class CouplingSchemeError : public precice::Error {
public:
  CouplingSchemeError(const std::string &what_arg)
      : precice::Error(what_arg){};
};

class ExportError : public precice::Error {
public:
  ExportError(const std::string &what_arg)
      : precice::Error(what_arg){};
};

class MappingError : public precice::Error {
public:
  MappingError(const std::string &what_arg)
      : precice::Error(what_arg){};
};

class InterfaceError : public precice::Error {
public:
  InterfaceError(const std::string &what_arg)
      : precice::Error(what_arg){};
};

class ConfigurationError : public precice::Error {
public:
  ConfigurationError(const std::string &what_arg)
      : precice::Error(what_arg){};
};

} // namespace precice
