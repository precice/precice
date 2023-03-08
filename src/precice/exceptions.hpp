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

class CommunicationError : public precice::Error {
public:
  CommunicationError(const std::string &what_arg)
      : precice::Error(what_arg){};
};

class IOError : public precice::Error {
public:
  IOError(const std::string &what_arg)
      : precice::Error(what_arg){};
};

class MappingError : public precice::Error {
public:
  MappingError(const std::string &what_arg)
      : precice::Error(what_arg){};
};

class MeshError : public precice::Error {
public:
  MeshError(const std::string &what_arg)
      : precice::Error(what_arg){};
};

class M2NError : public precice::Error {
public:
  M2NError(const std::string &what_arg)
      : precice::Error(what_arg){};
};

class APIError : public precice::Error {
public:
  APIError(const std::string &what_arg)
      : precice::Error(what_arg){};
};

class ConfigurationError : public precice::Error {
public:
  ConfigurationError(const std::string &what_arg)
      : precice::Error(what_arg){};
};

class PartitionError : public precice::Error {
public:
  PartitionError(const std::string &what_arg)
      : precice::Error(what_arg){};
};

} // namespace precice
