#pragma once
#ifndef PRECICE_NO_PYTHON

#include <string>
#include "action/Action.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"

struct _object;
using PyObject = _object;

namespace precice {
namespace action {

/// Action whose implementation is given in a Python file.
class PythonAction : public Action {
public:
  PythonAction(
      Timing               timing,
      const std::string &  modulePath,
      const std::string &  moduleName,
      const mesh::PtrMesh &mesh,
      int                  targetDataID,
      int                  sourceDataID);

  virtual ~PythonAction();

  virtual void performAction(
      double time,
      double dt,
      double computedPartFullDt,
      double fullDt);

private:
  logging::Logger _log{"action::PythonAction"};

  std::string _modulePath;

  std::string _moduleName;

  mesh::PtrData _targetData;

  mesh::PtrData _sourceData;

  int _numberArguments = 2;

  bool _isInitialized = false;

  PyObject *_moduleNameObject = nullptr;

  PyObject *_module = nullptr;

  PyObject *_sourceValues = nullptr;

  PyObject *_targetValues = nullptr;

  PyObject *_performAction = nullptr;

  PyObject *_vertexCallback = nullptr;

  PyObject *_postAction = nullptr;

  void initialize();

  int makeNumPyArraysAvailable();
};

} // namespace action
} // namespace precice

#endif
