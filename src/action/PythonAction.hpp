#pragma once
#ifndef PRECICE_NO_PYTHON

#include <string>
#include "action/Action.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"

struct _object;
using PyObject = _object;

namespace precice::action {

/// Action whose implementation is given in a Python file.
class PythonAction : public Action {
public:
  PythonAction(
      Timing               timing,
      std::string          modulePath,
      std::string          moduleName,
      const mesh::PtrMesh &mesh,
      int                  targetDataID,
      int                  sourceDataID);

  ~PythonAction() override;

  void performAction() final override;

private:
  logging::Logger _log{"action::PythonAction"};

  std::string _modulePath;

  std::string _moduleName;

  mesh::PtrData _targetData;

  mesh::PtrData _sourceData;

  int _numberArguments = 1;

  bool _isInitialized = false;

  PyObject *_moduleNameObject = nullptr;

  PyObject *_module = nullptr;

  PyObject *_sourceValues = nullptr;

  PyObject *_targetValues = nullptr;

  PyObject *_performAction = nullptr;

  void initialize();

  int makeNumPyArraysAvailable();
};

} // namespace precice::action

#endif
