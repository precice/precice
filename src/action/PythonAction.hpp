// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_ACTION_PYTHONACTION_HPP_
#define PRECICE_ACTION_PYTHONACTION_HPP_

#include "action/Action.hpp"
#include "mesh/SharedPointer.hpp"
#include "tarch/logging/Log.h"
#include <string>

struct _object;
typedef _object PyObject;

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace action {

/**
 * @brief Action whose implementation is given in a Python file.
 */
class PythonAction : public Action
{
public:

  PythonAction (
    Timing               timing,
    const std::string&   modulePath,
    const std::string&   moduleName,
    const mesh::PtrMesh& mesh,
    int                  targetDataID,
    int                  sourceDataID );

  virtual ~PythonAction();

  virtual void performAction (
    double time,
    double dt,
    double computedPartFullDt,
    double fullDt );

private:

  static tarch::logging::Log _log;

  std::string _modulePath;

  std::string _moduleName;

  mesh::PtrData _targetData;

  mesh::PtrData _sourceData;

  int _numberArguments;

  bool _isInitialized;

  PyObject* _moduleNameObject;

  PyObject* _module;

  PyObject* _sourceValues;

  PyObject* _targetValues;

  PyObject* _performAction;

  PyObject* _vertexCallback;

  PyObject* _postAction;

  void initialize();

  int makeNumPyArraysAvailable();
};

}} // namespace precice, action

#endif /* PRECICE_ACTION_PYTHONACTION_HPP_ */
