// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_NO_PYTHON
#include <Python.h>
#include <arrayobject.h>
#endif
#include "PythonAction.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Data.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace action {

tarch::logging::Log PythonAction:: _log ( "precice::action::PythonAction" );

PythonAction:: PythonAction
(
  Timing               timing,
  const std::string&   modulePath,
  const std::string&   moduleName,
  const mesh::PtrMesh& mesh,
  int                  targetDataID,
  int                  sourceDataID )
:
  Action(timing, mesh),
  _modulePath(modulePath),
  _moduleName(moduleName),
  _targetData(),
  _sourceData(),
  _numberArguments(1),
  _isInitialized(false),
  _moduleNameObject(NULL),
  _module(NULL),
  _sourceValues(NULL),
  _targetValues(NULL),
  _performAction(NULL),
  _vertexCallback(NULL),
  _postAction(NULL)
{
  if (targetDataID != -1){
    _targetData = getMesh()->data(targetDataID);
    _numberArguments++;
  }
  if (sourceDataID != -1){
    _sourceData = getMesh()->data(sourceDataID);
    _numberArguments++;
  }
}

PythonAction:: ~PythonAction()
{
# ifndef PRECICE_NO_PYTHON
  if (_module != NULL){
    assertion(_moduleNameObject != NULL);
    assertion(_module != NULL);
    Py_DECREF(_moduleNameObject);
    Py_DECREF(_module);
    Py_Finalize();
  }
# endif // not PRECICE_NO_PYTHON
}

void PythonAction:: performAction
(
  double time,
  double dt,
  double computedPartFullDt,
  double fullDt)
{
# ifndef PRECICE_NO_PYTHON
  preciceTrace4("performAction()", time, dt, computedPartFullDt, fullDt);

  if (not _isInitialized) initialize();

  PyObject* dataArgs = PyTuple_New(_numberArguments);
  if (_performAction != NULL){
    PyObject* pythonTime = PyFloat_FromDouble(time);
    PyTuple_SetItem(dataArgs, 0, pythonTime);
    if (_sourceData.use_count() > 0){
      npy_intp sourceDim[] = { _sourceData->values().size() };
      double* sourceValues = tarch::la::raw(_sourceData->values());
      //assertion(_sourceValues == NULL);
      _sourceValues =
          PyArray_SimpleNewFromData(1, sourceDim, NPY_DOUBLE, sourceValues);
      preciceCheck(_sourceValues != NULL, "PythonAction()",
                   "Creating python source values failed!");
      PyTuple_SetItem(dataArgs, 1, _sourceValues);
    }
    if (_targetData.use_count() > 0){
      npy_intp targetDim[] = { _targetData->values().size() };
      double* targetValues = tarch::la::raw(_targetData->values());
      //assertion(_targetValues == NULL);
      _targetValues =
          PyArray_SimpleNewFromData(1, targetDim, NPY_DOUBLE, targetValues);
      preciceCheck(_targetValues != NULL, "PythonAction()",
                   "Creating python target values failed!");
      int argumentIndex = _sourceData.use_count() > 0 ? 2 : 1;
      PyTuple_SetItem(dataArgs, argumentIndex, _targetValues);
    }
    PyObject_CallObject(_performAction, dataArgs);
    if(PyErr_Occurred()){
      PyErr_Print ();
      preciceError("performAction()", "Error occurred during call of function "
                     << "performAction() python module \"" << _moduleName << "\"!");
    }
  }

  if (_vertexCallback != NULL){
    PyObject* vertexArgs = PyTuple_New(3);
    mesh::PtrMesh mesh = getMesh();
    utils::DynVector coords(mesh->getDimensions());
    utils::DynVector normal(mesh->getDimensions());
    foreach (mesh::Vertex& vertex, mesh->vertices()){
      npy_intp vdim[] = { mesh->getDimensions() };
      int id = vertex.getID();
      coords = vertex.getCoords();
      normal = vertex.getNormal();
      using tarch::la::raw;
      PyObject* pythonID = PyInt_FromLong(id);
      PyObject* pythonCoords =
          PyArray_SimpleNewFromData(1, vdim, NPY_DOUBLE, raw(coords));
      PyObject* pythonNormal =
          PyArray_SimpleNewFromData(1, vdim, NPY_DOUBLE, raw(normal));
      preciceCheck(pythonID != NULL, "performAction()",
                     "Creating python ID failed!");
      preciceCheck(pythonCoords != NULL, "performAction()",
                     "Creating python coords failed!");
      preciceCheck(pythonNormal != NULL, "performAction()",
                     "Creating python normal failed!");
      PyTuple_SetItem(vertexArgs, 0, pythonID);
      PyTuple_SetItem(vertexArgs, 1, pythonCoords);
      PyTuple_SetItem(vertexArgs, 2, pythonNormal);
      PyObject_CallObject(_vertexCallback, vertexArgs);
      if (PyErr_Occurred()){
        PyErr_Print ();
        preciceError("performAction()", "Error occurred during call of function "
                       << "vertexCallback() python module \"" << _moduleName << "\"!");
      }
    }
    Py_DECREF(vertexArgs);
  }

  if (_postAction != NULL){
    PyObject* postActionArgs = PyTuple_New(0);
    PyObject_CallObject(_postAction, postActionArgs);
    if(PyErr_Occurred()){
      PyErr_Print ();
      preciceError("performAction()", "Error occurred during call of function "
                   << "postAction() in python module \"" << _moduleName << "\"!");
    }
    Py_DECREF(postActionArgs);
  }

  Py_DECREF(dataArgs);

  //if (_sourceValues != NULL){
  //  Py_DECREF(_sourceValues);
  //  _sourceValues = NULL;
  //}
  //if (_targetValues != NULL){
  //  Py_DECREF(_targetValues);
  //  _targetValues = NULL;
  //}

  //Py_Initialize();
  //makeNumPyArraysAvailable();
  //// Append execution path to find module to import
  //PyRun_SimpleString("import sys");
  //std::string appendPathCommand("sys.path.append('" + _modulePath + "')");
  //PyRun_SimpleString(appendPathCommand.c_str());
  //PyObject* moduleName = PyString_FromString(_moduleName.c_str());
  //PyObject* module = PyImport_Import(moduleName);
  //if (module == NULL){
  //  PyErr_Print();
  //  preciceError("performAction()", "Could not load python module \""
  //                 << _moduleName << "\" at path \"" << _modulePath << "\"!");
  //}

//  PyObject* dataArgs = PyTuple_New(_numberArguments);
//  if (_sourceData.use_count() > 0){
//    npy_intp sourceDim[] = { _sourceData->values().size() };
//    double* sourceValues = tarch::la::raw(_sourceData->values());
//    PyObject* pythonSourceValues =
//        PyArray_SimpleNewFromData(1, sourceDim, NPY_DOUBLE, sourceValues);
//    preciceCheck(pythonSourceValues != NULL, "performAction()",
//                   "Creating python source values failed!");
//    PyTuple_SetItem(dataArgs, 0, pythonSourceValues);
//  }
//  if (_targetData.use_count() > 0){
//    npy_intp targetDim[] = { _targetData->values().size() };
//    double* targetValues = tarch::la::raw(_targetData->values());
//    PyObject* pythonTargetValues =
//        PyArray_SimpleNewFromData(1, targetDim, NPY_DOUBLE, targetValues);
//    preciceCheck(pythonTargetValues != NULL, "performAction()",
//                   "Creating python target values failed!");
//    int argumentIndex = _sourceData.use_count() > 0 ? 1 : 0;
//    PyTuple_SetItem(dataArgs, argumentIndex, pythonTargetValues);
//  }

//  PyObject* function = PyObject_GetAttrString(module, "performAction");
//  if ((function != NULL) && PyCallable_Check(function)){
//    PyObject_CallObject(function, dataArgs);
//    if(PyErr_Occurred()){
//      PyErr_Print ();
//      preciceError("performAction()", "Error occurred during call of python module \""
//                     << _moduleName << "\"!");
//    }
//  }
//  else {
//    preciceWarning("performAction()", "No function void performAction() in python module \""
//                     << _moduleName << "\" found.");
//  }

//  function = PyObject_GetAttrString(module, "vertexCallback");
//  if ((function != NULL) && PyCallable_Check(function)){
//    PyObject* vertexArgs = PyTuple_New(3);
//    mesh::PtrMesh mesh = getMesh();
//    utils::DynVector coords(mesh->getDimensions());
//    utils::DynVector normal(mesh->getDimensions());
//    foreach (mesh::Vertex& vertex, mesh->vertices()){
//      npy_intp vdim[] = { mesh->getDimensions() };
//      int id = vertex.getID();
//      coords = vertex.getCoords();
//      normal = vertex.getNormal();
//      using tarch::la::raw;
//      PyObject * pythonID = PyInt_FromLong(id);
//      PyObject * pythonCoords =
//          PyArray_SimpleNewFromData(1, vdim, NPY_DOUBLE, raw(coords));
//      PyObject * pythonNormal =
//          PyArray_SimpleNewFromData(1, vdim, NPY_DOUBLE, raw(normal));
//      preciceCheck(pythonID != NULL, "performAction()",
//                     "Creating python ID failed!");
//      preciceCheck(pythonCoords != NULL, "performAction()",
//                     "Creating python coords failed!");
//      preciceCheck(pythonNormal != NULL, "performAction()",
//                     "Creating python normal failed!");
//      PyTuple_SetItem(vertexArgs, 0, pythonID);
//      PyTuple_SetItem(vertexArgs, 1, pythonCoords);
//      PyTuple_SetItem(vertexArgs, 2, pythonNormal);
//      PyObject_CallObject(function, vertexArgs);
//      if (PyErr_Occurred()){
//        PyErr_Print ();
//        preciceError("performAction()", "Error occurred during call of function "
//                       << "vertexCallback() python module \"" << _moduleName << "\"!");
//      }
//    }
//    Py_DECREF(vertexArgs);
//  }
//  else {
//    preciceWarning("performAction()", "No function void vertexCallback() in python module \""
//                     << _moduleName << "\" found.");
//  }

//  function = PyObject_GetAttrString(module, "postAction");
//  if ((function != NULL) && PyCallable_Check(function)){
//    PyObject* postActionArgs = PyTuple_New(0);
//    PyObject_CallObject(function, postActionArgs);
//    if(PyErr_Occurred()){
//      PyErr_Print ();
//      preciceError("performAction()", "Error occurred during call of python module \""
//                     << _moduleName << "\"!");
//    }
//    Py_DECREF(postActionArgs);
//  }
//  else {
//    preciceWarning("performAction()", "No function void postAction() in python module \""
//                     << _moduleName << "\" found.");
//  }

//  Py_DECREF(dataArgs);
//  Py_DECREF(moduleName);
//  Py_DECREF(module);
//  Py_Finalize();
# endif // not PRECICE_NO_PYTHON
}

void PythonAction:: initialize()
{
# ifndef PRECICE_NO_PYTHON
  assertion(not _isInitialized);
  // Initialize Python
  Py_Initialize();
  makeNumPyArraysAvailable();
  // Append execution path to find module to import
  PyRun_SimpleString("import sys");
  std::string appendPathCommand("sys.path.append('" + _modulePath + "')");
  PyRun_SimpleString(appendPathCommand.c_str());
  _moduleNameObject = PyString_FromString(_moduleName.c_str());
  _module = PyImport_Import(_moduleNameObject);
  if (_module == NULL){
    PyErr_Print();
    preciceError("PythonAction()", "Could not load python module \""
                   << _moduleName << "\" at path \"" << _modulePath << "\"!");
  }

  // Construct method performAction
  _performAction = PyObject_GetAttrString(_module, "performAction");
  if (PyErr_Occurred()){
    PyErr_Clear();
    preciceWarning("PythonAction()", "No function void performAction() in python module \""
                   << _moduleName << "\" found.");
    _performAction = NULL;
  }
//  bool valid = _performAction != NULL;
//  if (valid) valid = PyCallable_Check(_performAction);
//  if (not valid){
//  }

  // Construct method vertexCallback
  _vertexCallback = PyObject_GetAttrString(_module, "vertexCallback");
  if (PyErr_Occurred()){
    PyErr_Clear();
    preciceWarning("PythonAction()", "No function void vertexCallback() in python module \""
                   << _moduleName << "\" found.");
    _vertexCallback = NULL;
  }

  // Construct function postAction
  _postAction = PyObject_GetAttrString(_module, "postAction");
  if (PyErr_Occurred()){
    PyErr_Clear();
    preciceWarning("PythonAction()", "No function void postAction() in python module \""
                   << _moduleName << "\" found.");
    _postAction = NULL;
  }
# endif // not PRECICE_NO_PYTHON
}

int PythonAction:: makeNumPyArraysAvailable()
{
# ifndef PRECICE_NO_PYTHON
  static bool importedAlready = false;
  if (importedAlready) return 0;
  import_array1(-1); // this macro is defined be NumPy and must be included
  importedAlready = true;
# endif
  return 1;
}

}} // namespace precice, action
