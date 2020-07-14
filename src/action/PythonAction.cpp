#ifndef PRECICE_NO_PYTHON
#include "PythonAction.hpp"
#include <Eigen/Core>
#include <Python.h>
#include <memory>
#include <numpy/arrayobject.h>
#include <ostream>
#include "logging/LogMacros.hpp"
#include "mesh/Data.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace action {

PythonAction::PythonAction(
    Timing               timing,
    const std::string &  modulePath,
    const std::string &  moduleName,
    const mesh::PtrMesh &mesh,
    int                  targetDataID,
    int                  sourceDataID)
    : Action(timing, mesh),
      _modulePath(modulePath),
      _moduleName(moduleName)
{
  if (targetDataID != -1) {
    _targetData = getMesh()->data(targetDataID);
    _numberArguments++;
  }
  if (sourceDataID != -1) {
    _sourceData = getMesh()->data(sourceDataID);
    _numberArguments++;
  }
}

PythonAction::~PythonAction()
{
  if (_module != nullptr) {
    PRECICE_ASSERT(_moduleNameObject != nullptr);
    PRECICE_ASSERT(_module != nullptr);
    Py_DECREF(_moduleNameObject);
    Py_DECREF(_module);
    Py_Finalize();
  }
}

void PythonAction::performAction(double time,
                                 double dt,
                                 double computedPartFullDt,
                                 double fullDt)
{
  PRECICE_TRACE(time, dt, computedPartFullDt, fullDt);

  if (not _isInitialized)
    initialize();

  PyObject *dataArgs = PyTuple_New(_numberArguments);
  if (_performAction != nullptr) {
    PyObject *pythonTime = PyFloat_FromDouble(time);
    PyObject *pythonDt   = PyFloat_FromDouble(fullDt);
    PyTuple_SetItem(dataArgs, 0, pythonTime);
    PyTuple_SetItem(dataArgs, 1, pythonDt);
    if (_sourceData) {
      npy_intp sourceDim[]  = {_sourceData->values().size()};
      double * sourceValues = _sourceData->values().data();
      //PRECICE_ASSERT(_sourceValues == NULL);
      _sourceValues = PyArray_SimpleNewFromData(1, sourceDim, NPY_DOUBLE, sourceValues);
      PRECICE_CHECK(_sourceValues != nullptr, "Creating python source values failed!");
      PyTuple_SetItem(dataArgs, 2, _sourceValues);
    }
    if (_targetData) {
      npy_intp targetDim[]  = {_targetData->values().size()};
      double * targetValues = _targetData->values().data();
      //PRECICE_ASSERT(_targetValues == NULL);
      _targetValues =
          PyArray_SimpleNewFromData(1, targetDim, NPY_DOUBLE, targetValues);
      PRECICE_CHECK(_targetValues != nullptr, "Creating python target values failed!");
      int argumentIndex = _sourceData ? 3 : 2;
      PyTuple_SetItem(dataArgs, argumentIndex, _targetValues);
    }
    PyObject_CallObject(_performAction, dataArgs);
    if (PyErr_Occurred()) {
      PyErr_Print();
      PRECICE_ERROR("Error occurred during call of function "
                    << "performAction() python module \"" << _moduleName << "\"!");
    }
  }

  if (_vertexCallback != nullptr) {
    PyObject *      vertexArgs = PyTuple_New(3);
    mesh::PtrMesh   mesh       = getMesh();
    Eigen::VectorXd coords(mesh->getDimensions());
    Eigen::VectorXd normal(mesh->getDimensions());
    for (mesh::Vertex &vertex : mesh->vertices()) {
      npy_intp vdim[]        = {mesh->getDimensions()};
      int      id            = vertex.getID();
      coords                 = vertex.getCoords();
      normal                 = vertex.getNormal();
      PyObject *pythonID     = PyLong_FromLong(id);
      PyObject *pythonCoords = PyArray_SimpleNewFromData(1, vdim, NPY_DOUBLE, coords.data());
      PyObject *pythonNormal = PyArray_SimpleNewFromData(1, vdim, NPY_DOUBLE, coords.data());
      PRECICE_CHECK(pythonID != nullptr, "Creating python ID failed!");
      PRECICE_CHECK(pythonCoords != nullptr, "Creating python coords failed!");
      PRECICE_CHECK(pythonNormal != nullptr, "Creating python normal failed!");
      PyTuple_SetItem(vertexArgs, 0, pythonID);
      PyTuple_SetItem(vertexArgs, 1, pythonCoords);
      PyTuple_SetItem(vertexArgs, 2, pythonNormal);
      PyObject_CallObject(_vertexCallback, vertexArgs);
      if (PyErr_Occurred()) {
        PyErr_Print();
        PRECICE_ERROR("Error occurred during call of function "
                      << "vertexCallback() python module \"" << _moduleName << "\"!");
      }
    }
    Py_DECREF(vertexArgs);
  }

  if (_postAction != nullptr) {
    PyObject *postActionArgs = PyTuple_New(0);
    PyObject_CallObject(_postAction, postActionArgs);
    if (PyErr_Occurred()) {
      PyErr_Print();
      PRECICE_ERROR("Error occurred during call of function "
                    << "postAction() in python module \"" << _moduleName << "\"!");
    }
    Py_DECREF(postActionArgs);
  }

  Py_DECREF(dataArgs);
}

void PythonAction::initialize()
{
  PRECICE_ASSERT(not _isInitialized);
  // Initialize Python
  Py_Initialize();
  makeNumPyArraysAvailable();
  // Append execution path to find module to import
  PyRun_SimpleString("import sys");
  std::string appendPathCommand("sys.path.append('" + _modulePath + "')");
  PyRun_SimpleString(appendPathCommand.c_str());
  _moduleNameObject = PyUnicode_FromString(_moduleName.c_str());
  _module           = PyImport_Import(_moduleNameObject);
  if (_module == nullptr) {
    PyErr_Print();
    PRECICE_ERROR("Could not load python module \"" << _moduleName << "\" at path \"" << _modulePath << "\"!");
  }

  // Construct method performAction
  _performAction = PyObject_GetAttrString(_module, "performAction");
  if (PyErr_Occurred()) {
    PyErr_Clear();
    PRECICE_WARN("No function void performAction() in python module \"" << _moduleName << "\" found.");
    _performAction = nullptr;
  }
  //  bool valid = _performAction != NULL;
  //  if (valid) valid = PyCallable_Check(_performAction);
  //  if (not valid){
  //  }

  // Construct method vertexCallback
  _vertexCallback = PyObject_GetAttrString(_module, "vertexCallback");
  if (PyErr_Occurred()) {
    PyErr_Clear();
    PRECICE_WARN("No function void vertexCallback() in python module \"" << _moduleName << "\" found.");
    _vertexCallback = nullptr;
  }

  // Construct function postAction
  _postAction = PyObject_GetAttrString(_module, "postAction");
  if (PyErr_Occurred()) {
    PyErr_Clear();
    PRECICE_WARN("No function void postAction() in python module \"" << _moduleName << "\" found.");
    _postAction = nullptr;
  }
}

int PythonAction::makeNumPyArraysAvailable()
{
  static bool importedAlready = false;
  if (importedAlready)
    return 0;
  import_array1(-1); // this macro is defined be NumPy and must be included
  importedAlready = true;
  return 1;
}

} // namespace action
} // namespace precice

#endif
