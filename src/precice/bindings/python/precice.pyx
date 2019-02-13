from cpython       cimport array
from libcpp        cimport bool
from libc.stdlib   cimport free
from libc.stdlib   cimport malloc
from libcpp.set    cimport set
from libcpp.memory cimport shared_ptr
from libcpp.string cimport string
from libcpp.memory cimport unique_ptr
from libcpp.vector cimport vector

from cpython.version cimport PY_MAJOR_VERSION  # important for determining python version in order to properly normalize string input. See http://docs.cython.org/en/latest/src/tutorial/strings.html#general-notes-about-c-strings and https://github.com/precice/precice/issues/68 . 

cdef bytes convert(s):
    """
    source code from http://docs.cython.org/en/latest/src/tutorial/strings.html#general-notes-about-c-strings
    """
    if type(s) is bytes:
        return s
    elif type(s) is str:
        return s.encode()
    else:
        raise TypeError("Could not convert.")

include "constants.pyx"

cdef extern from "./src/precice/SolverInterface.hpp"  namespace "precice":
   cdef cppclass SolverInterface:
      SolverInterface (const string&, int, int) except +

      void configure (const string&)

      double initialize ()

      void initializeData ()

      double advance (double computedTimestepLength)

      void finalize()

      int getDimensions() const

      bool isCouplingOngoing()

      bool isReadDataAvailable()

      bool isWriteDataRequired (double computedTimestepLength)

      bool isTimestepComplete()

      bool isActionRequired (const string& action)

      void fulfilledAction (const string& action)

      bool hasMesh (const string& meshName ) const

      int getMeshID (const string& meshName)

      set[int] getMeshIDs ()

      bool hasData (const string& dataName, int meshID) const

      int getDataID (const string& dataName, int meshID)

      void setMeshVertices (int meshID, int size, double* positions, int* ids)

      int getMeshVertexSize(int meshID)

      void getMeshVertexIDsFromPositions (int meshID, int size, double* positions, int* ids)

      int setMeshEdge (int meshID, int firstVertexID, int secondVertexID)

      void setMeshTriangle (int meshID, int firstEdgeID, int secondEdgeID, int thirdEdgeID)

      void setMeshTriangleWithEdges (int meshID, int firstVertexID, int secondVertexID, int thirdVertexID)

      void setMeshQuad (int meshID, int firstEdgeID, int secondEdgeID, int thirdEdgeID, int fourthEdgeID)

      void setMeshQuadWithEdges (int meshID, int firstVertexID, int secondVertexID, int thirdVertexID, int fourthVertexID)

      void mapReadDataTo (int toMeshID)

      void mapWriteDataFrom (int fromMeshID)

      void writeBlockVectorData (int dataID, int size, int* valueIndices, double* values)

      void writeVectorData (int dataID, int valueIndex, const double* value)

      void writeBlockScalarData (int dataID, int size, int* valueIndices, double* values)

      void writeScalarData (int dataID, int valueIndex, double value)

      void readBlockVectorData (int dataID, int size, int* valueIndices, double* values)

      void readVectorData (int dataID, int valueIndex, double* value)

      void readBlockScalarData (int dataID, int size, int* valueIndices, double* values)

      void readScalarData (int dataID, int valueIndex, double& value)

#      MeshHandle getMeshHandle (const string& meshName)


cdef class SolverInterface:
   cdef SolverInterface *thisptr # hold a C++ instance being wrapped
   
   # constructor

   def __cinit__ (self, solverName, int solverProcessIndex, int solverProcessSize):
      self.thisptr = new SolverInterface (convert(solverName), solverProcessIndex, solverProcessSize)      
      pass

   # destructor
   def __dealloc__ (self):
      del self.thisptr

   # configure
   def configure (self, configurationFileName):
      self.thisptr.configure (convert(configurationFileName))

   # initialize
   def initialize (self):
      return self.thisptr.initialize ()

   # initialize data
   def initializeData (self):
      self.thisptr.initializeData ()

   # advance in time
   def advance (self, double computedTimestepLength):
      return self.thisptr.advance (computedTimestepLength)

   # finalize preCICE
   def finalize (self):
      self.thisptr.finalize ()

   # get dimensions
   def getDimensions (self):
      return self.thisptr.getDimensions ()

   # check if coupling is going on
   def isCouplingOngoing (self):
      return self.thisptr.isCouplingOngoing ()

   # check if data is available to be read
   def isReadDataAvailable (self):
      return self.thisptr.isReadDataAvailable ()

   # check if write data is needed
   def isWriteDataRequired (self, double computedTimestepLength):
      return self.thisptr.isWriteDataRequired (computedTimestepLength)

   # check if time-step is complete
   def isTimestepComplete (self):
      return self.thisptr.isTimestepComplete ()

   # check if action is needed
   def isActionRequired (self, action):
      return self.thisptr.isActionRequired (action)

   # notify of action being fulfilled
   def fulfilledAction (self, action):
      self.thisptr.fulfilledAction (action)

   # hasMesh
   def hasMesh(self, meshName):
      return self.thisptr.hasMesh (convert(meshName))

   # get mesh ID
   def getMeshID (self, meshName):
      return self.thisptr.getMeshID (convert(meshName))

   # get mesh IDs
   def getMeshIDs (self):
      return self.thisptr.getMeshIDs ()

   # hasData
   def hasData (self, dataName, meshID):
      return self.thisptr.hasData(convert(dataName), meshID)

   def getDataID (self, dataName, meshID):
      return self.thisptr.getDataID (convert(dataName), meshID)

   def setMeshVertices (self, meshID, size, positions, ids):
      cdef int* ids_
      cdef double* positions_
      ids_ = <int*> malloc(len(ids) * sizeof(int))
      positions_ = <double*> malloc(len(positions) * sizeof(double))

      if ids_ is NULL or positions_ is NULL:
         raise MemoryError()
      
      for i in xrange(len(ids)):
         ids_[i] = ids[i]
      for i in xrange(len(positions)):
         positions_[i] = positions[i]

      self.thisptr.setMeshVertices (meshID, size, positions_, ids_)

      for i in xrange(len(ids)):
         ids[i] = ids_[i]
      for i in xrange(len(positions)):
         positions[i] = positions_[i]

      free(ids_)
      free(positions_)

   def getMeshVertexSize (self, meshID):
      return self.thisptr.getMeshVertexSize(meshID)

   def getMeshVertexIDsFromPositions (self, meshID, size, positions, ids):
      cdef int* ids_
      cdef double* positions_

      ids_ = <int*> malloc(len(ids) * sizeof(int))
      positions_ = <double*> malloc(len(positions) * sizeof(double))

      if ids_ is NULL or positions_ is NULL:
         raise MemoryError()
      
      for i in xrange(len(ids)):
         ids_[i] = ids[i]
      for i in xrange(len(positions)):
         positions_[i] = positions[i]

      self.thisptr.getMeshVertexIDsFromPositions (meshID, size, positions_, ids_)

      for i in xrange(len(ids)):
         ids[i] = ids_[i]
      for i in xrange(len(positions)):
         positions[i] = positions_[i]

      free(ids_)
      free(positions_)

   def setMeshEdge (self, meshID, firstVertexID, secondVertexID):
      return self.thisptr.setMeshEdge (meshID, firstVertexID, secondVertexID)

   def setMeshTriangle (self, meshID, firstEdgeID, secondEdgeID, thirdEdgeID):
      self.thisptr.setMeshTriangle (meshID, firstEdgeID, secondEdgeID, thirdEdgeID)

   def setMeshTriangleWithEdges (self, meshID, firstVertexID, secondVertexID, thirdVertexID):
      self.thisptr.setMeshTriangleWithEdges (meshID, firstVertexID, secondVertexID, thirdVertexID)

   def setMeshQuad (self, int meshID, firstEdgeID, secondEdgeID, thirdEdgeID, fourthEdgeID):
      self.thisptr.setMeshQuad (meshID, firstEdgeID, secondEdgeID, thirdEdgeID, fourthEdgeID)

   def setMeshQuadWithEdges (self, meshID, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID):
      self.thisptr.setMeshQuadWithEdges (meshID, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID)

   def mapReadDataTo (self, toMeshID):
      self.thisptr.mapReadDataTo (toMeshID)

   def mapWriteDataFrom (self, fromMeshID):
      self.thisptr.mapWriteDataFrom (fromMeshID)

   def writeBlockVectorData (self, dataID, size, valueIndices, values):
      cdef int* valueIndices_
      cdef double* values_
      valueIndices_ = <int*> malloc(len(valueIndices) * sizeof(int))
      values_ = <double*> malloc(len(values) * sizeof(double))

      if valueIndices_ is NULL or values_ is NULL:
         raise MemoryError()
      
      for i in xrange(len(valueIndices)):
         valueIndices_[i] = valueIndices[i]
      for i in xrange(len(values)):
         values_[i] = values[i]

      self.thisptr.writeBlockVectorData (dataID, size, valueIndices_, values_)

      for i in xrange(len(valueIndices)):
         valueIndices[i] = valueIndices_[i]
      for i in xrange(len(values)):
         values[i] = values_[i]

      free(valueIndices_)
      free(values_)

   def writeVectorData (self, dataID, valueIndex, value):
      cdef double* value_
      value_ = <double*> malloc(len(value) * sizeof(double))

      if value_ is NULL:
         raise MemoryError()

      for i in xrange(len(value)):
         value_[i] = value[i]

      self.thisptr.writeVectorData (dataID, valueIndex, value_)

      for i in xrange(len(value)):
         value[i] = value_[i]

      free(value_)

   def writeBlockScalarData (self, dataID, size, valueIndices, values):
      cdef int* valueIndices_
      cdef double* values_
      valueIndices_ = <int*> malloc(len(valueIndices) * sizeof(int))
      values_ = <double*> malloc(len(values) * sizeof(double))

      if valueIndices_ is NULL or values_ is NULL:
         raise MemoryError()
      
      for i in xrange(len(valueIndices)):
         valueIndices_[i] = valueIndices[i]
      for i in xrange(len(values)):
         values_[i] = values[i]

      self.thisptr.writeBlockScalarData (dataID, size, valueIndices_, values_)

      for i in xrange(len(valueIndices)):
         valueIndices[i] = valueIndices_[i]
      for i in xrange(len(values)):
         values[i] = values_[i]

      free(valueIndices_)
      free(values_)

   def writeScalarData (self, dataID, valueIndex, value):
      self.thisptr.writeScalarData (dataID, valueIndex, value)

   def readBlockVectorData (self, dataID, size, valueIndices, values):
      cdef int* valueIndices_
      cdef double* values_
      valueIndices_ = <int*> malloc(len(valueIndices) * sizeof(int))
      values_ = <double*> malloc(len(values) * sizeof(double))

      if valueIndices_ is NULL or values_ is NULL:
         raise MemoryError()
      
      for i in xrange(len(valueIndices)):
         valueIndices_[i] = valueIndices[i]
      for i in xrange(len(values)):
         values_[i] = values[i]

      self.thisptr.readBlockVectorData (dataID, size, valueIndices_, values_)

      for i in xrange(len(valueIndices)):
         valueIndices[i] = valueIndices_[i]
      for i in xrange(len(values)):
         values[i] = values_[i]

      free(valueIndices_)
      free(values_)

   def readVectorData (self, dataID, valueIndex, value):
      cdef double* value_
      value_ = <double*> malloc(len(value) * sizeof(double))

      if value_ is NULL:
         raise MemoryError()

      for i in xrange(len(value)):
         value_[i] = value[i]

      self.thisptr.readVectorData (dataID, valueIndex, value_)

      for i in xrange(len(value)):
         value[i] = value_[i]

      free(value_)

   def readBlockScalarData (self, dataID, size, valueIndices, values):
      cdef int* valueIndices_
      cdef double* values_
      valueIndices_ = <int*> malloc(len(valueIndices) * sizeof(int))
      values_ = <double*> malloc(len(values) * sizeof(double))

      if valueIndices_ is NULL or values_ is NULL:
         raise MemoryError()
      
      for i in xrange(len(valueIndices)):
         valueIndices_[i] = valueIndices[i]
      for i in xrange(len(values)):
         values_[i] = values[i]

      self.thisptr.readBlockScalarData (dataID, size, valueIndices_, values_)

      for i in xrange(len(valueIndices)):
         valueIndices[i] = valueIndices_[i]
      for i in xrange(len(values)):
         values[i] = values_[i]

      free(valueIndices_)
      free(values_)

   def readScalarData (self, int dataID, int valueIndex, double& value):
      self.thisptr.readScalarData (dataID, valueIndex, value)

#   def getMeshHandle (self, meshName):
#      return self.thisptr.getMeshHandle (meshName)
   
