"""precice

The python module precice offers python language bindings to the C++ coupling library precice. Please refer to precice.org for further information.
"""

import numpy as np
cimport numpy as np
cimport cython
from libcpp        cimport bool
from libcpp.set    cimport set
from libcpp.string cimport string
from mpi4py import MPI

from libc.stdlib   cimport free
from libc.stdlib cimport malloc

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

cdef extern from "precice/SolverInterface.hpp"  namespace "precice":
    cdef cppclass SolverInterface:
        # construction and configuration

        SolverInterface (const string&, int, int) except +

        void configure (const string&)

        # steering methods

        double initialize ()

        void initializeData ()

        double advance (double computedTimestepLength)

        void finalize()

        # status queries

        int getDimensions() const

        bool isCouplingOngoing()

        bool isReadDataAvailable()

        bool isWriteDataRequired (double computedTimestepLength)

        bool isTimestepComplete()

        bool hasToEvaluateSurrogateModel ()

        bool hasToEvaluateFineModel ()

        # action methods

        bool isActionRequired (const string& action)

        void fulfilledAction (const string& action)

        # mesh access

        bool hasMesh (const string& meshName ) const

        int getMeshID (const string& meshName)

        set[int] getMeshIDs ()

        # MeshHandle getMeshHandle (const string& meshName)

        int setMeshVertex (int meshID, const double* position)

        int getMeshVertexSize (int meshID)

        void setMeshVertices (int meshID, int size, const double* positions, int* ids)

        void getMeshVertices (int meshID, int size, const int* ids, double* positions)

        void getMeshVertexIDsFromPositions (int meshID, int size, double* positions, int* ids)

        int setMeshEdge (int meshID, int firstVertexID, int secondVertexID)

        void setMeshTriangle (int meshID, int firstEdgeID, int secondEdgeID, int thirdEdgeID)

        void setMeshTriangleWithEdges (int meshID, int firstVertexID, int secondVertexID, int thirdVertexID)

        void setMeshQuad (int meshID, int firstEdgeID, int secondEdgeID, int thirdEdgeID, int fourthEdgeID)

        void setMeshQuadWithEdges (int meshID, int firstVertexID, int secondVertexID, int thirdVertexID, int fourthVertexID)

        # data access

        bool hasData (const string& dataName, int meshID) const

        int getDataID (const string& dataName, int meshID)

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


cdef extern from "precice/SolverInterface.hpp"  namespace "precice::constants":
    const string& actionWriteInitialData()
    const string& actionWriteIterationCheckpoint()
    const string& actionReadIterationCheckpoint()


cdef class Interface:
    cdef SolverInterface *thisptr # hold a C++ instance being wrapped

    # construction and configuration
    # constructor

    def __cinit__ (self, solver_name, int solver_process_index, int solver_process_size):
        self.thisptr = new SolverInterface (convert(solver_name), solver_process_index, solver_process_size)
        pass

    # destructor
    def __dealloc__ (self):
        del self.thisptr

    # configure
    def configure (self, configuration_file_name):
        self.thisptr.configure (convert(configuration_file_name))

    # steering methods
    # initialize
    def initialize (self):
        return self.thisptr.initialize ()

    # initialize data
    def initialize_data (self):
        self.thisptr.initializeData ()

    # advance in time
    def advance (self, double computed_timestep_length):
        return self.thisptr.advance (computed_timestep_length)

    # finalize preCICE
    def finalize (self):
        self.thisptr.finalize ()

    # status queries
    # get dimensions
    def get_dimensions (self):
        return self.thisptr.getDimensions ()

    # check if coupling is going on
    def is_coupling_ongoing (self):
        return self.thisptr.isCouplingOngoing ()

    # check if data is available to be read
    def is_read_data_available (self):
        return self.thisptr.isReadDataAvailable ()

    # check if write data is needed
    def is_write_data_required (self, double computed_timestep_length):
        return self.thisptr.isWriteDataRequired (computed_timestep_length)

    # check if time-step is complete
    def is_timestep_complete (self):
        return self.thisptr.isTimestepComplete ()

    # returns whether the solver has to evaluate the surrogate model representation
    def has_to_evaluate_surrogate_model (self):
        return self.thisptr.hasToEvaluateSurrogateModel ()

    # checks if the solver has to evaluate the fine model representation
    def has_to_evaluate_fine_model (self):
        return self.thisptr.hasToEvaluateFineModel ()

    # action methods
    # check if action is needed
    def is_action_required (self, action):
        return self.thisptr.isActionRequired (action)

    # notify of action being fulfilled
    def fulfilled_action (self, action):
        self.thisptr.fulfilledAction (action)

    # mesh access
    # hasMesh
    def has_mesh(self, mesh_name):
        return self.thisptr.hasMesh (convert(mesh_name))

    # get mesh ID
    def get_mesh_id (self, mesh_name):
        return self.thisptr.getMeshID (convert(mesh_name))

    # get mesh IDs
    def get_mesh_ids (self):
        return self.thisptr.getMeshIDs ()

    # returns a handle to a created mesh
    def get_mesh_handle(self, mesh_name):
        raise Exception("The API method get_mesh_handle is not yet available for the Python bindings.")

    # creates a mesh vertex
    def set_mesh_vertex(self, mesh_id, np.ndarray[np.double_t, ndim=1] position):
        vertex_id = self.thisptr.setMeshVertex(mesh_id, &position[0])
        return vertex_id

    # returns the number of vertices of a mesh
    def get_mesh_vertex_size (self, mesh_id):
        return self.thisptr.getMeshVertexSize(mesh_id)

    # creates multiple mesh vertices
    def set_mesh_vertices (self, mesh_id, np.ndarray[np.double_t, ndim=1] positions):
        size = positions.size/self.get_dimensions()
        assert(size.is_integer())
        size = int(size)
        cdef np.ndarray[int] ids = np.empty(int(size), dtype=np.int32)
        self.thisptr.setMeshVertices (mesh_id, size, &positions[0], &ids[0])
        return np.array(ids)

    # get vertex positions for multiple vertex ids from a given mesh
    def get_mesh_vertices(self, mesh_id, ids):
        cdef np.ndarray[int, ndim=1] ids_ = np.array(ids, dtype=np.int32)
        size = ids_.size
        cdef np.ndarray[double] positions = np.empty(size * self.get_dimensions(), dtype=np.double)
        self.thisptr.getMeshVertices (mesh_id, size, &ids_[0], &positions[0])
        return np.array(positions)

    # gets mesh vertex IDs from positions
    def get_mesh_vertex_ids_from_positions (self, mesh_id, np.ndarray[np.double_t, ndim=1] positions):
        size = positions.size/self.get_dimensions()
        assert(size.is_integer())
        cdef np.ndarray[int] ids = np.empty(int(size), dtype=np.int32)

        self.thisptr.getMeshVertexIDsFromPositions (mesh_id, size, &positions[0], &ids[0])

        return np.array(ids)

    # sets mesh edge from vertex IDs, returns edge ID
    def set_mesh_edge (self, mesh_id, first_vertex_id, second_vertex_id):
        return self.thisptr.setMeshEdge (mesh_id, first_vertex_id, second_vertex_id)

    # sets mesh triangle from edges
    def set_mesh_triangle (self, mesh_id, first_edge_id, second_edge_id, third_edge_id):
        self.thisptr.setMeshTriangle (mesh_id, first_edge_id, second_edge_id, third_edge_id)

    # sets mesh triangle from vertices
    def set_mesh_triangle_with_edges (self, mesh_id, first_vertex_id, second_vertex_id, third_vertex_id):
        self.thisptr.setMeshTriangleWithEdges (mesh_id, first_vertex_id, second_vertex_id, third_vertex_id)

    # sets mesh quad from edges
    def set_mesh_quad (self, int mesh_id, first_edge_id, second_edge_id, third_edge_id, fourth_edge_id):
        self.thisptr.setMeshQuad (mesh_id, first_edge_id, second_edge_id, third_edge_id, fourth_edge_id)

    # sets mesh quad from vertices
    def set_mesh_quad_with_edges (self, mesh_id, first_vertex_id, second_vertex_id, third_vertex_id, fourth_vertex_id):
        self.thisptr.setMeshQuadWithEdges (mesh_id, first_vertex_id, second_vertex_id, third_vertex_id, fourth_vertex_id)

    # data access
    # hasData
    def has_data (self, data_name, mesh_id):
        return self.thisptr.hasData(convert(data_name), mesh_id)

    def get_data_id (self, data_name, mesh_id):
        return self.thisptr.getDataID (convert(data_name), mesh_id)

    def map_read_data_to (self, to_mesh_id):
        self.thisptr.mapReadDataTo (to_mesh_id)

    def map_write_data_from (self, from_mesh_id):
        self.thisptr.mapWriteDataFrom (from_mesh_id)

    def write_block_vector_data (self, data_id, value_indices, np.ndarray[np.double_t, ndim=1] values):
        cdef np.ndarray[int, ndim=1] value_indices_ = np.array(value_indices, dtype=np.int32)
        size = value_indices_.size
        self.thisptr.writeBlockVectorData (data_id, size, &value_indices_[0], &values[0])

    def write_vector_data (self, data_id, value_index, np.ndarray[np.double_t, ndim=1] value):
        self.thisptr.writeVectorData (data_id, value_index, &value[0])

    def write_block_scalar_data (self, data_id, value_indices, np.ndarray[np.double_t, ndim=1] values):
        cdef np.ndarray[int, ndim=1] value_indices_ = np.array(value_indices, dtype=np.int32)
        size = value_indices_.size
        self.thisptr.writeBlockScalarData (data_id, size, &value_indices_[0], &values[0])

    def write_scalar_data (self, data_id, value_index, value):
        self.thisptr.writeScalarData (data_id, value_index, value)

    def read_block_vector_data (self, data_id, value_indices):
        cdef np.ndarray[int, ndim=1] value_indices_ = np.array(value_indices, dtype=np.int32)
        size = value_indices_.size
        cdef np.ndarray[np.double_t] values = np.empty(size * self.get_dimensions(), dtype=np.double)
        self.thisptr.readBlockVectorData (data_id, size, &value_indices_[0], &values[0])
        return values

    def read_vector_data (self, data_id, value_index):
        cdef np.ndarray[np.double_t] value = np.empty(self.get_dimensions(), dtype=np.double)
        self.thisptr.readVectorData (data_id, value_index, &value[0])
        return value

    def read_block_scalar_data (self, data_id, value_indices):
        cdef np.ndarray[int, ndim=1] value_indices_ = np.array(value_indices, dtype=np.int32)
        size = value_indices_.size
        cdef np.ndarray[np.double_t] values = np.empty(size, dtype=np.double)
        self.thisptr.readBlockScalarData (data_id, size, &value_indices_[0], &values[0])
        return values

    def read_scalar_data (self, data_id, value_index):
        cdef double value;
        self.thisptr.readScalarData (data_id, value_index, value)
        return value

def action_write_initial_data ():
    return actionWriteInitialData()
   
def action_write_iteration_checkpoint ():
    return actionWriteIterationCheckpoint()

def action_read_iteration_checkpoint ():
    return actionReadIterationCheckpoint()
