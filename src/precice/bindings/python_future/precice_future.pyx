# distutils: language = c++

"""precice

The python module precice offers python language bindings to the C++ coupling library precice. Please refer to precice.org for further information.
"""

import numpy as np
cimport numpy as np
cimport cython
from mpi4py import MPI


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


cdef class Interface:
    # construction and configuration
    # constructor

    def __cinit__ (self, solver_name, solver_process_index, solver_process_size):
        self.thisptr = new SolverInterface.SolverInterface (convert(solver_name), solver_process_index, solver_process_size)
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
    def set_mesh_vertex(self, mesh_id, position):
        if not isinstance(position, np.ndarray):
            position = np.asarray(position)
        dimensions = position.size
        assert(dimensions == self.get_dimensions())
        cdef np.ndarray[double, ndim=1] _position = np.ascontiguousarray(position, dtype=np.double)
        vertex_id = self.thisptr.setMeshVertex(mesh_id, <const double*>_position.data)
        return vertex_id

    # returns the number of vertices of a mesh
    def get_mesh_vertex_size (self, mesh_id):
        return self.thisptr.getMeshVertexSize(mesh_id)

    # creates multiple mesh vertices
    def set_mesh_vertices (self, mesh_id, positions):
        if not isinstance(positions, np.ndarray):
            positions = np.asarray(positions)
        size, dimensions = positions.shape
        assert(dimensions == self.get_dimensions())
        cdef np.ndarray[double, ndim=1] _positions = np.ascontiguousarray(positions.flatten(), dtype=np.double)
        cdef np.ndarray[int, ndim=1] _ids = np.empty(size, dtype=np.int32)
        self.thisptr.setMeshVertices (mesh_id, size, <const double*>_positions.data, <int*>_ids.data)
        return _ids

    # get vertex positions for multiple vertex ids from a given mesh
    def get_mesh_vertices(self, mesh_id, ids):
        cdef np.ndarray[int, ndim=1] _ids = np.ascontiguousarray(ids, dtype=np.int32)
        size = _ids.size
        cdef np.ndarray[double, ndim=1] _positions = np.empty(size * self.get_dimensions(), dtype=np.double)
        self.thisptr.getMeshVertices (mesh_id, size, <const int*>_ids.data, <double*>_positions.data)
        return _positions.reshape((size, self.get_dimensions()))

    # gets mesh vertex IDs from positions
    def get_mesh_vertex_ids_from_positions (self, mesh_id, positions):
        if not isinstance(positions, np.ndarray):
            positions = np.asarray(positions)
        size, dimensions = positions.shape
        assert(dimensions == self.get_dimensions())
        cdef np.ndarray[double, ndim=1] _positions = np.ascontiguousarray(positions.flatten(), dtype=np.double)
        cdef np.ndarray[int, ndim=1] _ids = np.empty(int(size), dtype=np.int32)
        self.thisptr.getMeshVertexIDsFromPositions (mesh_id, size, <const double*>_positions.data, <int*>_ids.data)
        return _ids

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
    def set_mesh_quad (self, mesh_id, first_edge_id, second_edge_id, third_edge_id, fourth_edge_id):
        self.thisptr.setMeshQuad (mesh_id, first_edge_id, second_edge_id, third_edge_id, fourth_edge_id)

    # sets mesh quad from vertices
    def set_mesh_quad_with_edges (self, mesh_id, first_vertex_id, second_vertex_id, third_vertex_id, fourth_vertex_id):
        self.thisptr.setMeshQuadWithEdges (mesh_id, first_vertex_id, second_vertex_id, third_vertex_id, fourth_vertex_id)

    # data access
    # hasData
    def has_data (self, str data_name, mesh_id):
        return self.thisptr.hasData(convert(data_name), mesh_id)

    def get_data_id (self, str data_name, mesh_id):
        return self.thisptr.getDataID (convert(data_name), mesh_id)

    def map_read_data_to (self, to_mesh_id):
        self.thisptr.mapReadDataTo (to_mesh_id)

    def map_write_data_from (self, from_mesh_id):
        self.thisptr.mapWriteDataFrom (from_mesh_id)

    def write_block_vector_data (self, data_id, value_indices, values):
        if not isinstance(values, np.ndarray):
            values = np.asarray(values)
        size, dimensions = values.shape
        assert(dimensions == self.get_dimensions())
        cdef np.ndarray[int, ndim=1] _value_indices = np.ascontiguousarray(value_indices, dtype=np.int32)
        cdef np.ndarray[double, ndim=1] _values = np.ascontiguousarray(values.flatten(), dtype=np.double)
        assert(size == _value_indices.size)
        size = value_indices.size
        self.thisptr.writeBlockVectorData (data_id, size, <const int*>_value_indices.data, <const double*>_values.data)

    def write_vector_data (self, data_id, value_index, value):
        if not isinstance(value, np.ndarray):
            value = np.asarray(value)
        dimensions = value.size
        assert(dimensions == self.get_dimensions())
        cdef np.ndarray[np.double_t, ndim=1] _value = np.ascontiguousarray(value, dtype=np.double)
        self.thisptr.writeVectorData (data_id, value_index, <const double*>_value.data)

    def write_block_scalar_data (self, data_id, value_indices, values):
        cdef np.ndarray[int, ndim=1] _value_indices = np.ascontiguousarray(value_indices, dtype=np.int32)
        cdef np.ndarray[double, ndim=1] _values = np.ascontiguousarray(values, dtype=np.double)
        assert(_values.size == _value_indices.size)
        size = value_indices.size
        self.thisptr.writeBlockScalarData (data_id, size, <const int*>_value_indices.data, <const double*>_values.data)

    def write_scalar_data (self, data_id, value_index, double value):
        self.thisptr.writeScalarData (data_id, value_index, value)

    def read_block_vector_data (self, data_id, value_indices):
        cdef np.ndarray[int, ndim=1] _value_indices = np.ascontiguousarray(value_indices, dtype=np.int32)
        size = _value_indices.size
        dimensions = self.get_dimensions()
        cdef np.ndarray[np.double_t, ndim=1] _values = np.empty(size * dimensions, dtype=np.double)
        self.thisptr.readBlockVectorData (data_id, size, <const int*>_value_indices.data, <double*>_values.data)
        return _values.reshape((size, dimensions))

    def read_vector_data (self, data_id, value_index):
        dimensions = self.get_dimensions()
        cdef np.ndarray[double, ndim=1] _value = np.empty(dimensions, dtype=np.double)
        self.thisptr.readVectorData (data_id, value_index, <double*>_value.data)
        return _value

    def read_block_scalar_data (self, data_id, value_indices):
        cdef np.ndarray[int, ndim=1] _value_indices = np.ascontiguousarray(value_indices, dtype=np.int32)
        size = _value_indices.size
        cdef np.ndarray[double, ndim=1] _values = np.empty(size, dtype=np.double)
        self.thisptr.readBlockScalarData (data_id, size, <const int*>_value_indices.data, <double*>_values.data)
        return _values

    def read_scalar_data (self, data_id, value_index):
        cdef double _value
        self.thisptr.readScalarData (data_id, value_index, _value)
        return _value

def action_write_initial_data ():
    return SolverInterface.actionWriteInitialData()
   
def action_write_iteration_checkpoint ():
    return SolverInterface.actionWriteIterationCheckpoint()

def action_read_iteration_checkpoint ():
    return SolverInterface.actionReadIterationCheckpoint()
