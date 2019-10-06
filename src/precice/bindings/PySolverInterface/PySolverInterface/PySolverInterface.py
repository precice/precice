from precice import Interface


class PySolverInterface:
    def __init__(self, solverName, solverProcessIndex, solverProcessSize):
        self.interface = Interface(solverName, solverProcessIndex, solverProcessSize)

    def configure(self, configurationFileName):
        self.interface.configure(configurationFileName)

    def initialize(self):
        return self.interface.initialize()

    def initializeData(self):
        self.interface.initialize_data()

    def advance(self, computedTimestepLength):
        return self.interface.advance(computedTimestepLength)

    def finalize(self):
        self.interface.finalize()

    def getDimensions(self):
        return self.interface.get_dimensions()

    def isCouplingOngoing(self):
        return self.interface.is_coupling_ongoing()

    def isReadDataAvailable(self):
        return self.interface.is_read_data_available()

    def isWriteDataRequired(self, computedTimestepLength):
        return self.interface.is_write_data_required(computedTimestepLength)

    def isTimestepComplete(self):
        return self.interface.is_timestep_complete()

    def isActionRequired(self, action):
        return self.interface.is_action_required(action)

    def fulfilledAction(self, action):
        self.interface.fulfilled_action(action)

    def hasMesh(self, meshName):
        return self.interface.has_mesh(meshName)

    def getMeshID(self, meshName):
        return self.interface.get_mesh_id(meshName)

    def getMeshIDs(self):
        return self.interface.get_mesh_ids()

    def hasData(self, dataName, meshID):
        return self.interface.has_data(dataName, meshID)

    def getDataID(self, dataName, meshID):
        return self.interface.get_data_id(dataName, meshID)

    def setMeshVertices(self, meshID, size, positions, ids):
        self.interface.set_mesh_vertices(meshID, size, positions, ids)

    def getMeshVertexSize(self, meshID):
        return self.interface.get_mesh_vertex_size(meshID)

    def getMeshVertexIDsFromPositions(self, meshID, size, positions, ids):
        self.interface.get_mesh_vertex_ids_from_positions(meshID, size, positions, ids)

    def setMeshEdge(self, meshID, firstVertexID, secondVertexID):
        return self.interface.set_mesh_edge(meshID, firstVertexID, secondVertexID)

    def setMeshTriangle(self, meshID, firstEdgeID, secondEdgeID, thirdEdgeID):
        self.interface.set_mesh_triangle(meshID, firstEdgeID, secondEdgeID, thirdEdgeID)

    def setMeshTriangleWithEdges(self, meshID, firstVertexID, secondVertexID, thirdVertexID):
        self.interface.set_mesh_triangle_with_edges(meshID, firstVertexID, secondVertexID, thirdVertexID)

    def setMeshQuad(self, meshID, firstEdgeID, secondEdgeID, thirdEdgeID, fourthEdgeID):
        self.interface.set_mesh_quad(meshID, firstEdgeID, secondEdgeID, thirdEdgeID, fourthEdgeID)

    def setMeshQuadWithEdges(self, meshID, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID):
        self.interface.set_mesh_quad_with_edges(meshID, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID)

    def mapReadDataTo(self, toMeshID):
        self.interface.map_read_data_to(toMeshID)

    def mapWriteDataFrom(self, fromMeshID):
        self.interface.map_write_data_from(fromMeshID)

    def writeBlockVectorData(self, dataID, size, valueIndices, values):
        self.interface.write_block_vector_data(dataID, size, valueIndices, values)

    def writeVectorData(self, dataID, valueIndex, value):
        self.interface.write_vector_data(dataID, valueIndex, value)

    def writeBlockScalarData(self, dataID, size, valueIndices, values):
        self.interface.write_block_scalar_data(dataID, size, valueIndices, values)

    def writeScalarData(self, dataID, valueIndex, value):
        self.interface.write_scalar_data(dataID, valueIndex, value)

    def readBlockVectorData(self, dataID, size, valueIndices, values):
        self.interface.read_block_vector_data(dataID, size, valueIndices, values)

    def readVectorData(self, dataID, valueIndex, value):
        self.interface.read_vector_data(dataID, valueIndex, value)

    def readBlockScalarData(self, dataID, size, valueIndices, values):
        self.interface.read_block_scalar_data(dataID, size, valueIndices, values)

    def readScalarData(self, dataID, valueIndex, value):
        self.interface.read_scalar_data(dataID, valueIndex, value)

