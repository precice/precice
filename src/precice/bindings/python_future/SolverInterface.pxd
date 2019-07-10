from libcpp        cimport bool
from libcpp.set    cimport set
from libcpp.string cimport string


cdef extern from "precice/SolverInterface.hpp" namespace "precice":
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

        void writeBlockVectorData (const int dataID, const int size, const int* valueIndices, const double* values)

        void writeVectorData (const int dataID, const int valueIndex, const double* value)

        void writeBlockScalarData (const int dataID, const int size, const int* valueIndices, const double* values)

        void writeScalarData (const int dataID, const int valueIndex, const double value)

        void readBlockVectorData (const int dataID, const int size, const int* valueIndices, double* values)

        void readVectorData (const int dataID, const int valueIndex, double* value)

        void readBlockScalarData (const int dataID, const int size, const int* valueIndices, double* values)

        void readScalarData (const int dataID, const int valueIndex, double& value)


cdef extern from "precice/SolverInterface.hpp"  namespace "precice::constants":
    const string& actionWriteInitialData()
    const string& actionWriteIterationCheckpoint()
    const string& actionReadIterationCheckpoint()
