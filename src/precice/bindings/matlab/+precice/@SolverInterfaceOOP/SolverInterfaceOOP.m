classdef SolverInterfaceOOP < precice.SolverInterface
    %SOLVERINTERFACE Matlab wrapper for the C++ SolverInterface
    %   This class will serve as a gateway for the precice Solver Interface
    %   in C++. It has access to a private mex Function written in the new
    %   MATLAB C++ API. Such a function is actually an instance of a C++ 
    %   class that is stored internally by MATLAB. Whenever they are 
    %   "called", a subroutine of the object (operator ()) is invoked which
    %   then has access to the members of the class.
    %   Hence, we can store the actual solver interface in this class and
    %   access it by invoking the mex function.
    
    properties(Access=private)
        % Possible changes:
        % - Allow the creation of multiple SolverInterfaces. For this, we
        % should replace the interface pointer in the Gateway by a list of 
        % interface pointers and add an InterfaceID as property of this
        % class.
        
        % mex host managing the out of process execution
        oMexHost;
        
        % bool set to true after oMexHost started running
        bMexHostRunning = 0;
    end
    
    methods
        %% Construction and configuration
        % Constructor
        function obj = SolverInterface(SolverName,configFileName,solverProcessIndex,solverProcessSize)
            %SOLVERINTERFACE Construct an instance of this class
            % Initialize the mex host
            obj.oMexHost = mexhost;
            obj.bMexHostRunning = true;
            
            if (solverProcessIndex > 0 || solverProcessSize > 1)
                error('Parallel runs are currently not supported with the MATLAB bindings.')    
            end
            if ischar(SolverName)
                SolverName = string(SolverName);
            end
            if ischar(configFileName)
                configFileName = string(configFileName);
            end
            feval(obj.oMexHost,"preciceGateway",uint8(0),SolverName,configFileName,int32(solverProcessIndex),int32(solverProcessSize));
        end
        
        % Destructor
        function delete(obj)
            % Delete the mex host
            feval(obj.oMexHost,"preciceGateway",uint8(1));
            delete(obj.oMexHost);
        end
        
        %% Steering methods
        % initialize
        function dt = initialize(obj)
            dt = feval(obj.oMexHost,"preciceGateway",uint8(10));
        end
        
        % initialize Data
        function initializeData(obj)
            feval(obj.oMexHost,"preciceGateway",uint8(11));
        end
        
        % advance
        function dt = advance(obj,dt)
            dt = feval(obj.oMexHost,"preciceGateway",uint8(12),dt);
        end
        
        % finalize
        function finalize(obj)
            feval(obj.oMexHost,"preciceGateway",uint8(13));
        end
        
        %% Status queries
        % getDimensions
        function dims = getDimensions(obj)
            dims = feval(obj.oMexHost,"preciceGateway",uint8(20));
        end
        
        % isCouplingOngoing
        function bool = isCouplingOngoing(obj)
            bool = feval(obj.oMexHost,"preciceGateway",uint8(21));
        end
        
        % isReadDataAvailable
        function bool = isReadDataAvailable(obj)
            bool = feval(obj.oMexHost,"preciceGateway",uint8(22));
        end
        
        % isWriteDataRequired
        function bool = isWriteDataRequired(obj,dt)
            bool = feval(obj.oMexHost,"preciceGateway",uint8(23),dt);
        end
        
        % isTimestepComplete
        function bool = isTimestepComplete(obj)
            bool = feval(obj.oMexHost,"preciceGateway",uint8(24));
        end
        
        % hasToEvaluateSurrogateModel
        function bool = hasToEvaluateSurrogateModel(obj)
            bool = feval(obj.oMexHost,"preciceGateway",uint8(25));
        end
        
        % hasToEvaluateFineModel
        function bool = hasToEvaluateFineModel(obj)
            bool = feval(obj.oMexHost,"preciceGateway",uint8(26));
        end
        
        %% Action Methods
        % isActionRequired
        function bool = isActionRequired(obj,action)
            if ischar(action)
                action = string(action);
            end
            bool = feval(obj.oMexHost,"preciceGateway",uint8(30),action);
        end
        
        % fulfilledAction
        function fulfilledAction(obj,action)
            if ischar(action)
                action = string(action);
            end
            feval(obj.oMexHost,"preciceGateway",uint8(31),action);
        end
        
        %% Mesh Access
        % hasMesh
        function bool = hasMesh(obj,meshName)
            if ischar(meshName)
                meshName = string(meshName);
            end
            bool = feval(obj.oMexHost,"preciceGateway",uint8(40),meshName);
        end
        
        % getMeshID
        function id = getMeshID(obj,meshName)
            if ischar(meshName)
                meshName = string(meshName);
            end
            id = feval(obj.oMexHost,"preciceGateway",uint8(41),meshName);
        end
        
        % getMeshIDs
        function ids = getMeshIDs(obj)
            ids = feval(obj.oMexHost,"preciceGateway",uint8(42));
        end
        
        % getMeshHandle not yet implemented
        
        % setMeshVertex
        function vertexId = setMeshVertex(obj,meshID,position)
            vertexId = feval(obj.oMexHost,"preciceGateway",uint8(44),int32(meshID),position);
        end
        
        % getMeshVertexSize
        function vertexId = getMeshVertexSize(obj,meshID)
            vertexId = feval(obj.oMexHost,"preciceGateway",uint8(45),int32(meshID));
        end
        
        % setMeshVertices
        function vertexIds = setMeshVertices(obj,meshID,positions)
            inSize = size(positions,2);
            vertexIds = feval(obj.oMexHost,"preciceGateway",uint8(46),int32(meshID),int32(inSize),positions);
        end
        
        % getMeshVertices
        function positions = getMeshVertices(obj,meshID,vertexIds)
            inSize = size(vertexIds,2);
            positions = feval(obj.oMexHost,"preciceGateway",uint8(47),int32(meshID),int32(inSize),vertexIds);
        end
        
        % getMeshVertexIDsFromPositions
        function vertexIds = getMeshVertexIDsFromPositions(obj,meshID,positions)
            inSize = size(positions,2);
            vertexIds = feval(obj.oMexHost,"preciceGateway",uint8(48),int32(meshID),int32(inSize),positions);
        end
        
        % setMeshEdge
        function edgeID = setMeshEdge(obj, meshID, firstVertexID, secondVertexID)
            edgeID = feval(obj.oMexHost,"preciceGateway",uint8(49),int32(meshID),int32(firstVertexID),int32(secondVertexID));
        end
        
        % setMeshTriangle
        function setMeshTriangle(obj, meshID, firstEdgeID, secondEdgeID, thirdEdgeID)
            feval(obj.oMexHost,"preciceGateway",uint8(50),int32(meshID),int32(firstEdgeID),int32(secondEdgeID),int32(thirdEdgeID));
        end
        
        % setMeshTriangleWithEdges
        function setMeshTriangleWithEdges(obj, meshID, firstVertexID, secondVertexID, thirdVertexID)
            feval(obj.oMexHost,"preciceGateway",uint8(51),int32(meshID),int32(firstVertexID),int32(secondVertexID),int32(thirdVertexID));
        end
        
        % setMeshQuad
        function setMeshQuad(obj, meshID, firstEdgeID, secondEdgeID, thirdEdgeID, fourthEdgeID)
            feval(obj.oMexHost,"preciceGateway",uint8(52),int32(meshID),int32(firstEdgeID),int32(secondEdgeID),int32(thirdEdgeID),int32(fourthEdgeID));
        end
        
        % setMeshQuadWithEdges
        function setMeshQuadWithEdges(obj, meshID, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID)
            feval(obj.oMexHost,"preciceGateway",uint8(53),int32(meshID),int32(firstVertexID),int32(secondVertexID),int32(thirdVertexID),int32(fourthVertexID));
        end
        
        %% Data Access
        % hasDataID
        function bool = hasData(obj,dataName,meshID)
            if ischar(dataName)
                dataName = string(dataName);
            end
            bool = feval(obj.oMexHost,"preciceGateway",uint8(60),dataName,int32(meshID));
        end
        
        % getDataID
        function id = getDataID(obj,dataName,meshID)
            if ischar(dataName)
                dataName = string(dataName);
            end
            id = feval(obj.oMexHost,"preciceGateway",uint8(61),dataName,int32(meshID));
        end
        
        % mapReadDataTo
        function mapReadDataTo(obj,meshID)
            feval(obj.oMexHost,"preciceGateway",uint8(62),int32(meshID));
        end
        
        % mapWriteDataFrom
        function mapWriteDataFrom(obj,meshID)
            feval(obj.oMexHost,"preciceGateway",uint8(63),int32(meshID));
        end
        
        % writeBlockVectorData
        function writeBlockVectorData(obj,dataID,valueIndices,values)
            if ~isa(valueIndices,'int32')
                warning('valueIndices should be allocated as int32 to prevent copying.');
                valueIndices = int32(valueIndices);
            end
            inSize = size(valueIndices, 2);
            feval(obj.oMexHost,"preciceGateway",uint8(64),int32(dataID),int32(inSize),valueIndices,values);
        end
        
        % writeVectorData
        function writeVectorData(obj,dataID,valueIndex,value)
            feval(obj.oMexHost,"preciceGateway",uint8(65),int32(dataID),int32(valueIndex),value);
        end
        
        % writeBlockScalarData
        function writeBlockScalarData(obj,dataID,valueIndices,values)
            if ~isa(valueIndices,'int32')
                warning('valueIndices should be allocated as int32 to prevent copying.');
                valueIndices = int32(valueIndices);
            end
            inSize = size(valueIndices, 2);
            feval(obj.oMexHost,"preciceGateway",uint8(66),int32(dataID),int32(inSize),valueIndices,values);
        end
        
        % writeScalarData
        function writeScalarData(obj,dataID,valueIndex,value)
            feval(obj.oMexHost,"preciceGateway",uint8(67),int32(dataID),int32(valueIndex),value);
        end
        
        % readBlockVectorData
        function values = readBlockVectorData(obj,dataID,valueIndices)
            if ~isa(valueIndices,'int32')
                warning('valueIndices should be allocated as int32 to prevent copying.');
                valueIndices = int32(valueIndices);
            end
            inSize = size(valueIndices, 2);
            values = feval(obj.oMexHost,"preciceGateway",uint8(68),int32(dataID),int32(inSize),valueIndices);
        end
        
        % readVectorData
        function value = readVectorData(obj,dataID,valueIndex)
            value = feval(obj.oMexHost,"preciceGateway",uint8(69),int32(dataID),int32(valueIndex));
        end
        
        % readBlockScalarData
        function values = readBlockScalarData(obj,dataID,valueIndices,transpose)
            if nargin<4
                transpose=false;
            end
            if ~isa(valueIndices,'int32')
                warning('valueIndices should be allocated as int32 to prevent copying.');
                valueIndices = int32(valueIndices);
            end
            inSize = size(valueIndices, 2);
            values = feval(obj.oMexHost,"preciceGateway",uint8(70),int32(dataID),int32(inSize),valueIndices,transpose);
        end
        
        % readScalarData
        function value = readScalarData(obj,dataID,valueIndex)
            value = feval(obj.oMexHost,"preciceGateway",uint8(71),int32(dataID),int32(valueIndex));
        end
    end
end
