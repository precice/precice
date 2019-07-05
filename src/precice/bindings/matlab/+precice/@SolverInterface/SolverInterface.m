classdef SolverInterface < handle
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
        function obj = SolverInterface(SolverName,OutOfProcess)
            %SOLVERINTERFACE Construct an instance of this class
            if nargin < 2
                OutOfProcess = false;
            end
            
            if OutOfProcess
                % Initialize the mex host
                obj.oMexHost = mexhost;
                obj.bMexHostRunning = true;
            end
            
            if ischar(SolverName)
                SolverName = string(SolverName);
            end
            if (obj.bMexHostRunning)
                feval(obj.oMexHost,"preciceGateway",uint8(0),SolverName);
            else
                preciceGateway(uint8(0),SolverName);
            end
        end
        
        % Destructor
        function delete(obj)
            if (obj.bMexHostRunning)
                feval(obj.oMexHost,"preciceGateway",uint8(1));
                delete(obj.oMexHost);
            end
            % Delete the mex host
            
        end
        
        % configure
        function configure(obj,configFileName)
            if ischar(configFileName)
                configFileName = string(configFileName);
            end
            if (obj.bMexHostRunning)
                feval(obj.oMexHost,"preciceGateway",uint8(2),configFileName);
            else
                preciceGateway(uint8(2),configFileName);
            end
        end
        
        %% Steering methods
        % initialize
        function dt = initialize(obj)
            if (obj.bMexHostRunning)
                dt = feval(obj.oMexHost,"preciceGateway",uint8(10));
            else
                dt = preciceGateway(uint8(10));
            end
        end
        
        % initialize Data
        function initializeData(obj)
            if (obj.bMexHostRunning)
                feval(obj.oMexHost,"preciceGateway",uint8(11));
            else
                preciceGateway(uint8(11));
            end
        end
        
        % advance
        function dt = advance(obj,dt)
            if (obj.bMexHostRunning)
                dt = feval(obj.oMexHost,"preciceGateway",uint8(12),dt);
            else
                dt = preciceGateway(uint8(12),dt);
            end
        end
        
        % finalize
        function finalize(obj)
            if (obj.bMexHostRunning)
                feval(obj.oMexHost,"preciceGateway",uint8(13));
            else
                preciceGateway(uint8(13));
            end
        end
        
        %% Status queries
        % getDimensions
        function dims = getDimensions(obj)
            if (obj.bMexHostRunning)
                dims = feval(obj.oMexHost,"preciceGateway",uint8(20));
            else
                dims = preciceGateway(uint8(20));
            end
        end
        
        % isCouplingOngoing
        function bool = isCouplingOngoing(obj)
            if (obj.bMexHostRunning)
                bool = feval(obj.oMexHost,"preciceGateway",uint8(21));
            else
                bool = preciceGateway(uint8(21));
            end
        end
        
        % isReadDataAvailable
        function bool = isReadDataAvailable(obj)
            if (obj.bMexHostRunning)
                bool = feval(obj.oMexHost,"preciceGateway",uint8(22));
            else
                bool = preciceGateway(uint8(22));
            end
        end
        
        % isWriteDataRequired
        function bool = isWriteDataRequired(obj,dt)
            if (obj.bMexHostRunning)
                bool = feval(obj.oMexHost,"preciceGateway",uint8(23),dt);
            else
                bool = preciceGateway(uint8(23),dt);
            end
        end
        
        % isTimestepComplete
        function bool = isTimestepComplete(obj)
            if (obj.bMexHostRunning)
                bool = feval(obj.oMexHost,"preciceGateway",uint8(24));
            else
                bool = preciceGateway(uint8(24));
            end
        end
        
        % hasToEvaluateSurrogateModel
        function bool = hasToEvaluateSurrogateModel(obj)
            if (obj.bMexHostRunning)
                bool = feval(obj.oMexHost,"preciceGateway",uint8(25));
            else
                bool = preciceGateway(uint8(25));
            end
        end
        
        % hasToEvaluateFineModel
        function bool = hasToEvaluateFineModel(obj)
            if (obj.bMexHostRunning)
                bool = feval(obj.oMexHost,"preciceGateway",uint8(26));
            else
                bool = preciceGateway(uint8(25));
            end
        end
        
        %% Action Methods
        % isActionRequired
        function bool = isActionRequired(obj,action)
            if ischar(action)
                action = string(action);
            end
            if (obj.bMexHostRunning)
                bool = feval(obj.oMexHost,"preciceGateway",uint8(30),action);
            else
                bool = preciceGateway(uint8(30),action);
            end
        end
        
        % fulfilledAction
        function fulfilledAction(obj,action)
            if ischar(action)
                action = string(action);
            end
            if (obj.bMexHostRunning)
                feval(obj.oMexHost,"preciceGateway",uint8(31),action);
            else
                preciceGateway(uint8(31),action);
            end
        end
        
        %% Mesh Access
        % hasMesh
        function bool = hasMesh(obj,meshName)
            if ischar(meshName)
                meshName = string(meshName);
            end
            if (obj.bMexHostRunning)
                bool = feval(obj.oMexHost,"preciceGateway",uint8(40),meshName);
            else
                bool = preciceGateway(uint8(40),meshName);
            end
        end
        
        % getMeshID
        function id = getMeshID(obj,meshName)
            if ischar(meshName)
                meshName = string(meshName);
            end
            if (obj.bMexHostRunning)
                id = feval(obj.oMexHost,"preciceGateway",uint8(41),meshName);
            else
                id = preciceGateway(uint8(41),meshName);
            end
        end
        
        % getMeshIDs
        function ids = getMeshIDs(obj)
            if (obj.bMexHostRunning)
                ids = feval(obj.oMexHost,"preciceGateway",uint8(42));
            else
                ids = preciceGateway(uint8(42));
            end
        end
        
        % getMeshHandle not yet implemented
        
        % setMeshVertex
        function vertexId = setMeshVertex(obj,meshID,position)
            if (obj.bMexHostRunning)
                vertexId = feval(obj.oMexHost,"preciceGateway",uint8(44),int32(meshID),position);
            else
                vertexId = preciceGateway(uint8(44),int32(meshID),position);
            end
        end
        
        % getMeshVertexSize
        function vertexId = getMeshVertexSize(obj,meshID)
            if (obj.bMexHostRunning)
                vertexId = feval(obj.oMexHost,"preciceGateway",uint8(45),int32(meshID));
            else
                vertexId = preciceGateway(uint8(45),int32(meshID));
            end
        end
        
        % setMeshVertices
        function vertexIds = setMeshVertices(obj,meshID,inSize,positions)
            if size(positions,2) ~= inSize
                error('Number of columns in position vector must match size!');
            end
            if (obj.bMexHostRunning)
                vertexIds = feval(obj.oMexHost,"preciceGateway",uint8(46),int32(meshID),int32(inSize),positions);
            else
                vertexIds = preciceGateway(uint8(46),int32(meshID),int32(inSize),positions);
            end
        end
        
        % getMeshVertices
        function positions = getMeshVertices(obj,meshID,inSize,vertexIds)
            if size(ids,2) ~= inSize
                error('Number of columns in position vector must match size!');
            end
            if (obj.bMexHostRunning)
                positions = feval(obj.oMexHost,"preciceGateway",uint8(47),int32(meshID),int32(inSize),vertexIds);
            else
                positions = preciceGateway(uint8(47),int32(meshID),int32(inSize),vertexIds);
            end
        end
        
        % getMeshVertexIDsFromPositions
        function vertexIds = getMeshVertexIDsFromPositions(obj,meshID,inSize,positions)
            if size(positions,2) ~= inSize
                error('Number of columns in position vector must match size!');
            end
            if (obj.bMexHostRunning)
                vertexIds = feval(obj.oMexHost,"preciceGateway",uint8(48),int32(meshID),int32(inSize),positions);
            else
                vertexIds = preciceGateway(uint8(48),int32(meshID),int32(inSize),positions);
            end
        end
        
        % setMeshEdge
        function edgeID = setMeshEdge(obj, meshID, firstVertexID, secondVertexID)
            if (obj.bMexHostRunning)
                edgeID = feval(obj.oMexHost,"preciceGateway",uint8(49),int32(meshID),int32(firstVertexID),int32(secondVertexID));
            else
                edgeID = preciceGateway(uint8(49),int32(meshID),int32(firstVertexID),int32(secondVertexID));
            end
        end
        
        % setMeshTriangle
        function setMeshTriangle(obj, meshID, firstEdgeID, secondEdgeID, thirdEdgeID)
            if (obj.bMexHostRunning)
                feval(obj.oMexHost,"preciceGateway",uint8(50),int32(meshID),int32(firstEdgeID),int32(secondEdgeID),int32(thirdEdgeID));
            else
                preciceGateway(uint8(50),int32(meshID),int32(firstEdgeID),int32(secondEdgeID),int32(thirdEdgeID));
            end
        end
        
        % setMeshTriangleWithEdges
        function setMeshTriangleWithEdges(obj, meshID, firstVertexID, secondVertexID, thirdVertexID)
            if (obj.bMexHostRunning)
                feval(obj.oMexHost,"preciceGateway",uint8(51),int32(meshID),int32(firstVertexID),int32(secondVertexID),int32(thirdVertexID));
            else
                preciceGateway(uint8(51),int32(meshID),int32(firstVertexID),int32(secondVertexID),int32(thirdVertexID));
            end
        end
        
        % setMeshQuad
        function setMeshQuad(obj, meshID, firstEdgeID, secondEdgeID, thirdEdgeID, fourthEdgeID)
            if (obj.bMexHostRunning)
                feval(obj.oMexHost,"preciceGateway",uint8(52),int32(meshID),int32(firstEdgeID),int32(secondEdgeID),int32(thirdEdgeID),int32(fourthEdgeID));
            else
                preciceGateway(uint8(52),int32(meshID),int32(firstEdgeID),int32(secondEdgeID),int32(thirdEdgeID),int32(fourthEdgeID));
            end
        end
        
        % setMeshQuadWithEdges
        function setMeshQuadWithEdges(obj, meshID, firstVertexID, secondVertexID, thirdVertexID, fourthVertexID)
            if (obj.bMexHostRunning)
                feval(obj.oMexHost,"preciceGateway",uint8(53),int32(meshID),int32(firstVertexID),int32(secondVertexID),int32(thirdVertexID),int32(fourthVertexID));
            else
                preciceGateway(uint8(53),int32(meshID),int32(firstVertexID),int32(secondVertexID),int32(thirdVertexID),int32(fourthVertexID));
            end
        end
        
        %% Data Access
        % hasDataID
        function bool = hasData(obj,dataName,meshID)
            if ischar(dataName)
                dataName = string(dataName);
            end
            if (obj.bMexHostRunning)
                bool = feval(obj.oMexHost,"preciceGateway",uint8(60),dataName,int32(meshID));
            else
                bool = preciceGateway(uint8(60),dataName,int32(meshID));
            end
        end
        
        % getDataID
        function id = getDataID(obj,dataName,meshID)
            if ischar(dataName)
                dataName = string(dataName);
            end
            if (obj.bMexHostRunning)
                id = feval(obj.oMexHost,"preciceGateway",uint8(61),dataName,int32(meshID));
            else
                id = preciceGateway(uint8(61),dataName,int32(meshID));
            end
        end
        
        % mapReadDataTo
        function mapReadDataTo(obj,meshID)
            if (obj.bMexHostRunning)
                feval(obj.oMexHost,"preciceGateway",uint8(62),int32(meshID));
            else
                preciceGateway(uint8(62),int32(meshID));
            end
        end
        
        % mapWriteDataFrom
        function mapWriteDataFrom(obj,meshID)
            if (obj.bMexHostRunning)
                feval(obj.oMexHost,"preciceGateway",uint8(63),int32(meshID));
            else
                preciceGateway(uint8(63),int32(meshID));
            end
        end
        
        % writeBlockVectorData
        function writeBlockVectorData(obj,dataID,inSize,valueIndices,values)
            if ~isa(valueIndices,'int32')
                warning('valueIndices should be allocated as int32 to prevent copying.');
                valueIndices = int32(valueIndices);
            end
            if (obj.bMexHostRunning)
                feval(obj.oMexHost,"preciceGateway",uint8(64),int32(dataID),int32(inSize),valueIndices,values);
            else
                preciceGateway(uint8(64),int32(dataID),int32(inSize),valueIndices,values);
            end
        end
        
        % writeVectorData
        function writeVectorData(obj,dataID,valueIndex,value)
            if (obj.bMexHostRunning)
                feval(obj.oMexHost,"preciceGateway",uint8(65),int32(dataID),int32(valueIndex),value);
            else
                preciceGateway(uint8(65),int32(dataID),int32(valueIndex),value);
            end
        end
        
        % writeBlockScalarData
        function writeBlockScalarData(obj,dataID,inSize,valueIndices,values)
            if ~isa(valueIndices,'int32')
                warning('valueIndices should be allocated as int32 to prevent copying.');
                valueIndices = int32(valueIndices);
            end
            if (obj.bMexHostRunning)
                feval(obj.oMexHost,"preciceGateway",uint8(66),int32(dataID),int32(inSize),valueIndices,values);
            else
                preciceGateway(uint8(66),int32(dataID),int32(inSize),valueIndices,values);
            end
        end
        
        % writeScalarData
        function writeScalarData(obj,dataID,valueIndex,value)
            if (obj.bMexHostRunning)
                feval(obj.oMexHost,"preciceGateway",uint8(67),int32(dataID),int32(valueIndex),value);
            else
                preciceGateway(uint8(67),int32(dataID),int32(valueIndex),value);
            end
        end
        
        % readBlockVectorData
        function values = readBlockVectorData(obj,dataID,inSize,valueIndices)
            if ~isa(valueIndices,'int32')
                warning('valueIndices should be allocated as int32 to prevent copying.');
                valueIndices = int32(valueIndices);
            end
            if (obj.bMexHostRunning)
                values = feval(obj.oMexHost,"preciceGateway",uint8(68),int32(dataID),int32(inSize),valueIndices);
            else
                values = preciceGateway(uint8(68),int32(dataID),int32(inSize),valueIndices);
            end
        end
        
        % readVectorData
        function value = readVectorData(obj,dataID,valueIndex)
            if (obj.bMexHostRunning)
                value = feval(obj.oMexHost,"preciceGateway",uint8(69),int32(dataID),int32(valueIndex));
            else
                value = preciceGateway(uint8(69),int32(dataID),int32(valueIndex));
            end
        end
        
        % readBlockScalarData
        function values = readBlockScalarData(obj,dataID,inSize,valueIndices,transpose)
            if nargin<5
                transpose=false;
            end
            if ~isa(valueIndices,'int32')
                warning('valueIndices should be allocated as int32 to prevent copying.');
                valueIndices = int32(valueIndices);
            end
            if (obj.bMexHostRunning)
                values = feval(obj.oMexHost,"preciceGateway",uint8(70),int32(dataID),int32(inSize),valueIndices,transpose);
            else
                values = preciceGateway(uint8(70),int32(dataID),int32(inSize),valueIndices,transpose);
            end
        end
        
        % readScalarData
        function value = readScalarData(obj,dataID,valueIndex)
            if (obj.bMexHostRunning)
                value = feval(obj.oMexHost,"preciceGateway",uint8(71),int32(dataID),int32(valueIndex));
            else
                value = preciceGateway(uint8(71),int32(dataID),int32(valueIndex));
            end
        end
    end
end