classdef SolverInterface < handle
    %SOLVERINTERFACE Matlab wrapper for the C++ SolverInterface
    %   This class will serve as a gateway for the precice Solver Interface
    %   in C++. It has access to a private mex Function written in the new
    %   MATLAB C++ API. Such a function is actually an instance of a C++ 
    %   class that is stored internally by MATLAB. Whenever they are 
    %   "called", a subroutine of the object (operator ()) is invoked which
    %   then has access to the of the class.
    %   Hence, we can store the actual solver interface in this class and
    %   access it by invoking the mex function.
    
    properties(Access=private)
        % This is currently just a dummy variable. It will be given meaning
        % later.
        interfaceID;
        
    end
    
    methods
        %% Construction and configuration
        % Constructor
        function obj = SolverInterface(SolverName)
            %SOLVERINTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.interfaceID = 0;
            if ischar(SolverName)
                SolverName = string(SolverName);
            end
            preciceGateway(uint8(0),SolverName);
        end
        
        % Destructor
        function delete(obj)
            obj.interfaceID = 0;
            preciceGateway(uint8(1));
        end
        
        % configure
        function configure(obj,configFileName)
            obj.interfaceID = 0;
            if ischar(configFileName)
                configFileName = string(configFileName);
            end
            preciceGateway(uint8(2),configFileName);
        end
        
        %% Steering methods
        % initialize
        function dt = initialize(obj)
            obj.interfaceID = 0;
            dt = preciceGateway(uint8(10));
        end
        
        % initialize Data
        function initializeData(obj)
            obj.interfaceID = 0;
            preciceGateway(uint8(11));
        end
        
        % advance
        function dt = advance(obj,dt)
            obj.interfaceID = 0;
            dt = preciceGateway(uint8(12),dt);
        end
        
        % finalize
        function finalize(obj)
            obj.interfaceID = 0;
            preciceGateway(uint8(13));
        end
        
        %% Status queries
        % getDimensions
        function dims = getDimensions(obj)
            obj.interfaceID = 0;
            dims = preciceGateway(uint8(20));
        end
        
        % isCouplingOngoing
        function bool = isCouplingOngoing(obj)
            obj.interfaceID = 0;
            bool = preciceGateway(uint8(21));
        end
        
        % isReadDataAvailable
        function bool = isReadDataAvailable(obj)
            obj.interfaceID = 0;
            bool = preciceGateway(uint8(22));
        end
        
        % isWriteDataRequired
        function bool = isWriteDataRequired(obj,dt)
            obj.interfaceID = 0;
            bool = preciceGateway(uint8(23),dt);
        end
        
        % isCouplingOngoing
        function bool = isTimestepComplete(obj)
            obj.interfaceID = 0;
            bool = preciceGateway(uint8(24));
        end
        
        % isCouplingOngoing
        function bool = hasToEvaluateSurrogateModel(obj)
            obj.interfaceID = 0;
            bool = preciceGateway(uint8(25));
        end
        
        % isCouplingOngoing
        function bool = hasToEvaluateFineModel(obj)
            obj.interfaceID = 0;
            bool = preciceGateway(uint8(26));
        end
        
        %% Action Methods
        % isActionRequired
        function bool = isActionRequired(obj,action)
            if ischar(action)
                action = string(action);
            end
            obj.interfaceID = 0;
            bool = preciceGateway(uint8(30),action);
        end
        
        % fulfilledAction
        function fulfilledAction(obj,action)
            if ischar(action)
                action = string(action);
            end
            obj.interfaceID = 0;
            preciceGateway(uint8(31),action);
        end
        
        %% Mesh Access
        % hasMesh
        function bool = hasMesh(obj,meshName)
            if ischar(meshName)
                meshName = string(meshName);
            end
            obj.interfaceID = 0;
            bool = preciceGateway(uint8(40),meshName);
        end
        
        % getMeshID
        function id = getMeshID(obj,meshName)
            if ischar(meshName)
                meshName = string(meshName);
            end
            obj.interfaceID = 0;
            id = preciceGateway(uint8(41),meshName);
        end
        
        % getMeshIDs
        function ids = getMeshIDs(obj)
            obj.interfaceID = 0;
            ids = preciceGateway(uint8(42));
        end
        
        % getMeshHandle not yet implemented
        
        % setMeshVertex
        function vertexId = setMeshVertex(obj,meshID,position)
            obj.interfaceID = 0;
            vertexId = preciceGateway(uint8(44),int32(meshID),position);
        end
        
        % getMeshVertexSize
        function vertexId = getMeshVertexSize(obj,meshID)
            obj.interfaceID = 0;
            vertexId = preciceGateway(uint8(45),int32(meshID));
        end
        
        % setMeshVertices
        function vertexIds = setMeshVertices(obj,meshID,inSize,positions)
            if size(positions,2) ~= inSize
                error('Number of columns in position vector must match size!');
            end
            obj.interfaceID = 0;
            vertexIds = preciceGateway(uint8(46),int32(meshID),uint64(inSize),positions);
        end
        
        % getMeshVertices
        function positions = getMeshVertices(obj,meshID,inSize,vertexIds)
            if size(ids,2) ~= inSize
                error('Number of columns in position vector must match size!');
            end
            obj.interfaceID = 0;
            positions = preciceGateway(uint8(47),int32(meshID),uint64(inSize),vertexIds);
        end
        
        % getMeshVertexIDsFromPositions
        function vertexIds = getMeshVertexIDsFromPositions(obj,meshID,inSize,positions)
            if size(positions,2) ~= inSize
                error('Number of columns in position vector must match size!');
            end
            obj.interfaceID = 0;
            vertexIds = preciceGateway(uint8(48),int32(meshID),uint64(inSize),positions);
        end
        
        %% Data Access
        % getMeshVertices
        function id = getDataID(obj,dataName,meshID)
            obj.interfaceID = 0;
            if ischar(dataName)
                dataName = string(dataName);
            end
            id = preciceGateway(uint8(61),dataName,int32(meshID));
        end
        
        function writeBlockScalarData(obj,dataID,inSize,valueIndices,values)
            obj.interfaceID = 0;
            if ~isa(valueIndices,'int32')
                warning('valueIndices should be allocated as int32 to prevent copying.');
                valueIndices = int32(valueIndices);
            end
            preciceGateway(uint8(66),int32(dataID),uint64(inSize),valueIndices,values);
        end
        
        function values = readBlockScalarData(obj,dataID,inSize,valueIndices)
            obj.interfaceID = 0;
            if ~isa(valueIndices,'int32')
                warning('valueIndices should be allocated as int32 to prevent copying.');
                valueIndices = int32(valueIndices);
            end
            values = preciceGateway(uint8(70),int32(dataID),uint64(inSize),valueIndices);
        end
    end
end