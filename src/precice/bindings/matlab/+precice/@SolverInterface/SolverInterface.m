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
        function obj = SolverInterface(SolverName)
            %SOLVERINTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.interfaceID = 0;
            if ischar(SolverName)
                SolverName = string(SolverName);
            end
            preciceGateway(uint8(0),SolverName);
        end
        
        function delete(obj)
            obj.interfaceID = 0;
            preciceGateway(uint8(1));
        end
        
        function configure(obj,configFileName)
            obj.interfaceID = 0;
            if ischar(configFileName)
                configFileName = string(configFileName);
            end
            preciceGateway(uint8(2),configFileName);
        end
        
        function dt = initialize(obj)
            obj.interfaceID = 0;
            dt = preciceGateway(uint8(10));
        end
        
        function dt = advance(obj,dt)
            obj.interfaceID = 0;
            dt = preciceGateway(uint8(12),dt);
        end
        
        function finalize(obj)
            obj.interfaceID = 0;
            preciceGateway(uint8(13));
        end
        
        function dims = getDimensions(obj)
            obj.interfaceID = 0;
            dims = preciceGateway(uint8(20));
        end
        
        function bool = isCouplingOngoing(obj)
            obj.interfaceID = 0;
            bool = preciceGateway(uint8(21));
        end
        
        function bool = isActionRequired(obj,action)
            if ischar(action)
                action = string(action);
            end
            obj.interfaceID = 0;
            bool = preciceGateway(uint8(30),action);
        end
        
        function fulfilledAction(obj,action)
            if ischar(action)
                action = string(action);
            end
            obj.interfaceID = 0;
            preciceGateway(uint8(31),action);
        end
        
        function id = getMeshID(obj,meshName)
            if ischar(meshName)
                meshName = string(meshName);
            end
            obj.interfaceID = 0;
            id = preciceGateway(uint8(41),meshName);
        end
        
        function vertexId = setMeshVertex(obj,meshID,position)
            obj.interfaceID = 0;
            vertexId = preciceGateway(uint8(44),int32(meshID),position);
        end
        
        function vertexIds = setMeshVertices(obj,meshID,inSize,positions)
            if size(positions,2) ~= inSize
                error('Number of columns in position vector must match size!');
            end
            obj.interfaceID = 0;
            vertexIds = preciceGateway(uint8(46),int32(meshID),uint64(inSize),positions);
        end
        
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