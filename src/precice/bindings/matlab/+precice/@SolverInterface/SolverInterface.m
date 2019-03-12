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
        interfaceID;
    end
    
    methods
        function obj = SolverInterface(SolverName)
            %SOLVERINTERFACE Construct an instance of this class
            %   Detailed explanation goes here
            obj.interfaceID = 0;
            preciceGateway(uint8(0),SolverName);
        end
        
        function delete(obj)
            obj.interfaceID = 0;
            preciceGateway(uint8(1));
        end
        
        function configure(obj,configFileName)
            obj.interfaceID = 0;
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
        
        function dims = getDimensions(~)
            dims = preciceGateway(uint8(20));
        end
        
        function bool = isCouplingOngoing(~)
            bool = preciceGateway(uint8(21));
        end
        
        function bool = isActionRequired(~,action)
            bool = preciceGateway(uint8(30),action);
        end
        
        function fulfilledAction(~,action)
            preciceGateway(uint8(31),action);
        end
        
        function id = getMeshID(~,meshName)
            id = preciceGateway(uint8(41),meshName);
        end
        
        function vertexId = setMeshVertex(~,meshID,position)
            vertexId = preciceGateway(uint8(44),meshID,position);
        end
    end
end