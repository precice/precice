classdef Constants
    %CONSTANTS Summary of this class goes here
    %   Detailed explanation goes here
    methods(Access=private)
        % Make the constructor private so the user does not initialize it.
        function obj = Constants()
        end
    end

    methods(Static)
        function result = nameConfiguration()
            persistent constantName;
            if isempty(constantName)
                constantName = preciceGateway(uint8(0));
            end
            result = constantName;
        end
        
        function result = dataDisplacements()
            persistent constantName;
            if isempty(constantName)
                constantName = preciceGateway(uint8(1));
            end
            result = constantName;
        end
        
        function result = dataForces()
            persistent constantName;
            if isempty(constantName)
                constantName = preciceGateway(uint8(2));
            end
            result = constantName;
        end
        
        function result = dataVelocities()
            persistent constantName;
            if isempty(constantName)
                constantName = preciceGateway(uint8(3));
            end
            result = constantName;
        end
        
        function result = actionWriteInitialData()
            persistent constantName;
            if isempty(constantName)
                constantName = preciceGateway(uint8(4));
            end
            result = constantName;
        end
        
        function result = actionWriteIterationCheckpoint()
            persistent constantName;
            if isempty(constantName)
                constantName = preciceGateway(uint8(5));
            end
            result = constantName;
        end
        
        function result = actionReadIterationCheckpoint()
            persistent constantName;
            if isempty(constantName)
                constantName = preciceGateway(uint8(6));
            end
            result = constantName;
        end
        
        function result = actionPlotOutput()
            persistent constantName;
            if isempty(constantName)
                constantName = preciceGateway(uint8(7));
            end
            result = constantName;
        end
        
        function result = exportVTK()
            persistent constantName;
            if isempty(constantName)
                constantName = preciceGateway(uint8(8));
            end
            result = constantName;
        end
        
        function result = exportAll()
            persistent constantName;
            if isempty(constantName)
                constantName = preciceGateway(uint8(9));
            end
            result = constantName;
        end
    end
end