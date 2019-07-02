classdef Constants
    %CONSTANTS Summary of this class goes here
    %   Detailed explanation goes here
    methods(Access=private)
        % Make the constructor private so the user does not initialize it.
        function obj = Constants()
        end
    end

    methods(Static)
        function result = actionWriteInitialData()
            persistent constantName;
            if isempty(constantName)
                constantName = preciceGateway(uint8(0));
            end
            result = constantName;
        end
        
        function result = actionWriteIterationCheckpoint()
            persistent constantName;
            if isempty(constantName)
                constantName = preciceGateway(uint8(1));
            end
            result = constantName;
        end
        
        function result = actionReadIterationCheckpoint()
            persistent constantName;
            if isempty(constantName)
                constantName = preciceGateway(uint8(2));
            end
            result = constantName;
        end
    end
end