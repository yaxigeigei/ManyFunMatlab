classdef ColorCode
    %PLOTKIT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        
    end
    
    methods(Static)
        function rainbowCode = Rainbow(numColors)
            rainbowBase = { [ 1 1 0 0 0 1 ]; [ 0 1 1 1 0 0 ]; [ 0 0 0 1 1 1 ] };
            rainbowCode = cellfun(@(x) interp1(1:6, x, linspace(1,6,numColors)), ...
                rainbowBase, 'UniformOutput', false);
            rainbowCode = cell2mat(rainbowCode)';
        end
    end
    
end

