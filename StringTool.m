classdef StringTool
    %MWUTILITIES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        function hit = MatchRegExps(string, keywords)
            % Handles string input
            if ischar(keywords)
                keywords = { keywords };
            end
            
            % Handles multidimensional cell array input
            keywords = keywords(:);
            
            % Matching
            hit = false;
            for i = 1 : length(keywords)
                if ~isempty(regexp(string, keywords{i}, 'once'))
                    hit = true;
                end
            end
        end
    end
    
end


