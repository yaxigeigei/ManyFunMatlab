classdef MData
    %MDATA Summary of this class goes here
    %   Detailed explanation goes here
    
    methods(Static)
        function s = CombineStructs(varargin)
            % Combine multiple structures into one. Field names must be unique across structures. 
            
            tbs = cellfun(@(x) struct2table(x, 'AsArray', true), varargin, 'Uni', false);
            s = table2struct(cat(2, tbs{:}));
        end
        
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

