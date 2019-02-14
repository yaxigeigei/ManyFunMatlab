classdef MUtil
    %MUtil Summary of this class goes here
    %   Detailed explanation goes here
    
    methods(Static)
        function s = CombineStructs(varargin)
            % Combine multiple structures into one. Field names must be unique across structures. 
            %
            %   s = MUtil.CombineStructs(varargin)
            % 
            
            tbs = cellfun(@(x) struct2table(x, 'AsArray', true), varargin, 'Uni', false);
            s = table2struct(cat(2, tbs{:}));
        end
        
        function dirTable = Dir2Table(varargin)
            % List folder contents. It is exactly the same function as dir(...) but output a table. 
            %
            %   dirTable = MUtil.DirTable( fileDir, fileName )
            %
            
            dirStruct = dir(varargin{:});
            dirTable = struct2table(dirStruct);
            
            for i = 1 : width(dirTable)
                if ischar(dirTable.(i))
                    dirTable.(i) = {dirTable.(i)};
                end
            end
        end
        
        function hit = MatchAnyRegExp(str, expressions)
            % Match any of the regular expressions
            % 
            %   hit = MUtil.MatchAnyRegExp(str, expressions)
            %
            
            expressions = cellstr(expressions);
            expressions = expressions(:);
            hit = cellfun(@(x) any(regexp(str, x, 'once')), expressions);
            hit = any(hit);
        end
    end
end

