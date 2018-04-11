classdef MFF
    %A group of functions that make manipulating files and folders easier
    %
    %   MFF.Dir2Table()
    %
    
    methods(Static)
        function dirTable = Dir2Table(varargin)
            % List folder contents. It is exactly the same function as dir(...) but output a table. 
            %
            %   dirTable = MFile.DirTable( fileDir, fileName )
            %
            
            dirStruct = dir(varargin{:});
            dirTable = struct2table(dirStruct);
            
            for i = 1 : width(dirTable)
                if ischar(dirTable.(i))
                    dirTable.(i) = {dirTable.(i)};
                end
            end
        end
        
        
    end
end

