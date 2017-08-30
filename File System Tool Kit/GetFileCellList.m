function cellsList = GetFileCellList( fileDir, fileName )
%Get a list of file paths as a cell array
%
% [ cellsList ] = GetFileCellList( fileDir, fileName )
%
%
%Input discription
%
% fileDir:
% The directory of file(s) where wildcard expression is supported
% 
% fileName:
% File name(s) in which wildcard expression is supported
% (e.g. '*sun*.mat' returns all .mat files whose names contain 'sun' )
%
%
%Output discription
%
% cellsList:
% A linear cell array containing path(s) of matched file(s)


cellsList = dir(fullfile(fileDir, fileName));
cellsList = struct2cell(cellsList);
cellsList = cellsList(1,:)';


end

