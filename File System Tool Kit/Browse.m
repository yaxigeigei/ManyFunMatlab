classdef Browse
    %BROWSE a convenient combination and extention of uigetfile/uigetdir and fileparts
    %   
    %   Browse.File()
    %   Browse.Files()
    %   Browse.Folder()
    %
    
    properties
    end
    
    methods(Static)
        function [ filePath, folderName, bareName, ext ] = File(folderName)
            % Browse to select a file and get path parts
            %
            %   [ filePath, folderName, bareName, ext ] = File()
            %   [ filePath, folderName, bareName, ext ] = File(folderName)
            %
            
            if nargin < 1
                try
                    load('lastimeDir');
                catch
                    folderName = pwd;
                end
            end
            
            [ fileName, folderName ] = uigetfile({'*.*'}, 'Select Data', [folderName '\']);
            
            if folderName
                filePath = fullfile(folderName, fileName);
                [ ~, bareName, ext ] = fileparts(fileName);
                save('lastimeDir.mat', 'folderName');
            else
                filePath = '';
                bareName = '';
                ext = '';
            end
        end
        
        function [ filePath, folderName, bareNames, exts ] = Files(folderName)
            % Browse to select multiple files and get path parts
            %
            %   [ filePath, folderName, bareNames, exts ] = Files()
            %   [ filePath, folderName, bareNames, exts ] = Files(folderName)
            %
            
            if nargin < 1
                try
                    load('lastimeDir');
                catch
                    folderName = pwd;
                end
            end
            
            [ fileName, folderName ] = uigetfile({'*.*'}, 'Select Data', [folderName '\'], 'MultiSelect', 'on');
            
            if folderName
                fileName = cellstr(fileName)';
                filePath = cellfun(@(x) fullfile(folderName, x), fileName, 'UniformOutput', false);
                [ ~, bareNames, exts ] = cellfun(@fileparts, fileName, 'UniformOutput', false);
                save('lastimeDir.mat', 'folderName');
            else
                filePath = {};
                bareNames = {};
                exts = {};
            end
        end
        
        function folderName = Folder(folderName)
            % Browse to select a folder and get path parts
            %
            %   [ filePath, folderName, bareName, ext ] = Folder()
            %   [ filePath, folderName, bareName, ext ] = Folder(folderName)
            %
            
            if nargin < 1
                try
                    load('lastimeDir');
                catch
                    folderName = pwd;
                end
            end
            
            folderName = uigetdir(folderName);
            
            if folderName
                save('lastimeDir.mat', 'folderName');
            end
        end
    end
    
end

