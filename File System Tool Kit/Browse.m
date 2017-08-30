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
        function [ filePath, folderPath, bareName, ext ] = File(defaultFolderPath, dialogTitle, filterSpec)
            % Browse to select a file and get path parts
            %
            %   [ filePath, folderPath, bareName, ext ] = Browse.File()
            %   [ filePath, folderPath, bareName, ext ] = Browse.File(defaultFolderPath)
            %   [ filePath, folderPath, bareName, ext ] = Browse.File(defaultFolderPath, dialogTitle)
            %   [ filePath, folderPath, bareName, ext ] = Browse.File(defaultFolderPath, dialogTitle, filterSpec)
            %
            
            % Handle user inputs
            if nargin < 3
                filterSpec = {'*.*'};
            end
            
            if nargin < 2
                dialogTitle = [];
            end
            
            if nargin < 1
                defaultFolderPath = [];
            end
            
            if isempty(dialogTitle)
                dialogTitle = 'Please select a file';
            end
            
            if isempty(defaultFolderPath)
                try
                    load('lastimeDir', 'defaultFolderPath');
                catch
                    defaultFolderPath = pwd;
                end
            end
            
            % Get file path info
            [ fileName, folderPath ] = uigetfile(filterSpec, dialogTitle, defaultFolderPath);
            
            % Check selection
            if folderPath
                % Decompose paths
                filePath = fullfile(folderPath, fileName);
                [ ~, bareName, ext ] = fileparts(fileName);
                
                % Remember new folder path
                defaultFolderPath = folderPath;
                save('lastimeDir.mat', 'defaultFolderPath');
            else
                filePath = '';
                bareName = '';
                ext = '';
            end
        end
        
        function [ filePath, folderPath, bareNames, exts ] = Files(defaultFolderPath, dialogTitle, filterSpec)
            % Browse to select multiple files and get path parts
            %
            %   [ filePath, folderPath, bareNames, exts ] = Browse.Files()
            %   [ filePath, folderPath, bareNames, exts ] = Browse.Files(defaultFolderPath)
            %   [ filePath, folderPath, bareNames, exts ] = Browse.Files(defaultFolderPath, dialogTitle)
            %   [ filePath, folderPath, bareNames, exts ] = Browse.Files(defaultFolderPath, dialogTitle, filterSpec)
            %
            
            % Handle user inputs
            if nargin < 3
                filterSpec = {'*.*'};
            end
            
            if nargin < 2
                dialogTitle = [];
            end
            
            if nargin < 1
                defaultFolderPath = [];
            end
            
            if isempty(dialogTitle)
                dialogTitle = 'Please select files';
            end
            
            if isempty(defaultFolderPath)
                try
                    load('lastimeDir', 'defaultFolderPath');
                catch
                    defaultFolderPath = pwd;
                end
            end
            
            % Get file path info
            [ fileName, folderPath ] = uigetfile(filterSpec, dialogTitle, defaultFolderPath, 'MultiSelect', 'on');
            
            if folderPath
                % Decompose paths
                fileName = cellstr(fileName)';
                filePath = cellfun(@(x) fullfile(folderPath, x), fileName, 'UniformOutput', false);
                [ ~, bareNames, exts ] = cellfun(@fileparts, fileName, 'UniformOutput', false);
                
                % Remember new folder path
                defaultFolderPath = folderPath;
                save('lastimeDir.mat', 'defaultFolderPath');
            else
                filePath = {};
                bareNames = {};
                exts = {};
            end
        end
        
        function folderPath = Folder(defaultFolderPath, dialogTitle)
            % Browse to select a folder and get path parts
            %
            %   folderPath = Browse.Folder()
            %   folderPath = Browse.Folder(defaultFolderPath)
            %   folderPath = Browse.Folder(defaultFolderPath, dialogTitle)
            %
            
            % Handle user inputs
            if nargin < 2
                dialogTitle = [];
            end
            
            if nargin < 1
                defaultFolderPath = [];
            end
            
            if isempty(dialogTitle)
                dialogTitle = 'Please select a folder';
            end
            
            if isempty(defaultFolderPath)
                try
                    load('lastimeDir', 'defaultFolderPath');
                catch
                    defaultFolderPath = pwd;
                end
            end
            
            % Get file path info
            folderPath = uigetdir(defaultFolderPath, dialogTitle);
            
            if folderPath
                % Remember new folder path
                defaultFolderPath = folderPath;
                save('lastimeDir.mat', 'defaultFolderPath');
            end
        end
    end
    
end




