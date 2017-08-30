classdef ManyFuncImgBaseClass
    %MANYFUNCIMGBASECLASS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Abstract)
    end
    
    methods(Static)
%         function img = Adjust(img)
%             %Adjust Summary of this function goes here
%             %   Detailed explanation goes here
%             
%             img = Img34.Array2Cell(img);
%             
%             maxValue = 0;
%             for i = 1 : length(img)
%                 maxValue = max(maxValue, max(img{i}(:)));
%             end
%             for i = 1 : length(img)
%                 img{i} = img{i} / maxValue;
%             end
%         end
        
        function Export(path, img, indexFormat)
            %Export image array as a TIFF file (will overwrite the file with the same name)
            %
            % TiffNow.Export(path, img, indexFormat)
            %
            %Input discription
            %
            % path: the path (directory + file name without suffix) to save
            % img: a two, three or four dimensional array, or a cell array of images(or stacks)
            % indexFormat: a string for formated indexing of individual stacks
            % (default is '_%04d', i.e. an underscore followed by four digits integer)
            
            if nargin < 3
                indexFormat = '_%04d';
            end
            
            img = Img34.Array2Cell(img);
            
            for i = 1 : length(img)
                a = img{i}(1);
                s = whos('a');
                bps = s.bytes * 8;
                
                if length(img) > 1
                    counting = sprintf(indexFormat, i);
                    tObj = Tiff([ path, counting '.tif' ], 'w');
                else
                    [ ~, ~, ext ] = fileparts(path);
                    if ~strcmp(ext, '.tif')
                        path = [ path, '.tif' ];
                    end
                    tObj = Tiff(path, 'w');
                end
                
                tagstruct.ImageLength = size(img{i}, 1);
                tagstruct.ImageWidth = size(img{i}, 2);
                tagstruct.BitsPerSample = bps;
                tagstruct.Compression = Tiff.Compression.None;
                tagstruct.SamplesPerPixel = 1;
                tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
                tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
                
                tObj.setTag(tagstruct);
                tObj.write(img{i}(:,:,1));
                
                for j = 2 : size(img{i}, 3)
                    tObj.writeDirectory();
                    tObj.setTag(tagstruct);
                    tObj.write(img{i}(:,:,j))
                end
                
                tObj.close();
            end
        end
        
        function img = Import(path, varargin)
            % Imports a TIFF file (conventional non-dynamic loading)
            %
            %   img = TiffNow.Import(path)
            %
            % Inputs
            %    path: The path of the TIFF file. File browsing UI will prompt if passing empty value.
            %   'range': Two-element vector indicating the begining and ending index of the range of interest. (default is the entire tif)
            %   'interleave': Two-element vector. The first specifies total number of channels; the second specifies the channel you want to keep. (default [1,1])
            %   'dataType': The data type of the image. (default 'uint8' for 8-bit TIFF, 'uint16' for 16-bit, 'double' for other depths)
            %
            % Output
            %   img: a two or three dimensional array
            
            p = inputParser;
            p.addParameter('range', [ 1 inf ], @isnumeric);
            p.addParameter('interleave', [ 1 1 ], @isnumeric);
            p.addParameter('dataType', '', @ischar);
            p.parse(varargin{:});
            frameBegin = p.Results.range(1);
            frameEnd = p.Results.range(2);
            channelTotal = p.Results.interleave(1);
            channelSelect = p.Results.interleave(2);
            datatype = p.Results.dataType;
            
            if nargin < 1 || isempty(path)
                path = Browse.File();
            end
            
            warning('off', 'MATLAB:imagesci:Tiff:libraryWarning');
            
            [ ~, ~, ext ] = fileparts(path);
            if ~strcmp(ext, '.tif')
                path = [ path, '.tif' ];
            end
            tObj = Tiff(path, 'r');
            imgWidth = tObj.getTag('ImageWidth');
            imgHeight = tObj.getTag('ImageLength');
            
            if isempty(datatype)
                pxlBytes = round(tObj.getTag('BitsPerSample') / 8);
                switch pxlBytes
                    case 1
                        datatype = 'uint8';
                    case 2
                        datatype = 'uint16';
                    otherwise
                        datatype = 'double';
                end
            end
            
            frameObj = 1;
            while ~tObj.lastDirectory()
                frameObj = frameObj + 1;
                tObj.nextDirectory();
            end
            
            frameEnd = floor(min([ frameObj / channelTotal, frameEnd ]));
            
            img = zeros(imgHeight, imgWidth, frameEnd-frameBegin+1, datatype);
            for i = frameBegin : frameEnd
                tObj.setDirectory((i-1) * channelTotal + channelSelect);
                img(:,:,i-frameBegin+1) = tObj.read();
            end
            
            tObj.close();
        end
        
        function stack = Import2Stack(dirPath, key)
            %Combine multiple TIFF files into a stack or stacks
            %
            % stack = TiffNow.Import2Stack(dirPath, key)
            %
            %Input discription
            %
            % dirPath (optional, default: prompting browser to select files (no 4-D grouping)):
            % The directory of TIFF files.
            %
            % key (optional, default: '*.tif'):
            % A string of keyword with wildcards to specify a group of files (e.g. 'ChanA_????_*.tif'). Files
            % within one groups must have the same size in every dimensions.
            % Also, you may use formated number (e.g. '*_%04d.tif') to incrementally specify different groups
            % of files and import them into respective stacks. Images across stacks are not required to have 
            % the same size.
            %
            %Output discription
            %
            % stacks: a three dimensional array or a cell array containing multiple 3D arrays, depending on the
            % grouping effect of the keyword
            
            if nargin < 1
                % Uses UI to select files
                [ filePaths, dirPath ] = Browse.Files();
                
                % Groups them into one stack
                list{1} = filePaths;
            else
                % Default keyword finds all TIFF files in the folder
                if nargin < 2
                    key = '*.tif';
                end
                
                % Initializes grouping
                list = cell(1);
                listEnd = false;
                i = 1;
                
                while(~listEnd)
                    % Applies indexing, if any, to the keyword
                    fileName = sprintf(key, i);
                    
                    % Finds all files satisfying the keyword
                    tempFilePaths = dir(fullfile(dirPath, fileName));
                    tempFilePaths = struct2cell(tempFilePaths);
                    tempFilePaths = tempFilePaths(1,:)';
                    
                    % Saves above file paths for current image stack
                    if ~isempty(tempFilePaths)
                        list{i} = tempFilePaths;
                    else
                        % Stops grouping if nothing more to do
                        listEnd = true;
                    end
                    
                    % Increments index for next stack
                    i = i + 1;
                    
                    % Stops grouping if the keyword is not incremental
                    if strcmp(fileName, sprintf(key, i));
                        listEnd = true;
                    end
                end
            end
            
            % Loading images to stack(s)
            stack = cell(length(list), 1);
            if ~isempty(list{1})
                for i = 1 : length(list)
                    % Gets the dimension info for preallocation
                    firstImg = TiffNow.Import(fullfile(dirPath, list{i}{1}));
                    [ imgHeight, imgWidth, imgSlices ] = size(firstImg);
                    stack{i} = zeros([ imgHeight, imgWidth, imgSlices * length(list{i}) ], 'like', firstImg);
                    
                    % Loads files into the current stack
                    stack{i}(:,:,1:imgSlices) = firstImg;
                    for j = 2 : length(list{i})
                        stack{i}(:, :, (j-1)*imgSlices+1 : j*imgSlices) = TiffNow.Import(fullfile(dirPath, list{i}{j}));
                    end
                end
            end
            
            if length(stack) == 1
                stack = stack{1};
            end
        end
        
        function proj = ProjCustom(img, hFnc)
            %ProjCustom Summary of this function goes here
            %   Detailed explanation goes here
            
            img = Img34.Array2Cell(img);
            
            proj = zeros(size(img{1},1), size(img{1},2), length(img), 'like', img{1});
            for i = 1 : length(img)
                proj(:,:,i) = hFnc(img{i});
            end
        end
        
        function proj = ProjMax(img)
            %ProjMax Summary of this function goes here
            %   Detailed explanation goes here
            
            img = Img34.Array2Cell(img);
            
            proj = zeros(size(img{1},1), size(img{1},2), length(img), 'like', img{1});
            for i = 1 : length(img)
                proj(:,:,i) = max(img{i}, [ ], 3);
            end
        end
        
        function proj = ProjMean(img)
            %ProjAverage Summary of this function goes here
            %   Detailed explanation goes here
            
            img = Img34.Array2Cell(img);
            
            proj = zeros(size(img{1},1), size(img{1},2), length(img), 'like', img{1});
            for i = 1 : length(img)
                proj(:,:,i) = mean(img{i}, 3);
            end
        end
        
        function [ comboStack, xProj, yProj, zProj ] = ProjXYZ(img, zFactor)
            %ProjCombo Summary of this function goes here
            %   Detailed explanation goes here
            
            if nargin < 2
                zFactor = 1;
            end
            
            img = Img34.Array2Cell(img);
            numStacks = length(img);
            [ imgHeight, imgWidth, imgSlices ] = size(img{1});
            
            xProj = zeros(imgHeight, imgSlices, numStacks, 'like', img{1});
            yProj = zeros(imgSlices, imgWidth, numStacks, 'like', img{1});
            zProj = zeros(imgHeight, imgWidth, numStacks, 'like', img{1});
            
            for i = 1 : numStacks
                xProj(:,:,i) = squeeze(max(img{i}, [ ], 2));
                yProj(:,:,i) = squeeze(max(img{i}, [ ], 1))';
                zProj(:,:,i) = max(img{i}, [ ], 3);
            end
            
            xProj = Img23.Scale(xProj, zFactor, 1, 1);
            yProj = Img23.Scale(yProj, 1, zFactor, 1);
            
            imgSlices = size(yProj,1);
            comboStack = [ yProj, ones(imgSlices, imgSlices, numStacks); zProj, xProj ];
        end
    end
    
end

