function stacks = TiffBatchProcess( varargin )
%Interactively import and preprocess TIFF image(s) according to your
%requirements, especially suitable for bath processing.
%
% stacks = TiffBatchImport()
% stacks = TiffBatchImport(path, ...)
%
%
%Input discription
%
% path(optional): the directory of the TIFF file(s)
% settings(optional): a construct containing pre-set options (then no question will be prompted anymore)
%
%
%Output discription
%
% stacks:
% A linear cell array in which each cell contains an image array from one imported
% TIFF file. 
%
% maxStack:
% An array of max-projections where each frame is a maximally projected
% image from one TIFF file.
%
% brightStack:
% Similar to 'maxStack' except that brightest slices, according to your selected
% region of interest, within each TIFF file are chosen for max-projection.


    p = inputParser;
    
    addParamValue(p, 'fileDir', '', @ischar);
    addParamValue(p, 'key', '', @ischar);
    addParamValue(p, 'ask', true, @islogical);
    addParamValue(p, 'channelTotal', 1, @isnumeric);
    addParamValue(p, 'channelSelect', 1, @isnumeric);
    addParamValue(p, 'frameStart', 1, @isnumeric);
    addParamValue(p, 'frameEnd', inf, @isnumeric);
    addParamValue(p, 'intraManualReg', false, @islogical);
    addParamValue(p, 'intraAutoReg', false, @islogical);
    addParamValue(p, 'interManualReg', false, @islogical);
    addParamValue(p, 'interAutoReg', false, @islogical);
    addParamValue(p, 'save', false, @islogical);
    
    parse(p, varargin{:});
    
    fileDir = p.Results.fileDir;
    key = p.Results.key;
    askCfg = p.Results.ask;
    channelTotal = p.Results.channelTotal;
    channelSelect = p.Results.channelSelect;
    frameStart = p.Results.frameStart;
    frameEnd = p.Results.frameEnd;
    intraAutoReg = p.Results.intraAutoReg;
    intraManualReg = p.Results.intraManualReg;
    interAutoReg = p.Results.interAutoReg;
    interManualReg = p.Results.interManualReg;
    saveEach = p.Results.save;
    



%% Configuration


    if strcmp(key, '')
        [ fileName, fileDir, ~ ] = uigetfile({'*.tif;*.tiff'}, 'Select Image or Images', [fileDir '\'], 'MultiSelect', 'on');
        if ischar(fileName)
            fileName = { fileName };
        end
    else
        fileList = dir(fullfile(fileDir, key));
        fileName = cell(1, length(fileList));
        for i = 1 : length(fileList)
            fileName{i} = fileList(i).name;
        end
    end
    

    if askCfg

        prompt = {'Total Number of Channels: ', 'Channel of Interest: ', 'Start Frame: ', 'End Frame: (inf means to the end)'};
        ansCfg = inputdlg(prompt, 'Import Configuration', 1, {'1', '1', '1', 'inf'});
        channelTotal = str2double(ansCfg{1});
        channelSelect = str2double(ansCfg{2});
        frameStart = str2double(ansCfg{3});
        frameEnd = str2double(ansCfg{4});
        
        
        intraAutoReg = false;
        intraManualReg = false;
        interAutoReg = false;
        interManualReg = false;
        
        options = {'Intra-stack manual registration', 'Intra-stack automatic registration', ...
            'Inter-stacks manual registration', 'Inter-stacks automatic registration'};
        [ registerAns, ~ ] = listdlg('ListString', options, 'SelectionMode', 'multiple', 'Name', 'Registration', 'ListSize', [350 150]);

        if ~isempty(find(registerAns == 1, 1))
            intraManualReg = true;
        end
        if ~isempty(find(registerAns == 2, 1))
            intraAutoReg = true;
        end
        if ~isempty(find(registerAns == 3, 1))
            interManualReg = true;
        end
        if ~isempty(find(registerAns == 4, 1))
            interAutoReg = true;
        end
        
        saveEach = false;
        
        options = {'Save individual processed stacks'};
        [ saveAns, ~ ] = listdlg('ListString', options, 'SelectionMode', 'multiple', 'Name', 'Save', 'ListSize', [350 90]);

        if ~isempty(find(saveAns == 1, 1))
            saveEach = true;
        end

    end



%% Import & Export


    fileNumber = length(fileName);
    stacks = cell([ fileNumber 1 ]);
    
    for i = 1 : fileNumber
        stacks{i} = TiffImport(fullfile(fileDir, fileName{1, i}), channelTotal, channelSelect, frameStart, frameEnd);
        if intraManualReg
            [ stacks{i}, ~, ~ ] = StackManualReg(stacks{i});
        end
        if intraAutoReg
            stacks{i} = StacksAutoReg(stacks{i});
        end
    end
    
    if interManualReg
        stacks = StacksManualReg(stacks);
    end
    
    if interAutoReg
        stacks = StacksAutoReg(stacks);
    end

    if saveEach
        for i = 1 : length(stacks)
            TiffExport(fullfile(fileDir, [ 'Processed_' fileName{i} ]), stacks{i});
        end
    end



end

