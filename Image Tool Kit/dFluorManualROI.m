function [ snROIs, bgROIs, im2show ] = dFluorManualROI( img, selectionShape, fileDir )
%Interactively select one or more region of interest (ROI) on an image or
%image stack, and optionally save individual ROI as .mat files for reuse
%
% [ snROI, bgROI, im2show ] = dFluorManualROI( img, selectionShape, fileDir )
%
%
%Inputs discription
%
% img (required):
% An two or three dimensional image array
% 
% selectionShape:
% 'rect'        Dragging editable rectangle
% 'ellipse'     Dragging editable ellipse
% 'freehand'    Arbitrary drawing (default)
% 
% fileDir:
% File directory where a .mat file of ROI infomation is saved (default is '', for no saving)
%
%
%Outpus discription
%
% snROIs:
% A cell array. One ROI for each row in which the first column contains
% coordinates of its contour points and the second column is its logic
% mask. (If you do not select signal ROI, an empty cell { } is returned)
%
% bgROIs:
% Similar to 'snROI'
%
% im2show:
% Two dimensional array of the image on which ROIs were selected

    
    if nargin < 3
        fileDir = '';
    end
    
    if nargin < 2
        selectionShape = 'freehand';
    end
    
    
    snROIs = { };
    bgROIs = { };
    
    figure('Name', 'ROI selection');
    im2show = imadjust(max(img, [ ], 3));
    imagesc(im2show);
    axis equal tight
    colormap gray
    
    
    
    reDraw = true;
    
    fileDirBackup = fileDir;
    [ fileName, fileDir, ~ ] = uigetfile({'*.mat'}, 'Select saved ROI profile(s) (otherwise just close the window)', [fileDir '\'], 'MultiSelect', 'on');
    if ~iscell(fileName)
        fileName = { fileName };
    end
    
    if fileDir ~= 0
        for i = 1 : length(fileName)
            load(fullfile(fileDir, fileName{i}));
            try
                RoiShow(bgROI, 'b');
                bgROIs = vertcat(bgROIs, bgROI);
            catch
            end
            try
                RoiShow(snROI, 'r', i);
                snROIs = vertcat(snROIs, snROI);
            catch
            end
        end
        ansRedraw = questdlg('Do you want to re-select ROIs', 'ROI', 'Yes', 'No', 'Yes');
        if strcmp(ansRedraw, 'No')
            reDraw = false;
        end
    end
    
    
    
    if reDraw
        
        imagesc(im2show);
        axis equal tight
        
        ansBgROI = questdlg('Do you need to select a background region? (double click on your selection to confirm)', 'Background', 'Yes', 'No', 'Yes');
        if strcmp(ansBgROI, 'Yes')
            [ bgRoiContour, bgRoiMask ] = RoiSelect(selectionShape, 'b', size(img));
            bgROIs = { bgRoiContour, bgRoiMask };
            text(max(bgRoiContour(:,1)), max(bgRoiContour(:,2)), num2str(1), 'Color', 'b');
            
            contSelect = true;
            while(contSelect)
                ansCont = questdlg('Do you need to select more background regions?', 'Multiple Background ROIs', 'Yes', 'No', 'No');
                if strcmp(ansCont, 'No')
                    contSelect = false;
                end
                if contSelect
                    [ bgRoiContour, bgRoiMask ] = RoiSelect(selectionShape, 'b', size(img));
                    bgROIs(size(bgROIs,1)+1, 1:2) = { bgRoiContour, bgRoiMask };
                    text(max(bgRoiContour(:,1)), max(bgRoiContour(:,2)), num2str(size(bgROIs,1)), 'Color', 'b');
                end
            end
        else
            bgROIs = { };
        end
        
        ansSnROI = questdlg('Do you need to select a signal region? (double click on your selection to confirm)', 'Signal', 'Yes', 'No', 'Yes');
        if strcmp(ansSnROI, 'Yes')
            [ snRoiContour, snRoiMask ] = RoiSelect(selectionShape, 'r', size(img));
            snROIs = { snRoiContour, snRoiMask };
            text(max(snRoiContour(:,1)), max(snRoiContour(:,2)), num2str(1), 'Color', 'r');
            
            contSelect = true;
            while(contSelect)
                ansCont = questdlg('Do you need to select more signal regions?', 'Multiple Signal ROIs', 'Yes', 'No', 'No');
                if strcmp(ansCont, 'No')
                    contSelect = false;
                end
                if contSelect
                    [ snRoiContour, snRoiMask ] = RoiSelect(selectionShape, 'r', size(img));
                    snROIs(size(snROIs,1)+1, 1:2) = { snRoiContour, snRoiMask };
                    text(max(snRoiContour(:,1)), max(snRoiContour(:,2)), num2str(size(snROIs,1)), 'Color', 'r');
                end
            end
        else
            snROIs = { };
        end
        
        
        if ~isempty(fileDirBackup)
            
            if fileDir == 0
                fileDir = fileDirBackup;
            end
            
            if ~isempty(bgROIs)
                for i = 1 : size(bgROIs, 1)
                    bgROI = bgROIs(i, :);
                    save(fullfile(fileDir, [ 'ROI_profile_bg' num2str(i) '.mat' ]), 'bgROI');
                end
            end
            
            if ~isempty(snROIs)
                for i = 1 : size(snROIs, 1)
                    snROI = snROIs(i, :);
                    save(fullfile(fileDir, [ 'ROI_profile_sn' num2str(i) '.mat' ]), 'snROI');
                end
            end
            
        end
        
    end

    
    
end

