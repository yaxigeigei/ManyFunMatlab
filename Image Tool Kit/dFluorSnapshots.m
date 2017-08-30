function [  ] = dFluorSnapshots(  )
%This function is typically used as a script


    %% Information

    fileDir = 'E:\神社\Lab\Data of Imaging\Two-photon\BDNF with 王良\20149304 pH-luorin high K+ B2 newborn\slice4';
    fileName = '';

    info_Date = '20140304';
    info_Discription = 'highK_bath_B2_newborn_brightest4';


    %% Import Image

    if exist(fullfile(fileDir, fileName), 'file') ~= 2
        [ fileName, fileDir, ~ ] = uigetfile({'*.tif;*.tiff'}, 'Select a Image File', [fileDir '\']);
    end
    img_Original = TiffImport(fullfile(fileDir, fileName));
    img_Original = im2double(img_Original);


    %% Generate Results

    [ ROI_signal, ROI_background, img_Projection ] = ...
        dFluorManualROI(img_Original, 'freehand', fileDir);   % 'rectangle', 'freehand', 'ellipse'
    close
    
    [ data_actualF, ~, ~ ] = dFluorCalc(img_Original, 'auto', ROI_signal);
    
    figure('Color', 'w');
    plot(data_actualF', 'LineWidth', 2);
    xlim([ 1 size(data_actualF,2) ]);
    ylim([ 0 1.2*max(max(data_actualF)) ]);
    
    figure('Name', 'ROI');
    imagesc(img_Projection);
    axis equal tight
    colormap gray
    RoiShow(ROI_signal, 'r', false);
    
    
    
    %% Save Results
    
    resultDir = 'E:\神社\Lab\Program\MATLAB\BDNF';
    resultName = [ info_Date '_' info_Discription ];
    
    save(fullfile(resultDir, resultName), 'img_*', 'info_*', 'ROI_*', 'data_*');
    
    
    
    
    
end

