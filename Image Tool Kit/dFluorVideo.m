function dFluorVideo()
%This function is typically used as a script
    
    
    fileDir = 'E:\神社\Lab\Data of Imaging\Two-photon\In-vivo Calcium with Motor Group\Data from 沈志明\cat';
    fileName = '';
    
    info_Date = '';
    info_Discription = 'cat-test';
    
    info_frameRate = 15.5;
    info_baselineWindow = [ 1 10 ];
    info_edgePixels = 16;
    

%% Import Image

    if exist(fullfile(fileDir, fileName), 'file') ~= 2
        [ fileName, fileDir, ~ ] = uigetfile({'*.tif;*.tiff'}, 'Select a Image File', [fileDir '\']);
    end
    img = TiffImport(fullfile(fileDir, fileName));
    img = im2single(img);
    img = ImageContrast(img);
    
    
%% Preprocessing

    img = FiltGaussian(img);
    [ img, ~ ] = dFluorBackgroundEven(img, round(size(img,1)/25));
    
    % Auto
    img(img < 2*mean(mean(mean(img)))) = 2*mean(mean(mean(img)));
    
    TiffExport(fullfile(fileDir,'cat-stack_cut'), img);
    
    
%% ROI Selection
    
    [ ROI_signal, ~, ~ ] = dFluorManualROI(img, 'freehand', fileDir);
    
    % Auto
    [ ~, img_fRatio ] = dFluorCalc(img, 'auto');
    img_FfRatio = FiltMedian(img_fRatio);
    [ ROI_signal, img_roiFMasks, img_roiCMasks ] = dFluorAutoROI(img_FfRatio, info_edgePixels, 40, 100);
    ImageShow(img_roiFMasks);
    ImageShow(img_roiCMasks);
    ImageShow(img, ROI_signal);
    
    
%% Calculating Results
    
    [ data_actualF, data_baselineF, data_fRatio ] = dFluorCalc(img, info_baselineWindow, ROI_signal);
    
    
%% Plot Video Results
    
    ImageShow(img, ROI_signal);
    dFluorPlotComboTraces(info_frameRate, data_actualF, data_baselineF, data_fRatio);
    dFluorPlotPopTraces(info_frameRate, data_fRatio);

    
%% Save Analysis Results

    resultDir = fileDir;
    save(fullfile(resultDir, 'maunual_result_package'), 'img*', 'ROI_*', 'info_*', 'data_*');
    


 
end

