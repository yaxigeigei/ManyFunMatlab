function [ varargout ] = dFluorCalc( img, varargin )
%Compute actual, baseline fluorescence intensities and deltaF/F of the
%input image series or of its ROI regions
% 
% [ baselineF, fRatio ] = dFluorCalc( img, baselineWin )
% [ actualF, baselineF, fRatio ] = dFluorCalc( img, baselineWin, ROI )
% 
% 
%Input discription
% 
% img:
% An three dimensional array in double
% 
% baselineWin:
% Range of frames taken as baseline fluorescence intensity. (e.g. [ 1 10 ])
% (default is 'auto' under which the program takes the average of lowest 15% intensity as baseline)
% 
% ROI (optional):
% A cell array. One ROI for each row in which the first column contains
% coordinates of its contour points and the second column is its logic mask.
%
%
%Output discription
%
% When input ROI is not given: [ baselineF, fRatio ]
% Return maps of baseline fluorescence intensities and deltaF/F of the
% whole image, respectively (resulting in arrays with the same size of 'img')
%
% When input ROI is available: [ actualF, baselineF, fRatio ]
% Return arrays of actual, baseline fluorescence intensities and
% deltaF/F, respectively. Each row corresponds to one ROI.

    
    
    baselineWin = varargin{1};
    [ imgHeight, imgWidth, imgFrames ] = size(img);
    
    
    % Get Blind dFluor Map
    if length(varargin) == 1
        
        baselineF = zeros(imgHeight, imgWidth, imgFrames, class(img));
        fRatio = zeros(imgHeight, imgWidth, imgFrames, class(img));
        
        for i = 1 : imgHeight
            for j = 1 : imgWidth
                [ ~, baselineF(i,j,:), fRatio(i,j,:) ] = dFluorComputeF(img(i,j,:), [ ], baselineWin);
            end
        end
        
        varargout{1} = baselineF;
        varargout{2} = fRatio;
        
    end


    % Get ROI based dFluor trace(s)
    
    if length(varargin) == 2
        
        ROI = varargin{2};
        
        actualF = zeros(size(ROI, 1), imgFrames, class(img));
        baselineF = zeros(size(ROI, 1), imgFrames, class(img));
        fRatio = zeros(size(ROI, 1), imgFrames, class(img));
        
        for i = 1 : size(ROI, 1)
            [ actualF(i,:), baselineF(i,:), fRatio(i,:) ] = dFluorComputeF(img, ROI{i,2}, baselineWin);
        end
        
        varargout{1} = actualF;
        varargout{2} = baselineF;
        varargout{3} = fRatio;
        
    end

    



        
end


%     difference = zeros(imgHeight, imgWidth, imgFrames);
%     for i = 2 : imgFrames
%         difference(:,:,i) = img(:,:,i) - img(:,:,1);
%     end