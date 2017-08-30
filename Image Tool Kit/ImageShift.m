function shift = ImageShift( moving, template, varargin )
%IMAGESHIFT Summary of this function goes here
%   Detailed explanation goes here

    
    
    p = inputParser;
    
    addParamValue(p, 'fineFactor', 1, @isnumeric);
    addParamValue(p, 'roughFactor', size(template,1)/64);
    
    parse(p, varargin{:});
    
    fineFactor = p.Results.fineFactor;
    roughFactor = p.Results.roughFactor;
    
    
    templateRough = ImageShrink(template, roughFactor);
    movingRough = ImageShrink(moving, roughFactor);
    
    templateFine = StackInterp(template, fineFactor, fineFactor);
    movingFine = StackInterp(moving, fineFactor, fineFactor);
    
    
    
    % 0.125X Screen
    msRough = XMeanSquare2(movingRough, templateRough, size(templateRough));
    [ ~, index ] = min(msRough(:));
    [ roughFocus(1), roughFocus(2) ] = ind2sub(size(msRough), index);
    
%     figure
%     subplot(1,3,1)
%     imagesc(msRough);
%     axis equal tight;

    
    
    % 1X Zoom
    nextFocus = roughFocus*roughFactor;
    msNorm = XMeanSquare2(moving, template, [ roughFactor*2, roughFactor*2 ], nextFocus);
    [ ~, index ] = min(msNorm(:));
    [ normFocus(1), normFocus(2) ] = ind2sub(size(msNorm), index);
    shiftNorm = normFocus - size(msNorm)/2; % Shift in 1X sub-image
    normFocus = nextFocus + shiftNorm; % Focus position in whole 1X image
    
    if fineFactor <= 1
        shift = normFocus - size(templateFine)/2; % Shift in whole 3X image
        shift = round(shift / fineFactor); % Shift in whole 1X image
    end
    
%     subplot(1,3,2)
%     imagesc(msNorm);
%     axis equal tight;
    

    
    if fineFactor > 1
        % 3X Zoom
        nextFocus = normFocus*fineFactor; % Focus position in whole 3X image
        msFine = XMeanSquare2(movingFine, templateFine, [ fineFactor*3, fineFactor*3 ], nextFocus);
        [ ~, index ] = min(msFine(:));
        [ fineFocus(1), fineFocus(2) ] = ind2sub(size(msFine), index);
        shiftFine = fineFocus - size(msFine)/2; % Shift in 3X sub-image
        fineFocus = nextFocus + shiftFine;
        
        shift = fineFocus - size(templateFine)/2; % Shift in whole 3X image
        shift = round(shift / fineFactor); % Shift in whole 1X image

%         subplot(1,3,3)
%         imagesc(msFine);
%         axis equal tight;
    end
    
    
    
    
end

