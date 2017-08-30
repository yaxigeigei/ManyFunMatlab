function [ msArray ] = XMeanSquare2( fixed, moving, outputSize, focus )
%MEANSQUARE2 Summary of this function goes here
%   Detailed explanation goes here
    
    
    
    [ imHeight, imWidth ] = size(fixed);
    
    if nargin < 3
        outputSize = [ 2*imHeight-1, 2*imWidth-1 ];
    end
    
    msArray = zeros(outputSize);
    [ msHeight, msWidth ] = size(msArray);
    
    
    if nargin < 4
        focus = [ ceil(imHeight/2), ceil(imWidth/2) ];
    end
    
    si = focus(1) - ceil(imHeight/2);
    sj = focus(2) - ceil(imWidth/2);
    
    iCenter = ceil(msHeight/2) - si;
    jCenter = ceil(msWidth/2) - sj;
    
    
    
    for j = 1 : msWidth
        
        dj = j - jCenter;
        
        jStartF = max(1, 1 + dj);
        jEndF = min(imWidth, imWidth + dj);
        
        jStartM = max(1, 1 - dj);
        jEndM = min(imWidth, imWidth - dj);
        
        for i = 1 : msHeight
            
            di = i - iCenter;
            
            iStartF = max(1, 1 + di);
            iEndF = min(imHeight, imHeight + di);
            
            iStartM = max(1, 1 - di);
            iEndM = min(imHeight, imHeight - di);
            
            fSub = fixed(iStartF : iEndF, jStartF : jEndF);
            mSub = moving(iStartM : iEndM, jStartM : jEndM);
            msArray(i,j) = MeanSquare2Norm(mSub, fSub);
            
        end
        
    end



end