function avgArray = AvgRunning( dataArray, binSize )
%AVGRUNNING Summary of this function goes here
%   Detailed explanation goes here

leftOffset = floor((binSize - 1) / 2);
rightOffset = ceil((binSize - 1) / 2);

avgArray = zeros(size(dataArray));

for i = 1 : size(dataArray, 1)
    
    beginIdx = i - leftOffset;
    if beginIdx < 1
        beginIdx = 1;
    end
    
    endIdx = i + rightOffset;
    if endIdx > size(dataArray, 1)
        endIdx = size(dataArray, 1);
    end
    
    avgArray(i,:) = nanmean(dataArray(beginIdx:endIdx,:));
    
end



end

