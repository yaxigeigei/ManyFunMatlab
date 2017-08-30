function ms = XMeanSquareV2( v1, v2 )
%MEANSQUARE2 Summary of this function goes here
%   Detailed explanation goes here
    
    
    [ r1, r2 ] = XIndexing(v1, v2);
    
    ms = zeros(size(r1,1), 1);
    
    for i = 1 : length(ms)
        ms(i) = MeanSquare2Norm(v1(r1(i,1):r1(i,2)), v2(r2(i,1):r2(i,2)));
    end
    
    
    
end





