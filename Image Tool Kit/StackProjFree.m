function imProj = StackProjFree( stack, thetas )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here



% Get Stack Information
[ imh, imw, iml ] = size(stack);
numPxl = numel(stack);

% Create Full Space Coordinate Matrix
fullSpace = zeros([ size(stack), 3 ]);
[ fullSpace(:,:,:,2), fullSpace(:,:,:,1), fullSpace(:,:,:,3) ] = meshgrid(1:imw, 1:imh, 1:iml);
fullSpace = permute(fullSpace, [ 4, 1, 2, 3 ]);
fullSpace = reshape(fullSpace, [ 3, numPxl ]);

% Preallocate Output Matrix
nimh = ceil(sqrt(imh^2+iml^2));
nimw = ceil(sqrt(imw^2+iml^2));
imProj = zeros([ nimh, nimw, length(thetas) ], 'uint16');


for i = 1 : length(thetas)
    
    % Calculate Subspace and Transformation Matrix
    subBases = [ 0, 1, 0; cos(thetas(i)), 0, 0; sin(thetas(i)), 0, 0 ];
    tfMat = subBases*subBases';
    
    % Transformation
    subSpace = tfMat * double(fullSpace);
    subSpace = uint16(subSpace);

    % Adjust Image Position
    rShift = ceil((nimh - max(subSpace(1,:)) - 1)/2);
    cShift = ceil((nimw - max(subSpace(2,:)) - 1)/2);
    subSpace(1,:) = subSpace(1,:) + rShift;
    subSpace(2,:) = subSpace(2,:) + cShift;
    
    % Fast Projection
    subSpace = reshape(subSpace, [ 3, imh, imw, iml ]);
    layerProj = zeros([ nimh, nimw, iml ], 'uint16');
    for j = 1 : iml
        for k = 1 : imw
            layerProj(1+rShift:end-rShift, subSpace(2,1,k,j), j) = ...
                max(layerProj(1+rShift:end-rShift, subSpace(2,1,k,j), j), stack(:,k,j));
        end
    end
    imProj(:,:,i) = max(layerProj, [ ], 3);
    
%     % Universal Projection
%     for j = 1 : numPxl
%         imProj(subSpace(1,j), subSpace(2,j), i) = ...
%             max(imProj(subSpace(1,j), subSpace(2,j), i), stack(j));
%     end
    
end



end



