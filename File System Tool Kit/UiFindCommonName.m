function [ folderName, commonName, ext ] = UiFindCommonName( folderName )
%FINDFILES Summary of this function goes here
%   Detailed explanation goes here

if nargin < 1
    try
        load('lastimeDir');
    catch
        folderName = '';
    end
end

[ fileName, folderName ] = uigetfile({'*.*'}, 'Select Data', [folderName '\']);

if folderName
    [ ~, commonName, ext ] = fileparts(fileName);
    save('lastimeDir.mat', 'folderName');
else
    commonName = '';
    ext = '';
end


end

