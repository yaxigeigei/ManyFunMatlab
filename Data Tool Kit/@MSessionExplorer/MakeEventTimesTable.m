function [tb, preTb] = MakeEventTimesTable(et, varargin)
% Organize event times data into an 'eventTimes' table. 
% Specifically, each row of this table contains data from a unique period of time, or 
% an epoch (e.g. a trial). Each column contains timestamps of a given event. Epochs 
% that do not have any data are filled with NaN (or an empty object). 
% 
%   [tb, preTb] = MSessionExplorer.MakeEventTimesTable(et)
%   [tb, preTb] = MSessionExplorer.MakeEventTimesTable(..., 'DelimiterTimes', [])
%   [tb, preTb] = MSessionExplorer.MakeEventTimesTable(..., 'VariableNames', [])
%   [tb, preTb] = MSessionExplorer.MakeEventTimesTable(..., 'Verbose', true)
% 
% Inputs
%   et      1) A numeric vector of timestamps of an event. 
%           2) A vector of objects of MSessionExplorer.Event class or superclass. 
%           2) A 1-D cell array of 1) or 2). Each cell is treated as a different event. 
%           3) A 2-D cell array of 1) or 2) where columns are different events and rows 
%              are different epochs. 'DelimeterTimes' is not supported for this input. 
%   'DelimiterTimes'
%           A vector of time values indicating when to cut data into different epochs. 
%           Data before the first delimiter time is separately stored in preTb. 
%           Timestamps will be converted to relative times wrt the preceeding delimiter 
%           time. The default value is [] which performs no delimiting. 
%   'VariableNames'
%           A cell array of column names for the output tables. The default names are 
%           {'event1', 'event2', 'event3', ...}. 
%   'Verbose'
%           A logical value that controls the display of progress. Default is true. 
% Outputs
%   tb      A table that works in MSessionExplorer objects as an 'eventTimes' table. 
%   preTb   Similar to tb but includes data before the first delimiter time, if any. 

% Parse inputs
p = inputParser();
p.addRequired('et', @(x) isnumeric(x) || iscell(x));
p.addParameter('DelimiterTimes', [], @isnumeric);
p.addParameter('VariableNames', []);
p.addParameter('Verbose', true, @islogical);
p.parse(et, varargin{:});
delimiterTimes = p.Results.DelimiterTimes;
varNames = p.Results.VariableNames;
isVerbose = p.Results.Verbose;

% Unify to cell array
if ~iscell(et)
    et = {et};
end

% Ensure column vectors
for i = 1 : numel(et)
    if ~iscolumn(et{i})
        et{i} = et{i}(:);
    end
end

% Make variable names
if isempty(varNames)
    varNames = arrayfun(@(x) ['event' num2str(x)], 1:size(et,2), 'Uni', false);
end

if ~isempty(delimiterTimes)
    % Delimit event times by delimiter times
    if isVerbose
        fprintf('Delimit event times data using %d delimiter times\n', numel(delimiterTimes));
    end
    assert(isvector(et), 'Cannot delimit event times in cell array with more than one dimension');
    et = et(:)';
    
    for i = numel(et) : -1 : 1
        for j = numel(delimiterTimes) : -1 : 1
            mask = et{i} >= delimiterTimes(j);
            etEpoch{j,i} = et{i}(mask) - delimiterTimes(j);
            et{i}(mask) = [];
        end
        et{i}(isnan(et{i})) = [];
    end
    etPre = et;
    
elseif isVerbose
    % Each cell is an epoch
    if isVerbose
        disp('Delimiter is not specified. Each cell is treated as an epoch.');
    end
    etEpoch = et;
    etPre = cell(0, size(et,2));
end

% Fill empty cells
etEpoch = FillEmpty(etEpoch);
etPre = FillEmpty(etPre);

% Put data into tables
tb = cell2table(etEpoch, 'VariableNames', varNames);
preTb = cell2table(etPre, 'VariableNames', varNames);

end

function C = FillEmpty(C)
% Fill empty cells
% If all values in a cloumn are numeric, fill empty cells with NaN
% Otherwise fill with an object from the non-numeric class (constructed with no input)

for j = 1 : size(C,2)
    % Cache values
    indEpt = find(cellfun(@isempty, C(:,j)));
    if isempty(indEpt)
        continue;
    end
    isNum = cellfun(@isnumeric, C(:,j));
    objIdx = find(~isNum, 1);
    
    for i = indEpt'
        if isempty(objIdx)
            % Fill NaN
            C{i,j} = NaN;
        else
            % Fill an empty object
            className = class(C{objIdx,j});
            C{i,j} = eval(className);
        end
    end
end

end


