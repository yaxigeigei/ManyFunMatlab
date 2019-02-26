function result = MakeEventTimesTable(et, varargin)
% 

% Handle user inputs
p = inputParser();
p.addRequired('et', @(x) isnumeric(x) || iscell(x));
p.addParameter('DelimiterTimes', [], @isnumeric);
p.addParameter('VariableNames', []);
p.addParameter('FillEmpty', [], @isnumeric);
p.addParameter('Verbose', true, @islogical);
p.parse(et, varargin{:});
delimiterTimes = p.Results.DelimiterTimes;
varNames = p.Results.VariableNames;
emptyVal = p.Results.FillEmpty;
isVerbose = p.Results.Verbose;

if ~iscell(et)
    et = {et};
end

if isempty(varNames)
    varNames = arrayfun(@(x) ['event' num2str(x)], 1:size(et,2), 'Uni', false);
end

if ~isempty(delimiterTimes)
    % Delimit event times by delimiter times
    if isVerbose
        fprintf('Delimit event times using %d delimiter times\n', numel(delimiterTimes));
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
for i = 1 : numel(etEpoch)
    if isempty(etEpoch{i})
        etEpoch{i} = emptyVal;
    end
end
for i = 1 : numel(etPre)
    if isempty(etPre{i})
        etPre{i} = emptyVal;
    end
end

% Put data into tables
result.etTable = cell2table(etEpoch, 'VariableNames', varNames);
result.preDelimEtTable = cell2table(etPre, 'VariableNames', varNames);
result.etInfo.delimiterTimes = delimiterTimes;

end