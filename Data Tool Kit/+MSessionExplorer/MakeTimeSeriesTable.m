function result = MakeTimeSeriesTable(t, sig, varargin)
% 

% Handle user inputs
p = inputParser();
p.addRequired('t', @(x) isnumeric(x) || iscell(x));
p.addRequired('sig', @(x) isnumeric(x) || iscell(x));
p.addParameter('DelimiterTimes', [], @isnumeric);
p.addParameter('VariableNames', []);
p.addParameter('Verbose', true, @islogical);
p.parse(t, sig, varargin{:});
delimiterTimes = p.Results.DelimiterTimes;
varNames = p.Results.VariableNames;
isVerbose = p.Results.Verbose;

if ~iscell(t)
    t = {t};
end
t = t(:);

if ~iscell(sig)
    sig = {sig};
end

% Make variable names
if isempty(varNames)
    varNames = [{'time'}, arrayfun(@(x) ['series' num2str(x)], 1:size(sig,2), 'Uni', false)];
elseif ~strcmp('time', varNames{1})
    varNames = [{'time'}; varNames(:)];
end

if isempty(delimiterTimes)
    % Each cell is an epoch
    if isVerbose
        disp('Delimiter is not specified. Each cell is treated as an epoch.');
    end
    assert(numel(t) == size(sig,1), ...
        'There are %d vectors of timestamp and but %d epochs of signal', numel(t), size(sig,1));
    
    tEpoch = t;
    sigEpoch = sig;
    tPre = cell(0,1);
    sigPre = cell(0, size(sig,2));
    
else
    % Delimit signals by delimiter times
    if isVerbose
        fprintf('Delimit signals using %d delimiter times\n', numel(delimiterTimes));
    end
    assert(isvector(sig), 'Cannot delimit signals in cell array with more than one dimension');
    sig = sig(:)';
    assert(numel(t) == 1, 'Cannot delimit more than one series of timestamp');
    t = t{1};
    assert(all(diff(t) > 0), 'Timestamps must be monotonically increasing');
    
    % Group timestamps
    for i = numel(delimiterTimes) : -1 : 1
        mask = t >= delimiterTimes(i);
        tEpoch{i,1} = t(mask) - delimiterTimes(i);
        t(mask) = [];
    end
    tPre = t;
    
    % Group data
    L = cellfun(@numel, [tPre; tEpoch]);
    sig = cellfun(@(x) mat2cell(x, L), sig, 'Uni', false);
    sig = cat(2, sig{:});
    sigPre = sig(1:numel(tPre),:);
    sigEpoch = sig(numel(tPre)+1:end,:);
end

% Put data into table
result.tsTable = cell2table([tEpoch, sigEpoch], 'VariableNames', varNames);
result.preDelimTsTable = cell2table([tPre, sigPre], 'VariableNames', varNames);
result.info.delimiterTimes = delimiterTimes;

end

