classdef MSessionExplorer < handle
    %MSessionExplorer is a data container that makes data munging easy
    
    properties(Constant)
        supportedTableTypes = ... % 'eventTimes', 'eventValues', 'timeSeries' (read-only)
            {'eventTimes', 'eventValues', 'timeSeries'};
    end
    properties(GetAccess = public, SetAccess = private)
        tot;                    % A table of data tables, reference times and related metadata (read-only)
    end
    properties
        userData;               % A struct where user stores arbitrary data
    end
    properties(Dependent)
        tableNames;             % A cell array of table names (read-only)
        isEventTimesTable;      % Whether each table is eventTimes table (read-only)
        isEventValuesTable;     % Whether each table is eventValues table (read-only)
        isTimesSeriesTable;     % Whether each table is timeSeries table (read-only)
        numEpochs;              % The number of epochs (read-only)
    end
    properties
        epochInd;               % Indices of epoch
        isVerbose = true;       % Whether or not to display progress and some warnings
    end
    properties(Access = private)
        originalTrialInd;       % for backward compatibility
    end
    properties(Dependent, Hidden)
        numTrials;              % for backward compatibility
    end
    
    % These methods are exposed to users
    methods
        % Constructor
        function this = MSessionExplorer()
            % Constructor of MSessionExplorer
            % 
            %   se = MSessionExplorer()
            
            % Initialize table of tables
            totHeaders = {'tableType', 'tableData', 'referenceTime'};
            this.tot = cell2table(cell(0, length(totHeaders)), 'VariableNames', totHeaders);
        end
        
        % Property IO
        function val = get.tableNames(this)
            val = this.tot.Properties.RowNames;
        end
        function val = get.isEventTimesTable(this)
            val = strcmp(this.tot.tableType, this.supportedTableTypes{1})';
        end
        function val = get.isEventValuesTable(this)
            val = strcmp(this.tot.tableType, this.supportedTableTypes{2})';
        end
        function val = get.isTimesSeriesTable(this)
            val = strcmp(this.tot.tableType, this.supportedTableTypes{3})';
        end
        function val = get.numEpochs(this)
            if ~isempty(this.tot)
                val = height(this.tot.tableData{1});
            else
                val = 0;
            end
        end
        function set.epochInd(this, val)
            assert(all(val > 0 & mod(val,1) == 0), 'Epoch indices must be positive integers');
            assert(numel(val) == numel(unique(val)), 'Epoch indices must be unique numbers');
            assert(numel(val) == this.numEpochs, 'The number of indices must match the number of epochs');
            this.epochInd = val(:);
        end
        function val = get.epochInd(this)
            if isempty(this.epochInd) && this.numEpochs > 0 && ~isempty(this.originalTrialInd)
                warning(['originalTrialInd property will be removed in a future version. ' ...
                    'Use epochInd instead and consider saving the updated object.']);
                this.epochInd = this.originalTrialInd;
            end
            val = this.epochInd;
        end
        function set.isVerbose(this, val)
            assert(islogical(val) && isscalar(val), 'isVerbose must be a logical scalar');
            this.isVerbose = val;
        end
        function val = get.numTrials(this)
            warning('numTrials property will be removed in a future version. Use numEpochs instead.');
            val = this.numEpochs;
        end
        
        % Content Management
        function Preview(this, varargin)
            % Print a summary of the main content of this object or the begining of specific table(s)
            %   
            %   Preview()
            %   Preview(tbName1, tbName2, tbName3, ...)
            % 
            % Inputs
            %   tbNameN         Name of a data table to display. 
            %                   If not specified, this method will display tot and userData. 
            
            if nargin < 2
                disp('tot');
                disp(this.tot);
                disp('userData');
                disp(this.userData);
            else
                this.IValidateTableNames(varargin, true);
                for i = 1 : numel(varargin)
                    head(this.tot{varargin{i}, 'tableData'}{1})
                end
            end
        end
        
        function se = Duplicate(this, varargin)
            % Make a deep copy of the current object
            % 
            %   se = Duplicate()
            %   se = Duplicate(tbNames)
            %   se = Duplicate(tbNames, isUserData)
            % 
            % Inputs
            %   tbNames             A string or cell array of table name(s) specifying which table(s) to include.
            %                       An empty array (default) indicates all tables. 
            %   isUserData          A logical variable indicating whether or not to copy userData. Default true. 
            % Output
            %   se                  The new MSessionExplorer object
            
            % Handle user input
            p = inputParser();
            p.addOptional('tbNames', []);
            p.addOptional('isUserData', true, @islogical);
            p.parse(varargin{:});
            tbNames = p.Results.tbNames;
            isUserData = p.Results.isUserData;
            
            if isempty(tbNames)
                tbNames = this.tableNames;
            end
            this.IValidateTableNames(tbNames, true);
            
            % Find table indices
            tbNames = cellstr(tbNames);
            tbInd = cellfun(@(x) find(strcmp(x, this.tableNames)), tbNames);
            
            % Copying
            se = MSessionExplorer();
            for i = tbInd(:)'
                se.SetTable(this.tableNames{i}, this.tot.tableData{i}, this.tot.tableType{i}, this.tot.referenceTime{i});
            end
            if isUserData
                se.userData = this.userData;
            end
            se.epochInd = this.epochInd;
        end
        
        function se = Merge(this, varargin)
            % Combine multiple MSessionExplorer objects into one by vertically concatenating data
            % 
            %   se = Merge(se1, se2, se3, ...)
            % 
            % Inputs
            %   se1, se2, se3, ...      Arbitrary number of MSessionExplorer objects. 
            % Output
            %   se                      The merged MSessionExplorer object
            
            % SEs to merge
            seArray = [this; cat(1, varargin{:})];
            
            % Output SE
            se = MSessionExplorer();
            
            % Concatenate and set each table
            for i = 1 : numel(this.tableNames)
                tbName = this.tableNames{i};
                tbData = arrayfun(@(x) x.GetTable(tbName), seArray, 'Uni', false);
                tbData = cat(1, tbData{:});
                refTime = arrayfun(@(x) x.GetReferenceTime(tbName), seArray, 'Uni', false);
                refTime = cat(1, refTime{:});
                se.SetTable(tbName, tbData, this.tot.tableType{i}, refTime);
            end
            
            % Add incremented epoch indices
            cumNumEp = cumsum(arrayfun(@(x) x.numEpochs, seArray));
            epInd = arrayfun(@(x) x.epochInd, seArray, 'Uni', false);
            epInd(2:end) = cellfun(@(x,y) x + y, num2cell(cumNumEp(1:end-1)), epInd(2:end), 'Uni', false);
            se.epochInd = cat(1, epInd{:});
            
            % Add user data
            se.userData = arrayfun(@(x) x.userData, seArray, 'Uni', false);
            try
                se.userData = cat(1, se.userData{:});
            catch
                warning('Not all userData have the same fields thus were stored into a cell array');
            end
        end
        
        function s = ToStruct(this)
            % Convert the current MSessionExplorer object and all tables to structures
            % 
            %   s = ToStruct()
            % 
            % Output
            %   s       The output structure with the following fields
            %           tableName       a cell array of table names
            %           tableType       a cell array of table types
            %           referenceTime   a cell array of referece time
            %           tableData       a cell array of structs whose fields are each table's variables 
            %                           (output from MATLAB table2struct function)
            %           userData        same as userData property
            
            s.tableName = this.tableNames;
            s.tableType = this.tot.tableType;
            s.referenceTime = this.tot.referenceTime;
            s.tableData = cellfun(@(x) table2struct(x), this.tot.tableData, 'Uni', false);
            s.userData = this.userData;
        end
        
        function SetTable(this, tbName, tb, varargin)
            % Add or update a data table
            % 
            %   SetTable(tbName, tb)
            %   SetTable(tbName, tb, tableType)
            %   SetTable(tbName, tb, tableType, referenceTime)
            %
            % Inputs:
            %   tbName              A string of the name of table to add or update. 
            %   tb                  The table variable.
            %   tableType           Either 'eventTimes', 'eventValues', or 'timeSeries' indicating 
            %                       the type of data table. It is required when adding a new table 
            %                       but is ignored when updating an existing table.  
            %   referenceTime       A numeric vector where each element stores the absolute time
            %                       of the zero relative time of each epoch. Default is empty. 
            
            % Handle user input
            p = inputParser();
            p.addRequired('tbName', @ischar);
            p.addRequired('tb', @istable);
            p.addOptional('tableType', [], @(x) any(strcmp(x, this.supportedTableTypes)));
            p.addOptional('referenceTime', [], @isnumeric);
            p.parse(tbName, tb, varargin{:});
            tbType = p.Results.tableType;
            refTimes = p.Results.referenceTime;
            
            % Check epoch number conflict
            if ~isempty(this.tableNames)
                assert(size(tb,1) == this.numEpochs, ...
                    'The new table cannot be added since it has %d rows whereas existing table has %d.', ...
                    size(tb,1), this.numEpochs);
            end
            
            % Set table
            if ismember(tbName, this.tableNames)
                % Replace an existing table
                this.tot{tbName, 'tableData'}{1} = tb;
                if ~isempty(tbType)
                    this.tot{tbName, 'tableType'}{1} = tbType;
                end
            else
                isFirst = isempty(this.tot);
                
                % Add a new table
                assert(~isempty(tbType), 'Table type is required for adding a new table');
                totHeaders = this.tot.Properties.VariableNames;
                totRow = cell2table({tbType, tb, []}, 'VariableNames', totHeaders, 'RowNames', {tbName});
                this.tot = [this.tot; totRow];
                
                % Initialize epoch indices
                if isFirst
                    this.epochInd = (1 : height(tb))';
                end
            end
            
            % Set reference time
            if ~isempty(refTimes)
                this.SetReferenceTime(refTimes, tbName);
            end
        end
        
        function varargout = GetTable(this, varargin)
            % Return specific data table(s)
            % 
            %   [tb1, tb2, tb3, ...] = GetTable(tbName1, tbName2, tbName3, ...)
            %
            % Input
            %   tbNameN         The name of table to return. 
            % Output
            %   tbN             The table data. 
            
            this.IValidateTableNames(varargin, true);
            for i = numel(varargin) : -1 : 1
                varargout{i} = this.tot{varargin{i}, 'tableData'}{1};
            end
        end
        
        function RemoveTable(this, varargin)
            % Remove specific data tables
            % 
            %   RemoveTable(tbName1, tbName2, tbName3, ...)
            %
            % Input
            %   tbNameN         A string or cell array of strings indicating the name of table(s) to remove.
            
            this.IValidateTableNames(varargin, true);
            this.tot(varargin,:) = [];
            if isempty(this.tot)
                this.epochInd = [];
            end
        end
        
        function SetReferenceTime(this, rt, tbNames)
            % Set reference times to table(s)
            % 
            %   SetReferenceTime(rt)
            %   SetReferenceTime(rt, tbNames)
            %
            % Inputs
            %   rt              A numeric vector where each element stores the absolute time of the 
            %                   zero time in each epoch. Or use [] to clear existing reference times. 
            %   tbNames         A string or cell array of the name(s) of table which reference time is 
            %                   set to. The default is empty and rt is set to all eligible tables.
            
            rt = rt(:);
            if ~all(diff(rt) > 0)
                warning('The reference time is not monotonically increasing');
            end
            
            if nargin < 3
                tbNames = this.tableNames(~this.isEventValuesTable);
            end
            this.IValidateTableNames(tbNames, true);
            
            tbNames = cellstr(tbNames);
            for i = 1 : numel(tbNames)
                assert(this.numEpochs == numel(rt) || isempty(rt), ...
                    'The reference time has %d elements which does not match the %d epochs', ...
                    numel(rt), this.numEpochs);
                
                if ~this.isEventValuesTable(strcmp(tbNames{i}, this.tableNames))
                    this.tot{tbNames{i}, 'referenceTime'}{1} = rt;
                else
                    warning('''%s'' is an eventValue table and reference time is not applicable', tbNames{i});
                end
            end
        end
        
        function varargout = GetReferenceTime(this, varargin)
            % Return reference times
            % 
            %   rt = GetReferenceTime()
            %   [rt1, rt2, rt3, ...] = GetReferenceTime(tbName1, tbName2, tbName3, ...)
            %
            % Input
            %   tbNameN         The name of a table which the reference time is associated with. If not 
            %                   specified, the first availble reference time will be returned. 
            % Output
            %   rtN             A vector of reference times. 
            
            if nargin > 1
                this.IValidateTableNames(varargin, true);
            else
                isRt = ~cellfun(@isempty, this.tot.referenceTime);
                assert(any(isRt), 'No table has reference time');
                varargin = this.tableNames(find(isRt,1));
            end
            
            for i = numel(varargin) : -1 : 1
                varargout{i} = this.tot{varargin{i}, 'referenceTime'}{1};
            end
        end
        
        function RemoveEpochs(this, ind2rm)
            % Remove specific epochs across all tables
            % 
            %   RemoveEpochs(ind2rm)
            %
            % Input
            %   ind2rm          Integer or logical indices of epochs to remove. 
            
            for k = 1 : height(this.tot)
                this.tot.tableData{k}(ind2rm,:) = [];
                if ~isempty(this.tot.referenceTime{k})
                    this.tot.referenceTime{k}(ind2rm) = [];
                end
            end
            this.epochInd(ind2rm) = [];
        end
        
        function colData = GetColumn(this, tbName, colNames)
            % Return specific columns from a data table
            %
            %   colData = GetColumn(tbName, colNames)
            %
            % Inputs
            %   tbName          The name of a table that contains the columns.
            %   colNames        Column name(s) as a string or a cell array of strings. 
            % Output
            %   colData         Requsted data. 
            
            this.IValidateTableName(tbName, true);
            colData = this.tot{tbName, 'tableData'}{1}{:,colNames};
        end
        
        function SetColumn(this, tbName, colNames, colData)
            % Add, update or delete column(s) in a data table
            %
            %   SetColumn(tbName, colNames, colData)
            %
            % Inputs
            %   tbName          The name of a table that contains the column. 
            %   colNames        Column name(s) as a string or a cell array of strings. 
            %   colData         Data to add or update. Use [] to delete. 
            
            tb = this.GetTable(tbName);
            colNames = cellstr(colNames);
            for i = 1 : numel(colNames)
                if isempty(colData)
                    % Remove variable
                    tb.(colNames{i}) = [];
                else
                    % Set values
                    tb.(colNames{i}) = colData(:,i);
                end
            end
            this.SetTable(tbName, tb);
        end
        
        % Data Operations
        function AlignTime(this, refEvent, refSourceTbName)
            % Align epochs by changing the origin of time in each epoch to the new reference event time
            % 
            %   AlignTime(refEvent)
            %   AlignTime(refEvent, refSourceTbName)
            %
            % Inputs
            %   refEvent            An event name in a eventTimes table or a numeric vector of reference times 
            %                       as zero time after alignment. Each element of the vector matches respective 
            %                       row in tables. 
            %   refSourceTbName     If refEvent is a string, then you must specify the name of eventTimes table 
            %                       where this event name should be found. This avoids ambiguity of the same 
            %                       variable name found in multiple tables. The default is empty. 
            
            % Check eligible tables
            indAlignable = find(this.isEventTimesTable | this.isTimesSeriesTable);
            assert(~isempty(indAlignable), 'Requires at least one eventTimes or timeSeries table to operate on');
            
            % Get reference event times
            if ischar(refEvent)
                assert(nargin == 3, 'Requires refSourceTableName to indicate where ''%s'' is in', refEvent);
                this.IValidateTableName(refSourceTbName, true);
                refEvent = this.tot{refSourceTbName, 'tableData'}{1}.(refEvent);
            end
            
            % Validate reference event times
            assert(isnumeric(refEvent) && isvector(refEvent), 'Reference times must be a numeric vector');
            assert(numel(refEvent) == this.numEpochs, ...
                'The number of reference times (%d) does not match the number of rows (%d) in data table.', ...
                numel(refEvent), this.numEpochs);
            refEvent = refEvent(:);
            
            % Align times
            for k = indAlignable
                if this.isEventTimesTable(k)
                    % For eventTimes table
                    tb = this.tot.tableData{k};
                    for i = 1 : size(tb, 2)
                        if isnumeric(tb{:,i})
                            % Numeric vector
                            tb{:,i} = tb{:,i} - refEvent;
                        elseif iscell(tb{:,i})
                            % Cell vector of numeric vectors
                            tb{:,i} = cellfun(@(x,r) x-r, tb{:,i}, num2cell(refEvent), 'Uni', false);
                        else
                            warning(fprintf('Failed to align data in the ''%s'' column of ''%s'' table\n', ...
                                tb.Properties.VariableNames{i}, this.tableName{k}));
                        end
                    end
                    this.tot.tableData{k} = tb;
                    
                elseif this.isTimesSeriesTable(k)
                    % For timeSeries table
                    this.tot.tableData{k}.time = cellfun(@(x,r) x-r, ...
                        this.tot.tableData{k}.time, num2cell(refEvent), 'Uni', false);
                end
                
                % Change referenceTime
                if ~isempty(this.tot.referenceTime{k})
                    this.tot.referenceTime{k} = this.tot.referenceTime{k} + refEvent;
                end
            end
        end
        
        function SortEpochs(this, ind)
            % Sort epochs across all tables
            % 
            %   SortEpochs()
            %   SortEpochs(ind)
            %
            % Input
            %   ind             Indices of the new order. If not specified, epochs are sorted back to the 
            %                   original order based on the epochInd property. 
            
            % Handle user inputs
            if nargin < 2
                [~, ind] = sort(this.epochInd);
            end
            ind = unique(ind(:), 'stable');
            assert(numel(ind) == this.numEpochs, ...
                'The number of unique indices (%d) does not equal to the number of epochs (%d)', ...
                numel(ind), this.numEpochs);
            
            % Process original indices
            this.epochInd = this.epochInd(ind);
            
            % Sort epochs
            for k = 1 : height(this.tot)
                this.tot.tableData{k} = this.tot.tableData{k}(ind,:);
                if ~isempty(this.tot.referenceTime{k})
                    this.tot.referenceTime{k} = this.tot.referenceTime{k}(ind);
                end
            end
        end
        
        function SliceSession(this, tSlice, refType)
            % Slice session into a different set of epochs
            % 
            %   SliceSession(this, tSlice, refType)
            %
            % Inputs
            %   tSlice          A vector of times at which slicing occurs. Note that data before the first slicing 
            %                   time will be irreversibly discarded. 
            %                   If refType is 'relative', tSlice can be a scalar which is used by all epochs, or 
            %                   a vector whose length matches the number epoch. 
            %                   If refType is 'absolute', tSlice can be a vector of arbitrary length. 
            %   refType         'absolute' or 'relative', indicating whether tSlice represents absolute times in 
            %                   the session or is relative to epoch reference times. If using 'absolute', existing 
            %                   epoch sorting, time alignment and all eventValues table will be lost due to an 
            %                   unknown relationship between current and new epochs. New data will sort epochs in 
            %                   temporal order and align epoch time to respective tSlice. 
            
            warning('This method is in beta version and should not be used for formal analysis');
            
            % Validate inputs
            indTable = find(this.isEventTimesTable | this.isTimesSeriesTable);
            assert(~isempty(indTable), 'Requires at least one eventTimes or timeSeries table to operate on');
            assert(isnumeric(tSlice) && isvector(tSlice), 'Slicing times must be a numeric vector');
            assert(ismember(refType, {'absolute', 'relative'}), ...
                'refType must be ''absolute'' or ''relative'' but instead was %s', refType);
            
            % Restore epoch order
            [~, indBack] = sort(this.epochInd);
            this.SortEpochs();
            
            % Loop through tables
            tSlice = tSlice(:);
            for k = indTable
                % Get reference times
                tRef = this.tot.referenceTime{k};
                assert(~isempty(tRef), 'To reslice, a table must have associated reference times');
                
                % Convert slicing times
                if strcmp(refType, 'absolute')
                    if numel(tSlice) ~= numel(this.numEpochs)
                        assert(~any(this.isEventValuesTable), ...
                            'The number of epoch cannot be changed due to the presence of eventValues table(s)');
                    end
                    tDelim = tSlice;
                else
                    if isscalar(tSlice)
                        tSlice = repmat(tSlice, [this.numEpochs 1]);
                    end
                    tDelim = tSlice + tRef;
                end
                assert(all(diff(tDelim) > 0), 'Slicing times must be monotonically increasing');
                
                % Slice table
                tb = this.tot.tableData{k};
                vect = cell(1,width(tb));
                if this.isEventTimesTable(k)
                    % For eventTimes table
                    for i = 1 : width(tb)
                        vect{i} = this.ICatColumn(tb.(i), tRef);
                    end
                    this.tot.tableData{k} = MSessionExplorer.MakeEventTimesTable(vect, ...
                        'DelimiterTimes', tDelim, ...
                        'VariableNames', tb.Properties.VariableNames, ...
                        'Verbose', this.isVerbose);
                    
                elseif this.isTimesSeriesTable(k)
                    % For timeSeries table
                    vect{1} = this.ICatColumn(tb.time, tRef);
                    for i = 2 : width(tb)
                        vect{i} = cat(1, tb.(i){:});
                    end
                    this.tot.tableData{k} = MSessionExplorer.MakeTimeSeriesTable(vect(1), vect(2:end), ...
                        'DelimiterTimes', tDelim, ...
                        'VariableNames', tb.Properties.VariableNames, ...
                        'Verbose', this.isVerbose);
                end
                this.tot.referenceTime{k} = tDelim;
            end
            
            if strcmp(refType, 'absolute')
                % Remove any eventValues table
                if any(this.isEventValuesTable)
                    warning('All eventValues table will be removed when using ''%s'' times', refType);
                    this.tot(this.isEventValuesTable,:) = [];
                end
                % Reset epochInd
                this.epochInd = (1 : this.numEpochs)';
            else
                % Sort epoch order back
                this.AlignTime(-tSlice);
                this.SortEpochs(indBack);
            end
        end
        
        function tbOut = SliceTimeSeries(this, tbIn, tWin, varargin)
            % Return slices of time series data using time windows
            % 
            %   tbOut = SliceTimeSeries(tbIn, tWin)
            %   tbOut = SliceTimeSeries(tbIn, tWin, rowInd)
            %   tbOut = SliceTimeSeries(tbIn, tWin, rowInd, colInd)
            %   tbOut = SliceTimeSeries(..., 'Fill', 'none')
            %   tbOut = SliceTimeSeries(..., 'ReferenceTime', [])
            % 
            % Inputs
            %   tbIn                A table of time series data or the name of a timeSeries table in the object.
            %   tWin                1) An n-by-2 matrix where each row has the begin and end time of a window. 
            %                       When n equals 1, this window is applied to every rows. When n equals the height 
            %                       of the table or the number of selected rows, each window is applied to respective 
            %                       row. Inf values will be substuted by min or max time available, respectively. 
            %                       2) An n-element cell array where each element is a m-by-2 matrix of time windows. 
            %                       All m windows will be applied to a single corresponding row, resulting in m rows 
            %                       in tbOut. Options for cell array length, n, are the same as described above. 
            %   rowInd              Integer or logical indices of rows to operate on and return. The default value is 
            %                       empty indicating all rows. 
            %   colInd              Integer or logical indices of columns to operate on and return. It can also be 
            %                       a cell array of column names of the input table. The default is empty indicating 
            %                       all columns. 
            %   'Fill'              When tWin exceeds data in an epoch, 'none' (default) does not fill anything in 
            %                       exceeded parts; 'bleed' will look for data from neighboring epochs up to the 
            %                       entire session. 
            %   'ReferenceTime'     Reference time to use with the 'bleed' option. The default value is empty and 
            %                       the method will look for reference time associated with the specified table. 
            % Output
            %   tbOut               The output table with inquired time series data.
            
            % Handle user inputs
            p = inputParser();
            p.addRequired('tbIn', @(x) ischar(x) || istable(x));
            p.addRequired('tWin', @(x) isnumeric(x) || iscell(x));
            p.addOptional('rowInd', [], @(x) isnumeric(x) || islogical(x));
            p.addOptional('colInd', [], @(x) isnumeric(x) || islogical(x) || iscellstr(x));
            p.addParameter('Fill', 'none', @(x) ismember(x, {'none', 'bleed', 'nan', 'nanbleed'}));
            p.addParameter('ReferenceTime', [], @(x) isnumeric(x) && isvector(x));
            
            p.parse(tbIn, tWin, varargin{:});
            rowInd = p.Results.rowInd(:);
            colInd = p.Results.colInd(:);
            fillMethod = p.Results.Fill;
            tRef = p.Results.ReferenceTime;
            
            if ~istable(tbIn)
                this.IValidateTableName(tbIn, true);
                assert(ismember(tbIn, this.tableNames(this.isTimesSeriesTable)), ...
                    '''%s'' is not a timeSeries table', tbIn);
                tRef = this.GetReferenceTime(tbIn);
                tbIn = this.GetTable(tbIn);
            end
            
            % Validate and standardize row and column indices
            [rowInd, colInd] = this.IValidateTableIndexing(rowInd, colInd, tbIn);
            if ~ismember(1, colInd)
                colInd = [1 colInd]; % make sure the time column is always included
            end
            tbIn = tbIn(:,colInd);
            
            % Validate and standardize time windows
            tWin = this.IValidateTimeWindows(tWin, height(tbIn), rowInd);
            
            % Replace Inf with ends of timestamp
            tBound = cellfun(@(x) x([1 end])', tbIn.time, 'Uni', false);
            for i = 1 : numel(tWin)
                bb = repmat(tBound{i}, size(tWin{i},1), 1);
                isWinInf = isinf(tWin{i});
                tWin{i}(isWinInf) = bb(isWinInf);
            end
            
            % Slicing
            switch fillMethod
                case 'none'
                    winCells = this.ISliceTsFillNone(tbIn, tWin, rowInd);
                case 'nan'
                    warning('''Fill'', ''nan'' may not be supported in a future version');
                    winCells = this.ISliceTsFillNaN(tbIn, tWin, rowInd);
                case 'bleed'
                    assert(~isempty(tRef), 'Requires ''ReferenceTime'' to use ''%s''', fillMethod);
                    winCells = this.ISliceTsFillBleed(tbIn, tWin, rowInd, tRef);
                case 'nanbleed'
                    warning('''Fill'', ''nanbleed'' may not be supported in a future version');
                    assert(~isempty(tRef), 'Requires ''ReferenceTime'' to use ''%s''', fillMethod);
                    tPan = cellfun(@(x) [min(x(:)) max(x(:))], tWin, 'Uni', false);
                    winCells = this.ISliceTsFillBleed(tbIn, tPan, rowInd, tRef);
                    winCells = cat(1, winCells{:});
                    tbIn = cell2table(winCells, 'VariableNames', tbIn.Properties.VariableNames);
                    winCells = this.ISliceTsFillNaN(tbIn, tWin(rowInd), 1:numel(rowInd));
            end
            winCells = cat(1, winCells{:});
            
            % Make table
            tbOut = table('Size', size(winCells), 'VariableTypes', repmat({'cell'}, [1 width(tbIn)]), ...
                'VariableNames', tbIn.Properties.VariableNames);
            tbOut{:,:} = winCells;
        end
        
        function tbOut = SliceEventTimes(this, tbIn, tWin, varargin)
            % Return slices of event time data using time windows
            % 
            %   tbOut = SliceEventTimes(tbIn, tWin)
            %   tbOut = SliceEventTimes(tbIn, tWin, rowInd)
            %   tbOut = SliceEventTimes(tbIn, tWin, rowInd, colInd)
            %   tbOut = SliceEventTimes(..., 'Fill', 'none')
            %   tbOut = SliceEventTimes(..., 'ReferenceTime', [])
            % 
            % Inputs
            %   tbIn                A table of event times data or the name of a eventTimes table in the object. 
            %   tWin                1) An n-by-2 matrix where each row has the begin and end time of a window. 
            %                       When n equals 1, this window is applied to every rows. When n equals the height 
            %                       of the table or the number of selected rows, each window is applied to respective 
            %                       row. 
            %                       2) An n-element cell array where each element is a m-by-2 matrix of time windows. 
            %                       n must equal to the height of the table or the number of selected rows. All m 
            %                       windows in a m-by-2 matrix are applied to a corresponding row. 
            %   rowInd              Integer or logical indices of rows to operate on and return. The default value is 
            %                       empty indicating all rows. 
            %   colInd              Integer or logical indices of columns to operate on and return. It can also be 
            %                       a cell array of column names of the input table. The default is empty indicating 
            %                       all columns. 
            %   'Fill'              When tWin exceeds an epoch, 'none' (default) ignores exceeded parts whereas 'bleed' 
            %                       looks for events from neighboring epochs up to the entire session. 
            %   'ReferenceTime'     Reference time to use with the 'bleed' option. The default value is empty and 
            %                       the method will look for reference time associated with the specified table. 
            % Output
            %   tbOut               The output table with inquired event times data.
            
            % Handle user inputs
            p = inputParser();
            p.addRequired('tbIn', @(x) ischar(x) || istable(x));
            p.addRequired('tWin', @(x) isnumeric(x) || iscell(x));
            p.addOptional('rowInd', [], @(x) isnumeric(x) || islogical(x));
            p.addOptional('colInd', [], @(x) isnumeric(x) || islogical(x) || iscellstr(x));
            p.addParameter('Fill', 'none', @(x) any(strcmpi(x, {'none', 'bleed'})));
            p.addParameter('ReferenceTime', [], @(x) isnumeric(x) && isvector(x));
            
            p.parse(tbIn, tWin, varargin{:});
            rowInd = p.Results.rowInd;
            colInd = p.Results.colInd;
            fillMethod = lower(p.Results.Fill);
            tRef = p.Results.ReferenceTime;
            
            if ~istable(tbIn)
                this.IValidateTableName(tbIn, true);
                assert(ismember(tbIn, this.tableNames(this.isEventTimesTable)), ...
                    '''%s'' is not a eventTimes table', tbIn);
                tRef = this.GetReferenceTime(tbIn);
                tbIn = this.GetTable(tbIn);
            end
            
            % Validate and standardize row and column indices
            [rowInd, colInd] = this.IValidateTableIndexing(rowInd, colInd, tbIn);
            tbIn = tbIn(:,colInd);
            
            % Validate and standardize time windows
            tWin = this.IValidateTimeWindows(tWin, height(tbIn), rowInd);
            
            % Slicing
            switch fillMethod
                case 'none'
                    winCells = this.ISliceEtFillNone(tbIn, tWin, rowInd);
                case 'bleed'
                    assert(~isempty(tRef), 'Requires ''ReferenceTime'' to use ''%s''', fillMethod);
                    winCells = this.ISliceEtFillBleed(tbIn, tWin, rowInd, tRef);
            end
            winCells = cat(1, winCells{:});
            
            % Make table
            for k = width(tbIn) : -1 : 1
                dtypes{k} = class(tbIn{1,k});
            end
            tbOut = table('Size', size(winCells), 'VariableTypes', dtypes, ...
                'VariableNames', tbIn.Properties.VariableNames);
            for k = 1 : width(tbOut)
                try
                    tbOut{:,k} = cell2mat(winCells(:,k));
                catch
                    tbOut.(k) = winCells(:,k);
                end
            end
        end
        
        function tbOut = ResampleTimeSeries(this, tbIn, tEdges, varargin)
            % Resample timeSeries table by interpolation
            % 
            %   tbOut = ResampleTimeSeries(tbIn, tEdges)
            %   tbOut = ResampleTimeSeries(tbIn, tEdges, rowInd)
            %   tbOut = ResampleTimeSeries(tbIn, tEdges, rowInd, colInd)
            % 
            % Inputs
            %   tbIn            A table of event times data or a name of a eventTimes table in the current object.
            %   tEdges          Edges of time bins. This can be one numeric vector that defines edges for every 
            %                   rows, or a cell array of vectors where each applies to a corresponding row. The 
            %                   number of element in cell array must equal the height of the table or the number 
            %                   of selected rows. 
            %   rowInd          Integer or logical indices of rows to operate on and return. The default value is 
            %                   empty indicating all rows. 
            %   colInd          Integer or logical indices of columns to operate on and return. It can also be 
            %                   a cell array of column names of the input table. The default is empty indicating 
            %                   all columns. 
            % Output
            %   tbOut           The output table of time series data where each value is the number of occurance. 
            
            warning('This method is in beta version and should not be used for formal analysis');
            
            % Handle user inputs
            p = inputParser();
            p.addRequired('tbIn', @(x) ischar(x) || istable(x));
            p.addRequired('tEdges', @(x) iscell(x) || isnumeric(x));
            p.addOptional('rowInd', [], @(x) isnumeric(x) || islogical(x));
            p.addOptional('colInd', [], @(x) isnumeric(x) || islogical(x) || iscellstr(x));
            
            p.parse(tbIn, tEdges, varargin{:});
            rowInd = p.Results.rowInd;
            colInd = p.Results.colInd;
            
            if ~istable(tbIn)
                this.IValidateTableName(tbIn, true);
                assert(ismember(tbIn, this.tableNames(this.isTimesSeriesTable)), ...
                    '''%s'' is not a timeSeries table', tbIn);
                tbIn = this.GetTable(tbIn);
            end
            
            % Validate and standardize row and column indices
            [rowInd, colInd] = this.IValidateTableIndexing(rowInd, colInd, tbIn);
            if ~ismember(1, colInd)
                colInd = [1 colInd]; % make sure the time column is always included
            end
            
            % Validate and standardize time edges
            tEdges = this.IValidateTimeEdges(tEdges, height(tbIn), rowInd);
            
            % Interpolate time series
            tbOut = tbIn(rowInd, colInd);
            tEdges = tEdges(rowInd);
            for i = 1 : height(tbOut)
                t = tbOut.time{i};
                tq = tEdges{i}(1:end-1) + diff(tEdges{i})/2;
                tbOut.time{i} = tq;
                tbOut{i,2:end} = cellfun(@(v) interp1(t, v, tq), tbOut{i,2:end}, 'Uni', false);
            end
        end
        
        function tbOut = ResampleEventTimes(this, tbIn, tEdges, varargin)
            % Resample an eventTimes table to a timeSeries table by binning events
            % 
            %   tbOut = ResampleEventTimes(tbIn, tEdges)
            %   tbOut = ResampleEventTimes(tbIn, tEdges, rowInd)
            %   tbOut = ResampleEventTimes(tbIn, tEdges, rowInd, colInd)
            %   
            % Inputs
            %   tbIn            A table of event times data or a name of a eventTimes table in the current object.
            %   tEdges          Edges of time bins. This can be one numeric vector that defines edges for every 
            %                   rows, or a cell array of vectors where each applies to a corresponding row. The 
            %                   number of element in cell array must equal the height of the table or the number 
            %                   of selected rows. 
            %   rowInd          Integer or logical indices of rows to operate on and return. The default value is 
            %                   empty indicating all rows. 
            %   colInd          Integer or logical indices of columns to operate on and return. It can also be 
            %                   a cell array of column names of the input table. The default is empty indicating 
            %                   all columns. 
            % Output
            %   tbOut           The output table of time series data where each value is the number of occurance. 
            
            warning('This method is in beta version and should not be used for formal analysis');
            
            % Handle user inputs
            p = inputParser();
            p.addRequired('tbIn', @(x) ischar(x) || istable(x));
            p.addRequired('tEdges', @(x) iscell(x) || isnumeric(x));
            p.addOptional('rowInd', [], @(x) isnumeric(x) || islogical(x));
            p.addOptional('colInd', [], @(x) isnumeric(x) || islogical(x) || iscellstr(x));
            
            p.parse(varargin{:});
            rowInd = p.Results.rowInd;
            colInd = p.Results.colInd;
            
            if ~istable(tbIn)
                this.IValidateTableName(tbIn, true);
                assert(ismember(tbIn, this.tableNames(this.isEventTimesTable)), ...
                    '''%s'' is not a eventTimes table', tbIn);
                tbIn = this.GetTable(tbIn);
            end
            
            assert(~ismember('time', tbIn.Properties.VariableNames), ...
                'The input table cannot contain variable named ''time'' since the output table requires it');
            
            % Validate and standardize row and column indices
            [rowInd, colInd] = this.IValidateTableIndexing(rowInd, colInd, tbIn);
            tbIn = tbIn(:,colInd);
            
            % Validate and standardize time edges
            tEdges = this.IValidateTimeEdges(tEdges, height(tbIn), rowInd);
            
            % Bin event times
            tbOut = tbIn(rowInd, colInd);
            tEdges = tEdges(rowInd);
            for k = 1 : width(tbOut)
                et = tbOut{:,k};
                if isnumeric(et)
                    et = num2cell(et);
                end
                tbOut.(k) = cellfun(@(x,y) histcounts(x,y)', et, tEdges, 'Uni', false);
            end
            
            % Add time column
            tbOut.time = cellfun(@(x) x(1:end-1)+diff(x)/2, tEdges, 'Uni', false);
            tbOut = tbOut(:,[end, 1:end-1]);
        end
    end
    
    % These methods are used internally
    methods(Access = protected)
        function val = IValidateTableName(this, tbName, isAssert)
            % The table name must be a string
            val = ischar(tbName);
            assert(isAssert && val, 'A table name must be a character array rather than %s.', class(tbName));
            if ~val
                return;
            end
            % The table must exist
            val = ismember(tbName, this.tableNames);
            assert(isAssert && val, 'A table named ''%s'' does not exist', tbName);
        end
        
        function val = IValidateTableNames(this, tbNames, isAssert)
            if iscell(tbNames)
                % A cell array of table names
                val = cellfun(@(x) this.IValidateTableName(x, isAssert), tbNames);
            else
                % One single table name
                val = this.IValidateTableName(tbNames, isAssert);
            end
        end
        
        function [rowInd, colInd] = IValidateTableIndexing(~, rowInd, colInd, tb)
            % Validate and convert row and column indices to integers 
            
            if isempty(rowInd)
                % include all rows
                rowInd = 1 : height(tb);
            elseif islogical(rowInd)
                % convert to numerical indices
                rowInd = find(rowInd);
            end
            assert(isnumeric(rowInd), 'Cannot interpret row indexing');
            
            if isempty(colInd)
                % include all columns
                colInd = 1 : width(tb);
            elseif islogical(colInd)
                % convert to numerical indices
                colInd = find(colInd);
            elseif iscellstr(colInd) 
                % find column ind based on variable names
                varNames = tb.Properties.VariableNames;
                colNames = colInd;
                for i = 1 : numel(colInd)
                    colInd{i} = find(strcmp(colNames{i}, varNames));
                    assert(~isempty(colInd{i}), '%s is not a valid column name in the table.', colNames{i});
                end
                colInd = cell2mat(colInd);
            end
            assert(isnumeric(colInd), 'Cannot interpret column indexing');
            
            % Reshape as row vectors
            rowInd = rowInd(:)';
            colInd = colInd(:)';
        end
        
        function winOut = IValidateTimeWindows(~, winIn, tbHeight, rowInd)
            % Validate and convert time windows to a cell array of window matrices
            
            if isnumeric(winIn)
                winIn = num2cell(winIn, 2);
            end
            winIn = winIn(:);
            
            winOut = num2cell(NaN(tbHeight,2), 2);
            if numel(winIn) == 1
                % Propagate value to all epochs
                winOut(rowInd) = repmat(winIn, [numel(rowInd) 1]);
            elseif numel(winIn) == numel(rowInd)
                % Add windows to selected rows
                winOut(rowInd) = winIn;
            elseif numel(winIn) == tbHeight
                % Full size windows
                winOut = winIn;
            else
                error('Incorrect number of time windows');
            end
            
            for i = 1 : numel(winOut)
                % Verify array size
                assert(~isempty(winOut{i}), 'Time window cannot be empty. Consider using [NaN NaN] instead');
                assert(size(winOut{i},2) == 2, 'Time window array must have 2 elements in each row');
                % Verify time increment (ignore NaN windows)
                dt = diff(winOut{i}, 1, 2);
                assert(all(isnan(dt) | dt > 0), 'The window end time must be greater than the begin time');
            end
        end
        
        function edgeOut = IValidateTimeEdges(~, edgeIn, tbHeight, rowInd)
            % Verify time edges and turn it into a standard format
            
            edgeOut = cell(tbHeight, 1);
            if isnumeric(edgeIn)
                % Propagate edges to all rows
                edgeOut(rowInd) = repmat({edgeIn}, [numel(rowInd) 1]);
            elseif numel(edgeIn) == numel(rowInd)
                % Convert edges array for full table
                edgeOut(rowInd) = edgeIn;
            elseif numel(edgeOut) == tbHeight
                % Use full edges
                edgeOut = edgeIn;
            else
                error('Incorrect size of time edges array');
            end
            edgeOut = cellfun(@(x) x(:), edgeOut, 'Uni', false);
            
            for i = rowInd
                x = edgeOut{i};
                assert(isnumeric(x) && isvector(x), 'Each time edges must be a numeric vector');
                assert(all(diff(x) > 0), 'Time edges must be increasing numbers');
            end
        end
        
        function vect = ICatColumn(~, col, tRef)
            if nargin < 3
                tRef = zeros(size(col));
            end
            if iscell(col)
                % Cell vector of numeric array
                colCells = cellfun(@(x,r) x+r, col, num2cell(tRef), 'Uni', false);
                vect = cat(1, colCells{:});
            else
                % Numeric vector
                vect = col + tRef;
            end
        end
        
        function winData = ISliceTsFillNone(~, tbIn, tWin, rowInd)
            % Find time series for each epoch
            winData = cell(size(tWin));
            for i = rowInd
                t = tbIn.time{i};
                w = tWin{i};
                winData{i} = cell(size(w,1), width(tbIn));
                
                % Find time series for each window
                for j = 1 : size(w,1)
                    isInWin = w(j,1) <= t & t < w(j,2);
                    winData{i}(j,:) = cellfun(@(x) x(isInWin,:), tbIn{i,:}, 'Uni', false);
                end
            end
        end
        
        function winData = ISliceTsFillNaN(~, tbIn, tWin, rowInd)
            % Find time series for each epoch
            winData = cell(size(tWin));
            for i = rowInd
                t = tbIn.time{i};
                dtPre = diff(t(1:2));
                dtPost = diff(t(end-1:end));
                w = tWin{i};
                winData{i} = cell(size(w,1), width(tbIn));
                
                % Find time series for each window
                for j = 1 : size(w,1)
                    isInWin = w(j,1) <= t & t < w(j,2);
                    
                    % Make time stamps in exceeded parts
                    tPre = flip(t(1)-dtPre : -dtPre : w(j,1))';
                    tPre(tPre >= w(j,2)) = [];
                    
                    tPost = (t(end)+dtPost : dtPost : w(j,2))';
                    tPost(tPost < w(j,1)) = [];
                    
                    % Indexing and filling
                    winData{i}{j,1} = [tPre; t(isInWin,:); tPost];
                    winData{i}(2:end) = cellfun( ...
                        @(x) [NaN(numel(tPre), size(x,2)); x(isInWin,:); NaN(numel(tPost), size(x,2))], ...
                        tbIn{i,2:end}, 'Uni', false);
                end
            end
        end
        
        function winData = ISliceTsFillBleed(~, tbIn, tWin, rowInd, tRef)
            % Cache variables
            tAbs = cellfun(@(x,r) x + r, tbIn.time, num2cell(tRef), 'Uni', false);
            tAbsBegin = cellfun(@(x) x(1), tAbs);
            tAbsEnd = cellfun(@(x) x(end), tAbs);
            
            % Find time series for each epoch
            winData = cell(size(tWin));
            for i = rowInd
                wAbs = tWin{i} + tRef(i);
                winData{i} = cell(size(wAbs,1), width(tbIn));
                
                % Find time series for each window
                for j = 1 : size(wAbs,1)
                    % Find relevant epochs
                    isWinBeforeEpoch = all(tAbsBegin - wAbs(j,:) >= 0, 2);
                    isWinAfterEpoch = all(wAbs(j,:) - tAbsEnd > 0, 2);
                    isWinOverlapEpoch = ~(isWinBeforeEpoch | isWinAfterEpoch);
                    srcRowInd = find(isWinOverlapEpoch);
                    
                    % Sort relevant epochs in temporal order
                    [~, ord] = sort(tRef(srcRowInd));
                    srcRowInd = srcRowInd(ord);
                    
                    % Get masks
                    isInWin = cellfun(@(x) wAbs(j,1) <= x & x < wAbs(j,2), tAbs(srcRowInd), 'Uni', false);
                    
                    % Find times
                    tCells = cellfun(@(x,y) x(y)-tRef(i), tAbs(srcRowInd), isInWin, 'Uni', false);
                    winData{i}{j,1} = vertcat(tCells{:});
                    
                    % Find data
                    for k = 2 : width(tbIn)
                        dataCells = cellfun(@(x,y) x(y,:), tbIn{srcRowInd,k}, isInWin, 'Uni', false);
                        winData{i}{j,k} = vertcat(dataCells{:});
                    end
                end
            end
        end
        
        function winData = ISliceEtFillNone(~, tbIn, tWin, rowInd)
            % Find event times for each epoch
            winData = cell(size(tWin));
            for i = rowInd
                w = tWin{i};
                winData{i} = cell(size(w,1), width(tbIn));
                
                % Find event times for each variable
                for j = 1 : width(tbIn)
                    % Get event times
                    t = tbIn{i,j};
                    if iscell(t)
                        t = t{1};
                    end
                    
                    % Find event times in each window
                    for k = 1 : size(w,1)
                        tHit = t(t >= w(k,1) & t < w(k,2));
                        if isempty(tHit)
                            tHit = NaN;
                        end
                        winData{i}{k,j} = tHit;
                    end
                end
            end
        end
        
        function winData = ISliceEtFillBleed(~, tbIn, tWin, rowInd, tRef)
            % Cache variables
            tAbs = cell(size(tbIn));
            for i = 1 : width(tbIn)
                if isnumeric(tbIn{:,i})
                    tAbs(:,i) = num2cell(double(tbIn{:,i}) + tRef);
                elseif iscell(tbIn{:,i})
                    tAbs(:,i) = cellfun(@(x,r) double(x)+r, tbIn{:,i}, num2cell(tRef), 'Uni', false);
                else
                    error('Cannot convert relative time with this data format');
                end
            end
            for i = 1 : numel(tAbs)
                if isempty(tAbs{i})
                    tAbs{i} = NaN;
                end
            end
            tAbsBegin = cellfun(@(x) x(1), tAbs);
            tAbsEnd = cellfun(@(x) x(end), tAbs);
            
            % Find event times for each epoch
            winData = cell(size(tWin));
            for i = rowInd
                wAbs = tWin{i} + tRef(i);
                winData{i} = cell(size(wAbs,1), width(tbIn));
                
                % Find event times for each variable
                for j = 1 : width(tbIn)
                    % Find event times in each window
                    for k = 1 : size(wAbs,1)
                        % Collect source epochs
                        isWinBeforeEpoch = all(tAbsBegin(:,j) - wAbs(k,:) >= 0, 2);
                        isWinAfterEpoch = all(wAbs(k,:) - tAbsEnd(:,j) > 0, 2);
                        isWinOverlapEpoch = ~(isWinBeforeEpoch | isWinAfterEpoch);
                        srcRowInd = find(isWinOverlapEpoch);
                        
                        % Sort source epochs in temporal order
                        [~, ord] = sort(tRef(srcRowInd));
                        srcRowInd = srcRowInd(ord);
                        
                        % Get masks
                        isInWin = cellfun(@(x) wAbs(k,1) <= x & x < wAbs(k,2), tAbs(srcRowInd,j), 'Uni', false);
                        
                        % Find times
                        tCells = cellfun(@(x,y) x(y)-tRef(i), tAbs(srcRowInd,j), isInWin, 'Uni', false);
                        t = vertcat(tCells{:});
                        if isempty(t)
                            t = NaN;
                        end
                        winData{i}{k,j} = t;
                    end
                end
            end
        end
    end
    
    methods(Hidden)
        % These methods are for backward compatibility
        function SortTrials(this, ind)
            warning('SortTrials method will be removed in a future version. Use SortEpochs instead.');
            this.SortEpochs(ind);
        end
        function RemoveTrials(this, ind2rm)
            warning('RemoveTrials method will be removed in a future version. Use RemoveEpochs instead.');
            this.RemoveEpochs(ind2rm);
        end
        
        % Hide methods from handle superclass
        function listener(~)
        end
        function addlistener(~)
        end
        function notify(~)
        end
        function findobj(~)
        end
        function findprop(~)
        end
    end
end




