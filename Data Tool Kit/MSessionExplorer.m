classdef MSessionExplorer < handle
    %MSESSIONEXPLORER Summary of this class goes here
    %   Detailed explanation goes here
    
    % Constant properties that can never be changed
    properties(Constant)
        supportedTableTypes = {'eventTimes', 'eventValues', 'timeSeries'};
    end
    
    % User have total control of these properties
    properties
        % User data
        userData;
    end
    
    % These properties are read-only by user but can be changed by methods
    % it guarentees that data across different containers are always consistent
    properties(GetAccess = public, SetAccess = private)
        % Table of data tables
        tot;
        originalTrialInd;
    end
    
    properties(Dependent)
        tableNames;
        isEventTimesTable;
        isEventValuesTable;
        isTimesSeriesTable;
        numTrials;
    end
    
    methods
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
        function val = get.numTrials(this)
            if ~isempty(this.tot)
                val = height(this.tot.tableData{1});
            else
                val = 0;
            end
        end
    end
    
    % These methods operate on data in this object
    methods
        function this = MSessionExplorer()
            % Constructor of MSessionExplorer
            
            % Initialize table of tables
            totHeaders = {'tableType', 'tableData', 'referenceTime'};
            this.tot = cell2table(cell(0, length(totHeaders)), 'VariableNames', totHeaders);
        end
        
        function Preview(this)
            % Print a summary of the content of this object
            %   
            %   Preview()
            %
            
            disp(this);
            disp('tot');
            disp(this.tot);
            disp('userData');
            disp(this.userData);
        end
        
        function se = Duplicate(this, varargin)
            % Duplicate the current MSessionExplorer object
            % 
            %   se = Duplicate()
            %   se = Duplicate(tableInclude)
            %   se = Duplicate(tableInclude, isUserData)
            % 
            % Inputs:
            %   tableInclude        A string or cell array of table name(s) specifying which table(s) to include.
            %                       An empty array indicates to include all table. 
            %   isUserData          A logical variable indicating whether or not to copy userData.
            % Output:
            %   se                  The new MSessionExplorer object
            % 
            
            % Handle user input
            p = inputParser();
            p.addOptional('tableInclude', []);
            p.addOptional('isUserData', true, @islogical);
            
            p.parse(varargin{:});
            tbInc = p.Results.tableInclude;
            isUserData = p.Results.isUserData;
            
            if isempty(tbInc)
                tbInc = this.tableNames;
            end
            tbInc = cellstr(tbInc);
            
            if ~all(cellfun(@(x) any(strcmp(x, this.tableNames)), tbInc))
                error('Not all table in tableInclude can be found.');
            end
            
            indTableInc = cellfun(@(x) find(strcmp(x, this.tableNames)), tbInc);
            indTableInc = indTableInc(:)';
            
            
            % Copying
            se = MSessionExplorer();
            
            for i = indTableInc
                se.SetTable(this.tableNames{i}, this.tot.tableData{i}, this.tot.tableType{i}, this.tot.referenceTime{i});
            end
            
            if isUserData
                se.userData = this.userData;
            end
        end
        
        function se = Merge(this, varargin)
            % Merge multiple MSessionExplorer objects into one
            % 
            %   se = Merge(se1, se2, se3, ...)
            % 
            % Inputs:
            %   se1, se2, se3, ...      Arbitrary number of MSessionExplorer objects
            % Output:
            %   se                      The merged MSessionExplorer object
            % 
            
            % Check user input
            seArray = [this; cat(1, varargin{:})];
            
            % Merging
            se = MSessionExplorer();
            
            for i = 1 : numel(this.tableNames)
                tbData = arrayfun(@(x) x.tot.tableData{i}, seArray, 'Uni', false);
                tbData = cat(1, tbData{:});
                refTime = arrayfun(@(x) x.tot.referenceTime{i}, seArray, 'Uni', false);
                refTime = cat(1, refTime{:});
                se.SetTable(this.tableNames{i}, tbData, this.tot.tableType{i}, refTime);
            end
            
            se.userData = this.userData;
        end
        
        function s = ToStruct(this)
            % Convert the current MSessionExplorer object to a structure
            % 
            %   s = ToStruct(this)
            % 
            % Output:
            %   s       The output structure with the following fields
            %           tableName       a cell array of table names
            %           tableType       a cell array of table types
            %           referenceTime   a cell array of referece time
            %           tableData       a cell array of structs whose fields are each table's
            %                           variables (i.e. output of table2struct function)
            %           userData        same as userData property
            %   
            
            s.tableName = this.tableNames;
            s.tableType = this.tot.tableType;
            s.referenceTime = this.tot.referenceTime;
            s.tableData = cellfun(@(x) table2struct(x), this.tot.tableData, 'Uni', false);
            s.userData = this.userData;
        end
        
        function SetTable(this, varargin)
            % Add or update a data table
            % 
            %   SetTable(tableName, dataTable)
            %   SetTable(tableName, dataTable, tableType)
            %   SetTable(tableName, dataTable, tableType, referenceTime)
            %
            % Inputs:
            %   tableName           A string of the name of table to add or update. 
            %   dataTable           The table variable.
            %   tableType           Either 'eventTimes', 'eventValues', or 'timeSeries' indicating 
            %                       the type of data table. It is required when adding a new table 
            %                       but is ignored when updating an existing table.  
            %   referenceTime       A numeric vector where each element stores the absolute time
            %                       of the zero relative time of each trial. Default is empty. 
            % 
            
            % Handle user input
            p = inputParser();
            p.addRequired('tableName', @ischar);
            p.addRequired('dataTable', @istable);
            p.addOptional('tableType', '', @(x) any(strcmp(x, this.supportedTableTypes)));
            p.addOptional('referenceTime', [], @isnumeric);
            
            p.parse(varargin{:});
            tbName = p.Results.tableName;
            tb = p.Results.dataTable;
            tbType = p.Results.tableType;
            refTimes = p.Results.referenceTime;
            
            % Check trial number conflict
            if ~isempty(this.tableNames)
                rowNums = cellfun(@(x) size(x,1), this.tot.tableData);
                if ~all(size(tb,1) == rowNums)
                    error('The new table cannot be added because it has different number of rows than existing tables.');
                end
            end
            
            % Set table
            if any(strcmp(tbName, this.tableNames))
                % Replace an existing table
                this.tot{tbName, 'tableData'}{1} = tb;
                if ~isempty(tbType)
                    this.tot{tbName, 'tableType'}{1} = tbType;
                end
                if ~isempty(refTimes)
                    this.tot{tbName, 'referenceTime'}{1} = refTimes;
                end
            else
                % Add a new table
                if isempty(tbType)
                    error('Table type is required for adding a new table.');
                end
                totHeaders = this.tot.Properties.VariableNames;
                totRow = cell2table({tbType, tb, refTimes}, 'VariableNames', totHeaders, 'RowNames', {tbName});
                this.tot = [this.tot; totRow];
            end
            
            if isempty(this.originalTrialInd)
                this.originalTrialInd = (1 : this.numTrials)';
            end
        end
        
        function SetReferenceTime(this, rt, tbName)
            % Set reference times to table(s)
            % 
            %   SetReferenceTime(rt)
            %   SetReferenceTime(rt, tbName)
            %
            % Inputs:
            %   rt              A numeric vector where each element stores the absolute time
            %                   of the zero relative time of each trial. 
            %   tbName          A string or cell array of the name(s) of table which reference time
            %                   is set to. The default is empty and rt is set to all applicable
            %                   tables. 
            % 
            
            if isempty(this.tot)
                warning('There is no table to set.');
                return;
            end
            
            if nargin < 3
                tbName = this.tableNames(~this.isEventValuesTable);
            end
            
            tbName = cellstr(tbName);
            
            for i = 1 : length(tbName)
                if height(this.tot{tbName{i},'tableData'}{1}) ~= length(rt)
                    error('Number of rows in table does not match the number of reference time.');
                end
                
                if ~this.isEventValuesTable(strcmp(tbName{i}, this.tableNames))
                    this.tot{tbName{i}, 'referenceTime'}{1} = rt;
                else
                    warning([tbName{i} ' is an eventValue table and reference time is not applicable.']);
                end
            end
        end
        
        function tb = GetTable(this, tbName)
            % Get table data from the table of tables
            % 
            %   tb = GetTable(tbName)
            %
            % Input:
            %   tbName          The name of table to return. 
            % Output:
            %   tb              The table data. 
            % 
            
            if isempty(this.tot)
                warning('There is no table to get.');
                return;
            end
            
            tb = this.tot{tbName, 'tableData'}{1};
        end
        
        function rt = GetReferenceTime(this, tbName)
            % Get a table from the table of tables
            % 
            %   rt = GetReferenceTime(tbName)
            %
            % Input:
            %   tbName          The name of a table which the reference time is associated with.
            % Output:
            %   rt              A vector of reference times. 
            % 
            
            if isempty(this.tot)
                warning('There is no table to get.');
                return;
            end
            
            rt = this.tot{tbName, 'referenceTime'}{1};
        end
        
        function RemoveTable(this, tbName)
            % Remove specified data table
            % 
            %   RemoveTable(tbName)
            %
            % Input:
            %   tbName          A string or cell array of strings indicating the name of table(s) to remove.
            % 
            
            if isempty(this.tot)
                warning('There is no table to remove.');
                return;
            end
            
            this.tot(tbName,:) = [];
        end
        
        function RemoveTrials(this, ind2rm)
            % Remove trials by indices in all tables
            % 
            %   RemoveTrials(ind2rm)
            %
            % Input:
            %   ind2rm          Integer or logical indices of trials to remove. 
            % 
            
            % Checking
            if isempty(this.tot)
                warning('There is no table to operate on.');
                return;
            end
            if (islogical(ind2rm) && all(ind2rm)) || length(unique(ind2rm)) == this.numTrials
                error('Removing all trials is not supported.');
            end
            
            % Remove trials
            for k = 1 : height(this.tot)
                this.tot.tableData{k}(ind2rm,:) = [];
                if ~isempty(this.tot.referenceTime{k})
                    this.tot.referenceTime{k}(ind2rm) = [];
                end
            end
            this.originalTrialInd(ind2rm) = [];
        end
        
        function SortTrials(this, ind)
            % Sort trials across all tables
            % 
            %   SortTrials(ind)
            %
            % Input:
            %   ind             Indices of the new order. 
            % 
            
            % Handle user inputs
            if length(unique(ind)) ~= this.numTrials
                error('The length of indices must equal the number of trials and each index must be unique.');
            end
            
            ind = ind(:);
            
            % Process original indices
            this.originalTrialInd = this.originalTrialInd(ind);
            
            % Sort trials
            for k = 1 : height(this.tot)
                this.tot.tableData{k} = this.tot.tableData{k}(ind,:);
                if ~isempty(this.tot.referenceTime{k})
                    this.tot.referenceTime{k} = this.tot.referenceTime{k}(ind);
                end
            end
        end
        
        function AlignTime(this, refEvent, refSourceTableName, tScale)
            % Shift time stamps in eventTimes and timeSeries tables wrt the reference event or time values. 
            % Optionally also scale time linearly after alignment. 
            % 
            %   AlignTime(refEvent)
            %   AlignTime(refEvent, refSourceTableName)
            %   AlignTime(refEvent, refSourceTableName, tScale)
            %
            % Input:
            %   refEvent                An event name in a eventTimes table or a numeric vector of reference times 
            %                           as zero time after alignment. Each element of the vector matches respective 
            %                           row in tables. 
            %   refSourceTableNames     If refEvent is a string, then you must specify the name of eventTimes table 
            %                           where this event names should be found. This avoids ambiguity of the same 
            %                           variable name found in multiple tables. The default is empty. 
            %   tScale                  A factor that linearly scales time after alingment (default 1, no scaling)
            % 
            
            % Handle user inputs
            if nargin < 4
                tScale = 1;
            end
            
            indAlignable = find(this.isEventTimesTable | this.isTimesSeriesTable);
            if isempty(indAlignable)
                error('You need to have at least one eventTimes or timeSeries table to operate on.');
            end
            
            if ischar(refEvent)
                % Find column names from specified tables
                if nargin < 3
                    error('You need to provide refSourceTableName indicating where refEvent is in.');
                end
                if ~ischar(refSourceTableName)
                    error('refSourceTableName must be a string of table name.');
                end
                refEvent = this.tot{refSourceTableName, 'tableData'}{1}.(refEvent);
            end
            if ~isnumeric(refEvent) || ~isvector(refEvent)
                error('Reference times must be a numeric vector.');
            end
            if length(refEvent) ~= size(this.tot.tableData{indAlignable(1)}, 1)
                error('The length of reference times does not match number of rows in data table.');
            end
            refEvent = refEvent(:);
            
            % Align times
            for k = indAlignable
                
                if this.isEventTimesTable(k)
                    % For eventTimes table
                    tb = this.tot.tableData{k};
                    
                    for i = 1 : size(tb, 2)
                        if isnumeric(tb{:,i})
                            % Numeric vector
                            tb{:,i} = (tb{:,i} - refEvent) * tScale;
                        elseif iscell(tb{:,i})
                            % Cell vector of numeric vectors
                            tb{:,i} = cellfun(@(x,r) (x-r)*tScale, tb{:,i}, num2cell(refEvent), 'Uni', false);
                        else
                            warning(fprintf('Failed to align data in %s column of %s table\n', ...
                                tb.Properties.VariableNames{i}, this.tableName{k}));
                        end
                    end
                    
                    this.tot.tableData{k} = tb;
                    
                elseif this.isTimesSeriesTable(k)
                    % For timeSeries table
                    this.tot.tableData{k}.time = ...
                        cellfun(@(x,r) (x-r)*tScale, this.tot.tableData{k}.time, num2cell(refEvent), 'Uni', false);
                end
                
                % Change referenceTime
                if ~isempty(this.tot.referenceTime{k})
                    this.tot.referenceTime{k} = this.tot.referenceTime{k} + refEvent;
                end
                
            end
        end
        
        function tbOut = SliceTimeSeries(this, varargin)
            % Return inquired time series
            % Note that the output table may contain incomplete and/or duplicating information and should not be used 
            % to replace the original data table stored in the object unless you are sure.
            % 
            %   SliceTimeSeries(table)
            %   SliceTimeSeries(table, timeWindow)
            %   SliceTimeSeries(table, timeWindow, rowInd)
            %   SliceTimeSeries(table, timeWindow, rowInd, colInd)
            %   SliceTimeSeries(..., 'Fill', 'none')
            %   SliceTimeSeries(..., 'ReferenceTime', [])
            %   SliceTimeSeries(..., 'Func', [])
            %   SliceTimeSeries(..., 'FuncScope', 'cell')
            % 
            % Inputs:
            %   table               A table of time series data or a name of a timeSeries table in the current object.
            %   timeWindow          An n-by-2 matrix of row vectors indicating the begin and end times. 
            %                       When n equals 1, this window is applied to every rows. When n equals the height of
            %                       the table or the number of selected rows, each windows is applied to respective row.
            %                       Inf values will be substuted by min or max time of each trial, respectively. 
            %                       The default of this argument is empty in which case no slicing is done but one can 
            %                       still apply user function. 
            %   rowInd              Integer or logical indices of rows to operate on and return. The default value is 
            %                       empty indicating all rows. 
            %   colInd              Integer or logical indices of columns to operate on and return. The default value is 
            %                       empty indicating all columns. 
            %   'Fill'              When timeWindow exceeds data in a trial, 'none' (default) does not fill anything in 
            %                       exceeded parts; 'nan' fills with NaNs (with the same sampling frequency); 'bleed' 
            %                       will look for data from neighboring trials up to the entire session but fill nothing 
            %                       where samples don't exist; 'nanbleed' is similar to 'bleed' but fills NaNs to time 
            %                       without samples. 
            %   'ReferenceTime'     Reference time to use with the 'bleed' option. The default value is empty and 
            %                       the method will look for reference time associated with the specified table. 
            %   'Func'              Function handle for manipulating (e.g. filtering) each time series. The returned 
            %                       vector must have the same size as the input vector. This is applied before 
            %                       resampling. 
            %   'FuncScope'         If set to 'cell' (default), 'func' operates in each cell in the table; if set to 
            %                       'column', 'func' operates on each column of the table as a continuous vector. 
            % Output:
            %   tbOut               The output table with inquired time series data.
            %
            
            % Handle user inputs
            p = inputParser();
            p.addRequired('table', @(x) ischar(x) || istable(x));
            p.addOptional('timeWindow', [], @isnumeric);
            p.addOptional('rowInd', [], @(x) isnumeric(x) || islogical(x));
            p.addOptional('colInd', [], @(x) isnumeric(x) || islogical(x) || iscell(x));
            p.addParameter('Fill', 'none', @(x) any(strcmpi(x, {'none', 'nan', 'bleed', 'nanbleed'})));
            p.addParameter('ReferenceTime', [], @(x) isnumeric(x) && isvector(x));
            p.addParameter('Func', [], @(x) isempty(x) || isa(x, 'function_handle'));
            p.addParameter('FuncScope', 'cell', @(x) any(strcmpi(x, {'cell', 'column'})));
            
            p.parse(varargin{:});
            tbIn = p.Results.table;
            tWin = p.Results.timeWindow;
            rowInd = p.Results.rowInd(:);
            colInd = p.Results.colInd(:);
            fillMethod = lower(p.Results.Fill);
            tRef = p.Results.ReferenceTime;
            funcHandle = p.Results.Func;
            funcScope = lower(p.Results.FuncScope);
            
            if ischar(tbIn)
                if ~ismember(tbIn, this.tableNames(this.isTimesSeriesTable))
                    error('Cannot find timeSeries table ''%s''.', tbIn);
                end
                if isempty(tRef)
                    tRef = this.tot{tbIn,'referenceTime'}{1};
                end
                tbIn = this.GetTable(tbIn);
            end
            
            [rowInd, colInd] = this.UnifyTableIndexing(rowInd, colInd, tbIn, true);
            tbIn = tbIn(:,colInd);
            
            % Apply time window
            if ~isempty(tWin)
                tWin = this.UnifyTimeWindow(tWin, height(tbIn), rowInd);
                
                % Substitute infinity to trial boundary
                tBound = cellfun(@(x) x([1 end])', tbIn.time, 'Uni', false);
                tBound = cell2mat(tBound);
                isWinInf = isinf(tWin);
                tWin(isWinInf) = tBound(isWinInf);
                
                % Find time series for each trial
                switch fillMethod
                    case 'none'
                        tbOut = this.SliceTsFillNone(tbIn, tWin, rowInd);
                    case 'nan'
                        tbOut = this.SliceTsFillNaN(tbIn, tWin, rowInd);
                    case 'bleed'
                        tbOut = this.SliceTsFillBleed(tbIn, tWin, rowInd, tRef);
                    case 'nanbleed'
                        tbOut = this.SliceTsFillBleed(tbIn, tWin, rowInd, tRef);
                        tbOut = this.SliceTsFillNaN(tbOut, tWin(rowInd,:), 1:numel(rowInd));
                end
            else
                tbOut = tbIn(rowInd,:);
            end
            
            % Apply user function
            if ~isempty(funcHandle)
                switch funcScope
                    case 'cell'
                        % Process each cell individually
                        tbOut{:,2:end} = cellfun(funcHandle, tbOut{:,2:end}, 'Uni', false);
                    case 'column'
                        % Process each column as a whole
                        cellLengths = cellfun(@length, tbOut{:,1});
                        for i = 2 : width(tbOut)
                            colVect = funcHandle(cell2mat(tbOut{:,i}));
                            tbOut{:,i} = mat2cell(colVect, cellLengths);
                        end
                end
            end
        end
        
        function tbOut = SliceEventTimes(this, varargin)
            % Retrieve events with given time window(s). 
            % (The output table may contain incomplete and/or duplicating information and should not be used 
            % to replace the original event table stored in the class.)
            % 
            %   SliceEventTimes(table, timeWindow)
            %   SliceEventTimes(table, timeWindow, rowInd)
            %   SliceEventTimes(table, timeWindow, rowInd, colInd)
            %   SliceEventTimes(..., 'Fill', 'none')
            %   SliceEventTimes(..., 'ReferenceTime', [])
            % 
            % Inputs:
            %   table               A table of event time data or a name of eventTimes table in this object.
            %   timeWindow          An n-by-2 matrix of row vectors indicating the begin and end times. When n equals 
            %                       1, this window is applied to every rows. When n equals the height of the table, 
            %                       n windows are applied to corresponding rows, respectively. 
            %   rowInd              Integer or logical indices of rows to operate on and return. The default value 
            %                       is empty indicating all rows. 
            %   colInd              Integer or logical indices of columns to operate on and return. The default value 
            %                       is empty indicating all columns. 
            %   'Fill'              When timeWindow exceeds a trial, 'none' ignores exceeded parts whereas 'bleed' 
            %                       (default) looks for events from neighboring trials up to the entire session. 
            %                       Only use 'bleed' when trial reference times are valid. 
            %   'ReferenceTime'     Reference time to use with the 'bleed' option. The default value is empty and 
            %                       the method will look for reference time associated with the specified table. 
            % Output:
            %   tbOut               The output table with inquired time series data.
            %
            
            % Handle user inputs
            p = inputParser();
            p.addRequired('table', @(x) ischar(x) || istable(x));
            p.addRequired('timeWindow', @isnumeric);
            p.addOptional('rowInd', [], @(x) isnumeric(x) || islogical(x));
            p.addOptional('colInd', [], @(x) isnumeric(x) || islogical(x));
            p.addParameter('Fill', 'none', @(x) any(strcmpi(x, {'none', 'bleed'})));
            p.addParameter('ReferenceTime', [], @(x) isnumeric(x) && isvector(x));
            
            p.parse(varargin{:});
            tbIn = p.Results.table;
            tWin = p.Results.timeWindow;
            rowInd = p.Results.rowInd;
            colInd = p.Results.colInd;
            fillMethod = lower(p.Results.Fill);
            tRef = p.Results.ReferenceTime;
            
            if ischar(tbIn)
                if ~ismember(tbIn, this.tableNames(this.isEventTimesTable))
                    error('Cannot find eventTimes table ''%s''.', tbIn);
                end
                if isempty(tRef)
                    tRef = this.tot{tbIn,'referenceTime'}{1};
                end
                tbIn = this.GetTable(tbIn);
            end
            if isempty(rowInd)
                rowInd = 1 : height(tbIn);
            end
            if islogical(rowInd)
                rowInd = find(rowInd);
            end
            if isempty(colInd)
                colInd = 1 : width(tbIn);
            end
            if islogical(colInd)
                colInd = find(colInd);
            end
            
            tWin = this.UnifyTimeWindow(tWin, height(tbIn), rowInd);
            
            switch fillMethod
                case 'none'
                    tbCells = this.SliceEtFillNone(tbIn, tWin, rowInd, colInd);
                case 'bleed'
                    tbCells = this.SliceEtFillBleed(tbIn, tWin, rowInd, colInd, tRef);
            end
            
            % Try to convert cell array to numeric array
            tbOut = tbIn(rowInd, colInd);
            for k = 1 : width(tbOut)
                try
                    tbOut{:,k} = cell2mat(tbCells(:,k));
                catch
                    tbOut.(k) = tbCells(:,k);
                end
            end
        end
        
        function tbOut = ResampleTimeSeries(this, varargin)
            
            % Handle user inputs
            p = inputParser();
            p.addRequired('table', @(x) ischar(x) || istable(x));
            p.addRequired('timeData', @(x) iscell(x) || isnumeric(x));
            p.addRequired('timeOpt', @isscalar);
            p.addOptional('rowInd', [], @(x) isnumeric(x) || islogical(x));
            p.addOptional('colInd', [], @(x) isnumeric(x) || islogical(x) || iscell(x));
            p.addParameter('Func', [], @(x) isa(x, 'function_handle'));
            p.addParameter('FuncScope', 'cell', @(x) any(strcmpi(x, {'cell', 'column'})));
            p.parse(varargin{:});
            tbIn = p.Results.table;
            tData = p.Results.timeData;
            tOpt = p.Results.timeOpt;
            rowInd = p.Results.rowInd;
            colInd = p.Results.colInd;
            funcHandle = p.Results.Func;
            funcScope = p.Results.FuncScope;
            
            if ischar(tbIn)
                if ~ismember(tbIn, this.tableNames(this.isTimesSeriesTable))
                    error('Cannot find timeSeries table ''%s''.', tbIn);
                end
                tbIn = this.GetTable(tbIn);
            end
            
            % Get table indices
            [rowInd, colInd] = this.UnifyTableIndexing(rowInd, colInd, tbIn, true);
            
            % Get time edges
            tEdges = this.UnifyTimeEdges(tData, tOpt, height(tbIn), rowInd);
            if tOpt
                tEdges = cellfun(@(x) [x(:); x(end)*2-x(end-1)], tEdges, 'Uni', false);
            end
            
            % Bin event times
            tbOut = this.RsTs(tbIn(rowInd, colInd), tEdges(rowInd));
            
            % Apply user function
            tbOut = this.SliceTimeSeries(tbOut, 'Func', funcHandle, 'FuncScope', funcScope);
        end
        
        function tbOut = ResampleEventTimes(this, varargin)
            % Resample an eventTimes table to a timeSeries table by binning events. 
            % 
            %   ResampleEventTimes(table, timeWindow, binSize)
            %   ResampleEventTimes(table, timeEdges, [])
            %   ResampleEventTimes(..., rowInd)
            %   ResampleEventTimes(..., rowInd, colInd)
            %   ResampleEventTimes(..., 'Func', [])
            %   ResampleEventTimes(..., 'FuncScope', 'cell')
            %   
            % Inputs:
            %   table               A table of event times data or a name of a eventTimes table in the current object.
            %   timeWindow          An n-by-2 matrix of row vectors indicating the begin and end times. 
            %                       When n equals 1, this window is applied to every rows. When n equals the height of
            %                       the table or the number of selected rows, each window is used on respective row.
            %   binSize             Size of time bin. 
            %   timeEdges           Edges of time bins. This can be one numeric vector that defines edges for every 
            %                       rows, or a cell array of vectors where each applies to a corresponding row. The 
            %                       number of element in cell array must equal the height of the table or the number 
            %                       of selected rows. 
            %   rowInd              Integer or logical indices of rows to operate on and return. The default value is 
            %                       empty indicating all rows. 
            %   colInd              Integer or logical indices of columns to operate on and return. It can also be 
            %                       a cell array of column names of the input table. The default is empty indicating 
            %                       all columns. 
            %   'Func'              Function handle for manipulating (e.g. filtering) each time series. The returned 
            %                       vector must have the same size as the input vector. This is applied before 
            %                       resampling. 
            %   'FuncScope'         If set to 'cell' (default), 'func' operates in each cell in the table; if set to 
            %                       'column', 'func' operates on each column of the table as a continuous vector. 
            % Output:
            %   tbOut               The output table of time series data where each value is the number of occurance. 
            %
            
            % Handle user inputs
            p = inputParser();
            p.addRequired('table', @(x) ischar(x) || istable(x));
            p.addRequired('timeData', @(x) iscell(x) || isnumeric(x));
            p.addRequired('timeOpt', @isscalar);
            p.addOptional('rowInd', [], @(x) isnumeric(x) || islogical(x));
            p.addOptional('colInd', [], @(x) isnumeric(x) || islogical(x) || iscell(x));
            p.addParameter('Func', [], @(x) isa(x, 'function_handle'));
            p.addParameter('FuncScope', 'cell', @(x) any(strcmpi(x, {'cell', 'column'})));
            
            p.parse(varargin{:});
            tbIn = p.Results.table;
            tData = p.Results.timeData;
            tOpt = p.Results.timeOpt;
            rowInd = p.Results.rowInd;
            colInd = p.Results.colInd;
            funcHandle = p.Results.Func;
            funcScope = p.Results.FuncScope;
            
            if ischar(tbIn)
                if ~ismember(tbIn, this.tableNames(this.isEventTimesTable))
                    error('Cannot find eventTimes table ''%s''.', tbIn);
                end
                tbIn = this.GetTable(tbIn);
            end
            
            % Get table indices
            [rowInd, colInd] = this.UnifyTableIndexing(rowInd, colInd, tbIn);
            
            % Get time edges
            tEdges = this.UnifyTimeEdges(tData, tOpt, height(tbIn), rowInd);
            if tOpt
                tEdges = cellfun(@(x) [x; x(end)*2-x(end-1)], tEdges, 'Uni', false);
            end
            
            % Bin event times
            tbOut = this.Et2Ts(tbIn(rowInd, colInd), tEdges(rowInd));
            
            % Apply user function
            tbOut = this.SliceTimeSeries(tbOut, 'Func', funcHandle, 'FuncScope', funcScope);
        end
    end
    
    % These methods are only used internally
    methods(Access = private)
        function tWin = UnifyTimeWindow(~, tWin, tbHeight, rowInd)
            % Verify time window(s) and turn it into a standard format
            if nargin < 4
                rowInd = 1 : tbHeight;
            end
            % Verify timeWindow array size
            if size(tWin,1) ~= tbHeight && size(tWin,1) ~= numel(rowInd) && size(tWin,2) ~= 2
                error('The array of timeWindow has incorrect size.');
            end
            % Verify time increment (ignore NaN windows)
            if ~all(tWin(:,2)-tWin(:,1) > 0 | any(isnan(tWin),2))
                error('The window end time must be greater than the begin time.');
            end
            % Propagate timeWindow value for all trials
            if numel(tWin) == 2
                tWin = repmat(tWin(:)', numel(rowInd), 1);
            end
            % Make timeWindow full size
            if size(tWin,1) < tbHeight
                tWinFull = NaN(tbHeight, 2);
                tWinFull(rowInd,:) = tWin;
                tWin = tWinFull;
            end
        end
        
        function tEdges = UnifyTimeEdges(this, tData, tOpt, tbHeight, rowInd)
            % Verify time window(s) and turn it into a standard format
            if nargin < 4
                rowInd = 1 : tbHeight;
            end
            % Standardize
            if isnumeric(tData) && ~iscolumn(tData)
                % Convert time window and sample time to time edges
                tWin = this.UnifyTimeWindow(tData, tbHeight, rowInd);
                tEdges = cellfun(@(x) x(1):tOpt:x(2), num2cell(tWin,2), 'Uni', false);
            else
                % Use time edges
                tEdges = tData;
                if isnumeric(tEdges) && isvector(tEdges)
                    % convert numeric vector to cell
                    tEdges = {tEdges};
                elseif ~iscell(tEdges)
                    error('Time edges must be a numeric vector or a cell array of numeric vector(s).');
                end
                if numel(tEdges) == 1
                    % propagate edges to all rows
                    tEdges = repmat(tEdges, [tbHeight, 1]);
                elseif numel(tEdges) == numel(rowInd)
                    % convert edges array for full table
                    tEdgesFull = cell(tbHeight, 1);
                    tEdgesFull(rowInd) = tEdges;
                    tEdges = tEdgesFull;
                elseif numel(tEdges) ~= tbHeight
                    error('The number of time edges vectors cannot be matched with the number of rows.');
                end
            end
            tEdges = cellfun(@(x) x(:), tEdges, 'Uni', false);
        end
        
        function [rowInd, colInd] = UnifyTableIndexing(~, rowInd, colInd, tb, isTs)
            
            if isempty(rowInd)
                % include all rows
                rowInd = 1 : height(tb);
            elseif islogical(rowInd)
                % convert to numerical indices
                rowInd = find(rowInd);
            end
            if ~isnumeric(rowInd)
                error('Cannot interpret row indices.');
            end
            
            if isempty(colInd)
                % include all columns
                colInd = 1 : width(tb);
            elseif islogical(colInd)
                % convert to numerical indices
                colInd = find(colInd);
            elseif iscell(colInd)
                % find column ind based on variable names
                varNames = tb.Properties.VariableNames;
                colNames = colInd;
                for i = 1 : numel(colInd)
                    colInd{i} = find(strcmp(colNames{i}, varNames));
                    if isempty(colInd{i})
                        error('%s is not a valid column name in the table.', colNames{i});
                    end
                end
                colInd = cell2mat(colInd);
            end
            if ~isnumeric(colInd)
                error('Cannot interpret column indices.');
            end
            
            rowInd = rowInd(:)';
            colInd = colInd(:)';
            
            if nargin < 5
                isTs = false;
            end
            if isTs && ~ismember(1, colInd)
                colInd = [1 colInd];
            end
        end
        
        function tbOut = SliceTsFillNone(~, tbIn, tWin, rowInd)
            
            tbOut = tbIn;
            
            % Find time series for each trial
            for i = rowInd(:)'
                % Find mask
                isInWin = tWin(i,1) <= tbIn.time{i} & tbIn.time{i} <= tWin(i,2);
                
                % Direct indexing
                tbOut{i,:} = cellfun(@(x) x(isInWin,:), tbIn{i,:}, 'Uni', false);
            end
            
            tbOut = tbOut(rowInd,:);
        end
        
        function tbOut = SliceTsFillNaN(~, tbIn, tWin, rowInd)
            
            tbOut = tbIn;
            
            % Find time series for each trial
            for i = rowInd(:)'
                % Find mask
                t = tbIn.time{i};
                isInWin = tWin(i,1) <= t & t <= tWin(i,2);
                
                % Make time stamps in exceeded parts
                dtPre = diff(t(1:2));
                tPre = flip(t(1)-dtPre : -dtPre : tWin(i,1))';
                tPre(tPre > tWin(i,2)) = [];
                
                dtPost = diff(t(end-1:end));
                tPost = (t(end)+dtPost : dtPost : tWin(i,2))';
                tPost(tPost < tWin(i,1)) = [];
                
                % Indexing and filling
                tbOut.time{i} = [tPre; t(isInWin,:); tPost];
                tbOut{i,2:end} = cellfun( ...
                    @(x) [NaN(length(tPre), size(x,2)); x(isInWin,:); NaN(length(tPost), size(x,2))], ...
                    tbIn{i,2:end}, 'Uni', false);
            end
            
            tbOut = tbOut(rowInd,:);
        end
        
        function tbOut = SliceTsFillBleed(~, tbIn, tWin, rowInd, tRef)
            
            tbOut = tbIn;
            
            % Cache variables
            tWinAbs = tWin + tRef;
            tAbs = cellfun(@(x,r) x + r, tbIn.time, num2cell(tRef), 'Uni', false);
            tAbsBegin = cellfun(@(x) x(1), tAbs);
            tAbsEnd = cellfun(@(x) x(end), tAbs);
            
            % Find time series for each trial
            for i = rowInd(:)'
                % Collect source trials
                isWinBeforeTrial = all(tWinAbs(i,:) - tAbsBegin < 0, 2);
                isWinAfterTrial = all(tWinAbs(i,:) - tAbsEnd > 0, 2);
                isWinOverlapTrial = ~(isWinBeforeTrial | isWinAfterTrial);
                srcRowInd = find(isWinOverlapTrial);
                
                % Sort source trials in temporal order
                [~, ord] = sort(tRef(srcRowInd));
                srcRowInd = srcRowInd(ord);
                
                % Find masks
                isInWin = cellfun(@(x) tWinAbs(i,1) <= x & x <= tWinAbs(i,2), tAbs(srcRowInd), 'Uni', false);
                
                % Find times
                tbOut.time{i} = cell2mat(cellfun(@(x,y) x(y)-tRef(i), tAbs(srcRowInd), isInWin, 'Uni', false));
                
                % Find data
                for j = 2 : width(tbIn)
                    dataCells = cellfun(@(x,y) x(y,:), tbIn{srcRowInd,j}, isInWin, 'Uni', false);
                    tbOut{i,j}{1} = vertcat(dataCells{:});
                end
            end
            
            tbOut = tbOut(rowInd,:);
        end
        
        function tbCells = SliceEtFillNone(~, tbIn, tWin, rowInd, colInd)
            
            tbCells = cell(size(tbIn));
            for j = colInd(:)'
                for i = rowInd(:)'
                    % Get event times
                    t = tbIn{i,j};
                    if iscell(t)
                        t = t{1};
                    end
                    
                    % Find events in window
                    tHit = t(t >= tWin(i,1) & t <= tWin(i,2));
                    if isempty(tHit)
                        tHit = NaN;
                    end
                    tbCells{i,j} = tHit;
                end
            end
            tbCells = tbCells(rowInd, colInd);
        end
        
        function tbCells = SliceEtFillBleed(~, tbIn, tWin, rowInd, colInd, tRef)
            
            % Convert relative time to absolute time
            tWinAbs = tWin + tRef;
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
            
            % Find events in each window from all events
            tbCells = cell(size(tbIn));
            for j = colInd(:)'
                for i = rowInd(:)'
                    % Collect source trials
                    isWinBeforeTrial = all(tWinAbs(i,:) - tAbsBegin(:,j) < 0, 2);
                    isWinAfterTrial = all(tWinAbs(i,:) - tAbsEnd(:,j) > 0, 2);
                    isWinOverlapTrial = ~(isWinBeforeTrial | isWinAfterTrial);
                    srcRowInd = find(isWinOverlapTrial);
                    
                    % Sort source trials in temporal order
                    [~, ord] = sort(tRef(srcRowInd));
                    srcRowInd = srcRowInd(ord);
                    
                    % Find masks
                    isInWin = cellfun(@(x) tWinAbs(i,1) <= x & x <= tWinAbs(i,2), tAbs(srcRowInd,j), 'Uni', false);
                    
                    % Find times
                    tCells = cellfun(@(x,y) x(y)-tRef(i), tAbs(srcRowInd,j), isInWin, 'Uni', false);
                    t = vertcat(tCells{:});
                    if isempty(t)
                        t = NaN;
                    end
                    tbCells{i,j} = t;
                end
            end
            tbCells = tbCells(rowInd, colInd);
        end
        
        function tb = Et2Ts(~, tb, tEdges)
            % Bin event times
            for k = 1 : width(tb)
                et = tb{:,k};
                if isnumeric(et)
                    et = num2cell(et);
                end
                tb.(k) = cellfun(@(x,y) histcounts(x,y)', et, tEdges, 'Uni', false);
                
                if strcmp(tb.Properties.VariableNames{k}, 'time')
                    warning(['The column named ''time'' is changed to ''time_event'' in the output timeSeries table' ...
                        ' which requires a dedicated ''time'' column for timestamps. ']);
                    tb.Properties.VariableNames{k} = 'time_event';
                end
            end
            
            % Add time column
            tb.time = cellfun(@(x) x(1:end-1)+diff(x)/2, tEdges, 'Uni', false);
            tb = tb(:,[end, 1:end-1]);
        end
        
        function tb = RsTs(~, tb, tEdges)
            for i = 1 : height(tb)
                t = tb.time{i};
                tNew = tEdges{i}(1:end-1) + diff(tEdges{i})/2;
                tb.time{i} = tNew;
                tb{i,2:end} = cellfun(@(x) interp1(t, x, tNew), tb{i,2:end}, 'Uni', false);
            end
        end
        
        function tb = ApplyFunc(~, tb, func, funcScope)
            switch funcScope
                case 'cell'
                    % Process each cell individually
                    tb{:,2:end} = cellfun(func, tb{:,2:end}, 'Uni', false);
                case 'column'
                    % Process each column as a whole
                    cellLengths = cellfun(@length, tb{:,1});
                    for i = 2 : width(tb)
                        colVect = func(cell2mat(tb{:,i}));
                        tb{:,i} = mat2cell(colVect, cellLengths);
                    end
            end
        end
    end
    
    % Hide methods from handle class
    methods(Hidden)
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







