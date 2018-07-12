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
        
        function Show(this)
            % Print a summary of content of this object
            
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
            %   se = Duplicate(..., 'tableInclude', all_tables)
            %   se = Duplicate(..., 'userData', true)
            % 
            % Input:
            %   'tableInclude'          A vector of table name(s) specifying which table(s) to include
            %   'userData'              Whether or not to copy userData
            % Output:
            %   se                      The new MSessionExplorer object
            % 
            
            % Handle user input
            p = inputParser();
            p.addParameter('tableInclude', this.tableNames);
            p.addParameter('userData', true, @islogical);
            
            p.parse(varargin{:});
            tbInc = cellstr(p.Results.tableInclude);
            isUserData = p.Results.userData;
            
            if ~all(cellfun(@(x) any(strcmp(x, this.tableNames)), tbInc))
                error('Not all table in tableInclude can be found.');
            end
            
            indTableInc = cellfun(@(x) find(strcmp(x, this.tableNames)), tbInc);
            indTableInc = indTableInc(:)';
            
            
            % Copying
            se = MSessionExplorer();
            
            for i = indTableInc
                se.AddTable(this.tableNames{i}, this.tot.tableType{i}, this.tot.tableData{i}, this.tot.referenceTime{i});
            end
            
            if isUserData
                se.userData = this.userData;
            end
        end
        
        function s = Save(this, matPath, varargin)
            % 
            
            s.tableName = this.tableNames;
            s.tableType = this.tot.tableType;
            s.referenceTime = this.tot.referenceTime;
            s.tableData = cellfun(@(x) table2struct(x), this.tot.tableData, 'Uni', false);
            s.userData = this.userData;
            
            save(matPath, '-struct', 's');
        end
        
        function SetTable(this, varargin)
            % Add or update a data table
            % 
            %   SetTable(tableName, dataTable)
            %   SetTable(tableName, dataTable, tableType)
            %   SetTable(tableName, dataTable, tableType, referenceTimes)
            %
            % Inputs:
            %   tableName           
            %   dataTable           
            %   tableType           
            %   referenceTimes      
            % 
            
            % Handle user input
            p = inputParser();
            p.addRequired('tableName', @ischar);
            p.addRequired('dataTable', @istable);
            p.addOptional('tableType', '', @(x) any(strcmp(x, this.supportedTableTypes)));
            p.addOptional('referenceTimes', [], @isnumeric);
            
            p.parse(varargin{:});
            tbName = p.Results.tableName;
            tb = p.Results.dataTable;
            tbType = p.Results.tableType;
            refTimes = p.Results.referenceTimes;
            
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
            % Get a table from the table of tables
            
            if isempty(this.tot)
                warning('There is no table to get.');
                return;
            end
            
            tb = this.tot{tbName, 'tableData'}{1};
        end
        
        function rt = GetReferenceTime(this, tbName)
            % Get a table from the table of tables
            
            if isempty(this.tot)
                warning('There is no table to get.');
                return;
            end
            
            rt = this.tot{tbName, 'referenceTime'}{1};
        end
        
        function RemoveTable(this, tbName)
            % Remove specified data table
            
            if isempty(this.tot)
                warning('There is no table to remove.');
                return;
            end
            
            this.tot(tbName,:) = [];
        end
        
        function AlignTime(this, refEvent, refSourceTableName)
            % Shift time stamps in eventTimes and timeSeries tables wrt the provided reference event name or time values
            % 
            %   AlignTime(refEvent)
            %   AlignTime(refEvent, refSourceTableName)
            %
            % Input:
            %   refEvent                The event name in a eventTimes table or a numeric vector of reference times
            %                           as zero time after alignment. Each element of the vector matches respective 
            %                           row in tables. 
            %   refSourceTableNames     If refEvent is a string, then you must specify the name of eventTimes table 
            %                           where this event names should be found. This avoids ambiguity of the same 
            %                           variable name found in multiple tables. 
            % 
            
            % Handle user input
            indAlignableTables = find(this.isEventTimesTable | this.isTimesSeriesTable);
            
            if isempty(indAlignableTables)
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
            
            if length(refEvent) ~= size(this.tot.tableData{indAlignableTables(1)}, 1)
                error('The length of reference times does not match number of rows in data table.');
            end
            
            refEvent = refEvent(:);
            
            
            % Align times
            for k = indAlignableTables
                
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
                            warning(fprintf('Failed to align data in %s column of %s table\n', ...
                                tb.Properties.VariableNames{i}, this.tableName{k}));
                        end
                    end
                    
                    this.tot.tableData{k} = tb;
                    
                elseif this.isTimesSeriesTable(k)
                    % For timeSeries table
                    this.tot.tableData{k}.time = ...
                        cellfun(@(x,r) x-r, this.tot.tableData{k}.time, num2cell(refEvent), 'Uni', false);
                end
                
                % Change referenceTime
                if ~isempty(this.tot.referenceTime{k})
                    this.tot.referenceTime{k} = this.tot.referenceTime{k} + refEvent;
                end
                
                % Print progress
                fprintf('%s table aligned\n', this.tableNames{k});
            end
        end
        
        function SortTrials(this, ind)
            % Sort trials across all tables
            % 
            %   SortTrials(ind)
            %
            % Input:
            %   ind                 Indices of the new order. 
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
        
        function tbOut = GetTimeSeries(this, varargin)
            % Return inquired time series
            % (The output table may contain incomplete and/or duplicating information and should not be used 
            % to replace the original data table stored in the class.)
            % 
            %   GetTimeSeries(table)
            %   GetTimeSeries(table, timeWindow)
            %   GetTimeSeries(table, timeWindow, rowInd)
            %   GetTimeSeries(..., 'fill', 'none')
            %   GetTimeSeries(..., 'func', [])
            %   GetTimeSeries(..., 'funcScope', 'cell')
            %   GetTimeSeries(..., 'resampleFactor', 1)
            % 
            % Inputs:
            %   table               A table of time series data or a name of a timeSeries table in the current object.
            %   timeWindow          n-by-2 matrix of row vectors indicating the begin and end times.
            %                       When n equals 1, this window is applied to every rows. When n equals the height of
            %                       the table or the number of selected rows, each windows is applied to respective row.
            %   rowInd              Indices of rows to operate on and return. Default is all rows. 
            %   'fill'              When timeWindow exceeds data in a trial, 'none'(default) does not fill anything in 
            %                       exceeded parts; 'nan' fills with NaNs (with the same sampling frequency); 'bleed' 
            %                       will look for data from neighboring trials up to the entire session. To use 'bleed', 
            %                       referenceTimes must be available. Time without samples will be filled with NaN. 
            %   'func'              Function handle for manipulating (e.g. filtering) each time series. The returned 
            %                       vector must have the same size as the input vector. This is applied before any 
            %                       resampling operation. 
            %   'funcScope'         If set to 'cell'(default), 'func' operates on individual cells in the table; if
            %                       set to 'column', 'func' operates on each column of the table as a continuous vector. 
            %   'resampleFactor'    Resample each time series by interpolation (after applying 'func'). Value greater 
            %                       than 1 upsamples and less than 1 downsamples. 
            % Output:
            %   tbOut               The output table with inquired time series data.
            %
            
            % Handle user inputs
            p = inputParser();
            p.addRequired('table', @(x) ischar(x) || istable(x));
            p.addOptional('timeWindow', [], @isnumeric);
            p.addOptional('rowInd', [], @(x) isnumeric(x) || islogical(x));
            p.addParameter('fill', 'none', @(x) any(strcmpi(x, {'none', 'nan', 'bleed'})));
            p.addParameter('func', [], @(x) isa(x, 'function_handle'));
            p.addParameter('funcScope', 'cell', @(x) any(strcmpi(x, {'cell', 'column'})));
            p.addParameter('resampleFactor', 1, @isscalar);
            
            p.parse(varargin{:});
            tbIn = p.Results.table;
            tWin = p.Results.timeWindow;
            rowInd = p.Results.rowInd(:);
            fillMethod = p.Results.fill;
            funcHandle = p.Results.func;
            funcScope = p.Results.funcScope;
            resampleFactor = p.Results.resampleFactor;
            
            if ~istable(tbIn)
                if ~any(strcmp(tbIn, this.tableNames(this.isTimesSeriesTable)))
                    error('Cannot find timeSeries table ''%s''.', tbIn);
                end
                
                if strcmpi(fillMethod, 'bleed')
                    tRef = this.tot{tbIn,'referenceTime'}{1};
                end
                
                tbIn = this.GetTable(tbIn);
            end
            
            if isempty(rowInd)
                rowInd = (1 : height(tbIn))';
            end
            
            if islogical(rowInd)
                rowInd = find(rowInd);
            end
            
            tbOut = tbIn;
            
            
            % Apply time window
            if ~isempty(tWin)
                
                tWin = this.UnifyTimeWindow(tWin, tbIn, rowInd);
                
                % Substitute infinity to trial boundary
                tBound = cellfun(@(x) x([1 end])', tbIn.time, 'Uni', false);
                tBound = cell2mat(tBound);
                isWinInf = isinf(tWin);
                tWin(isWinInf) = tBound(isWinInf);
                
                % Cache variables for 'bleed'
                if strcmpi(fillMethod, 'bleed')
                    tWinAbs = tWin + tRef;
                    tAbs = cellfun(@(x,r) x + r, tbIn.time, num2cell(tRef), 'Uni', false);
                    tAbsBegin = cellfun(@(x) x(1), tAbs);
                    tAbsEnd = cellfun(@(x) x(end), tAbs);
                end
                
                % Find time series for each trial
                for i = rowInd'
                    if strcmpi(fillMethod, 'none')
                        % Find mask
                        isInWin = tWin(i,1) < tbIn.time{i} & tbIn.time{i} < tWin(i,2);
                        
                        % Direct indexing
                        tbOut{i,:} = cellfun(@(x) x(isInWin,:), tbIn{i,:}, 'Uni', false);
                        
                    elseif strcmpi(fillMethod, 'nan')
                        % Find mask
                        isInWin = tWin(i,1) < tbIn.time{i} & tbIn.time{i} < tWin(i,2);
                        
                        % Make time stamp vectors in exceeded parts
                        tData = tbIn.time{i};
                        tPerSample = diff(tData(1:2));
                        
                        tPreData = flip(tData(1)-tPerSample : -tPerSample : tWin(i,1))';
                        tPreData(tPreData > tWin(i,2)) = [];
                        
                        tPostData = (tData(end)+tPerSample : tPerSample : tWin(i,2))';
                        tPostData(tPostData < tWin(i,1)) = [];
                        
                        % Indexing and filling
                        tbOut.time{i} = [tPreData; tData(isInWin,:); tPostData];
                        tbOut{i,2:end} = cellfun( ...
                            @(x) [NaN(length(tPreData), size(x,2)); x(isInWin,:); NaN(length(tPostData), size(x,2))], ...
                            tbIn{i,2:end}, 'Uni', false);
                        
                    elseif strcmpi(fillMethod, 'bleed')
                        % Collect source trials (in temporal order)
                        srcRowInd = find(~(all(tAbsBegin - tWinAbs(i,:) > 0, 2) | all(tAbsEnd - tWinAbs(i,:) < 0, 2)));
                        [~, ord] = sort(tRef(srcRowInd));
                        srcRowInd = srcRowInd(ord);
                        
                        % Find mask(s)
                        isInWin = cellfun(@(x) tWinAbs(i,1) < x & x < tWinAbs(i,2), tAbs(srcRowInd), 'Uni', false);
                        
                        % Find times
                        tbOut.time{i} = cell2mat(cellfun(@(x,y) x(y)-tRef(i), tAbs(srcRowInd), isInWin, 'Uni', false));
                        
                        % Find data
                        for j = 2 : width(tbIn)
                            tbOut{i,j}{1} = cell2mat(cellfun(@(x,y) x(y,:), tbIn{srcRowInd,j}, isInWin, 'Uni', false));
                        end
                    end
                end
                
                tbOut = tbOut(rowInd,:);
            end
            
            
            % Apply user function
            if ~isempty(funcHandle)
                if strcmpi(funcScope, 'cell')
                    % Process each cell individually
                    tbOut{:,2:end} = cellfun(funcHandle, tbOut{:,2:end}, 'Uni', false);
                    
                elseif strcmpi(funcScope, 'column')
                    % Process each column as a whole
                    cellLengths = cellfun(@length, tbOut{:,1});
                    
                    for i = 2 : width(tbOut)
                        colVect = funcHandle(cell2mat(tbOut{:,i}));
                        tbOut{:,i} = mat2cell(colVect, cellLengths);
                    end
                end
            end
            
            
            % Resampling
            if resampleFactor ~= 1
                for i = 1 : height(tbOut)
                    
                    tPerSample = diff(tbOut.time{i}(1:2));
                    tPerSampleNew = tPerSample / resampleFactor;
                    
                    if resampleFactor > 1
                        tResample = (tbOut.time{i}(1) : tPerSampleNew : tbOut.time{i}(end))';
                        tbOut{i,:} = cellfun(@(x) interp1(tbOut.time{i}, x, tResample, 'linear'), tbOut{i,:}, 'Uni', false);
                    else
                        tbOut{i,2:end} = cellfun(@(x) decimate(double(x), 1/resampleFactor), tbOut{i,2:end}, 'Uni', false);
                        tbOut.time{i} = linspace(tbOut.time{i}(1)+tPerSampleNew/2, tbOut.time{i}(end), length(tbOut{i,2}{1}))';
                    end
                end
            end
            
        end
        
        function tbOut = GetEventTimes(this, varargin)
            % Retrieve events with given time window(s). 
            % (The output table may contain incomplete and/or duplicating information and should not be used 
            % to replace the original event table stored in the class.)
            % 
            %   GetEventTimes(table, timeWindow)
            %   GetEventTimes(table, timeWindow, rowInd)
            %   GetEventTimes(..., 'fill', 'bleed')
            % 
            % Inputs:
            %   table           A table of event time data or a name of eventTimes table in this object.
            %   timeWindow      n-by-2 matrix of row vectors indicating the begin and end times. When n equals 1, 
            %                   this window is applied to every rows. When n equals the height of the table, 
            %                   n windows are applied to corresponding rows, respectively. 
            %   'fill'          When timeWindow exceeds a trial, 'none' ignores exceeded parts whereas 'bleed' 
            %                   (default) looks for events from neighboring trials up to the entore session. 
            %                   Only use 'bleed' when trial reference times are valid. 
            % Output:
            %   tbOut           The output table with inquired time series data.
            %
            
            % Handle user inputs
            p = inputParser();
            p.addRequired('table', @(x) ischar(x) || istable(x));
            p.addRequired('timeWindow', @isnumeric);
            p.addOptional('rowInd', [], @(x) isnumeric(x) || islogical(x));
            p.addOptional('colInd', [], @(x) isnumeric(x) || islogical(x));
            p.addParameter('fill', 'none', @(x) any(strcmpi(x, {'none', 'bleed'})));
            
            p.parse(varargin{:});
            tbIn = p.Results.table;
            tWin = p.Results.timeWindow;
            rowInd = p.Results.rowInd;
            colInd = p.Results.colInd;
            fillMethod = p.Results.fill;
            
            if ~istable(tbIn)
                % Verify table name
                if ~any(strcmp(tbIn, this.tableNames(this.isEventTimesTable)))
                    error('Cannot find eventTimes table ''%s''.', tbIn);
                else
                    tRef = this.tot{tbIn,'referenceTime'}{1};
                    tbIn = this.GetTable(tbIn);
                end
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
            
            isBleed = strcmpi(fillMethod, 'bleed');
            
            tWin = this.UnifyTimeWindow(tWin, tbIn, rowInd);
            
            
            % Convert relative time to absolute time
            eventAbsTimes = cell(1, width(tbIn));
            for i = 1 : width(tbIn)
                if isnumeric(tbIn{:,i})
                    eventAbsTimes{i} = tbIn{:,i} + tRef;
                elseif iscell(tbIn{:,i})
                    eventAbsTimes{i} = cell2mat(cellfun(@(x,r) x(:) + r, tbIn{:,i}, num2cell(tRef), 'Uni', false));
                else
                    error('Cannot convert relative time with this data format');
                end
            end
            
            tWin = tWin + tRef;
            
            
            % Find events in each window from all events
            tbCells = cell(size(tbIn));
            for j = colInd
                % All events
                evts = eventAbsTimes{j};
                
                for i = rowInd
                    % Find events in time window from all events in j for trial i
                    evtsHit = evts(evts >= tWin(i,1) & evts <= tWin(i,2));
                    if isempty(evtsHit)
                        evtsHit = NaN;
                    end
                    
                    % Convert absolute time to relative time as in input table
                    evtsHit = evtsHit - tRef(i);
                    
                    % Sort times in case rows in table are not in order
                    evtsHit = sort(evtsHit);
                    
                    tbCells{i,j} = evtsHit;
                end
            end
            tbCells = tbCells(rowInd, colInd);
            tbOut = tbIn(rowInd, colInd);
            
            
            % Try to convert cell array to numeric array
            for k = 1 : width(tbOut)
                try
                    tbOut{:,k} = cell2mat(tbCells(:,k));
                catch
                    tbOut.(k) = tbCells(:,k);
                end
            end
            
        end
        
    end
    
    % These methods are only used internally
    methods(Access = private)
        function tWin = UnifyTimeWindow(~, tWin, tb, rowInd)
            % Verify time window(s) and turn it into a standard format
            %
            %   tWin            A n-by-2 matrix of row vectors indicating the start and end times.
            %                   When n equals 1, this window is applied to every row. When n equals the height of
            %                   the table or the number of selected rows, each window is applied to respective row.
            %   roiInd          Indices of rows to apply window(s). Default is all rows. 
            %
            
            [tbHeight, ~] = size(tb);
            
            if nargin < 4
                rowInd = 1 : tbHeight;
            end
            
            if islogical(rowInd)
                rowInd = find(rowInd);
            end
            
            % Verify timeWindow array size
            if size(tWin,1) ~= tbHeight && size(tWin,1) ~= length(rowInd) && size(tWin,2) ~= 2
                error('The array of timeWindow has incorrect size.');
            end
            
            % Verify time increment (ignore NaN windows)
            if ~all(tWin(:,2)-tWin(:,1) > 0 | any(isnan(tWin),2))
                error('The window end time must be greater than the begin time.');
            end
            
            % Propagate timeWindow value for all trials
            if numel(tWin) == 2
                tWin = repmat(tWin(:)', length(rowInd), 1);
            end
            
            % Make timeWindow full size
            if size(tWin,1) < tbHeight
                tWinFull = NaN(tbHeight, 2);
                tWinFull(rowInd,:) = tWin;
                tWin = tWinFull;
            end
        end
    end
    
    % Hide methods from handle class
    methods(Hidden)
        function addlistener(this)
        end
        function notify(this)
        end
        function findobj(this)
        end
        function findprop(this)
        end
    end
end









