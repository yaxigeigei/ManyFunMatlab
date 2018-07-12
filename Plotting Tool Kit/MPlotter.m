classdef MPlotter < handle
    %MPlotter Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        layoutTable;
        plotTable;
        notes;
    end
    
    properties(Access = private)
        gui;
    end
    
    methods
        function this = MPlotter()
            % Constructor of MPlotter
            
            % Initialize layout table
            columnNames = {'figureNumber', 'figureObj', 'subplot', 'axesObj', 'functionHandle', 'variableName'};
            exampleRows = { ...
                1, [], '3,1,1', [], [], ''; ...
                1, [], '3,1,2:3', [], [], ''};
            this.layoutTable = cell2table(exampleRows, 'VariableNames', columnNames);
            
            % Initialize plot table
            columnNames = {'axesIdx', 'functionHandle', 'variableName', 'updateOption'};
            exampleRows = { ...
                1, @MPlotter.PlotTimeIndicator, '', 'time'; ...
                2, @MPlotter.PlotTimeIndicator, '', 'trial'};
            this.plotTable = cell2table(exampleRows, 'VariableNames', columnNames);
            
            % Open GUI
            this.GUI();
        end
        
        function GUI(this)
            % Open GUI to add, modify and control plots
            
            % Close existing GUI
            try
                close(this.gui.fig);
                this.gui = [];
            catch
            end
            
            % Window
            wh = 150; ww = 500;
            uh = 20; uw = 60;
            s = 5;
            
            this.gui.fig = figure(...
                'Name', 'MSessionPlotter', ...
                'NumberTitle', 'off', ...
                'IntegerHandle', 'off', ...
                'MenuBar', 'none', ...
                'Resize', 'off', ...
                'WindowKeyPressFcn', @this.KeyPress, ...
                'WindowKeyReleaseFcn', @this.KeyRelease, ...
                'CloseRequestFcn', @this.CloseRequest, ...
                'Interruptible', 'off', ...
                'BusyAction', 'cancel');
            
            figPos = get(this.gui.fig, 'Position');
            set(this.gui.fig, 'Position', [figPos(1:2) ww wh]);
            
            % Trial
            x = s; y = s;
            uicontrol(this.gui.fig, ...
                'Style', 'text', ...
                'String', 'Trial #  ', ...
                'FontSize', 9, ...
                'HorizontalAlignment', 'right', ...
                'Units', 'pixel', ...
                'Position', [x y-2 45 uh]);
            x = x + 45;
            
            this.gui.trialEdit = uicontrol(this.gui.fig, ...
                'Style', 'edit', ...
                'String', '1', ...
                'FontSize', 9, ...
                'HorizontalAlignment', 'left', ...
                'Units', 'pixel', ...
                'Position', [x y 45 uh], ...
                'KeyReleaseFcn', @this.TrialEditChange, ...
                'Interruptible', 'off', ...
                'BusyAction', 'cancel');
            x = x + 45;
            
            % Time
            x = x + 10;
            uicontrol(this.gui.fig, ...
                'Style', 'text', ...
                'String', 'Time  ', ...
                'FontSize', 9, ...
                'HorizontalAlignment', 'right', ...
                'Units', 'pixel', ...
                'Position', [x y-2 40 uh]);
            x = x + 40;
            
            this.gui.timeEdit = uicontrol(this.gui.fig, ...
                'Style', 'edit', ...
                'String', '0', ...
                'FontSize', 9, ...
                'HorizontalAlignment', 'left', ...
                'Units', 'pixel', ...
                'Position', [x y 60 uh], ...
                'KeyReleaseFcn', @this.TimeEditChange, ...
                'Interruptible', 'off', ...
                'BusyAction', 'cancel');
            x = x + 60;
            
            % Limits
            x = x + 15;
            uicontrol(this.gui.fig, ...
                'Style', 'text', ...
                'String', 'Limits  ', ...
                'FontSize', 9, ...
                'HorizontalAlignment', 'right', ...
                'Units', 'pixel', ...
                'Position', [x y-2 40 uh]);
            x = x + 40;
            
            this.gui.limitEdit1 = uicontrol(this.gui.fig, ...
                'Style', 'edit', ...
                'String', '-Inf', ...
                'FontSize', 9, ...
                'HorizontalAlignment', 'left', ...
                'Units', 'pixel', ...
                'Position', [x y 60 uh], ...
                'KeyReleaseFcn', @this.LimitEditChange, ...
                'Interruptible', 'off', ...
                'BusyAction', 'cancel');
            x = x + 60;
            
            this.gui.limitEdit2 = uicontrol(this.gui.fig, ...
                'Style', 'edit', ...
                'String', 'Inf', ...
                'FontSize', 9, ...
                'HorizontalAlignment', 'left', ...
                'Units', 'pixel', ...
                'Position', [x y 60 uh], ...
                'KeyReleaseFcn', @this.LimitEditChange, ...
                'Interruptible', 'off', ...
                'BusyAction', 'cancel');
            x = x + 60;
            
            % Refresh
            x = x + 15;
            uicontrol(this.gui.fig, ...
                'Style', 'pushbutton', ...
                'String', 'Refresh all', ...
                'FontSize', 9, ...
                'Units', 'pixel', ...
                'Position', [x y-1 100 uh+2], ...
                'Callback', @this.UpdateRoutine, ...
                'Interruptible', 'off', ...
                'BusyAction', 'cancel');
            x = x + 80;
            
            % Note
            x = 0;
            y = y + uh + s;
            
            this.gui.noteEdit = uicontrol(this.gui.fig, ...
                'Style', 'edit', ...
                'Max', Inf, ...
                'FontSize', 9, ...
                'HorizontalAlignment', 'left', ...
                'Units', 'pixel', ...
                'Position', [x y ww wh-uh-y]);
            y = wh-uh;
            
            this.gui.noteMenu = uicontrol(this.gui.fig, ...
                'Style', 'popupmenu', ...
                'String', {'Write note in the text box below'}, ...
                'FontSize', 9, ...
                'HorizontalAlignment', 'left', ...
                'Units', 'pixel', ...
                'Position', [x y ww uh]);
            
            % Update GUI
            this.UpdateRoutine();
            
            % Show instruction
%             eval('help MImgBaseClass.Viewer');
        end
    end
    
    methods(Static)
        function PlotTimeIndicator(ax, inputVar)
            
            tr = ax.UserData.trialNum;
            t = ax.UserData.time;
            tLims = ax.UserData.timeLimits;
            
            if isfield(ax.UserData, 'indicatorObj') && ~isempty(ax.UserData.indicatorObj) && ishandle(ax.UserData.indicatorObj)
                ax.UserData.indicatorObj.XData = [t t]';
                ax.UserData.indicatorObj.YData = ax.YLim';
            else
                ax.UserData.indicatorObj = plot(ax, [t t]', ax.YLim', 'LineWidth', 2, 'Color', [0 0 0 .2]);
                hold(ax, 'on');
            end
            
            ax.XLim = tLims;
        end
    end
    
    methods(Access = private)
        function TrialEditChange(this, src, eventdata)
            if strcmp(eventdata.Key, 'return')
                this.UpdateRoutine('trial', str2double(src.String));
            end
        end
        
        function TimeEditChange(this, src, eventdata)
            if strcmp(eventdata.Key, 'return')
                this.UpdateRoutine('time', str2double(src.String));
            end
        end
        
        function LimitEditChange(this, src, eventdata)
            if strcmp(eventdata.Key, 'return')
                this.UpdateRoutine('time');
            end
        end
        
        function KeyPress(this, src, eventdata)
            
            uiTrialNum = str2double(this.gui.trialEdit.String);
            uiTime = str2double(this.gui.timeEdit.String);
            
            dTrial = 1;
            dTime = 0.001;
            if any(strcmp(eventdata.Modifier, 'control'))
                dTrial = dTrial * 2;
                dTime = dTime * 5;
            end
            if any(strcmp(eventdata.Modifier, 'shift'))
                dTrial = dTrial * 10;
                dTime = dTime * 20;
            end
            
            switch eventdata.Key
                case 'rightarrow'
                    this.UpdateRoutine('time', uiTime + dTime);
                case 'leftarrow'
                    this.UpdateRoutine('time', uiTime - dTime);
                case 'downarrow'
                    this.UpdateRoutine('trial', uiTrialNum + dTrial);
                case 'uparrow'
                    this.UpdateRoutine('trial', uiTrialNum - dTrial);
            end
            
            pause(0.07)
        end
        
        function KeyRelease(this, src, eventdata)
            val = str2double(eventdata.Key(end));
            if ~isnan(val)
                
            end
        end
        
        function CloseRequest(this, src, eventdata)
            for i = 1 : height(this.layoutTable)
                try
                    delete(this.layoutTable.figureObj{i});
                catch
                end
            end
            delete(gcf);
        end
        
        function IntializeLayout(this)
            % Create figures and axes according to layoutTable
            
            isNewFig = false;
            
            for i = 1 : height(this.layoutTable)
                
                % Unpack variables
                figNum = this.layoutTable.figureNumber(i);
                figObj = this.layoutTable.figureObj{i};
                sp = this.layoutTable.subplot{i};
                sp = eval(['{' sp '}']);
                axesObj = this.layoutTable.axesObj{i};
                funcHandle = this.layoutTable.functionHandle{i};
                varName = this.layoutTable.variableName{i};
                
                % Make figure if not exsiting
                if isempty(figObj) || ~ishandle(figObj) || ~isvalid(figObj)
                    this.layoutTable.figureObj{i} = figure(figNum);
                    isNewFig = true;
                end
                
                % Make axes if not exsiting
                if isempty(axesObj) || ~ishandle(axesObj) || ~isvalid(axesObj)
                    this.layoutTable.axesObj{i} = subplot(sp{:}, 'Parent', figNum);
                    
                    if ~isempty(funcHandle) && isa(funcHandle, 'function_handle')
                        try
                            funcHandle(figObj, axesObj, evalin('base', varName));
                        catch
                        end
                    end
                end
            end
            
            % Trun focus back to main window
            if isNewFig
                figure(this.gui.fig);
            end
            
        end
        
        function UpdateRoutine(this, updateType, newVal)
            
            if nargin < 2
                updateType = 'trial';
            end
            
            % Apply new value to UI
            if nargin > 2
                switch updateType
                    case 'time'
                        this.gui.timeEdit.String = num2str(newVal);
                    case 'trial'
                        this.gui.trialEdit.String = num2str(newVal);
                end
            end
            
            % Get variables
            trialNum = str2double(this.gui.trialEdit.String);
            timePt = str2double(this.gui.timeEdit.String);
            timeLims(1) = str2double(this.gui.limitEdit1.String);
            timeLims(2) = str2double(this.gui.limitEdit2.String);
            
            % Make figures and axes if not present
            this.IntializeLayout();
            
            % Iterate through plots
            for i = 1 : height(this.plotTable)
                
                % Unpack plotting variables
                axIdx = this.plotTable.axesIdx(i);
                axObj = this.layoutTable.axesObj{axIdx};
                axObj.UserData.trialNum = trialNum;
                axObj.UserData.time = timePt;
                
                if ~isnan(timeLims(1))
                    axObj.UserData.timeLimits(1) = timeLims(1);
                else
                    this.gui.limitEdit1.String = num2str(axObj.UserData.timeLimits(1));
                end
                
                if ~isnan(timeLims(2))
                    axObj.UserData.timeLimits(2) = timeLims(2);
                else
                    this.gui.limitEdit2.String = num2str(axObj.UserData.timeLimits(2));
                end
                
                varName = this.plotTable.variableName{i};
                funcHandle = this.plotTable.functionHandle{i};
                updateOption = this.plotTable.updateOption{i};
                
                % Check the need for update
                switch updateOption
                    case 'time'
                        % always update
                    case 'trial'
                        if strcmp(updateType, {'time'})
                            continue;
                        end
                    case 'manual'
                        if any(strcmp(updateType, {'time', 'trial'}))
                            continue;
                        end
                    otherwise
                        continue;
                end
                
                % Run function
                if isempty(varName)
                    funcHandle(axObj);
                else
                    funcHandle(axObj, evalin('base', varName));
                end
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


