classdef MPlot
    %MPlot A collection of functions useful for plotting
    %   
    
    methods(Static)
        function Blocks(xRange, yRange, color, varargin)
            % Plot rectangle blocks based on X anf Y ranges (based on third party function drawShadedRectangle)
            % 
            %   MPlot.Blocks(xRange, yRange)
            %   MPlot.Blocks(xRange, yRange, color)
            %   MPlot.Blocks(..., 'gradientOrientation', 'vertical')
            %
            % Inputs:
            %   xRange      Boundary coordinates of block(s) along X-axis in a n-by-2 matrix where n is 
            %               the number of block(s) to plot. Alternatively, you can provide logical vector 
            %               for X-axis where region(s) containing block(s) are set to true. If multiple 
            %               Y-ranges are provided, 
            %   yRange      Same cnvention as xRange but for Y-axis. 
            %   color       A 1-by-3 matrix of RGB color for uniform coloring or a 3-by-3 matrix of row 
            %               vectors for a gradient of colors. This parameter applies to all blocks. The 
            %               default color is uniform gray ([.9 .9 .9]).
            %
            
            % Handles user inputs
            if nargin < 3
                % Default color is gray
                color = [.9 .9 .9];
            end
            
            % Convert logical mask to boundaries
            if islogical(xRange)
                xRange = MMath.Logical2Bounds(xRange);
            end
            if islogical(yRange)
                yRange = MMath.Logical2Bounds(yRange);
            end
            
            % Apply yRange to all xRanges if necessary
            if numel(yRange) == 2
                yRange = repmat(yRange, size(xRange,1), 1);
            end
            if numel(xRange) == 2
                xRange = repmat(xRange, size(yRange,1), 1);
            end
            
            % Plots blocks
            xx = [xRange, flip(xRange,2)]';
            yy = yRange(:,[1 1 2 2])';
            patch(xx, yy, color, 'EdgeColor', 'none', varargin{:});
        end
        
        function h = Circle(x, y, r, c)
            % Plot a shape of circle
            % 
            %   h = MPlot.Circle(x, y, r, c)
            %
            % Inputs:
            %   x       X-coordinate of the center
            %   y       Y-coordinate of the center
            %   r       Radius of the circle
            %   c       Color of the circle
            % Output:
            %   h       Object handle of the circle shape. By nature, it is a rectangle with rounded corners. 
            %
            
            d = r*2;
            px = x-r;
            py = y-r;
            h = rectangle('Position', [ px py d d ], 'Curvature', [1,1], 'FaceColor', c, 'LineStyle', 'none');
            daspect([ 1 1 1 ]);
        end
        
        function cc = Color2Str(cc)
            cc = mat2str(cc);
            cc = strrep(cc(2:end), ' ', ',');
            cc = strsplit(cc, ';');
        end
        
        function ErrorShade(varargin)
            % Plot error as shading
            % 
            %   MPlot.ErrorShade(y, err)
            %   MPlot.ErrorShade(x, y, err)
            %   MPlot.ErrorShade(x, y, errPos, errNeg)
            %   MPlot.ErrorShade(..., 'IsRelative', true)
            %   MPlot.ErrorShade(..., 'Orientation', 'vertical')
            %   MPlot.ErrorShade(..., 'Color', 'k')
            %   MPlot.ErrorShade(..., 'Alpha', 0.3)
            %
            % Inputs:
            %   x               X-coordinates. Default is indices of elements in y. 
            %   y               Y-coordinates. 
            %   err             Errors, which apply to both sides. 
            %   errPos          Errors on the positive side. 
            %   errNeg          Errors on the negative side.
            %   'IsRelative'    Logical variable indicate whether (default) or not error inputs are relative to y. 
            %   'Orientation'   Orientation along which err is applied. 'vertical' (default) for Y-axis 
            %                   and 'horizontal' for X-axis.
            %   'Color'         Color of the shade. Default is black. 
            %   'Alpha'         Transparancy of the shade. Default 0.3. 
            %
            
            % Handles user inputs
            p = inputParser();
            p.addRequired('arg1');
            p.addRequired('arg2');
            p.addOptional('arg3', []);
            p.addOptional('arg4', []);
            p.addParameter('IsRelative', true, @islogical);
            p.addParameter('Color', 'k');
            p.addParameter('Alpha', 0.3, @isnumeric);
            p.addParameter('Orientation', 'vertical', @(x) any(strcmp(x, {'vertical', 'horizontal'})));
            
            p.parse(varargin{:});
            arg1 = p.Results.arg1;
            arg2 = p.Results.arg2;
            arg3 = p.Results.arg3;
            arg4 = p.Results.arg4;
            if ~isempty(arg4)
                x = arg1;
                y = arg2;
                errPos = arg3;
                errNeg = arg4;
            elseif ~isempty(arg3)
                x = arg1;
                y = arg2;
                errPos = arg3;
                errNeg = arg3;
            else
                y = arg1;
                x = cumsum(ones(size(y)));
                errPos = arg2;
                errNeg = arg2;
            end
            isRelative = p.Results.IsRelative;
            color = p.Results.Color;
            faceAlpha = p.Results.Alpha;
            ori = p.Results.Orientation;
            
            if isvector(y)
                x = x(:);
                y = y(:);
                errPos = errPos(:);
                errNeg = errNeg(:);
            end
            
            % Ploting
            for k = 1 : size(y,2)
                x = [x(:,k); flip(x(:,k))];
                if isRelative
                    err = [y(:,k)+errPos(:,k); flip(y(:,k)-errNeg(:,k))];
                else
                    err = [errPos(:,k); flip(errNeg(:,k))];
                end
                if strcmp(ori, 'vertical')
                    patch(x, err, color, 'FaceAlpha', faceAlpha, 'LineStyle', 'none');
                else
                    patch(err, x, color, 'FaceAlpha', faceAlpha, 'LineStyle', 'none');
                end
            end
        end
        
        function Paperize(varargin)
            % Make axes comply with conventions of publication
            %
            %   MPlot.Paperize(h)
            %   MPlot.Paperize(..., 'FontSize', 6)
            %   MPlot.Paperize(..., 'ColumnsWide', [])
            %   MPlot.Paperize(..., 'ColumnsHigh', [])
            %   MPlot.Paperize(..., 'AspectRatio', [])
            %   MPlot.Paperize(..., 'JournalStyle', 'Cell')
            %
            % Inputs:
            %   h           Array of Axes or Figure handle(s). Default is the current figure.
            %               If h is empty, all existing Axes will be operated on. 
            %   fontSize    Font size. Default 6. 
            %
            
            p = inputParser();
            p.addOptional('h', gcf, @ishandle);
            p.addParameter('FontSize', 6, @isscalar);
            p.addParameter('ColumnsWide', [], @isscalar);
            p.addParameter('ColumnsHigh', [], @isscalar);
            p.addParameter('AspectRatio', [], @isscalar);
            p.addParameter('JournalStyle', 'cell', @(x) any(strcmpi(x, {'nature', 'cell'})));
            
            p.parse(varargin{:});
            h = p.Results.h;
            fontSize = p.Results.FontSize;
            colsWide = p.Results.ColumnsWide;
            colsHigh = p.Results.ColumnsHigh;
            aRatio = p.Results.AspectRatio;
            journalStyle = lower(p.Results.JournalStyle);
            
            switch journalStyle
                case 'nature'
                    widthSet = [8.9 12 18.3];
                case 'cell'
                    widthSet = [8.5 11.4 17.4];
            end
            
            % Resolve figure width
            if ~isempty(colsWide)
                % Calculate width by fold of cols
                figWidth = widthSet(1) * colsWide;
                
                % Overwrite if at specific #cols
                colOpts = [1 1.5 2];
                optIdx = colsWide == colOpts;
                if any(optIdx)
                    figWidth = widthSet(colsWide == colOpts);
                end
            end
            
            % 
            if isempty(h)
                h = findobj('Type', 'Axes');
            end
            
            for i = 1 : numel(h)
                if isa(h(i), 'matlab.ui.Figure')
                    if ~isempty(colsWide)
                        h(i).Color = 'w';
                        h(i).Units = 'centimeter';
                        h(i).Position(3) = figWidth;
                        if ~isempty(colsHigh)
                            h(i).Position(4) = figWidth / colsWide * colsHigh;
                        elseif ~isempty(aRatio)
                            h(i).Position(4) = figWidth * aRatio;
                        end
                    end
%                     tightfig(h(i));
                    ax = findobj(h, 'Type', 'Axes');
                else
                    ax = h(i);
                end
                
                for j = 1 : numel(ax)
                    set(ax(j), 'TickDir', 'out', 'FontSize', fontSize, 'FontName', 'arial');
                end
            end
        end
        
        function varargout = PlotPointAsLine(x, y, d, varargin)
            % Plot points as lines with specified length. It also accepts arguments of Line properties 
            % like 'line' function. 
            % 
            %   MPlot.PlotPointAsLine(x, y, d)
            %   MPlot.PlotPointAsLine(..., 'Orientation', 'vertical')
            %   MPlot.PlotPointAsLine(..., 'LinePropertyName', value)
            %   hh = MPlot.PlotPointAsLine(...)
            % 
            % Inputs
            %   x                   x coordinates of line centers. 
            %   y                   y coordinates of line centers. 
            %   d                   Length of lines. It can be a scalar or an array for each line. 
            %   'Orientation'       The orientation of lines. Default is 'vertical'.
            %   Any PropertyName-Value pair for MATLAB built-in 'line' function. 
            % 
            % Output
            %   hh                  Handles of Line objects
            % 
            
            p = inputParser();
            p.KeepUnmatched = true;
            p.addParameter('Orientation', 'vertical', @(x) any(strcmpi(x, {'vertical', 'horizontal'})));
            
            p.parse(varargin{:});
            orientationStr = p.Results.Orientation;
            
            varargin = p.Unmatched;
            
            x = x(:);
            y = y(:);
            d = d(:);
            
            if strcmpi(orientationStr, 'horizontal')
                xx = [x-d/2, x+d/2, NaN(size(x))]';
                yy = [y, y, NaN(size(y))]';
            else
                xx = [x, x, NaN(size(x))]';
                yy = [y-d/2, y+d/2, NaN(size(y))]';
            end
            
            xx = xx(:);
            yy = yy(:);
            
            hh = plot(xx, yy, varargin);
            
            if nargout == 1
                varargout{1} = hh;
            end
        end
        
        function varargout = PlotTraceLadder(varargin)
            % Plot traces as a ladder
            % 
            %   MPlot.PlotTraceLadder(yy)
            %   MPlot.PlotTraceLadder(xx, yy)
            %   MPlot.PlotTraceLadder(xx, yy, yPos)
            %   MPlot.PlotTraceLadder(..., 'ColorArray', [])
            %   MPlot.PlotTraceLadder(..., 'TraceFun', [])
            %   MPlot.PlotTraceLadder(..., 'LinePropertyName', value)
            %   hh = MPlot.PlotTraceLadder(...)
            % 
            % Inputs:
            %   xx                  1) A vector of x coordinates for all traces in yy. 
            %                       2) A matrix of x coordinates for every point in yy.
            %   yy                  Original y coordinates for each trace in column. 
            %   yPos                Y position of each trace's zero after shifting them into a ladder. 
            %   'ColorArray'        1) An n-by-3 array of RGB colors. n is the number of traces. 
            %                       2) An n-by-4 array of RGBA colors. 'A'(alpha) controls transparency. 
            %                       3) An n-element char vector of colors. (e.g. 'k', 'r', 'm')
            %   'TraceFun'          A function handle that process each trace in yy. This function must not change
            %                       the number of elements in each trace. 
            %   Any PropertyName-Value pair for MATLAB built-in 'line' function. 
            % 
            % Output:
            %   hh                  Handles of plotted line objects. 
            %
            
            % Handle inputs
            p = inputParser();
            p.KeepUnmatched = true;
            p.PartialMatching = false;
            
            p.addRequired('arg1');
            p.addOptional('arg2', []);
            p.addOptional('yPos', []);
            p.addParameter('ColorArray', [], @isnumeric);
            p.addParameter('TraceFunc', [], @(x) isa(x,'function_handle'));
            
            p.parse(varargin{:});
            arg1 = p.Results.arg1;
            arg2 = p.Results.arg2;
            yPos = p.Results.yPos;
            colorArray = p.Results.ColorArray;
            traceFunc = p.Results.TraceFunc;
            
            varargin = p.Unmatched;
            
            if ~isempty(arg2)
                xx = arg1;
                yy = arg2;
            else
                yy = arg1;
                xx = [];
            end
            
            if isvector(yy)
                yy = yy(:);
            end
            if isempty(xx)
                xx = 1 : size(yy,1);
            end
            if isvector(xx)
                xx = repmat(xx(:), 1, size(yy,2));
            end
            if ~isempty(traceFunc)
                for i = 1 : size(yy,2)
                    yy(:,i) = traceFunc(yy(:,i));
                end
            end
            if isempty(yPos)
                yPos = cumsum(-min(yy) + [0, max(yy(:,1:end-1))]);
            end
            yPos = yPos(:)';
            
            % Plot traces
            hh = plot(xx, yy+yPos, varargin);
            
            if ~isempty(colorArray)
                if isvector(colorArray)
                    colorArray = colorArray(:);
                end
                
                for i = 1 : length(hh)
                    hh(i).Color = colorArray(i,:);
                end
            end
            
            if nargout == 1
                varargout{1} = hh;
            end
        end
        
        function cc = Rainbow(numColors)
            % Returns a color map based on rainbow spectrum
            %
            %   cc = MPlot.Rainbow(numColors)
            %
            % Input:
            %   numColors       The number of colors to return
            % Output:
            %   cc              A series of color in numColors-by-3 matrix ranging from red to violet. 
            %
            
            cBase = { [ 1 1 0 0 0 1 ]; [ 0 1 1 1 0 0 ]; [ 0 0 0 1 1 1 ] };
            cc = cellfun(@(x) interp1(1:6, x, linspace(1,6,numColors)), cBase, 'UniformOutput', false);
            cc = cell2mat(cc)';
        end
        
        function SaveFigure(h, filePath)
            % Save figure as a PNG file
            %
            %   MPlot.SaveFigure(h, filePath)
            %
            
            f = getframe(h);
            imwrite(frame2im(f), [ filePath '.png' ]);
        end
        
        function barLength = ScaleBar(data, heightRatio)
            % Return the length of a scale bar with appropriate length
            %
            %   barLength = MPlot.ScaleBar(data, heightRatio)
            %
            % Inputs:
            %   data                Numeric array. 
            %   heightRatio         The desired ratio of the length of scale bar to the span of data. 
            %                       Default ratio is 0.5.
            % Output:
            %   barLength           The length of scale bar in the same unit as data. 
            %
            
            if nargin < 2
                heightRatio = 0.5;
            end
            
            if isvector(data)
                data = data(:);
            end
            
            barLength = (max(data) - min(data)) * heightRatio;
            barLength = round(barLength, 1, 'significant');
        end
        
        function ax = SetAxes(figHandle, gridSize, panelCoor, panelSize, panelSpacing, gridSpacing)
            % Set an axes with hard-coded positions
            % 
            %   MPlot.SetAxes(figHandle, gridSize, panelCoor, panelSize, panelSpacing, gridSpacing)
            %
            
            if numel(panelSize) == 1
                panelSize = panelSize([1 1]);
            end
            
            if numel(panelSpacing) == 1
                panelSpacing = panelSpacing([1 1]);
            end
            
            if numel(gridSpacing) == 1
                gridSpacing = gridSpacing([1 1]);
            end
            
            % Make sure the figure has the right size
            figSize = gridSpacing*2 + panelSize.*flip(gridSize) + panelSpacing.*flip(gridSize-1);
            
            figPos = get(figHandle, 'Position');
            set(gcf, 'Position', [figPos(1:2), figSize]);
            
            % Compute panel parameters
            panelPos(1) = gridSpacing(1) + (panelCoor(2)-1) * (panelSize(1) + panelSpacing(1));
            panelPos(2) = figSize(2) - gridSpacing(2) - (panelCoor(1)-1) * (panelSize(2) + panelSpacing(2)) - panelSize(2);
            
            ax = axes;
            set(ax, 'Unit', 'pixel', 'Position', [panelPos, panelSize]);
        end
        
        function Violin(nn, bb, hh, varargin)
            % Plot histograms as violins
            % 
            %   MPlot.Violin(nn, bb, hh)
            %   MPlot.Violin(nn, bb, hh, s)
            %   MPlot.Violin(..., 'Color', 'k')
            %   MPlot.Violin(..., 'Alpha', 0.3)
            %   MPlot.Violin(..., 'Orientation', 'vertical')
            %   MPlot.Violin(..., 'Alignment', 'center')
            %
            % Inputs:
            %   nn              
            %   bb              
            %   hh              
            %   s               
            %   'Color'         Color of the shade. Default is black. 
            %   'Alpha'         Transparancy of the shade. Default 0.3. 
            %   'Orientation'   Orientation of violins. Default is 'vertical'.
            %   'Alignment'     
            %
            
            % Handles user inputs
            p = inputParser();
            p.addOptional('s', 1, @isscalar);
            p.addParameter('Quantiles', [], @isnumeric);
            p.addParameter('Color', 'k');
            p.addParameter('Alpha', 1, @isnumeric);
            p.addParameter('Orientation', 'vertical', @(x) any(strcmpi(x, {'vertical', 'horizontal'})));
            p.addParameter('Alignment', 'center', @(x) any(strcmpi(x, {'low', 'high', 'center'})));
            
            p.parse(varargin{:});
            s = p.Results.s;
            color = p.Results.Color;
            faceAlpha = p.Results.Alpha;
            qt = p.Results.Quantiles;
            ori = p.Results.Orientation;
            alignType = p.Results.Alignment;
            
            % Ploting
            hold on;
            for k = 1 : size(hh,2)
                b = [bb(:,k); flip(bb(:,k))];
                h = hh(:,k)/s;
                
                switch alignType
                    case 'center'
                        h = nn(k) + [-h/2; flip(h/2)];
                    case 'low'
                        h = nn(k) + [zeros(size(h)); flip(h)] - nanmax(h)/2;
                    case 'high'
                        h = nn(k) + [-h; zeros(size(h))] + nanmax(h)/2;
                end
                
                if strcmp(ori, 'vertical')
                    patch(h, b, color, 'FaceAlpha', faceAlpha, 'LineStyle', 'none');
                else
                    patch(b, h, color, 'FaceAlpha', faceAlpha, 'LineStyle', 'none');
                end
            end
            
            function computeMarks(d, qtVect)
                
            end
        end
    end
    
end




