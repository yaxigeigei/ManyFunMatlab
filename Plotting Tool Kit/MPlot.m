classdef MPlot
    %MPLOT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        function h = BarSurf(v)
            h = bar3(v, 1);
            for k = 1 : length(h)
                zdata = get(h(k), 'Zdata');
                set(h(k), 'CData', zdata, 'FaceColor', 'interp');
            end
        end
        
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
            %               default color is uniform gray ([.8 .8 .8]). 
            %   'gradientOrientation'
            %               The orientation of gradient, 'horizontal' or 'vertical' (default).
            %
            
            % Handles user inputs
            p = inputParser();
            p.addParameter('gradientOrientation', 'vertical', @ischar);
            p.parse(varargin{:});
            gradOrient = p.Results.gradientOrientation;
            
            if nargin < 3
                % Default color is gray
                color = [ .8 .8 .8 ];
            end
            if isvector(color)
                % No gradient by default
                color = repmat(color, 3, 1);
            end
            
            % Convert logical mask to boundaries
            if islogical(xRange)
                xRange = MMath.Logical2Boundaries(xRange);
            end
            if islogical(yRange)
                yRange = MMath.Logical2Boundaries(yRange);
            end
            
            % Apply yRange to all xRanges if necessary
            if numel(yRange) == 2
                yRange = repmat(yRange, size(xRange,1), 1);
            end
            if numel(xRange) == 2
                xRange = repmat(xRange, size(yRange,1), 1);
            end
            
            % Plots blocks
            for i = 1 : size(xRange,1)
                drawShadedRectangle(xRange(i,:), yRange(i,:), color(1,:), color(2,:), color(3,:), gradOrient);
            end
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
        
        function rainbowCode = ColorCodeRainbow(numColors)
            % Returns a color map based on rainbow spectrum
            %
            %   rainbowCode = MPlot.ColorCodeRainbow(numColors)
            %
            % Input:
            %   numColors       The number of colors to return
            % Output:
            %   rainbowCode     A series of color in numColors-by-3 matrix ranging from red to violet. 
            %
            
            rainbowBase = { [ 1 1 0 0 0 1 ]; [ 0 1 1 1 0 0 ]; [ 0 0 0 1 1 1 ] };
            rainbowCode = cellfun(@(x) interp1(1:6, x, linspace(1,6,numColors)), rainbowBase, 'UniformOutput', false);
            rainbowCode = cell2mat(rainbowCode)';
        end
        
        function ErrorShade(y, err, varargin)
            % Plot error as shading
            % 
            %   MPlot.ErrorShade(y, err)
            %   MPlot.ErrorShade(y, err, x)
            %   MPlot.ErrorShade(..., 'orientation', 'vertical')
            %   MPlot.ErrorShade(..., 'color', 'k')
            %   MPlot.ErrorShade(..., 'alpha', 0.3)
            %
            % Inputs:
            %   y               Y-coordinates.
            %   err             Errors, one sided. 
            %   x               X-coordinates. Default is indices of member in y. 
            %   'orientation'   Orientation along which err is applied. 'vertical' (default) for Y-axis 
            %                   and 'horizontal' for X-axis.
            %   'color'         Color of the shade. Default is black. 
            %   'alpha'         Transparancy of the shade. Default 0.3. 
            %
            
            % Handles user inputs
            p = inputParser();
            p.addOptional('x', []);
            p.addParameter('color', 'k');
            p.addParameter('alpha', 0.3, @isnumeric);
            p.addParameter('orientation', 'vertical', @(x) any(strcmp(x, {'vertical', 'horizontal'})));
            
            p.parse(varargin{:});
            x = p.Results.x;
            color = p.Results.color;
            faceAlpha = p.Results.alpha;
            ori = p.Results.orientation;
            
            if isempty(x)
                x = 1 : length(y);
            end
            
            x = x(:);
            y = y(:);
            err = err(:);
            
            % Ploting
            if strcmp(ori, 'vertical')
                patch([ x; flip(x) ], [ y+err; flip(y-err) ], color, ...
                    'FaceAlpha', faceAlpha, ...
                    'LineStyle', 'none');
            else
                patch([ y+err; flip(y-err) ], [ x; flip(x) ], color, ...
                    'FaceAlpha', faceAlpha, ...
                    'LineStyle', 'none');
            end
        end
        
        function Figure(n, w, h, varargin)
            % Create or modify a new figure window
            % 
            %   Figure(h, varargin)
            % 
            
%             % Handle user input
%             p = inputParser();
%             p.addOptional('x', []);
%             p.addParameter('color', 'k');
%             p.addParameter('alpha', 0.3, @isnumeric);
%             p.addParameter('orientation', 'vertical', @(x) any(strcmp(x, {'vertical', 'horizontal'})));
%             
%             p.parse(varargin{:});
%             x = p.Results.x;
%             color = p.Results.color;
%             faceAlpha = p.Results.alpha;
%             ori = p.Results.orientation;
            
            figure(n);
            figPos = get(gcf, 'Position');
            
            if ~isempty(w)
                figPos(3) = w;
            end
            
            if ~isempty(h)
                figPos(4) = h;
            end
            
            set(gcf, 'Color', 'w', 'Position', figPos);
            
        end
        
        function Paperize(varargin)
            % Make an axes comply with conventions of publication
            %
            %   MPlot.Paperize(h)
            %   MPlot.Paperize(h, fontSize)
            %
            % Inputs:
            %   h           Axes handle. Default is the current axes.
            %   fontSize    Font size. Default 6. 
            %
            
            p = inputParser();
            p.addOptional('h', gca, @ishandle);
            p.addOptional('fontSize', 6, @isscalar);
            
            p.parse(varargin{:});
            h = p.Results.h;
            fontSize = p.Results.fontSize;
            
            set(h, 'TickDir', 'out', 'FontSize', fontSize, 'FontName', 'arial');
            
        end
        
        function Pcolor(x, y, c)
            % Pseudocolor plot similar to MATLAB pcolor() but avoids "missing" edges
            %
            %   MPlot.Pcolor(x, y, c)
            %
            
            x = [x(:); x(end) + diff(x(end-1:end))];
            y = [y(:); y(end) + diff(y(end-1:end))];
            c = [c, c(:,end)];
            c = [c; c(end,:)];
            
            pcolor(x, y, c);
            
        end
        
        function SaveFigure(h, filePath)
            % Save figure as a PNG file
            %
            %   MPlot.SaveFigure(h, filePath)
            %
            
            f = getframe(h);
            imwrite(frame2im(f), [ filePath '.png' ]);
        end
        
        function ax = SetAxes(figHandle, gridSize, panelCoor, panelSize, panelSpacing, gridSpacing)
            % Set an axes with hard-coded positions
            % 
            %   SetAxes(figHandle, gridSize, panelCoor, panelSize, panelSpacing, gridSpacing)
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
        
        function valBar = ScaleBar(data, percentHeight)
            % Return the length of a scale bar with appropriate length
            %
            %   valBar = MPlot.ScaleBar(data, percentHeight)
            %
            % Inputs:
            %   data                Vector of data. 
            %   percentHeight       Desiered height of the scale bar wrt the percentage of the span of data. 
            % Output:
            %   valBar              The length of scale bar in the same unit as data. 
            %
            
            if nargin < 2
                percentHeight = 0.5;
            end
            
            if isvector(data)
                data = data(:);
            end
            
            valBarRaw = (max(data) - min(data)) * percentHeight;
            valBar = round(valBarRaw, 1, 'significant');
        end
    end
    
end

