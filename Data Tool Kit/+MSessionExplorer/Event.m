classdef Event
    %Event Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Dependent)
        t;
        T;
        V;
        isValid;
    end
    
    properties(Hidden)
        i_t = NaN;
        i_T = struct();
        i_V;
    end
    
    methods
        function this = Event(t, T, V)
            % MSessionExplorer.Event constructs an instance of this class
            %
            %   obj = MSessionExplorer.Event()
            %   obj = MSessionExplorer.Event(t)
            %   obj = MSessionExplorer.Event(t, T)
            %   obj = MSessionExplorer.Event(t, T, V)
            %
            if nargin > 0
                this.t = t;
            end
            if nargin > 1
                this.T = T;
            end
            if nargin > 2
                this.V = V;
            end
        end
        
        % Setters and Getters
        function this = set.t(this, val)
            assert(isscalar(val) && isnumeric(val), 't must be a numeric scalar');
            this.i_t = val;
        end
        function val = get.t(this)
            val = this.i_t;
        end
        function this = set.T(this, val)
            if isempty(val)
                val = struct();
            end
            assert(isstruct(val), 'eventTimes must be a struct or empty but instead was %s', class(val));
            structfun(@(x) assert(isnumeric(x), 'Fields in eventTimes must be of numeric types'), val);
            this.i_T = val;
        end
        function val = get.T(this)
            val = this.i_T;
        end
        function this = set.V(this, val)
            this.i_V = val;
        end
        function val = get.V(this)
            val = this.i_V;
        end
        function val = get.isValid(this)
            val = ~isnan(this.i_t);
        end
        
        % Function Overloading
        function val = double(this)
            val = arrayfun(@(x) x.t, this);
        end
        function [this, ind] = sort(this)
            [~, ind] = sort(double(this));
            this = this(ind);
        end
        function val = plus(obj1, obj2)
            if isa(obj2, 'double')
                val = obj1.IOffsetTime(obj2);
            else
                val = obj2.IOffsetTime(obj1);
            end
        end
        function val = minus(obj1, obj2)
            if isa(obj2, 'double')
                val = obj1.IOffsetTime(-obj2);
            else
                val = obj2.IOffsetTime(-obj1);
            end
        end
        function val = times(obj1, obj2)
            if isa(obj2, 'double')
                val = obj1.IScaleTime(obj2);
            else
                val = obj2.IScaleTime(obj1);
            end
        end
        function val = mtimes(obj1, obj2)
            val = times(obj1, obj2);
        end
        function val = lt(obj1, obj2)
            val = double(obj1) < double(obj2);
        end
        function val = gt(obj1, obj2)
            val = double(obj1) > double(obj2);
        end
        function val = le(obj1, obj2)
            val = double(obj1) <= double(obj2);
        end
        function val = ge(obj1, obj2)
            val = double(obj1) >= double(obj2);
        end
        function val = eq(obj1, obj2)
            val = double(obj1) == double(obj2);
        end
        function val = ne(obj1, obj2)
            val = double(obj1) ~= double(obj2);
        end
    end
    
    methods(Access = private)
        function this = IOffsetTime(this, val)
            for i = 1 : numel(this)
                this(i).t = this(i).t + val;
                this(i).i_T = structfun(@(x) x+val, this(i).i_T);
            end
        end
        function this = IScaleTime(this, val)
            for i = 1 : numel(this)
                this(i).t = this(i).t * val;
                this(i).i_T = structfun(@(x) x*val, this(i).i_T);
            end
        end
    end
end

