classdef Event
    %Event Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Dependent)
        t;
        T;
    end
    properties
        V;
    end
    properties(Access = private)
        t_ = NaN;
        T_ = struct();
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
            this.t_ = val;
        end
        function val = get.t(this)
            val = this.t_;
        end
        function this = set.T(this, val)
            if isempty(val)
                val = struct();
            end
            assert(isstruct(val), 'T must be a struct or empty array but instead was %s', class(val));
            isNum = structfun(@isnumeric, val);
            assert(all(isNum), 'Fields of T must be numeric array');
            this.T_ = val;
        end
        function val = get.T(this)
            val = this.T_;
        end
        function val = Tfield(this, fieldName, varargin)
            val = arrayfun(@(x) x.T_.(fieldName), this, varargin{:});
        end
        function val = Vfield(this, fieldName, varargin)
            val = arrayfun(@(x) x.V.(fieldName), this, varargin{:});
        end
        
        % Function Overloading
        function val = double(this)
            val = double(this.IGet_t());
        end
        function val = single(this)
            val = single(this.IGet_t());
        end
        function val = isnan(this)
            val = isnan(this.IGet_t());
        end
        function [this, ind] = sort(this, varargin)
            [~, ind] = sort(this.IGet_t(), varargin{:});
            this = this(ind);
        end
        
        % Operator Overloading
        function obj = plus(obj1, obj2)
            if isa(obj2, 'double')
                obj = IOffsetTime(obj1, obj2);
            else
                obj = IOffsetTime(obj2, obj1);
            end
        end
        function obj = minus(obj1, obj2)
            if isa(obj2, 'double')
                obj = IOffsetTime(obj1, -obj2);
            else
                obj = IOffsetTime(obj2, -obj1);
            end
        end
        function obj = times(obj1, obj2)
            if isa(obj2, 'double')
                obj = IScaleTime(obj1, obj2);
            else
                obj = IScaleTime(obj2, obj1);
            end
        end
        function obj = mtimes(obj1, obj2)
            obj = times(obj1, obj2);
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
        % Utilities
        function val = IGet_t(this)
            val = zeros(size(this));
            for i = 1 : numel(this)
                val(i) = this(i).t_;
            end
        end
        function this = IOffsetTime(this, val)
            for i = 1 : numel(this)
                this(i).t_ = this(i).t_ + val;
                this(i).T_ = structfun(@(x) x + val, this(i).T_, 'Uni', false);
            end
        end
        function this = IScaleTime(this, val)
            for i = 1 : numel(this)
                this(i).t_ = this(i).t * val;
                this(i).T_ = structfun(@(x) x * val, this(i).T_, 'Uni', false);
            end
        end
    end
end

