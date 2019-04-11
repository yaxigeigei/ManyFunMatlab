classdef MMath
    % MMath is a collection of functions useful for doing math and manipulating data
    
    methods(Static)
        function [ci, btDist] = BootCI(nboot, bootfun, vect, varargin)
            %Bootstrap confidence interval similar to the built-in bootci function but faster
            %
            %   [ci, btDist] = BootCI(nboot, bootfun, vect)
            %   [ci, btDist] = BootCI(nboot, bootfun, vect, 'alpha', 0.05)
            %
            % Inputs:
            %   nboot           Number of resampling.
            %   bootfun         Function handle for computation.
            %   vect            Data to be sent into bootfun.
            % Output:
            %   ci              Confidence interval.
            %   btDist          The distribution of computed values from bootstrap.
            
            % Handles user inputs
            p = inputParser();
            p.addParameter('alpha', 0.05, @(x) isscalar(x) && x<1 && x>0);
            p.parse(varargin{:});
            a = p.Results.alpha * 100;
            
            btDist = bootstrp(nboot, bootfun, vect);
            
            ci(1,1) = prctile(btDist, a/2);
            ci(2,1) = prctile(btDist, 100-a/2);
        end
        
        function val = Bound(val, valRange)
            %Bounds input values to specified range or allowed members
            %
            %   val = MMath.Bound(val, valRange)
            %
            % Inputs:
            %   val         Numeric scalar, vector or array
            %   valRange    Either (1) a tuple of boundaries (e.g. [ minVal, maxVal ])
            %               or (2) a numeric vector or array of all allowed values (number of elements must > 2) and
            %               each raw input value will be bound to the closest allowed value.
            % Output:
            %   val         Same size as the input val with bound values
            
            if numel(valRange) > 2
                for i = 1 : numel(val)
                    dists = abs(valRange(:) - val(i));
                    [ ~, idx ] = min(dists);
                    val(i) = valRange(idx);
                end
            else
                val = max(val, min(valRange));
                val = min(val, max(valRange));
            end
        end
        
        function [ condDist, margDist ] = ConditionalDist(jointDist, givenWhich)
            %Compute conditional probability distribution from joint distribution
            %
            %   [ condDist, margDist ] = MMath.ConditionalDist(jointDist, givenWhich)
            %
            % Inputs:
            %   jointDist       multidimentional (>=2) array of joint probability
            %   givenWhich      dimentions to condition on (currently has to be ndims(jointDist)-1 dimensions)
            % Output:
            %   condDist        array of conditional distribution
            %   margDist        marginal distribution
            
            % Collapsing other dimensions
            numDims = ndims(jointDist);
            margWhich = setdiff(1:numDims, givenWhich);
            if isempty(margWhich)
                margDist = jointDist;
            else
                for i = 1 : length(margWhich)
                    margDist = nansum(jointDist, margWhich(i));
                end
            end
            
            % Normalization
            margSize = ones(1, numDims);
            margSize(1:ndims(margDist)) = size(margDist);
            condDist = jointDist ./ repmat(margDist, size(jointDist)./margSize);
            
            % Set NaN to zero
            condDist(isnan(condDist)) = 0;
        end
        
        function [x, D] = Decimate(x, r, n, D)
            %Similar to the MATLAB decimate function but with the following modifications
            % 1) Accepts data types other than double
            % 2) If x is an matrix, treats each column as a time series
            % 3) To be consistent with MATLAB downsample function, it uses the (n+1)-th 
            % 	 value in each bin of r samples
            % 4) By default, it uses 10th-order elliptic (IIR) lowpass filter rather than 
            %    cheby1 or FIRs. 
            %    - fpass is fs/r/2 (i.e. Nyquist freq. of the output), fstop is ~1.2*fpass
            %    - passband has ~1/1000th uneveness (ripple) in amplitude (or 0.01dB)
            %    - stopband is attenuated >10000 times in amplitude (or 80dB)
            % 5) User can provide custom digitalFilter object
            %   
            %   [x, D] = MMath.Decimate(x, r)
            %   [x, D] = MMath.Decimate(x, r, n)
            %   [x, D] = MMath.Decimate(x, r, n, D)
            
            if nargin < 4
                D = designfilt('lowpassiir', ...
                    'FilterOrder', 10, ...
                    'PassbandFrequency', 1/r, ...
                    'PassbandRipple', .01, ...
                    'StopbandAttenuation', 80);
            end
            if nargin < 3
                n = 0;
            end
            
            isRow = isrow(x);
            if isRow
                x = x';
            end
            
            % Filter one column at a time to save memory and conserve data type
            for i = 1 : size(x,2)
                x(:,i) = filtfilt(D, double(x(:,i)));
            end
            
            x = downsample(x, r, n);
            
            if isRow
                x = x';
            end
        end
        
        function H = Entropy(probDist)
            %Calculates the Shannon's entropy of a given probability distribution
            %
            %   H = MMath.Entropy(probDist)
            %
            % Inputs:
            %   probDist    A numeric array of probability values
            % Output:
            %   H           Entropy in bits
            
            probDist = probDist(:);
            probDist = probDist / nansum(probDist);
            H = -nansum(probDist .* log2(probDist));
        end
        
        function expect = Expectation(distMat, val, idxDim)
            %Calculate expected values in the dimension of interest given other variables
            %
            %   expect = MMath.Expectation(distMat, val, idxDim)
            %
            % Inputs:
            %   distMat     Joint probability distribution
            %   val         Values of each probabilities in the dimension of interest
            %   idxDim      Index of the dimension of interest
            % Output:
            %   expect      Expected values
            
            if nargin < 3
                idxDim = 1;
            end
            
            sizeDist = size(distMat);
            reshapeVect = ones(1, length(sizeDist));
            reshapeVect(idxDim) = length(val);
            val = reshape(val, reshapeVect);
            
            sizeDist(idxDim) = 1;
            expect = nansum(distMat .* repmat(val, sizeDist), idxDim);
        end
        
        function roiInd = Ind2Roi(ind, winRoi, valRange)
            %Converts each index to indices for indexing a corresponding region of interest (ROI)
            %e.g. indices [ 23; 56 ] with a window of [ -1 0 1 2 ] => ROIs of [ 22 23 24 25; 55 56 57 58 ]
            %
            %   roiInd = MMath.Ind2Roi(ind, winRoi, valRange)
            %
            % Inputs:
            %   ind         An index or a vector of indices
            %   winRoi      A vector of relative indices (not a boundary tuple) indicating the region of interest
            %   valRange    Usually (1) a tuple of indexing boundaries (e.g. [ minIdx, maxIdx ])
            %               or (2) a numeric vector or array of all allowed index values (number of elements must > 2) and
            %               each index will be bound to the closest allowed value.
            % Output:
            %   roiInd      An array where each row is a set of ROI indices for one input index
            %               ROIs causing edge issue are removed from the reuslt.
            
            if islogical(ind) || ~isempty(find(ind == 0, 1))
                if nargin < 3
                    valRange = [ 1 length(ind) ];
                end
                ind = find(ind);
            end
            
            roiInd = zeros(length(ind), length(winRoi));
            for i = 1 : length(winRoi)
                roiInd(:,i) = ind + winRoi(i);
            end
            
            roiInd = MMath.Bound(roiInd, valRange);
            
            validRoiMask = true(size(roiInd,1), 1);
            roiLength = size(roiInd, 2);
            for i = 1 : size(roiInd, 1)
                if unique(roiInd(i,:)) < roiLength
                    validRoiMask(i) = false;
                end
            end
            roiInd = roiInd(validRoiMask, :);
        end
        
        function mask = Ind2Logical(ind, vecLength)
            % Convert a vector of indices to logical mask
            %
            %   mask = MMath.Ind2Logical(ind, vecLength)
            %
            % Inputs:
            %   ind         A vector of indices
            %   vecLength   The length of the output vector
            % Output:
            %   mask        A logical vector where elements at 'ind' are true
            
            if isrow(ind)
                mask = false(1, vecLength);
            else
                mask = false(vecLength, 1);
            end
            
            mask(ind) = true;
        end
        
        function [ jointDist, axisLabels, binCoor ] = JointDist(randVars, numBins, axisLabels)
            % Calculates the joint probability distribution of multiple random variables
            %
            %   [ jointDist, axisLabels, binInd ] = MMath.JointDist(randVars, numBins)
            %
            % Inputs:
            %   randVars	   	Random variables in an array of column vectors. Each row is an observation.
            %   numBins         The number(s) of bins for input random variables. Use NaN to indicate categorical
            %                   variable. You may assign different bin numbers to different variables in a
            %                   vector otherwise they all use the same bin number.
            %   axisLabels      For categorical variables, you can provide predefined categories (e.g. to account
            %                   for unobserved values). It expects a 1-by-n cell array where n is the number of
            %                   variables in randVars. Each element in this cell array is a vector of numeric
            %                   categories along that dimension.
            % Output:
            %   jointDist       N-dimensional matrix of joint probability. The order of dimensions is the same
            %                   as the order of input variables.
            %   axisLabels      Labels of each axis in cell array.
            %   binCoor         Each row is a coordinate in the joint distribution sapce for the corresponding
            %                   observation (row) in randVars.
            
            % Handles user inputs
            if nargin < 3
                axisLabels = cell(1, size(randVars,2));
            end
            
            if length(numBins) ~= size(randVars,2)
                numBins = repmat(numBins(1), 1, size(randVars,2));
            end
            
            % Bins variables
            binCoor = zeros(size(randVars));
            for i = size(randVars,2) : -1 : 1
                if isnan(numBins(i))
                    % Catagorical
                    axisLabels{i} = union(axisLabels{i}, unique(randVars(:,i)));
                    binCoor(:,i) = arrayfun(@(x) find(axisLabels{i} == x), randVars(:,i));
                    numBins(i) = length(axisLabels{i});
                else
                    % Continuous
                    edges = linspace(nanmin(randVars(:,i)), nanmax(randVars(:,i)), numBins(i)+1);
                    axisLabels{i} = mean([ edges(1:end-1); edges(2:end) ])';
                    edges(end) = edges(end) + 1;
                    [ ~, binCoor(:,i) ] = histc(randVars(:,i), edges);
                end
            end
            
            % Multidimensional tally
            dimInputStr = 'ndgrid(';
            for i = 1 : length(numBins)
                dimInputStr = [ dimInputStr, '1:numBins(', num2str(i), '), ' ];
            end
            dimInputStr = [ dimInputStr(1:end-2), ')' ];
            
            dimOutput = cell(1, length(numBins));
            [ dimOutput{:} ] = eval(dimInputStr);
            dimOutput = cell2mat(cellfun(@(x) x(:), dimOutput, 'UniformOutput', false));
            
            jointDist = zeros(numBins(:)');
            for i = size(dimOutput,1) : -1 : 1
                subBinCoor = binCoor;
                for j = 1 : length(numBins)
                    hitVect = subBinCoor(:,j) == dimOutput(i,j);
                    subBinCoor = subBinCoor(hitVect, :);
                end
                jointDist(i) = sum(hitVect);
            end
            
            jointDist = jointDist / size(randVars,1);
        end
        
        function [H, pValue, KSstatistic] = KStest2CDF(ecdf1, ecdf2, varargin)
            %KSTEST2 Two-sample Kolmogorov-Smirnov goodness-of-fit hypothesis test.
            
            if nargin > 2
                [varargin{:}] = convertStringsToChars(varargin{:});
            end
            
            if nargin < 2
                error(message('stats:kstest2:TooFewInputs'));
            end
            
            % Parse optional inputs
            alpha = []; tail = [];
            if nargin >=3
                if isnumeric(varargin{1})
                    % Old syntax
                    alpha = varargin{1};
                    if nargin == 4
                        tail = varargin{2};
                    end
                else
                    % New syntax
                    params = {'alpha', 'tail'};
                    dflts =  { []     , []};
                    [alpha, tail] =...
                        internal.stats.parseArgs(params, dflts, varargin{:});
                end
            end
            
            % Ensure each sample is a VECTOR.
            ecdf1  =  ecdf1(:);
            ecdf2  =  ecdf2(:);
            if isempty(ecdf1)
                error(message('stats:kstest2:NotEnoughData', 'X1'));
            end
            if isempty(ecdf2)
                error(message('stats:kstest2:NotEnoughData', 'X2'));
            end
            
            % Ensure the significance level, ALPHA, is a scalar
            % between 0 and 1 and set default if necessary.
            
            if ~isempty(alpha)
                if ~isscalar(alpha) || (alpha <= 0 || alpha >= 1)
                    error(message('stats:kstest2:BadAlpha'));
                end
            else
                alpha  =  0.05;
            end
            
            % Ensure the type-of-test indicator, TAIL, is a string or scalar integer
            % from the allowable set, and set default if necessary.
            if ~isempty(tail)
                if ischar(tail)
                    try
                        [~,tail] = internal.stats.getParamVal(tail, ...
                            {'smaller','unequal','larger'},'Tail');
                    catch
                        error(message('stats:kstest2:BadTail'));
                    end
                    tail = tail - 2;
                elseif ~isscalar(tail) || ~((tail==-1) || (tail==0) || (tail==1))
                    error(message('stats:kstest2:BadTail'));
                end
            else
                tail  =  0;
            end
            
            % Compute the test statistic of interest.
            switch tail
                case  0      %  2-sided test: T = max|F1(x) - F2(x)|.
                    deltaCDF  =  abs(ecdf1 - ecdf2);
                case -1      %  1-sided test: T = max[F2(x) - F1(x)].
                    deltaCDF  =  ecdf2 - ecdf1;
                case  1      %  1-sided test: T = max[F1(x) - F2(x)].
                    deltaCDF  =  ecdf1 - ecdf2;
            end
            KSstatistic   =  max(deltaCDF);
            
            % Compute the asymptotic P-value approximation and accept or
            % reject the null hypothesis on the basis of the P-value.
            n1     =  length(ecdf1);
            n2     =  length(ecdf2);
            n      =  n1 * n2 /(n1 + n2);
            lambda =  max((sqrt(n) + 0.12 + 0.11/sqrt(n)) * KSstatistic , 0);
            if tail ~= 0        % 1-sided test.
                pValue  =  exp(-2 * lambda * lambda);
            else                % 2-sided test (default).
                %  Use the asymptotic Q-function to approximate the 2-sided P-value.
                j       =  (1:101)';
                pValue  =  2 * sum((-1).^(j-1).*exp(-2*lambda*lambda*j.^2));
                pValue  =  min(max(pValue, 0), 1);
            end
            
            H  =  (alpha >= pValue);
        end
        
        function boundariesInd = Logical2Bounds(vect)
            %Converts logical(boolean/binary) vector to tuples indicating the bondaries of 1 regions
            %e.g. [ 0 0 0 1 1 1 0 0 1 ] => bondaries [ 4 6; 9 9 ]
            %
            %   bondariesInd = MMath.Logical2Bounds(vect)
            %
            % Inputs:
            %   vect            A logical(boolean/binary) vector
            % Output:
            %   bondariesInd    Tuples (each row) indicating the bondaries of 1 regions
            
            if isempty(vect)
                boundariesInd = zeros(0,2);
                return;
            end
            
            vect = logical(vect(:));
            vect = [ zeros(1,size(vect,2)); vect; zeros(1,size(vect,2)) ];
            dVect = diff(vect);
            boundariesInd = [ find(dVect == 1), find(dVect == -1) - 1 ];
        end
        
        function [m, sd, se, ci] = MeanStats(A, dim, varargin)
            % Compute means, SDs, SEMs and bootstrap CIs from samples
            %
            %   [m, sd, se, ci] = MMath.MeanStats(A)
            %   [m, sd, se, ci] = MMath.MeanStats(A, dim)
            %   [m, sd, se, ci] = MMath.MeanStats(A, dim, isoutlierArg, ...)
            %
            % Inputs
            %   A           Numeric array of samples.
            %   dim         Dimension to operate along.
            %   isoutlierArg, ...
            %               One or more arguments for isoutlier function to remove outliers in A.
            % Outputs
            %   All otuputs has the same dimensionality as A with statistical values in the dim 
            %   dimension. 
            %   m           Mean values.
            %   sd          Standard deviations.
            %   se          Standard error of the mean.
            %   ci          95% bootstrap confidence intervals.
            
            if isempty(A)
                A = NaN;
            end
            
            if nargin < 2
                if isscalar(A)
                    dim = 1;
                else
                    dim = find(size(A) > 1, 1);
                end
            end
            
            if numel(varargin) > 0
                A(isoutlier(A, varargin{:}, dim)) = [];
            end
            
            m = nanmean(A, dim);
            sd = nanstd(A, 0, dim);
            se = sd ./ sqrt(size(A,dim));
            
            if size(A,dim) < 2
                ci = cat(dim, m, m);
            else
                % bootci can only sample along the first dimension, thus permuting A
                dimOrder = [dim setdiff(1:ndims(A), dim)];
                A = permute(A, dimOrder);
                ci = bootci(1000, @nanmean, A);
                
                % Restore original dimension order
                ci = permute(ci, [1 3:ndims(ci) 2]); % squeeze out the second (mean value) dimension
                ci = ipermute(ci, dimOrder);
            end
        end
        
        function [m, qt, ad] = MedianStats(A, dim)
            % Compute medians, 1st and 3rd quartiles and absolute deviations from samples
            %
            %   [m, qt, ad] = MMath.MedianStats(A)
            %   [m, qt, ad] = MMath.MedianStats(A, dim)
            %   [m, qt, ad] = MMath.MedianStats(A, dim, isoutlierArg, ...)
            %
            % Inputs
            %   A           Numeric array of samples.
            %   dim         Dimension to operate along. Default is the first non-singleton dimension. 
            % Outputs
            %   m           Median values.
            %   qt          1st and 3rd quartiles.
            %   ad          Median absolute deviation.
            
            if isempty(A)
                A = NaN;
            end
            
            if nargin < 2
                if isscalar(A)
                    dim = 1;
                else
                    dim = find(size(A) > 1, 1);
                end
            end
            
            m = nanmedian(A, dim);
            qt = prctile(A, [25 75], dim);
            ad = mad(A, 1, dim);
        end
        
        function I = MutualInfo(jointDist)
            %Calculates the mutual information of two random variables from their joint probability distribution
            %
            %   I = MMath.MutualInfo(jointDist)
            %
            % Inputs:
            %   jointDist   The joint probability distribution matrix of two random variables
            % Output:
            %   I           Mutual information in bits
            
            marginDist{1} = nansum(jointDist, 2);
            marginDist{2} = nansum(jointDist, 1);
            [ marginDist{2}, marginDist{1} ] = meshgrid(marginDist{2}, marginDist{1});
            I = jointDist .* log2(jointDist ./ (marginDist{1}.*marginDist{2}));
            I = nansum(I(:));
        end
        
        function [corrVal, pVal] = NanCorr(x, y)
            % Similar to MATLAB built-in function corr but ignores NaN values
            
            x = x(:);
            y = y(:);
            indNan = isnan(x) | isnan(y);
            [corrVal, pVal] = corr(x(~indNan), y(~indNan));
        end
        
        function y = SigMf(x, params)
            %SIGMF Sigmoid curve membership function.
            %   SIGMF(X, PARAMS) returns a matrix which is the sigmoid
            %   membership function evaluated at X. PARAMS is a 2-element vector
            %   that determines the shape and position of this membership function.
            %   Specifically, the formula for this membership function is:
            %
            %   SIGMF(X, [A, C]) = 1./(1 + EXP(-A*(X-C)))
            %
            %   For example:
            %
            %       x = (0:0.2:10)';
            %       y1 = sigmf(x, [-1 5]);
            %       y2 = sigmf(x, [-3 5]);
            %       y3 = sigmf(x, [4 5]);
            %       y4 = sigmf(x, [8 5]);
            %       subplot(211); plot(x, [y1 y2 y3 y4]);
            %       y1 = sigmf(x, [5 2]);
            %       y2 = sigmf(x, [5 4]);
            %       y3 = sigmf(x, [5 6]);
            %       y4 = sigmf(x, [5 8]);
            %       subplot(212); plot(x, [y1 y2 y3 y4]);
            %       set(gcf, 'name', 'sigmf', 'numbertitle', 'off');
            %
            %   See also DSIGMF, EVALMF, GAUSS2MF, GAUSSMF, GBELLMF, MF2MF, PIMF, PSIGMF, SMF,
            %   TRAPMF, TRIMF, ZMF.
            
            %       Roger Jang, 6-29-93, 4-17-93.
            %   Copyright 1994-2002 The MathWorks, Inc.
            
            if nargin ~= 2
                error('Two arguments are required by the sigmoidal MF.');
            elseif length(params) < 2
                error('The sigmoidal MF needs at least two parameters.');
            end
            
            a = params(1); c = params(2);
            y = 1./(1 + exp(-a*(x-c)));
        end
        
        function [se, sd] = StandardError(A, dim)
            % Calculates the standard error of the population mean
            %
            %   [se, sd] = MMath.StandardError(A)
            %   [se, sd] = MMath.StandardError(A, dim)
            %
            % Inputs
            %   A           Numeric array of samples.
            %   dim         Dimension to operate along.
            % Outputs
            %   se          Standard error of the population mean.
            %   sd          Standard deviation of the population mean.
            
            if isempty(A)
                A = NaN;
            end
            
            if nargin < 2
                if isscalar(A)
                    dim = 1;
                else
                    dim = find(size(A) > 1, 1);
                end
            end
            
            sd = nanstd(A, 0, dim);
            se = sd ./ sqrt(size(A,dim));
        end
    end
end

