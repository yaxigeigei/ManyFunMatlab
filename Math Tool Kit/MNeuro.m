classdef MNeuro
    
    methods(Static)
        function ccg = CCG(tEdges, varargin)
            % Compute all pairwise cross-correlograms and auto-correlograms
            %
            %   ccg = MNeuro.CCG(tEdges, spkTrain1, spkTrain2, spkTrain3, ...)
            %
            % Inputs
            %   tEdges          A vector of time bin, two-sided.
            %   spkTrainN       A vector of spike time.
            % Output
            %   ccg             A [numTrain, numTrain, numBin] array of histograms. Each element 
            %                   stores the number of pairwise intervals between two spike trains 
            %                   whose values are in a range defined by the time bin. 
            
            % Handle user inputs
            spkTrains = varargin(:);
            nTrains = numel(spkTrains);
            for i = 1 : nTrains
                spkTrains{i} = spkTrains{i}(:);
                spkTrains{i}(isnan(spkTrains{i})) = [];
            end
            
            % Create labels of train membership for each spike
            spkIds = cellfun(@(x,y) zeros(size(x))+y, spkTrains, num2cell(1:nTrains)', 'Uni', false);
            
            % Combine spike trains
            fullTrain = cell2mat(spkTrains);
            [fullTrain, tOrder] = sort(fullTrain, 'ascend');
            fullId = cell2mat(spkIds);
            fullId = fullId(tOrder);
            
            % Compute CCGs
            ccg = zeros(nTrains, nTrains, numel(tEdges)-1);
            isInWin = true;
            s = 1;
            while any(isInWin)
                % Compute difference of spike times for a given shift
                id1 = fullId(1+s:end);
                id2 = fullId(1:end-s);
                dt = fullTrain(1+s:end) - fullTrain(1:end-s);
                
                % Find time intervals that are in the window of interest
                isInWin = dt <= tEdges(end);
                id1 = id1(isInWin);
                id2 = id2(isInWin);
                dt = dt(isInWin);
                
                % Bin time intervals
                dt = discretize(dt, tEdges);
                
                % Tally occurances
                ind = sub2ind(size(ccg), id1, id2, dt);
                for i = 1 : numel(ind)
                    ccg(ind(i)) = ccg(ind(i)) + 1;
                end
                
                % Increment shift
                s = s + 1;
            end
            
%             % Zero center bins of auto-correlograms (Phy does it but don't know why)
%             for i = 1 : nTrains
%                 ccg(i,i,ceil(numel(tEdges)/2)) = 0;
%             end
            
            % Symmetrize matrices
            ccg = ccg + flip(permute(ccg, [2 1 3]), 3);
        end
        
        function [r, varargout] = Filter1(r, fs, methodOpt, varargin)
            % Filter activity or bin spike times using standard methods
            % 
            % Smooth trace with Gaussian kernel
            %   [r, ker] = MNeuro.Filter1(r, fs, 'gaussian', sigma)
            %   [r, ker] = MNeuro.Filter1(r, fs, 'gaussian', sigma, kerSize)
            % 
            % Smooth trace with exponential kernel
            %   [r, ker] = MNeuro.Filter1(r, fs, 'exponential', tau)
            %   [r, ker] = MNeuro.Filter1(r, fs, 'exponential', tau, kerSize)
            % 
            % Bin spike times
            %   [r, tEdges] = MNeuro.Filter1(t, [], 'bin', tEdges)
            %   [r, tEdges] = MNeuro.Filter1(t, [], 'bin', tWin, binSize)
            % 
            % Inputs
            %   t           A vector of spike time.
            %   r           A vector of spike rate.
            %   fs          Sampling rate of r.
            %   methodOpt   'gaussian', 'exponential', 'bin'.
            %   sigma       Standard deviation of the Gausssian kernel in second.
            %   tau         Time constant of the exponential kernal in second.
            %   kerSize     Length of the kernel in second. The default is 6 sigma for 'gaussian' 
            %               (i.e. 3 sigma each side) or 5 tau for 'exponential'. 
            %   tEdges      A vector of bin edges.
            %   tWin        A 2-element array indicating the start and end of binning.
            %   binSize     Width of time bin.
            % Outputs
            %   r           A vector of spike rate.
            %   ker         The kernel used for convolution.
            %   tEdges      A vector of bin edges.
            %
            % See also: smooth, conv
            
            if isvector(r)
                r = r(:);
            end
            
            switch lower(methodOpt)
                case 'bin'
                    % Bin spike times
                    assert(all(diff(r) > 0), 'Spike times in r must be monotonically increasing');
                    if numel(varargin) == 1
                        tEdges = varargin{1};
                    else
                        tWin = varargin{1};
                        binSize = varargin{2};
                        tEdges = tWin(1) : binSize : tWin(2);
                    end
                    r = histcounts(r, tEdges, 'Normalization', 'countdensity');
                    varargout{1} = tEdges;
                    
                case 'gaussian'
                    % Get parameters
                    sigma = varargin{1} * fs;
                    if numel(varargin) == 1
                        kerSize = 2 * 3 * sigma;
                    else
                        kerSize = varargin{2} * fs;
                    end
                    kerSize = ceil(kerSize);
                    
                    % Construct gaussian kernel
                    a = (kerSize-1) / (2*sigma);
                    ker = gausswin(kerSize, a);
                    ker = ker / sum(ker);
                    
                    % Filtering
                    for i = 1 : size(r,2)
                        r(:,i) = conv(r(:,i), ker, 'same');
                    end
                    varargout{1} = ker;
                    
                case 'exponential'
                    % Get parameters
                    tau = varargin{1} * fs;
                    if numel(varargin) == 1
                        kerSize = tau * 5;
                    else
                        kerSize = varargin{2} * fs;
                    end
                    kerSize = ceil(kerSize);
                    
                    % Construct exponential kernel
                    x = (1 : kerSize)';
                    y = exp(-x/tau);
                    ker = zeros(kerSize*2-1, 1);
                    ker(kerSize:end) = y;
                    ker = ker / sum(ker);
                    
                    % Filtering
                    for i = 1 : size(r,2)
                        r(:,i) = conv(r(:,i), ker, 'same');
                    end
                    varargout{1} = ker;
                    
                otherwise
                    error('%s (case-insensitive) is not a valid method', methodOpt);
            end
        end
        
        function t = JointTuning(varargin)
            % Calculates the joint probability distribution of multiple random variables
            %
            %   t = MNeuro.JointTuning(stimuli)
            %   t = MNeuro.JointTuning(stimuli, response)
            %   t = MNeuro.JointTuning(..., 'mask', logicals)
            %   t = MNeuro.JointTuning(..., 'numBins', numbers)
            %   t = MNeuro.JointTuning(..., 'rmMethod', options)
            %   t = MNeuro.JointTuning(..., 'rmParam', params)
            %   
            % Inputs:
            %   stimuli         Data of stimuli (each column is one stimulus).
            %   response        Array of responses corresponding to the stimuli.
            %   'mask'          Numeric or cell vector of binary indices used to mask selected data. (default is
            %                   no masking)
            %   'numBins'       Number of bins used to group the stimuli range (default is 30). Use NaN to
            %                   specify catagorical variable.
            %   'rmMethod'      A string (or a cell array of strings) specifying the method of outlier removal for 
            %                   all (or individual) stimulus variable (default 'percentile', other methods include
            %                   'std', 'whisker', 'cut'). See help of MMath.RemoveOutliers() for more details. 
            %   'rmParam'       Parameter(s) for outlier removal (default 0.5). You may provide an array of row 
            %                   vectors that specify parameters for individual stimulus variables. See the help
            %                   of MMath.RemoveOutliers() for more details. 
            % Output:
            %   t               A struct of various tuning information
            
            % Handle user input
            p = inputParser();
            p.addRequired('stimuli');
            p.addRequired('response');
            p.addParameter('mask', []);
            p.addParameter('numBins', 30, @isnumeric);
            p.addParameter('rmParam', .5, @isnumeric);
            p.addParameter('rmMethod', 'percentile');
            p.addParameter('tfParam', NaN, @isnumeric);
            p.parse(varargin{:});
            stim = MNeural.Vars2Mat(p.Results.stimuli);
            resp = MNeural.Vars2Mat(p.Results.response);
            dataMask = p.Results.mask;
            numBins = p.Results.numBins;
            rmParam = p.Results.rmParam;
            rmMethod = cellstr(p.Results.rmMethod);
            tfParam = p.Results.tfParam;
            
            % Propagate settings
            if length(numBins) < size(stim,2)
                numBins(end+1:size(stim,2),1) = numBins(1);
            end
            if length(rmMethod) < size(stim,2)
                rmMethod = repmat(rmMethod(1), size(stim,2), 1);
            end
            if size(rmParam,1) < size(stim,2)
                rmParam = repmat(rmParam(1,:), size(stim,2), 1);
            end
            if ~isempty(tfParam) && size(tfParam,1) < size(stim,2)
                tfParam = repmat(tfParam(1), size(stim,2), 1);
            end
            
            % Data masking
            if ~isempty(dataMask)
                if iscell(dataMask)
                    dataMask = cell2mat(dataMask);
                end
                dataMask = logical(dataMask);
                stim = stim(dataMask, :);
                resp = resp(dataMask, :);
            end
            
            % Tansform data
            for i = size(stim,2) : -1 : 1
                if isnan(tfParam(i))
                    f{i} = @(x) x;
                    fi{i} = @(y) y;
                else
                    k = 1 / prctile(stim(:,i), tfParam(i));
                    f{i} = @(x) 2./(1 + exp(-k*x)) - 1;
                    fi{i} = @(y) -1/k * (log(1-y) - log(1+y));
                    stim(:,i) = f{i}(stim(:,i));
                end
            end
            
            % Remove outliers
            indClean = true(size(stim,1), 1);
            for i = 1 : size(stim,2)
                [ ~, varMaskKeep ] = MMath.RemoveOutliers(stim(:,i), rmParam(i,:), rmMethod{i});
                indClean = indClean & varMaskKeep;
            end
            
            % Bin variables
            randVars = [stim, resp];
            randVars = randVars(indClean, :);
            if length(numBins) <= size(stim,2)
                numBins(size(stim,2)+1) = NaN;
            end
            [t.jointProb, t.binCenters, binCoor] = MMath.JointDist(randVars, numBins);
            
            t.condProb = MMath.ConditionalDist(t.jointProb, 1:size(stim,2));                % P(resp|stim)
            t.condExp = MMath.Expectation(t.condProb, t.binCenters{end}, size(stim,2)+1);   % E(resp|stim)
            
            binCoor = num2cell(binCoor, 1);
            indSample = sub2ind(size(t.condExp), binCoor{1:end-1});
            t.sampleVal = arrayfun(@(x) randVars(indSample==x,end), (1:numel(t.condExp))', 'Uni', false);
            t.sampleVal = reshape(t.sampleVal, size(t.condExp));
            t.sampleSE = cellfun(@MMath.StandardError, t.sampleVal);
            
            t.margProb = sum(t.jointProb, ndims(t.jointProb));                              % P(stim)
            t.sampleNum = round(t.margProb * sum(indClean));
            t.sampleLog = log10(t.sampleNum);
            t.sampleLog(t.sampleLog == -Inf) = NaN;
            
            t.condExpClean = t.condExp;
            t.condExpClean(t.sampleNum < 25) = NaN;
            
            % Transform back
            for i = size(stim,2) : -1 : 1
                t.binCenters{i} = fi{i}(t.binCenters{i});
            end
        end
        
        function [mm, ee, stats] = MeanEventRate(T, edges)
            % Compute mean event rates and related stats from event times across repetitions
            %
            %   [mm, ee, stats] = MNeuro.MeanEventRate(T, edges)
            %
            % Inputs
            %   T               A cell array or a table of event time vectors. Rows are repeats 
            %                   (e.g. trials); columns are different types of event (e.g. units). 
            %   edges           Edges of time bins in a numeric vector. 
            % Outputs
            %   mm              An array of mean event rates. Rows are time bins and columns are 
            %                   different types of event. 
            %   ee              Standard error of elements in mm. 
            %   stats           A table with the following variables.
            %     colNum          Column index of each event in T.
            %     pkIdx           Index of the time bin where each trace in mm peaks.
            %     pkVal           The value of event rate at pkIdx.
            %     pkProb          Probability of observing any event across repeats at pkIdx.
            %     AUC             Area under the curve (i.e. the integral of each trace in mm). 
            %     entropy         Shannon's entropy of each trace in mm.
            
            % Check and standardize inputs
            if istable(T)
                T = table2cell(T);
            end
            assert(iscell(T), 'T must be a cell array or table.');
            
            for i = size(T,2) : -1 : 1
                % Compute histogram for each repetition
                hh = cell(size(T,1), 1);
                for j = 1 : numel(hh)
                    hh{j} = histcounts(T{j,i}, edges, 'Normalization', 'countdensity');
                end
                hh = cell2mat(hh);
                
                % Compute mean and sem of event rate
                m = mean(hh);
                e = MMath.StandardError(hh);
                mm(:,i) = m;
                ee(:,i) = e;
                
                % Compute other stats
                [pkRate(i), pkBin(i)] = max(m);
                pkProb(i) = mean(hh(:,pkBin(i)) > 0);
                I(i) = MMath.Entropy(m/sum(m));
            end
            
            stats = table();
            stats.colNum = (1:size(T,2))';
            stats.pkIdx = pkBin';
            stats.pkVal = pkRate';
            stats.pkProb = pkProb';
            stats.AUC = sum(mm)';
            stats.entropy = I';
        end
        
        function [mm, ee, stats] = MeanTimeSeries(S)
            % Compute mean and related stats of multiple time series across repetitions
            %
            %   [mm, ee, stats] = MNeuro.MeanTimeSeries(S)
            %
            % Inputs
            %   S               A cell array or a table of time series vectors. Rows are repeats 
            %                   (e.g. trials); columns are different signals (e.g. neurons). 
            % Outputs
            %   mm              An array of mean time series. Rows are samples and columns are for 
            %                   different signals. 
            %   ee              Standard error of elements in mm. 
            %   stats           A table with the following variables.
            %     colNum          Column index of each signal in T.
            %     pkIdx           Index of the time bin where each trace in mm peaks.
            %     pkVal           The value at pkIdx.
            %     AUC             Area under the curve (i.e. the integral of each trace in mm). 
            %     entropy         Shannon's entropy of each trace in mm.
            
            % Check and standardize inputs
            if isnumeric(S)
                S = {S};
            end
            if istable(S)
                S = table2cell(S);
            end
            assert(iscell(S), 'S must be a numeric vector, a cell array or a table.');
            for i = 1 : numel(S)
                assert(isnumeric(S{i}) && isvector(S{i}), 'Individual series must be numeric vector.')
                if iscolumn(S{i})
                    S{i} = S{i}';
                end
            end
            
            for i = size(S,2) : -1 : 1
                % Concatenate rows
                s = cell2mat(S(:,i));
                
                % Compute mean and sem
                m = mean(s);
                e = MMath.StandardError(s);
                mm(:,i) = m;
                ee(:,i) = e;
                
                % Compute other stats
                [pkVal(i), pkIdx(i)] = max(m);
                I(i) = MMath.Entropy(m/sum(m));
            end
            
            stats = table();
            stats.colNum = (1:size(S,2))';
            stats.pkIdx = pkIdx';
            stats.pkVal = pkVal';
            stats.AUC = sum(mm)';
            stats.entropy = I';
        end
        
        function tunings = Tuning(varargin)
            %Computes tuning curves (and errors) for given stimuli
            %
            %   tunings = MNeuro.Tuning()
            %   tunings = MNeuro.Tuning(stimuli)
            %   tunings = MNeuro.Tuning(stimuli, response)
            %   tunings = MNeuro.Tuning(..., 'mask', logicals)
            %   tunings = MNeuro.Tuning(..., 'numBins', value)
            %   tunings = MNeuro.Tuning(..., 'rmMethod', options)
            %   tunings = MNeuro.Tuning(..., 'rmParam', params)
            %   
            % Inputs:
            %   stimuli         Data of stimuli (each column is one stimulus)
            %   response        Array of spike rates corresponding to the stimuli.
            %   'mask'          Numeric or cell array of binary indices used to mask data. (default is no masking)
            %   'numBins'       Number of bins used to group the stimuli range (default is 50)
            %   'rmMethod'      A string (or a cell array of strings) specifying the method of outlier removal for 
            %                   all (or individual) stimulus variable (default 'percentile', other methods include
            %                   'std', 'whisker', 'cut'). See help of MMath.RemoveOutliers() for more details. 
            %   'rmParam'       Parameter(s) for outlier removal (default NaN, no outlier removal). You may provide an array of row 
            %                   vectors that specify parameters for individual stimulus variables. See the help
            %                   of MMath.RemoveOutliers() for more details. 
            % Output:
            %   tunings         Cell array of tuning curves, each contains 4 column vectors - bin centers, 
            %                   mean firing rate, standard error of the mean, and number of samples - from 1 to 4.
            
            % Handles user inputs
            p = inputParser();
            p.addRequired('stimuli');
            p.addRequired('response');
            p.addParameter('mask', []);
            p.addParameter('numBins', 50, @isscalar);
            p.addParameter('rmMethod', 'percentile');
            p.addParameter('rmParam', NaN, @isnumeric);
            p.addParameter('tfParam', NaN, @isnumeric);
            p.parse(varargin{:});
            stim = p.Results.stimuli;
            resp = p.Results.response;
            dataMask = p.Results.mask;
            numBins = p.Results.numBins;
            rmParam = p.Results.rmParam;
            rmMethod = cellstr(p.Results.rmMethod);
            tfParam = p.Results.tfParam;
            
            if size(stim,1) ~= length(resp)
                error([ 'The length of stimuli does not match with that of the spike rates. ', ...
                    'If you are not using the full length of stimuli, ', ...
                    'the corresponding array of spike rates needs to be provided.' ]);
            end
            
            % Propagates settings
            if length(numBins) < size(stim,2)
                numBins(end+1:size(stim,2),1) = numBins(1);
            end
            if size(rmParam,1) < size(stim,2)
                rmParam = repmat(rmParam(1,:), size(stim,2), 1);
            end
            if length(rmMethod) < size(stim,2)
                rmMethod = repmat(rmMethod(1), size(stim,2), 1);
            end
            if ~isempty(tfParam) && size(tfParam,1) < size(stim,2)
                tfParam = repmat(tfParam(1), size(stim,2), 1);
            end
            
            % Data masking
            if ~isempty(dataMask)
                if iscell(dataMask)
                    dataMask = cell2mat(dataMask);
                end
                if size(stim,1) ~= length(dataMask)
                    error('The length of stimuli does not match with that of the mask.');
                end
                dataMask = logical(dataMask);
                stim = stim(dataMask, :);
                resp = resp(dataMask);
            end
            
            % Tansform data
            for i = size(stim,2) : -1 : 1
                if isnan(tfParam(i))
                    f{i} = @(x) x;
                    fi{i} = @(y) y;
                else
                    k = 1 / prctile(stim(:,i), tfParam(i));
                    f{i} = @(x) 2./(1 + exp(-k*x)) - 1;
                    fi{i} = @(y) -1/k * (log(1-y) - log(1+y));
                    stim(:,i) = f{i}(stim(:,i));
                end
            end
            
            % Computes for each stimulus
            for i = size(stim,2) : -1 : 1
                % Removes outliers in the stimulus vector and corresponding entries in response vector
                [ stimClean, indKeep ] = MMath.RemoveOutliers(stim(:,i), rmParam(i,:), rmMethod{i});
                respClean = resp(indKeep);
                
                % Bining response by stimulus value ranges
                edges = linspace(min(stimClean), max(stimClean), numBins(i)+1);
                centers = mean([ edges(1:end-1); edges(2:end) ]);
                [ ~, binInd ] = histc(stimClean, edges);
                for j = length(edges)-1 : -1 : 1
                    ratesVect = respClean(binInd == j);
                    numSample = numel(ratesVect);
                    if numSample == 0
                        meanSpikeRate = 0;
                        steSpikeRate = 0;
                    else
                        meanSpikeRate = nanmean(ratesVect);
                        steSpikeRate = MMath.StandardError(ratesVect);
                    end
                    binNumSample(j) = numSample;
                    binMeanSpikeRate(j) = meanSpikeRate;
                    binSteSpikeRate(j) = steSpikeRate;
                end
                
                % Transforms axes back
                centers = fi{i}(centers);
                
                tunings{i,1} = [ centers; binMeanSpikeRate; binSteSpikeRate; binNumSample ]';
            end
        end
        
        function c = ClusterContamination(p, F, tauR, tauC)
            % Estimate the contamination of a unit based on Hill et al. 2011
            % 
            %   c = MNeuro.ClusterContamination(p, F, tauR, tauC)
            % 
            % Inputs
            %   p       Fraction of ISI violation.
            %   F       Mean firing rate (spikes per second) of the unit.
            %   tauR    Length of the refractory period in second. The default is 0.0025s.
            %   tauC    The censored period (in second) following a spike during which spikes are 
            %           not detected by the recording system. The default is 0s (e.g. Kilosort).
            % Output
            %   c       Fraction of contaminated spikes. Note that the maximal possible value is 
            %           0.5.
            
            % Compute time window constant
            if nargin < 4
                tauC = 0;
            end
            if nargin < 3
                tauR = 0.0025;
            end
            T = 2 * (tauR - tauC);
            
            % Compute term inside the square root
            x = (4*p)./(F*T);
            x(x>1) = 1;
            
            % Compute contamination rate
            c = (1 - sqrt(1 - x)) / 2;
        end
    end
end


