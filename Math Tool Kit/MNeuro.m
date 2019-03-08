classdef MNeuro
    
    methods(Static)
        function ccg = CCG(tEdges, varargin)
            % Compute all pairwise cross-correlograms and auto-correlograms
            %
            %   ccg = CCG(tEdges, spkTrain1, spkTrain2, spkTrain3, ...)
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
            
            % Zero center bins of auto-correlograms (Phy does it but don't know why)
            for i = 1 : nTrains
                ccg(i,i,ceil(numel(tEdges)/2)) = 0;
            end
            
            % Symmetrize matrices
            ccg = ccg + flip(permute(ccg, [2 1 3]), 3);
        end
        
        function [r, varargout] = FilterSpikeRate(r, fs, methodOpt, varargin)
            % Filter spike rate or bin spike times using standard methods
            % 
            % Bin spike times
            %   [r, tEdges] = MNeural.Filter1(t, [], 'bin', tEdges)
            %   [r, tEdges] = MNeural.Filter1(t, [], 'bin', tWin, binSize)
            % 
            % Smooth spike rate with Gaussian kernel
            %   [r, ker] = MNeural.Filter1(r, fs, 'gaussian', sigma)
            %   [r, ker] = MNeural.Filter1(r, fs, 'gaussian', sigma, kerSize)
            % 
            % Smooth spike rate with exponential kernel
            %   [r, ker] = MNeural.Filter1(r, fs, 'exponential', tau)
            %   [r, ker] = MNeural.Filter1(r, fs, 'exponential', tau, kerSize)
            % 
            % Inputs
            %   t           A vector of spike time.
            %   r           A vector of spike rate.
            %   fs          Sampling rate of spike rate.
            %   methodOpt   'bin', 'gaussian' or 'exponential'.
            %   tEdges      A vector of bin edges.
            %   tWin        A 2-element array indicating the start and end of binning.
            %   binSize     Width of time bin.
            %   sigma       Standard deviation of the Gausssian kernel in second.
            %   tau         Time constant of the exponential kernal in second.
            %   kerSize     Length of the kernel in second. The default is 6 sigma for 'gaussian' 
            %               (i.e. 3 sigma each side) or 5 tau for 'exponential'. 
            % Outputs
            %   r           A vector of spike rate.
            %   tEdges      A vector of bin edges.
            %   ker         The kernel used for convolution.
            
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
            %   t = GetJointDist(stimuli)
            %   t = GetJointDist(stimuli, response)
            %   t = GetJointDist(..., 'mask', logicals)
            %   t = GetJointDist(..., 'numBins', numbers)
            %   t = GetJointDist(..., 'rmMethod', options)
            %   t = GetJointDist(..., 'rmParam', params)
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
        
        function [hh_mean, hh_sem, hh_stats] = MeanTimeHistogram(spk_times, t_edges, filterFunc)
            % Compute time histogram and related stats from spike times
            
            if istable(spk_times)
                spk_times = table2cell(spk_times);
            end
            assert(iscell(spk_times), '');
            if isvector(spk_times)
                spk_times = spk_times(:);
            end
            
            for i = size(spk_times,2) : -1 : 1
                % Compute histogram for each trial
                h_trials = cellfun(@(x) histcounts(x, t_edges), spk_times(:,i), 'Uni', false);
                h_trials = cell2mat(h_trials);
                
                % Optionally filter spike counts
                if nargin > 2
                    h_trials = filterFunc(h_trials')';
                end
                
                % Compute mean and sem of spike count
                h_mean = mean(h_trials);
                h_sem = MMath.StandardError(h_trials);
                hh_mean(:,i) = h_mean;
                hh_sem(:,i) = h_sem;
                
                % Compute other stats
                [pk_count(i), pk_bin(i)] = max(h_mean);
                pk_prob(i) = mean(h_trials(:,pk_bin(i)) > 0);
                I(i) = MMath.Entropy(h_mean/sum(h_mean));
            end
            
            hh_stats = table();
            hh_stats.unit_num = (1:size(spk_times,2))';
            hh_stats.pk_bin = pk_bin';
            hh_stats.pk_count = pk_count';
            hh_stats.pk_prob = pk_prob';
            hh_stats.spk_sum = sum(hh_mean)';
            hh_stats.spk_scs = mean(cumsum(hh_mean./pk_count))';
            hh_stats.entropy = I';
        end
        
        function tunings = Tuning(varargin)
            %Computes tuning curves (and errors) for given stimuli
            %
            %   tunings = Tuning()
            %   tunings = Tuning(stimuli)
            %   tunings = Tuning(stimuli, response)
            %   tunings = Tuning(..., 'mask', logicals)
            %   tunings = Tuning(..., 'numBins', value)
            %   tunings = Tuning(..., 'rmMethod', options)
            %   tunings = Tuning(..., 'rmParam', params)
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
            stim = MNeural.Vars2Mat(p.Results.stimuli);
            resp = MNeural.Vars2Mat(p.Results.response);
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
    end
    
    methods(Static, Access = private)
        function vars = Vars2Mat(vars)
            % Converts heterogeneous and/or nested variable specifications to homogeneous numeric array
            %
            %   vars = Vars2Mat(vars)
            %   
            % Inputs
            %   vars            Data of stimuli (each column is one stimulus).
            % Output:
            %   vars            Homogeneous numeric array of column vectors
            
            try
                if iscell(vars)
                    % Nested cell array
                    if isvector(vars)
                        for i = 1 : length(vars)
                            % De-nesting
                            if iscell(vars{i})
                                vars{i} = cell2mat(vars{i});
                            end
                        end
                    end
                    % Cell array
                    vars = cell2mat(vars);
                end
            catch
                error('The input cannot be resolved. Please check for the consistency of potential length of each column.');
            end
        end
    end
end


