classdef MKilosort
    %MKilosort Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        chanMapH3 = [ ... % from tip to base or anatomically basal to apical
            27, 25, 23, 20, 22, 31, 29, 18, ...
            30, 32, 15, 21, 19, 24, 26, 28, ...
            02, 04, 05, 07, 09, 11, 13, 17, ...
            16, 14, 12, 10, 08, 06, 03, 01, ...
            64, 62, 59, 57, 55, 53, 51, 47, ...
            50, 52, 54, 56, 58, 60, 61, 63, ...
            37, 39, 41, 46, 44, 33, 35, 48, ...
            36, 34, 49, 43, 45, 42, 40, 38 ];
        
        chanMapH2 = [ ... % Shank A to B, from tip to base or anatomically basal to apical, 
            64, 62, 59, 57, 55, 53, 51, 47, ...
            50, 52, 54, 56, 58, 60, 61, 63, ...
            37, 39, 41, 46, 44, 33, 35, 48, ...
            36, 34, 49, 43, 45, 42, 40, 38, ...
            27, 25, 23, 20, 22, 31, 29, 18, ...
            30, 32, 15, 21, 19, 24, 26, 28, ...
            02, 04, 05, 07, 09, 11, 13, 17, ...
            16, 14, 12, 10, 08, 06, 03, 01 ];
        
        chanMapP64 = [ ... % Shank A to D, from tip to base or anatomically basal to apical
            37, 39, 41, 46, 44, 33, 43, 35, 45, 48, 42, 36, 40, 34, 38, 49, ... % Shank A
            64, 62, 59, 57, 55, 53, 56, 51, 58, 47, 60, 50, 61, 52, 63, 54, ... % Shank B
            02, 04, 05, 07, 09, 11, 10, 13, 08, 17, 06, 16, 03, 14, 01, 12, ... % Shank C
            27, 25, 23, 20, 22, 31, 21, 29, 19, 18, 24, 30, 26, 32, 28, 15];    % Shank D
        
        chanMapTetrode32 = [ 9:24, 8:-1:1, 32:-1:25 ];
    end
    
    methods(Static)
        function Sort(datFilePath, chanMapName)
            % Run spike sorting on a binary file
            
            % Find .dat file
            if nargin < 1
                datFilePath = MBrowse.File([], 'Select a Kilosort .dat file');
            end
            if isempty(datFilePath)
                return;
            end
            ksDir = fileparts(datFilePath);
            
            % Choose channel map if not provided
            if nargin < 2
                chanMapList = properties(MKilosort);
                selected = listdlg('PromptString', 'Select a channel map', ...
                    'SelectionMode', 'single', ...
                    'ListString', chanMapList);
                if isempty(selected)
                    return;
                end
                chanMapName = chanMapList{selected};
            end
            
            % Load correct option structure and save the channel map
            switch chanMapName
                case 'chanMapH3'
                    ops = MKilosort.Config64(datFilePath);
                    MKilosort.SaveChanMapH3(ksDir);
                case 'chanMapH2'
                    ops = MKilosort.Config64(datFilePath);
                    MKilosort.SaveChanMapH2(ksDir);
                case 'chanMapP64'
                    ops = MKilosort.Config64(datFilePath);
                    MKilosort.SaveChanMapP64(ksDir);
                case 'chanMapTetrode32'
                    ops = MKilosort.Config32(datFilePath);
                    MKilosort.SaveChanMapTetrode32(ksDir);
                otherwise
                    error('%s is not a valid name of channel map.', chanMapName);
            end
            
            % Initialize GPU (will erase any existing GPU arrays)
            if ops.GPU
                gpuDevice(1);
            end
            
            % Preprocess data and extract spikes for initialization
            reset(parallel.gpu.GPUDevice.current);
            [rez, DATA, uproj] = preprocessData(ops);
            
            % Fit templates iteratively
            reset(parallel.gpu.GPUDevice.current);
            rez = fitTemplates(rez, DATA, uproj);
            
            % Extract final spike times (overlapping extraction)
            reset(parallel.gpu.GPUDevice.current);
            rez = fullMPMU(rez, DATA);
            
            % % AutoMerge. rez2Phy will use for clusters the new 5th column of st3 if you run this
            % rez = merge_posthoc2(rez);
            
            % save python results file for Phy
            rezToPhy(rez, ksDir);
            
            % % save matlab results file for future use (although you should really only be using the manually validated spike_clusters.npy file)
            % save(fullfile(fpath, 'rez.mat'), 'rez', '-v7.3');
            
            % remove temporary file
            delete(ops.fproc);
            reset(parallel.gpu.GPUDevice.current) % to free memory space!
            
        end
        
        function result = Output(ksDir)
            % Import results after reviewing in TemplateGUI
            
            if nargin < 1
                ksDir = MBrowse.Folder([], 'Select the Kilosort folder');
            end
            
            % Get parameters from params.py by running lines in MATLAB
            %   these include dat_path, dtype, n_channels_dat, sample_rate ...
            fid = fopen(fullfile(ksDir, 'params.py'));
            while ~feof(fid)
                try
                    l = fgetl(fid);
                    eval([l ';']);
                catch
                    %warning('''%s'' cannot be evaluated.', l);
                end
            end
            fclose(fid);
            
            % Map raw data to memory
            datPath = fullfile(ksDir, dat_path);
            mdat = memmapfile(datPath, 'Format', dtype);
            nSample = numel(mdat.Data) / n_channels_dat;
            mdat = memmapfile(datPath, 'Format', {dtype, [n_channels_dat nSample], 'v'});
            
            % Load sorting outputs
            chanMap = readNPY(fullfile(ksDir, 'channel_map.npy')) + 1;          % channel map
            spkTimeInd = readNPY(fullfile(ksDir, 'spike_times.npy'));           % sample indices of all spikes
            spkClusterIDs = readNPY(fullfile(ksDir, 'spike_clusters.npy'));     % cluster ID of all spikes
            spkTemplateIDs = readNPY(fullfile(ksDir, 'spike_templates.npy'));   % template ID of all spikes
            spkAmplitudes = readNPY(fullfile(ksDir, 'amplitudes.npy'));         % template amplitude of all spikes
            templates = readNPY(fullfile(ksDir, 'templates.npy'));              % [nTemp, time, chan]
            qualityTable = readtable(fullfile(ksDir, 'cluster_groups.csv'));    % contains ID of good clusters
            
            % Find IDs of good clusters
            isGoodCluster = ismember(qualityTable.group, 'good');
            goodClusterIDs = qualityTable.cluster_id(isGoodCluster)';
            
            % Prepare parameters for waveform extraction
            %   1) Phy uses bandpass between 500Hz and 14.25kHz(.475*sample_rate) to visulize waveform; 
            %      here uses 300Hz highpass to be consistent with Kilosort configuration
            %   2) 82 samples per waveform is consistent with Kilosort template size
            %   3) Consistent with Phy, using samples 3 times of the filter order as margins when filtering
            D = designfilt('highpassiir', ...
                'PassbandFrequency', 300, ...
                'StopbandFrequency', 250, ...
                'StopbandAttenuation', 60, ...
                'PassbandRipple', 0.1, ...
                'SampleRate', sample_rate, ...
                'DesignMethod', 'ellip'); % order of this filter is 8
            nHalfW = 41;
            nMarg = filtord(D) * 3;
            nSide = nHalfW + nMarg;
            
            % Loop through good clusters
            for k = numel(goodClusterIDs) : -1 : 1
                % Get sample indices and template IDs of the current cluster
                isUnit = spkClusterIDs == goodClusterIDs(k);
                stInd = spkTimeInd(isUnit);
                tpID = spkTemplateIDs(isUnit);
                tpAmp = spkAmplitudes(isUnit);
                
                % Compute mean template
                tps = arrayfun(@(i,a) templates(i+1,:,:) * a, tpID, tpAmp, 'Uni', false);
                tps = cat(1, tps{:});
                meanTp = squeeze(mean(tps,1));
                
                % Determine the primary channel
                tpPower = sum(meanTp.^2);
                [~, chanIdx] = max(tpPower);
                
                % Collect raw spike waveforms
                rawChanIdx = chanMap(chanIdx);
                W = zeros(nSide*2, numel(stInd));
                for i = 1 : numel(stInd)
                    W(:,i) = mdat.Data.v(rawChanIdx, stInd(i)-nSide+1 : stInd(i)+nSide);
                end
                
                % Highpass-filter waveforms
                W = single(filtfilt(D, W));
                W = W(nMarg+1:end-nMarg, :);
                
                % Scale back to original unit in Intan
                W = W * 0.195; % 0.195 is the resolution of Intan. previoely used for scaling up to integers.
                
                % Save reuslts
                uChanInd(k) = chanIdx;
                uSpkTimes{k} = double(stInd) / sample_rate;
                uSpkTemplateID{k} = tpID;
                uSpkTemplateAmp{k} = tpAmp;
                uMeanTemplates(:,:,k) = meanTp;
                uWaveforms{k} = W';
                uMeanWaveform(:,k) = mean(W,2);
            end
            [uChanInd, sortInd] = sort(uChanInd);
            
            % Output data
            result.spike_times = uSpkTimes(sortInd);
            result.spike_template_ids = uSpkTemplateID(sortInd);
            result.spike_template_amplitudes = uSpkTemplateAmp(sortInd);
            result.spike_waveforms = uWaveforms(sortInd);
            result.info.recording_time = nSample / sample_rate;
            result.info.sample_rate = sample_rate;
            result.info.templates = templates;
            result.info.unit_mean_template = uMeanTemplates(:,:,sortInd);
            result.info.unit_mean_waveform = uMeanWaveform(:,sortInd);
            result.info.unit_channel_ind = uChanInd;
        end
        
        function ops = Config64(datFilePath)
            % This configuration works with H3, H2, and P-64chan when running on GPU with 3GB+ memory
            
            ksDir = fileparts(datFilePath);
            
            ops.GPU                 = 1; % whether to run this code on an Nvidia GPU (much faster, mexGPUall first)
            ops.parfor              = 0; % whether to use parfor to accelerate some parts of the algorithm
            ops.verbose             = 1; % whether to print command line progress
            ops.showfigures         = 1; % whether to plot figures during optimization
            
            ops.datatype            = 'dat';  % binary ('dat', 'bin') or 'openEphys'
            ops.fbinary             = datFilePath; % will be created for 'openEphys'
            ops.fproc               = fullfile(ksDir, 'temp_wh.dat'); % residual from RAM of preprocessed data
            ops.root                = ksDir; % 'openEphys' only: where raw files ar
            
            % ops.fs                  = 30000;        % sampling rate		(omit if already in chanMap file)
            % ops.NchanTOT            = 32;           % total number of channels (omit if already in chanMap file)
            % ops.Nchan               = 32;           % number of active channels (omit if already in chanMap file)
            ops.Nfilt               = 128;  % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)
            ops.nNeighPC            = 12; % visualization only (Phy): number of channnels to mask the PCs, leave empty to skip (12)
            ops.nNeigh              = 16; % visualization only (Phy): number of neighboring templates to retain projections of (16)
            
            % options for channel whitening
            ops.whitening           = 'full'; % type of whitening (default 'full', for 'noSpikes' set options for spike detection below)
            ops.nSkipCov            = 1; % compute whitening matrix from every N-th batch (1)
            ops.whiteningRange      = 32; % how many channels to whiten together (Inf for whole probe whitening, should be fine if Nchan<=32)
            
            % define the channel map as a filename (string) or simply an array
            ops.chanMap             = fullfile(ksDir, 'chanMap.mat'); % make this file using createChannelMapFile.m
            ops.criterionNoiseChannels = 0.2; % fraction of "noise" templates allowed to span all channel groups (see createChannelMapFile for more info).
            % ops.chanMap = 1:ops.Nchan; % treated as linear probe if a chanMap file
            
            % other options for controlling the model and optimization
            ops.Nrank               = 3;    % matrix rank of spike template model (3)
            ops.nfullpasses         = 6;    % number of complete passes through data during optimization (6)
            ops.maxFR               = 20000;  % maximum number of spikes to extract per batch (20000)
            ops.fshigh              = 300;   % frequency for high pass filtering
            % ops.fslow               = 6000;   % frequency for low pass filtering (optional)
            ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
            ops.scaleproc           = 200;   % int16 scaling of whitened data
            ops.NT                  = 64*1024*2 + ops.ntbuff;% this is the batch size (try decreasing if out of memory)
            % for GPU should be multiple of 32 + ntbuff		(originally 128*1024+ ops.ntbuff)
            
            % the following options can improve/deteriorate results.
            % when multiple values are provided for an option, the first two are beginning and ending anneal values,
            % the third is the value used in the final pass.
            ops.Th                  = [4 10 10];    % threshold for detecting spikes on template-filtered data ([6 12 12])
            ops.lam                 = [5 20 20];    % large means amplitudes are forced around the mean ([10 30 30])
            ops.nannealpasses       = 4;            % should be less than nfullpasses (4)
            ops.momentum            = 1./[20 400];  % start with high momentum and anneal (1./[20 1000])
            ops.shuffle_clusters    = 1;            % allow merges and splits during optimization (1)
            ops.mergeT              = .1;           % upper threshold for merging (.1)
            ops.splitT              = .1;           % lower threshold for splitting (.1)
            
            % options for initializing spikes from data
            ops.initialize          = 'fromData';    %'fromData' or 'no'
            ops.spkTh               = -4;      % spike threshold in standard deviations (4)
            ops.loc_range           = [3 1];   % ranges to detect peaks; plus/minus in time and channel ([3 1])
            ops.long_range          = [30 6];  % ranges to detect isolated peaks ([30 6])
            ops.maskMaxChannels     = 5;       % how many channels to mask up/down ([5])
            ops.crit                = .65;     % upper criterion for discarding spike repeates (0.65)
            ops.nFiltMax            = 10000;   % maximum "unique" spikes to consider (10000)
            
            % load predefined principal components (visualization only (Phy): used for features)
            dd                      = load('PCspikes2.mat'); % you might want to recompute this from your own data
            ops.wPCA                = dd.Wi(:,1:7);   % PCs
            
            % options for posthoc merges (under construction)
            ops.fracse              = 0.1; % binning step along discriminant axis for posthoc merges (in units of sd)
            ops.epu                 = Inf;
            
            ops.ForceMaxRAMforDat   = 20e9; % maximum RAM the algorithm will try to use; on Windows it will autodetect.
            
        end
        
        function ops = Config32(datFilePath)
            % This configuration works with H3, H2, and P-64chan when running on GPU with 3GB+ memory
            
            ksDir = fileparts(datFilePath);
            
            ops.GPU                 = 1; % whether to run this code on an Nvidia GPU (much faster, mexGPUall first)
            ops.parfor              = 0; % whether to use parfor to accelerate some parts of the algorithm
            ops.verbose             = 1; % whether to print command line progress
            ops.showfigures         = 1; % whether to plot figures during optimization
            
            ops.datatype            = 'dat';  % binary ('dat', 'bin') or 'openEphys'
            ops.fbinary             = datFilePath; % will be created for 'openEphys'
            ops.fproc               = fullfile(ksDir, 'temp_wh.dat'); % residual from RAM of preprocessed data
            ops.root                = ksDir; % 'openEphys' only: where raw files ar
            
            % ops.fs                  = 30000;        % sampling rate		(omit if already in chanMap file)
            % ops.NchanTOT            = 32;           % total number of channels (omit if already in chanMap file)
            % ops.Nchan               = 32;           % number of active channels (omit if already in chanMap file)
            ops.Nfilt               = 64;  % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)
            ops.nNeighPC            = 12; % visualization only (Phy): number of channnels to mask the PCs, leave empty to skip (12)
            ops.nNeigh              = 16; % visualization only (Phy): number of neighboring templates to retain projections of (16)
            
            % options for channel whitening
            ops.whitening           = 'full'; % type of whitening (default 'full', for 'noSpikes' set options for spike detection below)
            ops.nSkipCov            = 1; % compute whitening matrix from every N-th batch (1)
            ops.whiteningRange      = 32; % how many channels to whiten together (Inf for whole probe whitening, should be fine if Nchan<=32)
            
            % define the channel map as a filename (string) or simply an array
            ops.chanMap             = fullfile(ksDir, 'chanMap.mat'); % make this file using createChannelMapFile.m
            ops.criterionNoiseChannels = 0.2; % fraction of "noise" templates allowed to span all channel groups (see createChannelMapFile for more info).
            % ops.chanMap = 1:ops.Nchan; % treated as linear probe if a chanMap file
            
            % other options for controlling the model and optimization
            ops.Nrank               = 3;    % matrix rank of spike template model (3)
            ops.nfullpasses         = 6;    % number of complete passes through data during optimization (6)
            ops.maxFR               = 20000;  % maximum number of spikes to extract per batch (20000)
            ops.fshigh              = 300;   % frequency for high pass filtering
            % ops.fslow               = 6000;   % frequency for low pass filtering (optional)
            ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
            ops.scaleproc           = 200;   % int16 scaling of whitened data
            ops.NT                  = 64*1024*2 + ops.ntbuff;% this is the batch size (try decreasing if out of memory)
            % for GPU should be multiple of 32 + ntbuff		(originally 128*1024+ ops.ntbuff)
            
            % the following options can improve/deteriorate results.
            % when multiple values are provided for an option, the first two are beginning and ending anneal values,
            % the third is the value used in the final pass.
            ops.Th                  = [4 10 10];    % threshold for detecting spikes on template-filtered data ([6 12 12])
            ops.lam                 = [5 20 20];    % large means amplitudes are forced around the mean ([10 30 30])
            ops.nannealpasses       = 4;            % should be less than nfullpasses (4)
            ops.momentum            = 1./[20 400];  % start with high momentum and anneal (1./[20 1000])
            ops.shuffle_clusters    = 1;            % allow merges and splits during optimization (1)
            ops.mergeT              = .1;           % upper threshold for merging (.1)
            ops.splitT              = .1;           % lower threshold for splitting (.1)
            
            % options for initializing spikes from data
            ops.initialize          = 'fromData';    %'fromData' or 'no'
            ops.spkTh               = -4;      % spike threshold in standard deviations (4)
            ops.loc_range           = [3 1];   % ranges to detect peaks; plus/minus in time and channel ([3 1])
            ops.long_range          = [30 6];  % ranges to detect isolated peaks ([30 6])
            ops.maskMaxChannels     = 5;       % how many channels to mask up/down ([5])
            ops.crit                = .65;     % upper criterion for discarding spike repeates (0.65)
            ops.nFiltMax            = 10000;   % maximum "unique" spikes to consider (10000)
            
            % load predefined principal components (visualization only (Phy): used for features)
            dd                      = load('PCspikes2.mat'); % you might want to recompute this from your own data
            ops.wPCA                = dd.Wi(:,1:7);   % PCs
            
            % options for posthoc merges (under construction)
            ops.fracse              = 0.1; % binning step along discriminant axis for posthoc merges (in units of sd)
            ops.epu                 = Inf;
            
            ops.ForceMaxRAMforDat   = 20e9; % maximum RAM the algorithm will try to use; on Windows it will autodetect.
            
        end
        
        function SaveChanMapH3(ksDir)
            % Create a channel map file
            
            % here I know a priori what order my channels are in.  So I just manually
            % make a list of channel indices (and give
            % an index to dead channels too). chanMap(1) is the row in the raw binary
            % file for the first channel. chanMap(1:2) = [33 34] in my case, which happen to
            % be dead channels.
            chanMap = MKilosort.chanMapH3;
            
            % chanMap order in 'F:\O'Connor lab data (other than
            % Intan)\electrophysiology\Silicon Probe\correspondence of probe
            % channels.xls' plus one (so that the chanNum will not have zero)
            
            % the first thing Kilosort does is reorder the data with data = data(chanMap, :).
            % Now we declare which channels are "connected" in this normal ordering,
            % meaning not dead or used for non-ephys data
            connected = true(1, 64);
            % connected(1:2) = 0;
            
            % now we define the horizontal (x) and vertical (y) coordinates of these
            % 34 channels. For dead or nonephys channels the values won't matter. Again
            % I will take this information from the specifications of the probe. These
            % are in um here, but the absolute scaling doesn't really matter in the
            % algorithm.
            xcoords = zeros(1,64);
            ycoords = (0:63)*20;
            
            % Often, multi-shank probes or tetrodes will be organized into groups of
            % channels that cannot possibly share spikes with the rest of the probe. This helps
            % the algorithm discard noisy templates shared across groups. In
            % this case, we set kcoords to indicate which group the channel belongs to.
            % In our case all channels are on the same shank in a single group so we
            % assign them all to group 1.
            kcoords = ones(1,64);
            
            % at this point in Kilosort we do data = data(connected, :), ycoords =
            % ycoords(connected), xcoords = xcoords(connected) and kcoords =
            % kcoords(connected) and no more channel map information is needed (in particular
            % no "adjacency graphs" like in KlustaKwik).
            % Now we can save our channel map for the eMouse.
            
            % would be good to also save the sampling frequency here
            fs = 30000;
            
            save(fullfile(ksDir, 'chanMap.mat'), 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'fs')
        end
        
        function SaveChanMapH2(ksDir)
            % Create a channel map file
            
            chanMap = MKilosort.chanMapH2;  % channel order
            connected = true(1, 64);    % dead channels are set to false
            
            % Probe geometry
            xcoords = [zeros(1,32), zeros(1,32)+250];
            ycoords = [(0:31)*25, (0:31)*25];
            
            % Channel grouping
            kcoords = [ones(1,32), ones(1,32)*2];
            
            fs = 30000;                 % sampling frequency
            
            save(fullfile(ksDir, 'chanMap.mat'), 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'fs')
        end
        
        function SaveChanMapP64(ksDir)
            % Create a channel map file
            
            chanMap = MKilosort.chanMapP64; % channel order
            connected = true(1, 64);    % dead channels are set to false
            
            % Probe geometry
            xshank = repmat([22.5 0], 1, 8);
            xcoords = [xshank, xshank+250, xshank+500, xshank+750];
            ycoords = repmat((0:15)*12.5, 1, 4);
            
            % Channel grouping
            kcoords = [ones(1,16), ones(1,16)*2, ones(1,16)*3, ones(1,16)*4];
            
            fs = 30000;                 % sampling frequency
            
            save(fullfile(ksDir, 'chanMap.mat'), 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'fs')
        end
        
        function SaveChanMapTetrode32(ksDir)
            % Create a channel map file
            
            chanMap = MKilosort.chanMapTetrode32; % channel order
            connected = true(1, 32);    % dead channels are set to false
            
            % Probe geometry
            xcoords = zeros(1,32);
            ycoords = ((0:31) + repelem(0:7, 4))*10;
            
            % Channel grouping
            kcoords = repelem(1:8, 4);
            
            fs = 30000;                 % sampling frequency
            
            save(fullfile(ksDir, 'chanMap.mat'), 'chanMap', 'connected', 'xcoords', 'ycoords', 'kcoords', 'fs')
        end
    end
end

