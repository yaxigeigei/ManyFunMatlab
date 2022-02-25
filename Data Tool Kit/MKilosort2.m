classdef MKilosort2
    %MKilosort2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        
    end
    
    methods(Static)
        function result = ImportResultsOld(ksDir)
            % Import and process results after sorting and reviewing in TemplateGUI. 
            % (~4300 spikes/second if running on SSD and >10 times slower if on HDD)
            % 
            %   result = ImportResults()
            %   result = ImportResults(ksDir)
            % 
            % Input
            %   ksDir                         The folder which saves Kilosort and Phy outputs. If not provided, 
            %                                 a folder selection window will be prompted
            % Outputs
            %   result                        A structure with the following fields.
            % 
            %     recording_time                The duration of recording in second
            %     sample_rate                   Recording sampling rate
            %     channel_map                   Indices that reorders recording channels
            %     cluster_info                  Information about each cluster
            %     templates                     A [nTemp,time,nChan] array of all templates used in spike sorting
            % 
            %     spike_times                   A [1,nUnit] cell array of [nSpk,1] spike times
            %     spike_template_ids            A [1,nUnit] cell array of [nSpk,1] spike template IDs
            %     spike_template_amplitudes     A [1,nUnit] cell array of [nSpk,1] spike template amplitudes
            %     spike_waveforms               A [1,nUnit] cell array of [nSpk,time] spike waveform extracted from 
            %                                   the primary channel of each unit. Scaled by 0.195 microvolt
            % 
            %     unit_channel_ind              Indices of primary channel (after mapping) for each unit
            %     unit_mean_template            A [time,nChan,nUnit] array of mean templates
            %     unit_mean_waveform            A [time,nUnit] array of mean waveform
            % 
            
            if nargin < 1
                ksDir = MBrowse.Folder([], 'Select the Kilosort output folder');
                if ~ksDir
                    result = [];
                    return;
                end
            end
            
            disp('Import and process outputs from Kilosort and TemplateGUI');
            
            % Load Kilosort and TemplateGUI outputs
            chanMap = readNPY(fullfile(ksDir, 'channel_map.npy')) + 1;          % channel map
            spkTimeInd = readNPY(fullfile(ksDir, 'spike_times.npy'));           % sample indices of all spikes
            spkTimeSec = readNPY(fullfile(ksDir, 'spike_times_in_sec_adj.npy')); % adjusted time of all spikes
            spkClusterIDs = readNPY(fullfile(ksDir, 'spike_clusters.npy'));     % cluster ID of all spikes
            spkTemplateIDs = readNPY(fullfile(ksDir, 'spike_templates.npy'));   % template ID of all spikes
            spkAmplitudes = readNPY(fullfile(ksDir, 'amplitudes.npy'));         % template amplitude of all spikes
            templates = readNPY(fullfile(ksDir, 'templates.npy'));              % [nTemp, time, chan]
            [clusInfo, clusGroup] = MKilosort2.ReadTsvFiles(ksDir, ...
                'cluster_info.tsv', ...     % Summary info and stats for all clusters
                'cluster_group.tsv');       % KS or manually labeled quality for each cluster
            
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
            [~, datName, datExt] = fileparts(dat_path);
            dat_path = fullfile(ksDir, [datName datExt]);
            mdat = memmapfile(dat_path, 'Format', dtype);
            nSample = numel(mdat.Data) / n_channels_dat;
            mdat = memmapfile(dat_path, 'Format', {dtype, [n_channels_dat nSample], 'V'});
            
            % Find IDs of good clusters
            isGoodCluster = ismember(clusGroup.(2), 'good'); % use number indexing since it can be 'KSLabel' or 'group'
            goodClusterIDs = clusGroup.cluster_id(isGoodCluster)';
            
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
                fprintf('Process unit #%i\n', numel(goodClusterIDs)-k+1);
                
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
                %   when a waveform is at the edge of recording, we pad the waveform to the required length 
                %   by replicating the last value
                rawChanIdx = chanMap(chanIdx);
                W = zeros(nSide*2, numel(stInd));
                for i = 1 : numel(stInd)
                    idxStart = stInd(i) - nSide + 1;
                    idxEnd = stInd(i) + nSide;
                    if idxStart < 1
                        idxStart = 1;
                        w = mdat.Data.V(rawChanIdx, idxStart : idxEnd);
                        wPad = repmat(w(1), 1, nSide*2-numel(w));
                        W(:,i) = [wPad w];
                    elseif idxEnd > nSample
                        idxEnd = nSample;
                        w = mdat.Data.V(rawChanIdx, idxStart : idxEnd);
                        wPad = repmat(w(end), 1, nSide*2-numel(w));
                        W(:,i) = [w wPad];
                    else
                        W(:,i) = mdat.Data.V(rawChanIdx, idxStart : idxEnd);
                    end
                end
                
                % Highpass-filter waveforms
                W = single(filtfilt(D, W));
                W = W(nMarg+1:end-nMarg, :);
                
                % Save reuslts
                uChanInd(k) = chanIdx;
                uSpkTimes{k} = spkTimeSec(isUnit);
                uSpkTemplateID{k} = tpID;
                uSpkTemplateAmp{k} = tpAmp;
                uMeanTemplates(:,:,k) = meanTp;
                uWaveforms{k} = W';
                uMeanWaveform(:,k) = mean(W,2);
            end
            [uChanInd, sortInd] = sort(uChanInd);
            
            % Output data
            result.recording_time = nSample / sample_rate;
            result.sample_rate = sample_rate;
            result.channel_map = chanMap;
            result.cluster_info = clusInfo;
            result.templates = templates;
            
            result.spike_times = uSpkTimes(sortInd);
            result.spike_template_ids = uSpkTemplateID(sortInd);
            result.spike_template_amplitudes = uSpkTemplateAmp(sortInd);
            result.spike_waveforms = uWaveforms(sortInd);
            
            result.unit_channel_ind = uChanInd;
            result.unit_mean_template = uMeanTemplates(:,:,sortInd);
            result.unit_mean_waveform = uMeanWaveform(:,sortInd);
        end
        
        function s = ImportResults(ksDir, varargin)
            % Import and process results after sorting and reviewing in Phy. 
            % (~8300 spikes/second if running on SSD and >10 times slower if on HDD)
            % 
            %   result = ImportResults()
            %   result = ImportResults(ksDir)
            % 
            % Input
            %   ksDir                         The folder which saves Kilosort and Phy outputs. If not provided, 
            %                                 a folder selection window will be prompted
            % Outputs
            %   s                             A structure with the following fields.
            % 
            %     channel_map                   Indices that reorders recording channels
            %     cluster_info                  Information about each cluster
            %     templates                     A [nTemp,time,nChan] array of all templates used in spike sorting
            % 
            %     spike_times                   A [1,nUnit] cell array of [nSpk,1] spike times
            %     spike_template_ids            A [1,nUnit] cell array of [nSpk,1] spike template IDs
            %     spike_template_amplitudes     A [1,nUnit] cell array of [nSpk,1] spike template amplitudes
            %     spike_waveforms               A [1,nUnit] cell array of [nSpk,time] spike waveform extracted from 
            %                                   the primary channel of each unit. Scaled by 0.195 microvolt
            % 
            %     unit_channel_ind              Indices of primary channel (after mapping) for each unit
            %     unit_mean_template            A [time,nChan,nUnit] array of mean templates
            %     unit_mean_waveform            A [time,nUnit] array of mean waveform
            % 
            
            if nargin < 1 || isempty(ksDir)
                ksDir = MBrowse.Folder([], 'Select the Kilosort output folder');
                if ~ksDir
                    s = [];
                    return;
                end
            end
            
            p = inputParser();
            p.addParameter('Groups', {'good', 'mua'}, @(x) ismember(x, {'good', 'mua', 'noise'}));
            p.parse(varargin{:});
            groups = lower(cellstr(p.Results.Groups));
            
            % Load channel configuration
            chanInfo = load('NP1_NHP_HalfCol_kilosortChanMap.mat');
            chanInfo = rmfield(chanInfo, 'name');
            chanTb = struct2table(chanInfo);
            chanTb = chanTb(chanTb.connected, :);
            [chanTb, I] = sortrows(chanTb, 'ycoords', 'descend');
            
            % Map binary dat file
            mdat = MKilosort2.MapDatFile(ksDir);
            
            % Load Kilosort and Phy outputs
%             chanMap = readNPY(fullfile(ksDir, 'channel_map.npy')) + 1;          % channel map
            spkTimeInd = readNPY(fullfile(ksDir, 'spike_times.npy'));           % sample indices of all spikes
%             spkTimeSec = readNPY(fullfile(ksDir, 'spike_times_in_sec_adj.npy')); % adjusted time of all spikes
            spkClusIds = readNPY(fullfile(ksDir, 'spike_clusters.npy'));        % cluster ID of all spikes
            spkTempIds = readNPY(fullfile(ksDir, 'spike_templates.npy'));       % template ID of all spikes
            spkTempAmp = readNPY(fullfile(ksDir, 'amplitudes.npy'));            % template amplitude of all spikes
            templates = readNPY(fullfile(ksDir, 'templates.npy'));              % [nTemp, time, chan]
            templates = permute(templates(:,:,I), [3 2 1]);                     % [chan, time, nTemp]
            pcFeats = readNPY(fullfile(ksDir, 'pc_features.npy'));
            pcInd = readNPY(fullfile(ksDir, 'pc_feature_ind.npy'));
            clusGroup = MKilosort2.ReadTsvFiles(ksDir, 'cluster_group.tsv');    % KS or manually labeled quality for each cluster
            
            % Find primary channels (those with highest power) from templates
            nTemp = size(templates, 3);
            primChan = zeros(nTemp, 1);
            for i = 1 : nTemp
                T = templates(:,:,i);
                tpPower = sum(T.^2, 2);
                [~, primChan(i)] = max(tpPower);
            end
            
            % Find clusters
            isGroup = ismember(clusGroup.(2), groups); % use number indexing since it can be 'KSLabel' or 'group'
            clusIds = clusGroup.cluster_id(isGroup);
            
            % Select spikes for the clusters of interest
            m = ismember(spkClusIds, clusIds);
            spkTimeInd = spkTimeInd(m);
            spkClusIds = spkClusIds(m);
            spkTempIds = spkTempIds(m);
            spkTempAmp = spkTempAmp(m);
            spkPrimChans = primChan(spkTempIds+1);
            pcFeats = pcFeats(m,:,:);
            
            % Chunking
            nSpk = sum(m);
            chunkSz = 1e3;
            nChunk = ceil(nSpk / chunkSz);
            WW = cell(nChunk,1);
            
            % Extract waveform
            disp('Extracting spike waveform');
            for i = 1 : nChunk
                a = (i-1)*chunkSz + 1;
                b = min(i*chunkSz, nSpk);
                WW{i} = MKilosort2.ReadSnippets(mdat, spkTimeInd(a:b), 82, spkPrimChans(a:b), 'ChannelOrder', I);
            end
            WW = cat(3, WW{:});
            
            % Group data by clusters
            G = findgroups(spkClusIds);
            spkTimeInd = splitapply(@(x) {x}, spkTimeInd, G);
%             spkTimeSec = splitapply(@(x) {x}, spkTimeSec, G);
            spkTempIds = splitapply(@(x) {x}, spkTempIds, G);
            spkTempAmp = splitapply(@(x) {x}, spkTempAmp, G);
            spkPrimChans = splitapply(@(x) {x}, spkPrimChans, G);
            for i = numel(clusIds) : -1 : 1
                % Compute mean template
                tps = arrayfun(@(i,a) templates(:,:,i+1) * a, spkTempIds{i}, spkTempAmp{i}, 'Uni', false);
                tps = cat(3, tps{:});
                meanTemp(:,:,i) = squeeze(mean(tps,3));
                
                % Split pc features
                m = spkClusIds == clusIds(i);
                pcFeatClus{i,1} = pcFeats(m,:,:);
                
                % Split waveform and compute mean stats
                W = WW(:,:,m);
                clusWaveform{i,1} = W;
                [meanWaveform(:,:,i), sdWaveform(:,:,i)] = MMath.MeanStats(double(W), 3);
            end
            
            % Output data
            s.channel_info = chanTb;
%             result.cluster_info = clusInfo;
            s.templates = templates;
            s.pc_features = pcFeatClus;
            s.pc_feature_ind = pcInd;
            s.cluster_ids = unique(G, 'stable');
            s.spike_times = spkTimeInd;
            s.spike_template_ids = spkTempIds;
            s.spike_template_amplitudes = spkTempAmp;
            s.spike_prim_chan = spkPrimChans;
            s.spike_waveform = clusWaveform;
            s.mean_template = meanTemp;
            s.mean_waveform = meanWaveform;
            s.sd_waveform = sdWaveform;
        end
        
        function varargout = ReadTsvFiles(ksDir, varargin)
            % 
            
            varargout = cell(size(varargin));
            
            for k = 1 : numel(varargin)
                % Check file
                fileName = varargin{k};
                filePath = fullfile(ksDir, fileName);
                if ~exist(filePath, 'file')
                    warning('%s does not exist', fileName);
                    continue
                end
                
                % Specify the pattern to read
                if strcmp(fileName, 'cluster_info.tsv')
                    formatStr = '%f %f %f %s %f %f %f %f %s %f %f';
                    headers = {'cluster_id', 'Amplitude', 'ContamPct', 'KSLabel', 'amp', 'ch', ...
                        'depth', 'fr', 'group', 'n_spikes', 'sh'};
                    
                elseif strcmp(fileName, 'cluster_group.tsv')
                    formatStr = '%f %s';
                    headers = {'cluster_id', 'group'};
                    
                elseif strcmp(fileName, 'cluster_KSLabel.tsv')
                    formatStr = '%f %s';
                    headers = {'cluster_id', 'KSLabel'};
                    
                else
                    warning('Reading %s is not supported', fileName);
                    continue
                end
                
                % Read file
                fid = fopen(filePath);
                C = textscan(fid, formatStr, 'HeaderLines', 1);
                fclose(fid);
                
                % Make a tbale
                tb = table;
                for i = 1 : numel(headers)
                    tb.(headers{i}) = C{i};
                end
                varargout{k} = tb;
            end
        end
        
        function mdat = MapDatFile(ksDir)
            % Make a memory map to temp_wh.dat
            
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
            
            % Map binary data to memory
            if ~exist(dat_path, 'file')
                [~, datName, datExt] = fileparts(dat_path);
                dat_path = fullfile(ksDir, [datName datExt]);
            end
            mdat = memmapfile(dat_path, 'Format', dtype);
            nSample = numel(mdat.Data) / n_channels_dat;
            mdat = memmapfile(dat_path, 'Format', {dtype, [n_channels_dat nSample], 'V'});
        end
        
        function sn = ReadSnippets(mmap, tmInd, tmWin, varargin)
            % Read snippets from mapped binary data array
            %
            %   sn = ReadSnippets(mmap, tmInd, tmWin)
            %   sn = ReadSnippets(mmap, tmInd, tmWin, chInd)
            %   sn = ReadSnippets(mmap, tmInd, tmWin, chInd, chWin)
            %   sn = ReadSnippets(mmap, tmInd, tmWin, ..., 'ChannelOrder', [])
            %   sn = ReadSnippets(mmap, tmInd, tmWin, ..., 'VoltScale', 1)
            %   sn = ReadSnippets(mmap, tmInd, tmWin, ..., 'Filter', [])
            % 
            % Inputs
            %   mmap                A memmapfile object to the binary data with the array size specified.
            %   tmInd               A vector of time indices.
            %   tmWin               The time window. The i-th snippet is read from tmInd(i)+tmWin(1) to tmInd(i)+tmWin(2), inclusively.
            %   chInd               A vector of channel indices. If scalar, the same channel will be extracted for all snippets.
            %                       The default is empty, reading all channels.
            %   chWin               The window of channels. The i-th snippet is read from tmInd(i)+tmWin(1) to tmInd(i)+tmWin(2), inclusively.
            %                       The default is [0 0], reading only the channel specified by chInd. chWin is ignored if chInd is empty.
            %   'ChannelOrder'      A vector of indices that reorder the channels in binary array before extracting snippets.
            %   'VoltScale'         A scaling factor applied to the data.
            %   'Filter'            A digitalFilter object used to filter the snippet data (by each channel in each snippet).
            % Output
            %   sn                  A channel-by-time-by-#snippets array of signal values.
            % 
            
            tmInd = double(tmInd(:));
            nSn = numel(tmInd);
            [nCh, nTm] = size(mmap.Data.V);
            
            p = inputParser();
            p.addOptional('chInd', [], @(x) isnumeric(x) && isvector(x));
            p.addOptional('chWin', [0 0], @(x) isnumeric(x) && numel(x)==2);
            p.addParameter('ChannelOrder', [], @(x) isnumeric(x) && isvector(x));
            p.addParameter('VoltScale', 1, @(x) isnumeric(x) && isvector(x));
            p.addParameter('Filter', [], @(x) isa(x, 'digitalFilter'));
            p.parse(varargin{:});
            chInd = double(p.Results.chInd(:));
            chWin = p.Results.chWin(:)';
            cOrder = p.Results.ChannelOrder;
            vScale = p.Results.VoltScale;
            filtObj = p.Results.Filter;
            
            if isempty(chInd)
                chInd = 1;
                chWin = [0 nCh-1];
            end
            if isscalar(chInd)
                chInd = ones(size(tmInd))*chInd;
            end
            
            if isempty(cOrder)
                cOrder = 1 : nCh;
            end
            
            % Read data
            tmWins = tmInd + tmWin;
            chWins = chInd + chWin;
            for i = numel(tmInd) : -1 : 1
                % Get channel indices
                ind1 = chWins(i,1) : chWins(i,2);
                isChOut = ind1 < 1 | ind1 > nCh;
                ind1 = min(max(ind1, 1), nCh);
                ind1 = cOrder(ind1);
                
                % Get sample indices
                ind2 = tmWins(i,1) : tmWins(i,2);
                ind2 = min(max(ind2, 1), nTm); % propagate end values to fill out of range samples
                
                % Read from mapped data
                sn(:,:,i) = mmap.Data.V(ind1, ind2);
                sn(isChOut,:,i) = 0; % reset values in out of range channels to zeros
            end
            
            if vScale ~= 1
                sn = double(sn) .* vScale(:); % apply voltage scaling
            end
            
            % Highpass filtering
            if ~isempty(filtObj)
                sn = double(sn);
                for i = 1 : nCh
                    for j = 1 : nSn
                        sn(i,:,j) = filtfilt(filtObj, sn(i,:,j));
                    end
                end
            end
        end
        
        function D = MakeHighpassFilter(fs)
            % Prepare parameters for filtering
            %   1) Phy uses bandpass between 500Hz and 14.25kHz(.475*sample_rate) to visulize waveform;
            %      here uses 300Hz highpass to be consistent with Kilosort configuration
            %   2) 82 samples per waveform is consistent with Kilosort template size
            %   3) Consistent with Phy, using samples 3 times of the filter order as margins when filtering
            if nargin < 1
                fs = 30e3;
            end
            D = designfilt('highpassiir', ...
                'PassbandFrequency', 300, ...
                'StopbandFrequency', 250, ...
                'StopbandAttenuation', 60, ...
                'PassbandRipple', 0.1, ...
                'SampleRate', fs, ...
                'DesignMethod', 'ellip'); % order of this filter is 8
        end
        
    end
end
