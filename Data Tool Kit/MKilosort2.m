classdef MKilosort2
    %MKilosort2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties(Constant)
        
    end
    
    methods(Static)
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
