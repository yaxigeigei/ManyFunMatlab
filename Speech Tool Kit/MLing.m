classdef MLing
    % MLing is...
    
    methods(Static)
        function [k1, k2, alignInfo] = FindAlignedWords(s1, s2)
            % Globally align two sentences at the level of letters, then find pairs of corresponding words.
            % 
            %   [k1, k2, alignInfo] = MLing.FindAlignedWords(s1, s2)
            % 
            % Inputs
            %   s1, s2          Two sequences of words (e.g. sentences). If an input is a vector of 
            %                   strings, elements will be joined together by spaces into a sentence.
            % Outputs
            %   k1, k2          Indices of matching words in s1 and s2, respectively. Namely, the 
            %                   k1(i)-th word in s1 matches the k2(i)-th word in s2.
            %   alignInfo       A struct storing the alignment result in the following fields.
            %                   s1, s2                  The two input sequences. Concatenated if any 
            %                                           of the input is a vector of words.
            %                   aa, steps, tscore       Outputs from MLing.SeqAlign function. See 
            %                                           its documentation for details.
            % 
            
            ss = {s1, s2};
            wdInd = cell(size(ss));
            for i = 1 : numel(ss)
                s = ss{i};
                
                % Find which word (by index) each character belongs to
                s = strtrim(s);
                if ~ischar(s) && numel(s) > 1
                    w = s;
                    s = strjoin(w, ' ');
                else
                    w = strsplit(s, ' ');
                end
                wLen = strlength(w) + 1; % plus one to include a trailing space
                wInd = repelem(1:numel(w), wLen);
                wInd(end) = []; % remove index of the trailing space of the last word
                
                ss{i} = s;
                wdInd{i} = wInd;
            end
            
            % Global alignment
            [aa, steps, tscore] = MLing.SeqAlign(char(ss{1}), char(ss{2}));
            
            % Find each character's word index in the aligned strings
            chInd = cumsum(steps, 2);
            chInd(~chInd) = 1; % avoid zero indexing
            wdInd = [wdInd{1}(chInd(1,:)); wdInd{2}(chInd(2,:))];
            
            % Check word pairs
            isGap = any(~steps, 1);
            wdIndNoGap = wdInd(:, ~isGap);
            prInd = unique(wdIndNoGap', 'rows', 'stable');
            
            % Output
            k1 = prInd(:,1);
            k2 = prInd(:,2);
            alignInfo.s1 = ss{1};
            alignInfo.s2 = ss{2};
            alignInfo.aa = aa;
            alignInfo.steps = steps;
            alignInfo.tscore = tscore;
        end
        
        function [k1, k2, alignInfo] = FindAlignedTokens(s1, s2)
            % Globally align two sequences of tokens, and find pairs of corresponding ones.
            % 
            %   [k1, k2, alignInfo] = MLing.FindAlignedTokens(s1, s2)
            % 
            % Inputs
            %   s1, s2          Two sequences of tokens (e.g. phonemes), each is a vector of string
            %                   objects or a cell vector of char strings.
            % Outputs
            %   k1, k2          Indices of matching tokens in s1 and s2, respectively. Namely, the 
            %                   k1(i)-th token in s1 matches the k2(i)-th token in s2.
            %   alignInfo       A struct storing the alignment result in the following fields.
            %                   s1, s2                  The two input sequences.
            %                   aa, steps, tscore       Outputs from MLing.SeqAlign function. See 
            %                                           its documentation for details.
            % 
            
            % Global alignment
            s1 = string(s1);
            s2 = string(s2);
            [aa, steps, tscore] = MLing.SeqAlign(s1, s2);
            
            % Get token indices
            tkInd = cumsum(steps, 2);
            
            % Find token pairs
            isGap = any(~steps, 1);
            tkInd = tkInd(:, ~isGap);
            
            % Output
            k1 = tkInd(1,:)';
            k2 = tkInd(2,:)';
            alignInfo.s1 = s1;
            alignInfo.s2 = s2;
            alignInfo.aa = aa;
            alignInfo.steps = steps;
            alignInfo.tscore = tscore;
        end
        
        function [wf, t] = ReadTimitWaveform(timitDir, id)
            % Read TIMIT audio waveform
            %
            %   tg = ReadTimitWaveform(timitDir, id)
            % 
            % Inputs
            %   srcDir      The source directory path of TIMIT feature/label files.
            %   id          A cell or string array of TIMIT stim IDs.
            % Output
            %   wf          A cell array of audio waveform.
            %   t           Timestamps associated with wf in seconds.
            % 
            
            id = cellstr(id);
            assert(exist(timitDir, 'dir'), "The source directory does not exist.");
            wf = cell(size(id));
            t = wf;
            fs = 16e3; % hardcode the TIMIT sampling frequency
            for i = 1 : numel(wf)
                p = fullfile(timitDir, id{i} + ".wav");
                wf{i} = audioread(p);
                t{i} = (1:length(wf{i}))' / fs;
            end
        end
        
        function tg = ReadTimitFeatures(srcDir, id, varargin)
            % Construct textgrid compatible struct from wrd, phn
            %
            %   tg = ReadTimitFeatures(timitDir, id)
            % 
            % Inputs
            %   srcDir      The source directory path of TIMIT feature/label files.
            %   id          A list of TIMIT stim IDs.
            % Output
            %   tg          An array of structs compatible to use with the TextGrid MATLAB package.
            %               Currently it only contains words and phones level.
            % 
            % TODO          Add syllable level to tg.
            %               Read quantity features such as formants.
            % 
            
            p = inputParser();
            p.addParameter('UniformOutput', true, @islogical);
            p.parse(varargin{:});
            isUni = p.Results.UniformOutput;
            
            id = cellstr(id);
            tg = cell(size(id));
            
            for i = 1 : numel(id)
                % Read wrd and phn files as tables
                wrdFile = fullfile(srcDir, id{i} + ".wrd");
                phnFile = fullfile(srcDir, id{i} + ".phn");
                if ~exist(wrdFile, 'file') || ~exist(wrdFile, 'file')
                    warning('Cannot find wrd and/or phn files with stim ID: ''%s''', id{i});
                    continue
                end
                wrdTb = readtable(wrdFile, 'FileType', 'text');
                phnTb = readtable(phnFile, 'FileType', 'text');
                
                % Remove silent phones
                phnTb(strcmp(phnTb.(3), 'h#'), :) = [];
                phnTb.(3) = upper(phnTb.(3));
                
                % Hardcode the audio sampling rate of TIMIT
                fs = 16000;
                
                % Construct TextGrid compatible struct
                w = struct;
                w.name = 'words';
                w.type = 'interval';
                w.T1 = wrdTb.(1) / fs;
                w.T2 = wrdTb.(2) / fs;
                w.Label = wrdTb.(3)';
                
                p = struct;
                p.name = 'phones';
                p.type = 'interval';
                p.T1 = phnTb.(1) / fs;
                p.T2 = phnTb.(2) / fs;
                p.Label = phnTb.(3)';
                
                s = struct;
                s.tier = {w, p};
                s.tmin = w.T1(1);
                s.tmax = w.T2(end);
                
                tg{i} = s;
            end
            
            % Convert to struct array or remain as cell array
            if isUni
                tg = cat(1, tg{:});
            end
        end
        
    end
end

