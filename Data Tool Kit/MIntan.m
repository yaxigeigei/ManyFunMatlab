classdef MIntan
    %MDATA Summary of this class goes here
    %   Detailed explanation goes here
    
    methods(Static)
        function ops = GetOptions(signalList)
            % Get option structures to use with MIntan.ReadPhdFiles
            % 
            %   ops = MIntan.GetOptions()
            %   ops = MIntan.GetOptions(signalList)
            % 
            % Input
            %   signalList          A string or a cell array of strings specifying the signal(s) of interest. 
            %                       Possible options are 'amplifier', 'aux_in', 'adc', 'dig_in'. The default
            %                       is {'amplifier', 'aux_in', 'adc', 'dig_in'}.
            % Output
            %   ops                 Structure(s) of options used to specify the processig in MIntan.ReadRhdFiles. 
            %     signalName        The signal of interest from each element in signalList.
            %     signalFunc        A handle to a function that operates the signal read from each file. The input
            %                       to this function is an array with #channel rows and #sample columns. The 
            %                       output should have the same size as input. (default @(x) x)
            %     downsampleFactor  Folds to downsample the signal (used by MATLAB downsample function under the hood). 
            %     isReturn          Whether or not this signal should be returned in the output of MIntan.ReadRhdFiles. 
            %                       The default is true. 
            %     varBaseName       Variable base name used to store processed signal in the output of MIntan.ReadRhdFiles.
            %                       Time stamps and data will be specified by '_time' and '_data' suffix to the base name. 
            %                       The default is simply the signalName. 
            %     binFilePath       If specified, the signal data will be saved as a binary file. The data type will 
            %                       inherit from the data array. However, you can use signalFunc to change data type. 
            %                       Time stamps are not saved. 
            % 
            
            % Handle user inputs
            if nargin < 1
                signalList = {'amplifier', 'aux_in', 'adc', 'dig_in'};
            end
            signalList = cellstr(signalList);
            
            % Create option structs
            for i = numel(signalList) : -1 : 1
                ops(i) = signalOptions(signalList{i});
            end
            
            function s = signalOptions(signalName)
                % Signal processing
                s.signalName = signalName;
                s.signalFunc = @(x) x;
                s.downsampleFactor = 1;
                
                % Returning data in function output
                s.isReturn = true;
                s.varBaseName = signalName;
                
                % Generation of binary files
                s.binFilePath = '';
            end
        end
        
        function result = ReadRhdFiles(rhdFilePaths, ops)
            % Read Intan RHD2000 files (a fancy wrapper of Intan.read_Intan_RHD2000_file_noUI)
            % 
            %   result = MIntan.ReadRhdFiles()
            %   result = MIntan.ReadRhdFiles(rhdFilePaths)
            %   result = MIntan.ReadRhdFiles(rhdFilePaths, ops)
            % 
            % Inputs
            %   rhdFilePaths        A string or a cell array of strings of RHD2000 file paths. If left empty, 
            %                       a file selection window will show up. 
            %   ops                 Structure(s) returned from MIntan.GetOptions method and optionally
            %                       modified by user to customize loading. See MIntan.GetOptions for details.
            %                       The default is extracting all data with no custom preprocessing. 
            % Output
            %   result              A structure with read/processed data and metadata. 
            %   
            
            % Handles user inputs
            if nargin < 2
                ops = MIntan.GetOptions();
            end
            if nargin < 1 || isempty(rhdFilePaths)
                rhdFilePaths = MBrowse.Files();
            end
            rhdFilePaths = cellstr(rhdFilePaths);
            
            % Preallocation in output structure
            dataCells = cell(numel(rhdFilePaths), numel(ops));
            timeCells = cell(numel(rhdFilePaths), numel(ops));
            
            for k = 1 : numel(rhdFilePaths)
                % Load data from a Intan file
                [notes, frequency_parameters, ...
                    amplifier_channels, amplifier_data, amplifier_time, ...
                    aux_in_channels, aux_in_data, aux_in_time, ...
                    adc_channels, adc_data, adc_time, ...
                    dig_in_channels, dig_in_data, dig_in_time] = ...
                    Intan.read_Intan_RHD2000_file_noUI(rhdFilePaths{k});
                
                amplifier_data = single(amplifier_data);
                aux_in_data = single(aux_in_data);
                adc_data = single(adc_data);
                dig_in_data = uint8(dig_in_data);
                
                % Process requested signals
                for i = 1 : numel(ops)
                    % Cache variables
                    sigName = ops(i).signalName;
                    assert(ismember(sigName, {'amplifier', 'adc', 'dig_in', 'aux_in'}), ...
                        '%s (case-sensitive) is not a valid signal name.', sigName);
                    sigData = eval([sigName '_data;']);
                    sigTime = eval([sigName '_time;']);
                    fprintf('Processing %s data\n', sigName);
                    
                    % Apply user's signal function
                    sigData = ops(i).signalFunc(sigData);
                    
                    % Downsampling
                    ds = ops(i).downsampleFactor;
                    if ds > 1
                        r = mod(numel(sigTime), ds);
                        if r ~= 0
                            warning('The last downsampled value correspond to only %i original samples', r);
                        end
                        sigData = downsample(sigData', ds)';
                        sigTime = downsample(sigTime', ds)';
                    end
                    
                    % Store processed data
                    if ops(i).isReturn
                        dataCells{k,i} = sigData';
                        timeCells{k,i} = sigTime';
                    end
                    
                    % Save processed data to binary file
                    bPath = ops(i).binFilePath;
                    if ~isempty(bPath)
                        bDir = fileparts(bPath);
                        if ~isempty(bDir) && ~exist(bDir, 'dir')
                            mkdir(bDir);
                        end
                        if k == 1 && exist(bPath, 'file')
                            warning('The existing binary file will be replaced. \n%s', bPath);
                            delete(bPath);
                        end
                        fid = fopen(bPath, 'a');
                        fwrite(fid, sigData, class(sigData));
                        fclose(fid);
                    end
                end
            end
            
            % Combine data chunks
            for i = 1 : numel(ops)
                if ops(i).isReturn
                    result.([ops(i).varBaseName '_data']) = cell2mat(dataCells(:,i));
                    result.([ops(i).varBaseName '_time']) = cell2mat(timeCells(:,i));
                end
            end
            
            % Add metadata
            result.info.files = rhdFilePaths;
            result.info.frequency_parameters = frequency_parameters;
            result.info.amplifier_channels = amplifier_channels;
            result.info.aux_in_channels = aux_in_channels;
            result.info.adc_channels = adc_channels;
            result.info.dig_in_channels = dig_in_channels;
            result.info.ops = ops;
        end
    end
end

