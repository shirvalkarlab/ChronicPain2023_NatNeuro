% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Omid Sani, Yuxiao Yang, Maryam Shanechi
%   Shanechi Lab, University of Southern California, 2017
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%doPreprocessing Perform preprocessing on ECoG data
%   Inputs:
%     - (1) data: raw ECoG data. columns are channels.
%     - (2) secTime: column vector representing time of samples.
%     - (3) settings (optional): A structure with additionalsettings. For a 
%               list of supported fields and default values, refer to the 
%               beginning of the source code for the function
%     - (4) Fs (optional, default: settings.loadFs): Sampling rate of 
%               "data". If not provided, will be set to settings.loadFs
%     - (5) chansToLoad (optional, default: []): channel numbers of the
%               channels in the data file. If empty will assume 1 to m
%               where m is the total number of columns in "data"
%     - (6) filters (optional, default: struct): structure with all filter
%               objects. Provide to continue causal filtering from previous
%               data segment.
%   Outputs:
%     - (1) data: preprocessed data
%     - (2) time: time vector for the preprocesed data
%     - (3) source: source structure array for the output data. Contains
%           channel numbers and reference info for each output column
%     - (4) settings: updated settings struct with final values
%     - (5) filters: updated filters struct with final state of the
%               filters. Provide this as input with the next data chunk to
%               continue causal filtering. 
%
%   Usage example:
%       [ data, time, source, settings, filters ] = doPreprocessing(secAllChanData, secTime, settings, Fs, chansToLoad, dsSavePath, filters );

function [ data, time, source, settings, filters ] = doPreprocessing( secAllChanData, secTime, settings, Fs, chansToLoad, filters )

if nargin < 3 || ~isstruct(settings), settings = struct; end

if ~isfield(settings, 'desiredFs'), settings.desiredFs = 250; end % [Hz] Downsamle to this before extracting features
if ~isfield(settings, 'lineFs'), settings.lineFs = 60; end % [Hz] Line noise frequency [60Hz in USA]
if ~isfield(settings, 'removeNoisy'), settings.removeNoisy = true; end % If true, removes channels and times that are marked as noisy
if ~isfield(settings, 'doHPFiltering'), settings.doHPFiltering = true; end % If true, does high-pass filtering
if ~isfield(settings, 'doBPFiltering'), settings.doBPFiltering = false; end % If true, does band-pass filtering
if ~isfield(settings, 'doLNFiltering'), settings.doLNFiltering = true; end % If true, does line-noise notch filtering
if ~isfield(settings, 'doLNHarmonicsFiltering'), settings.doLNHarmonicsFiltering = length(settings.lineFs)==1; end % If true, does line-noise notch filtering
if ~isfield(settings, 'causal'), settings.causal = false; end % If true, will apply filters causally
if ~isfield(settings, 'adaptRefToExcludeNoise'), settings.adaptRefToExcludeNoise = false; end % If true, will exclude channels from ref if they are noisy
if ~isfield(settings, 'noiseAreas'), settings.noiseAreas = []; end % If not empty, will consider this as the noise area

if ~isfield(settings, 'doArtifactRemoval'), settings.doArtifactRemoval = false; end % If true, will call the artifact removal code
if ~isfield(settings, 'artifactRemovalSettings'), settings.artifactRemovalSettings = struct; end % Structure with settings for artifact removal

if nargin < 4, Fs = settings.loadFs; end
if nargin < 5, chansToLoad = []; end
if nargin < 6, filters = struct; end

if isempty(chansToLoad), chansToLoad = 1:size(secAllChanData, 2); end

sourceStruct = struct('chans', num2cell(chansToLoad(:)'), 'ref', cell(1, length(chansToLoad)));

if settings.doArtifactRemoval
    [secAllChanData, secTime, settings.artifactRemovalSettings] = ...
            doArtifactRemoval(secAllChanData, secTime, settings.artifactRemovalSettings, Fs, sourceStruct);
end

% Downsample if required
% Initial decimation before preprocessing
decimationRatio = 2^round(log2(Fs/settings.desiredFs)); % Make sure decimation ratio is a power of 2
newFs = Fs/decimationRatio;

if isfield(settings, 'newFs')&&~isequal(newFs, settings.newFs)
%     fprintf('WARNING: Did not expect sampling rate of %.2fHz! Used to be different (%.2fHz) in sys id!\n', newFs, settings.newFs);
end
settings.newFs = newFs;

if ~isempty(filters), settings.usedFilters = filters; end

if (decimationRatio > 1) % Do initial decimation
    NBU = size(secAllChanData, 1);
    [ secAllChanData, secTime, settings, filters ] = doDownSampling( secAllChanData, secTime, settings, Fs, decimationRatio, filters );
    actualNewFs = Fs/NBU*size(secAllChanData, 1);
    settings.newFs = actualNewFs; % To fix any DS sample count rounding issues 
end

% Prepare filter if not already ready
[settings, filters] = initPreProcessingFilters(settings, filters, settings.newFs);

if isfield(filters, 'DFilt') && ~isempty(filters.DFilt) % General digital filters
    for i = 1:length(filters.DFilt)
        secAllChanData = filters.DFilt(i).apply(secAllChanData);
        if ~isequal(filters.DFilt(i).Fs, settings.newFs), warning('featExt:FsMismatch', 'Fs mismatch between data and digital filter #%d!\n', i); end
    end
end
if (settings.doHPFiltering)&&isfield(filters, 'HPFilt')
    % Do high pass filtering
    secAllChanData = filters.HPFilt.apply(secAllChanData);
    if ~isequal(filters.HPFilt.Fs, settings.newFs), warning('featExt:FsMismatch', 'Fs mismatch between data and high-pass filter!\n'); end
end
if (settings.doBPFiltering)&&isfield(filters, 'BPFilt')
    % Do band pass filtering
    secAllChanData = filters.BPFilt.apply(secAllChanData);
    if ~isequal(filters.BPFilt.Fs, settings.newFs), warning('featExt:FsMismatch', 'Fs mismatch between data and band-pass filter!\n'); end
end
if (settings.doLNFiltering)&&isfield(filters, 'BSFilt')
    % Remove line noise by Band Stop filtering
    for lnfi = 1:length(filters.BSFilt)
        secAllChanData = filters.BSFilt(lnfi).apply(secAllChanData);
        if ~isequal(filters.BSFilt(lnfi).Fs, settings.newFs), warning('featExt:FsMismatch', 'Fs mismatch between data and band-stop filter!\n'); end
    end
end

if settings.adaptRefToExcludeNoise && ~isempty(settings.noiseAreas) && ...
   (isfield(settings, 'refMap')&&~isempty(settings.refMap))
    % Make a copy and set noisy areas to NaN
    dataStruct = struct('data', secAllChanData, 'time', secTime);
    dataStruct = setFeatValuesWithinGivenAreasTo( dataStruct, settings.noiseAreas, NaN, sourceStruct, false );
    [secAllChanData, sourceStruct] = applyReferencingBasedOnRefMap(secAllChanData, settings.refMap, sourceStruct, dataStruct.data);
else
    
    if (isfield(settings, 'refMat')&&~isempty(settings.refMat))&&(size(secAllChanData, 2) == size(settings.refMat, 1))
        [secAllChanData, sourceStruct] = applyReferencingBasedOnRefMatrix(secAllChanData, settings.refMat, sourceStruct);
    elseif (isfield(settings, 'refMap')&&~isempty(settings.refMap))
        % Do the rereferencing
        %{
        % Old
        secDataStruct = struct('time', secTime, 'data', secAllChanData);
        [reRefData, newSourceStruct] = applyReferencingBasedOnRefMap(secDataStruct.data, settings.refMap, sourceStruct, secDataStruct.data);
        secAllChanData = reRefData;
        sourceStruct = newSourceStruct;
        %}
        fprintf('Building refMat from RefMap...\n');
        [refMat, newSourceStruct] = generateReferencingMatrixBasedOnRefMap(settings.refMap, sourceStruct);
        [secAllChanData, sourceStruct] = applyReferencingBasedOnRefMatrix(secAllChanData, refMat, sourceStruct);
        fprintf('Saving refMat for future use...\n');
        settings.refMat = refMat; % Save for future use
    end

end

data = secAllChanData;
time = secTime;
source = sourceStruct;


end

function [settings, filters] = initPreProcessingFilters(settings, filters, Fs)

if ~isfield(settings, 'usedFilters'), settings.usedFilters = struct; end

if isfield(settings, 'DFilterSpecs') && numel(settings.DFilterSpecs) > 0 % General digital filters
    for i = 1:length(settings.DFilterSpecs)
        if ~iscell(settings.DFilterSpecs)
            dFiltSpecs = settings.DFilterSpecs(i);
        else
            dFiltSpecs = settings.DFilterSpecs{i};
        end
        
        
        if ~isfield(settings.usedFilters, 'DFilt') || (length(settings.usedFilters.DFilt) < length(settings.DFilterSpecs)) || ...
           (~isfield(settings.usedFilters.DFilt(i), 'Fs') && ~isprop(settings.usedFilters.DFilt(i), 'Fs')) || ...
           (settings.usedFilters.DFilt(i).Fs ~= Fs)
            
            if strcmpi(dFiltSpecs.shape, 'HP')
                
                fprintf('Designing a high-pass filter in [ %s]Hz...\n', sprintf('%.1f ', dFiltSpecs.Fc));
                [b, a] = designFilterHPIIRButter(Fs, dFiltSpecs.Fc, dFiltSpecs.order);
                thisDFilt = struct('b', b, 'a', a, 'Fs', Fs);

            elseif strcmpi(dFiltSpecs.shape, 'BP')
                
                fprintf('Designing a band-pass filter in [ %s]Hz...\n', sprintf('%.1f ', dFiltSpecs.Fc));
                [b, a] = designFilterBPIIRButter(Fs, dFiltSpecs.Fc, dFiltSpecs.order);
                thisDFilt = struct('b', b, 'a', a, 'Fs', Fs);

            elseif strcmpi(dFiltSpecs.shape, 'BS')
                
                fprintf('Designing a band-stop filter in [ %s]Hz...\n', sprintf('%.1f ', dFiltSpecs.Fc));
                [b, a] = designFilterBSIIRButter(Fs, dFiltSpecs.Fc, dFiltSpecs.order); 
                thisDFilt = struct('b', b, 'a', a, 'Fs', Fs);
                
            elseif strcmpi(dFiltSpecs.shape, 'Notch')
                
                fprintf('Designing a notch filter in [ %s]Hz...\n', sprintf('%.1f ', dFiltSpecs.Fc));
                [b, a] = designFilterNotchIIR(Fs, dFiltSpecs.Fc, dFiltSpecs.Fc*dFiltSpecs.BW); 
                thisDFilt = struct('b', b, 'a', a, 'Fs', Fs);
                
            end     
        else
            thisDFilt = settings.usedFilters.DFilt(i);
        end
        if i == 1
            DFilt = thisDFilt;
        else
            DFilt = cat(1, DFilt, thisDFilt);
        end
    end
    
    % Convert to filter objects
    if ~isa(DFilt, 'dFilter')
        DFiltBU = DFilt;
        for lnfi = 1:length(DFiltBU)
            thisDFilt = dFilter(DFiltBU(lnfi).b, DFiltBU(lnfi).a, settings.causal, Fs);
            if lnfi == 1
                DFilt = thisDFilt;
            else
                DFilt = cat(1, DFilt, thisDFilt);
            end
        end
    end
    settings.usedFilters.DFilt = DFilt;
    filters = settings.usedFilters;
end
    
if (settings.doHPFiltering)
    if ~isfield(settings, 'HPFilterSpecs'), settings.HPFilterSpecs = struct('Fc', 1, 'type', 'IIR', 'order', 2); end % If high-pass filter specs
    
    if ~isfield(settings.usedFilters, 'HPFilt') || ...
       (~isfield(settings.usedFilters.HPFilt, 'Fs') && ~isprop(settings.usedFilters.HPFilt, 'Fs')) || ...
       (settings.usedFilters.HPFilt.Fs ~= Fs)
        % HPFilt0 = designFilter('HPIIR', Fs); % OLD function
        fprintf('Designing high-pass filter... (Fs = %.1fHz)\n', Fs);
        [b, a] = designFilterHPIIRButter(Fs, settings.HPFilterSpecs.Fc, settings.HPFilterSpecs.order);
        HPFilt = struct('b', b, 'a', a, 'Fs', Fs);
    else
        HPFilt = settings.usedFilters.HPFilt;
    end
    % Convert to filter objects
    if ~isa(HPFilt, 'dFilter')
        HPFilt = dFilter(HPFilt.b, HPFilt.a, settings.causal, Fs);
    end
    settings.usedFilters.HPFilt = HPFilt;
    filters = settings.usedFilters;
end
if (settings.doBPFiltering)
    if ~isfield(settings, 'BPFilterSpecs'), settings.BPFilterSpecs = struct('Fc', [1 Fs/2], 'type', 'IIR', 'order', 4); end % If band-pass filter specs
    settings.BPFilterSpecs.Fc(settings.BPFilterSpecs.Fc > Fs/2) = Fs/2;
    
    if ~isfield(settings.usedFilters, 'BPFilt') || ...
       (~isfield(settings.usedFilters.BPFilt, 'Fs') && ~isprop(settings.usedFilters.BPFilt, 'Fs')) || ...
       (settings.usedFilters.BPFilt.Fs ~= Fs)
        fprintf('Designing band-pass filter... (Fs = %.1fHz)\n', Fs);
        [b, a] = designFilterBPIIRButter(Fs, settings.BPFilterSpecs.Fc, settings.BPFilterSpecs.order);
        BPFilt = struct('b', b, 'a', a, 'Fs', Fs);
    else
        BPFilt = settings.usedFilters.BPFilt;
    end
    % Convert to filter objects
    if ~isa(BPFilt, 'dFilter')
        BPFilt = dFilter(BPFilt.b, BPFilt.a, settings.causal, Fs);
    end
    settings.usedFilters.BPFilt = BPFilt;
    filters = settings.usedFilters;
end
if (settings.doLNFiltering)
    if ~isfield(settings, 'BSFilterSpecs'), settings.BSFilterSpecs = struct('Fc', settings.lineFs(1)+[-1 1], 'type', 'IIR', 'order', 4); end % If band-pass filter specs

    lineFreqs = settings.lineFs;
    if (settings.doLNHarmonicsFiltering)
        lineFreqs = settings.lineFs*(1:floor((Fs/2)/settings.lineFs));
    end
    if (Fs < settings.lineFs)
        lineFreqs = [];
        settings.usedFilters.BSFilt = [];
    end

    if ~isfield(settings.usedFilters, 'BSFilt') || (length(settings.usedFilters.BSFilt)~=length(lineFreqs)) || ...
       (~isfield(settings.usedFilters.BSFilt(1), 'Fs') && ~isprop(settings.usedFilters.BSFilt(1), 'Fs')) || ...
       (settings.usedFilters.BSFilt(1).Fs ~= Fs)
        for lnfi = 1:length(lineFreqs)
            if ~isfield(settings.BSFilterSpecs, 'BW') % Band stop
                fprintf('Designing band-stop filter #%d/%d at %.1fHz... (Fs = %.1fHz)\n', lnfi, length(lineFreqs), lineFreqs(lnfi), Fs);
                % thisBSFilt0 = designFilter('BSIIR', Fs, lineFreqs(lnfi)); % OLD function
                [b, a] = designFilterBSIIRButter(Fs, lineFreqs(lnfi)+(settings.BSFilterSpecs.Fc-mean(settings.BSFilterSpecs.Fc)), settings.BSFilterSpecs.order); 
                thisBSFilt = struct('b', b, 'a', a, 'Fs', Fs);
            else % Notch
                fprintf('Designing Notch filter #%d/%d at %.1fHz... (Fs = %.1fHz)\n', lnfi, length(lineFreqs), lineFreqs(lnfi), Fs);
                [b, a] = designFilterNotchIIR(Fs, lineFreqs(lnfi), lineFreqs(lnfi)*settings.BSFilterSpecs.BW); 
                thisBSFilt = struct('b', b, 'a', a, 'Fs', Fs);
            end
            if lnfi == 1
                BSFilt = thisBSFilt;
            else
                BSFilt = cat(1, BSFilt, thisBSFilt);
            end
        end
    else
        BSFilt = settings.usedFilters.BSFilt;
    end
    % Convert to filter objects
    if ~isa(BSFilt, 'dFilter')
        BSFiltBU = BSFilt;
        for lnfi = 1:length(BSFiltBU)
            thisBSFilt = dFilter(BSFiltBU(lnfi).b, BSFiltBU(lnfi).a, settings.causal, Fs);
            if lnfi == 1
                BSFilt = thisBSFilt;
            else
                BSFilt = cat(1, BSFilt, thisBSFilt);
            end
        end
    end
    settings.usedFilters.BSFilt = BSFilt;
    filters = settings.usedFilters;
end
 
end