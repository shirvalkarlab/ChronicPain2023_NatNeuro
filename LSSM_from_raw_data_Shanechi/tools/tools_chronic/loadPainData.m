% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Omid Sani, Parima Ahmadipour, Maryam Shanechi
%   Shanechi Lab, University of Southern California, 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%loadPainData This function loads the LFP and pain scores from raw data files
%   Inputs: 
%     - (1) fileDir: the directory with the raw CP data
%     - (2) fileName: the name of the raw CP data file
%     - (3) settings: structre with additional settings. For options, see
%           below.
%   Outputs: 
%     - (1) allData: structure with all LFP and pain data 

function allData = loadPainData(dataDir, fileName, settings)

if nargin < 1, dataDir = ''; end
if nargin < 2, fileName = ''; end
if nargin < 3, settings = struct; end

if ~isfield(settings, 'discardFirstNSamples'), settings.discardFirstNSamples = []; end
if ~isfield(settings, 'maxSamples'), settings.maxSamples = []; end
if ~isfield(settings, 'hemi'), settings.hemi = {}; end  % List of hemispheres to load in bilateral
if ~isfield(settings, 'factor'), settings.factor = 'autoNRS'; end
if ~isfield(settings, 'refFilePath'), settings.refFilePath = ''; end

if isempty(dataDir), dataDir = ''; end
if isempty(fileName), fileName = ''; end

filePath = fullfile(dataDir, fileName);
fD = load(filePath);

if isfield(fD, 'LFP')
    dataFieldName = 'LFP';
elseif isfield(fD, 'LFPbl')
    dataFieldName = 'LFPbl';
else
    error('Could not find data field');
end

if isfield(fD, 'LFPmeta')
    metaFieldName = 'LFPmeta';
elseif isfield(fD, 'LFPmetabl')
    metaFieldName = 'LFPmetabl';
else
    error('Could not find data field');
end

chanCodes = fieldnames(fD.(dataFieldName));
chanCodesAll = chanCodes;
if ~isempty(settings.hemi)
    keepInds = contains(chanCodes, settings.hemi);
    chanCodes = chanCodes(keepInds);
end
fprintf('Keeping %d/%d electrodes: %s\n', numel(chanCodes), numel(chanCodesAll), sprintf('%s ', chanCodes{:}));

source = struct('chans', cell(1, numel(chanCodes)), 'ref', cell(1, numel(chanCodes)), 'label', cell(1, numel(chanCodes)));
for ci = 1:length(chanCodes)
    hemiMeta = fD.(metaFieldName);
    metaContactFieldName = [chanCodes{ci},'contact'];
    if contains('left', fieldnames(hemiMeta)) && contains(chanCodes{ci}, 'left')
        hemiMeta = fD.(metaFieldName).left;
        metaContactFieldName = [strrep(chanCodes{ci}, 'left', ''),'contact'];
    elseif contains('right', fieldnames(hemiMeta)) && contains(chanCodes{ci}, 'right')
        hemiMeta = fD.(metaFieldName).right;
        metaContactFieldName = [strrep(chanCodes{ci}, 'right', ''),'contact'];
    end
    [A] = sscanf(hemiMeta.(metaContactFieldName){1}, 'E%dE%d');
    source(ci).chans = A(1);
    refChans = A(2);
    source(ci).ref = struct('refChans', refChans, 'type', 'bipolar');
    source(ci).label = chanCodes{ci};
end

[ measureVals ] = generateMeasureStructPain( fD, settings.factor );

allRaw = [];

sessionCount = size(fD.(dataFieldName).(chanCodes{1}), 2);
for si = 1:sessionCount
    sessionData = [];
    for chi = 1:length(chanCodes)
        chanCode = chanCodes{chi};
        chanData = fD.(dataFieldName).(chanCode);
        sessionData = cat(2, sessionData, chanData(:, si));
    end
    
    hemiMeta = fD.(metaFieldName);
    if contains('left', fieldnames(hemiMeta))
        hemiMeta = fD.(metaFieldName).left;
    elseif contains('right', fieldnames(hemiMeta))
        hemiMeta = fD.(metaFieldName).right;
    end
    
    Fs = hemiMeta.fs(si); % Hz

    x = sessionData;
    try
        [x, ~, removedInd] = removeConstantSamplesAtTheEnd(x, [], [], 1e6 * eps);
    catch
    end
    
    t = linspace(0, length(x)/Fs, length(x))';

    if ~isempty(settings.discardFirstNSamples)
        t = t(settings.discardFirstNSamples:end);
        x = x(settings.discardFirstNSamples:end, :);
    end    
    if ~isempty(settings.maxSamples)
        if numel(x) > settings.maxSamples
            t = t(1:settings.maxSamples);
            x = x(1:settings.maxSamples, :);
        end
    end

    % Find absolute data time
    dSTimeNum = datenum(hemiMeta.time(si));
    tNum = t / (24 * 3600);
    tzone = 0 ;
    time = dSTimeNum + tNum;

    thisRaw = struct;
    thisRaw.data = x;

    thisRaw.time = time;
    thisRaw.Fs = Fs;
    thisRaw.tzone = tzone;
    thisRaw.source = source;

    
    if si == 1
        allRaw = thisRaw;
    else
        allRaw = cat(1, allRaw, thisRaw);
    end
end

allData = struct( ...
    'raw', allRaw, ...
    'source', source, ...
    'measureVals', measureVals, ...
    'subject', strrep(fileName(1:min(5, numel(fileName))), '_', '-'), ...
    'settings', settings ...
);

% Remove noisy segments
dStd = arrayfun(@(r)(mean(std(r.data, 'omitnan'), 'omitnan')), allData.raw);
% figure; histogram(dStd);
    
qTiles = quantile(dStd,[0.25 0.5 0.75]);
noiseSegmenst = find( dStd > (qTiles(3)+10*diff(qTiles([1 3]))) | dStd < (qTiles(1)-10*diff(qTiles([1 3]))) );

if contains(fileName, 'CP1') % Manually noted outliers for CP1
   noiseSegmenst = cat(1, noiseSegmenst(:),[99; 100; 101]); %[99; 100; 101]% Manually detected outlier recording
elseif contains(fileName, 'CP2') % Manually noted outliers for CP2
   noiseSegmenst = cat(1, noiseSegmenst(:), [137]); %[137]% Manually detected outlier recording  % 137 pain score outlier, no need to remove
elseif contains(fileName, 'CP3') % Manually noted outliers for CP3
   noiseSegmenst = cat(1, noiseSegmenst(:), []); % Manually detected outlier recording
elseif contains(fileName, 'CP4') % Manually noted outliers for CP4
    noiseSegmenst = cat(1, noiseSegmenst(:), [280; 306; 372; 435; 456]); % Manually detected outlier recording % 309: no problem at all, no need to remove, % 456: autoNRS extreme outlier
end
noiseSegmenst = noiseSegmenst( noiseSegmenst <= size(allData.raw) );
if ~isempty(noiseSegmenst)
    fprintf('%d segments will be discarded due to large signal noise: %s\n', numel(noiseSegmenst), sprintf('%d ', noiseSegmenst));
end

try
    if ~isempty(settings.refFilePath)
        if ~exist(settings.refFilePath, 'file')
            error('Could not find reference file here: %s\n', settings.refFilePath);
        end
        noiseSegmenst = [];
        fprintf('Loading reference data file from %s\n', settings.refFilePath);
        fD = jsondecode(fileread(settings.refFilePath));
        %{
        fNames = fieldnames(fD); fNames = fNames(contains(fNames, {'acc', 'ofc'}) & ~contains(fNames, {'var'}) & contains(fNames, {'gamma'})); 
        figure('Color', 'w'); 
        ax = subplot(4,1,1:2); hold(ax, 'on'); box(ax, 'off');
        allFeats = [];
        for fi = 1:numel(fNames), allFeats = cat(2, allFeats, [fD.(fNames{fi})]'); end
        h = plot(ax, allFeats); hm = plot(mean(allFeats, 2), 'k', 'LineWidth', 3);
        legend([h; hm], fNames{:}, 'Mean'); 
        title(ax, fileName(1:3)); ylabel(ax, 'Feature value');
%         xlabel(ax, 'Rec'); 
        ax = subplot(4,1,3); hold(ax, 'on'); box(ax, 'off');
        plot(ax, arrayfun(@(r)(size(r.data, 1)), allData.raw)/allData.raw(1).Fs, 'LineWidth', 2);
        ylabel(ax, 'Duration (seconds)');
%         title(ax, fileName(1:3)); xlabel(ax, 'Recording number'); 
        ax = subplot(4,1,4); hold(ax, 'on'); box(ax, 'off');
        P = [allData.measureVals.value]' >= median([allData.measureVals.value]);
%         P = [fD.nrs]' >= median([fD.nrs]);
        plot(ax, P, 'LineWidth', 2);
        xlabel(ax, 'Recording number'); ylabel(ax, 'Pain class');
        title(ax, fileName(1:3)); 
        ri = 50; %numel(allData.raw)-1;
        t = (1:size(allData.raw(ri).data, 1))/allData.raw(ri).Fs;
        figure; ax = axes; hold(ax, 'on'); box(ax, 'off');
        h = plot(allData.raw(ri).data); 
        legend(ax, h, {allData.raw(ri).source.label});
        ylabel(ax, 'Voltage'); xlabel(ax, 'Sample'); % xlabel('Time [s]');
        title(sprintf('%s - Recording #%d', fileName(1:3), ri));
        %}
        if strcmpi(settings.factor, 'relativeVAS')
            fieldName = 'relativevas';
        elseif strcmpi(settings.factor, 'painVAS')
            fieldName = 'vas';
        elseif strcmpi(settings.factor, 'autoNRS')
            fieldName = 'nrs';
        elseif strcmpi(settings.factor, 'MPQsum')
            fieldName = 'mpq';
        elseif strcmpi(settings.factor, 'unpleasantNRS')
            fieldName = 'unpnrs';
        elseif strcmpi(settings.factor, 'unpleasantVAS')
            fieldName = 'unpvas';
        end
        hasScore = arrayfun( @(v)(~isempty(v.(fieldName))), fD);
        usedDataTime = datenum(datetime({fD(hasScore).righttimes}'));
        usedDataVals = [fD(hasScore).(fieldName)]';

        allDataTime = [allData.measureVals.time]';
        
        if isempty(usedDataTime)
            toBeDiscardedInd = 1:numel(allDataTime);
        else
            closestTimeDist = arrayfun(@(thisT)( min(abs(thisT-usedDataTime)) ), allDataTime);
            closestTimeDistInSecs = closestTimeDist*24*3600;
            toBeDiscardedInd = find( closestTimeDistInSecs > 3660 );
        end
        
        fprintf('Discarding %d/%d data points because they were not in the reference file\n', numel(toBeDiscardedInd), numel(allDataTime));
        noiseSegmenst = cat(1, noiseSegmenst, toBeDiscardedInd);
    end
catch ME
    fprintErrorObj( ME );
end


allData.raw(noiseSegmenst) = [];
allData.measureVals(noiseSegmenst) = [];

%{
figure; ax = axes; hold(ax, 'on');
h1 = scatter(usedDataTime, usedDataVals,'r');
h2 = scatter([allData.measureVals.time]', [allData.measureVals.value]','b');
% legend([h1, h2], {'UCSF', 'USC'});
markTimeUnitsOnAxes(ax, [0 0 7 0 0 0], 'mm/dd/yy'); 
title(ax, allData.subject);
length(usedDataTime)
length([allData.measureVals.time])
length(noiseSegmenst)

%}

end