% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Omid Sani, Parima Ahmadipour, Maryam Shanechi 
%   Shanechi Lab, University of Southern California, 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%decodePain_chronicPain This function runs the cross-validated model 
% fitting with the chronic 
% pain data 
%   Inputs: 
%     - (1) fileDir: the directory with the raw CP data
%     - (2) fileName: the name of the raw CP data file
%     - (3) refFilePath: the path to final json files to keep the same datapoints
%     - (4) measureName: name of the pain measure
%     - (5) fitSettings: name of the pain additional modeling settings such as
%             number of permutations, etc
%   Outputs:
%     - (1) results: structure with all cross-validation results

function results = decodePain_chronicPain(fileDir, fileName, refFilePath, measureName, fitSettings)

if nargin < 5, fitSettings = struct; end
if ~isfield(fitSettings, 'useLSSM'), fitSettings.useLSSM = true; end % Set to true to fit LSSM model and use its states for prediction

allData = loadPainData(fileDir, fileName, struct( ...
    'factor', measureName, ...
    'hemi', {{}}, ...  % List of hemispheres to load (leave empty to load all)
    'refFilePath', refFilePath ...
));
%%
% Find non-NaN samples
isNonNaN = ~isnan([allData.measureVals.value]');
scoreVals=[allData.measureVals(isNonNaN).value]';

results = [];
if isempty(scoreVals) || numel(unique(scoreVals)) == 1
    fprintf('%s does not have %s scores.... skipping...\n', allData.subject, measureName);
    return
end
%%
% Frequency bands to extract power from 
bands  = [1 4; 4 8; 8 12; 12 30; 30 70; 70 100];
% bands = [1 4; 4 8; 8 12; 12 20; 20 30; 30 55; 65 100];  % Original USC bands

% Preprocessing settings
Fs = allData.raw(1).Fs;
settings = struct( ...
    'desiredFs', Fs, ...
    'causal', true, ...
    'doHPFiltering', false, ...
    'doBPFiltering', false, ...
    'doLNFiltering', false, ...
    'DFilterSpecs', [] ...
);
% Uncomment the following to add original USC preprocessing for all subjects:
% settings.DFilterSpecs = cat(1, settings.DFilterSpecs, ...
%     [struct('shape', 'BS', 'Fc', [ 66  70], 'type', 'IIR', 'order', 4); ...
%      struct('shape', 'BS', 'Fc', [ 81  85], 'type', 'IIR', 'order', 4); ...
%      struct('shape', 'BS', 'Fc', [103 107], 'type', 'IIR', 'order', 4); ...
%      struct('shape', 'BS', 'Fc', [103 107], 'type', 'IIR', 'order', 4); ...
%      struct('shape', 'BS', 'Fc', [126 130], 'type', 'IIR', 'order', 4); ...
%      struct('shape', 'HP', 'Fc', [1], 'type', 'IIR', 'order', 4); ...
%     ]);

settingsBU = settings;

sourceDependentSettings = getSubjectSpecificPreprocessingSettings(fileName);
% sourceDependentSettings = []; % Uncomment to use original USC preprocessing

keepSamples = 60 * allData.raw(1).Fs; % Keep first 30s (some recordings are longer in CP1)
NOkRawSamples = arrayfun(@(r)( size(r.data,1) ), allData.raw);

if contains(fileName, 'CP1')
    stTime = arrayfun( @(r)( r.time(1) ), allData.raw );
    badInds = find(stTime >= datenum([2018,4,3,16,0,0]));
    NOkRawSamples(badInds) = min(NOkRawSamples(badInds), floor(30*allData.raw(1).Fs));  % To do: only use 29.86s
end


regionsToInclude = {'rightofc', 'leftofc'};  % Specifify regions to keep: 'rightacc', 'rightofc', 'leftacc', 'leftofc'
chansToInclude = find(contains({allData.source.label}, regionsToInclude)); 
chansToInclude = 1:size(allData.raw(1).data, 2); % Uncomment to use regionsToInclude

allChanPredictors = [];
for ci = chansToInclude(:)' % Loop over channels (ofc, acc, etc)
    allRaw = zeros(max(NOkRawSamples), numel(allData.raw));
    for ri = 1:numel(allData.raw) % Collect recordings
        allRaw(1:NOkRawSamples(ri), ri) = allData.raw(ri).data(1:NOkRawSamples(ri), ci);
    end
    settings = settingsBU;
    % Add any channel specific settings
    if ~isempty(sourceDependentSettings)
        chanSettingsInd = find(ismember({sourceDependentSettings.label}, allData.source(ci).label));
        if ~isempty(chanSettingsInd)
            chanDFilterSpecs = sourceDependentSettings(chanSettingsInd).DFilterSpecs;
            settings.DFilterSpecs = cat(1, settings.DFilterSpecs, chanDFilterSpecs);
        end
    end
    % Do preprocessing
    [ dataPrepro, time, source, settings, filters ] = doPreprocessing(allRaw, [], settings, Fs, [], [] );
    % Extract power features
    featSettings = struct('windowSize', 1, 'method', 'pwelch', 'pwWindow', []);
    chanPredictors = getBandPowers(bands, dataPrepro, Fs, featSettings); % Compute powers
    NOkFeatSamples = floor(NOkRawSamples / (Fs*featSettings.windowSize));
    for ri = 1:numel(NOkFeatSamples)
        chanPredictors(ri, :, (NOkFeatSamples(ri)+1):end) = nan;
    end
    allChanPredictors = cat(2, allChanPredictors, chanPredictors);
end
% figure; plot(mean(allChanPredictors, 3, 'omitnan'));

%%
if fitSettings.useLSSM 
    allDataCat = reshape(permute(allChanPredictors, [3 1 2]), [size(allChanPredictors,1)*size(allChanPredictors, 3), size(allChanPredictors, 2)]);
    allDataCat(any(isnan(allDataCat), 2), :) = [];
    [allDataCat,mu,sigma] = zscore(allDataCat);
    horizon = 10;
    % Determine latent state dimension using AIC
    nxVals = 1:20;
    AIC = nan(size(nxVals));
    for nxi = 1:numel(nxVals) % Model order (i.e. state dimension)
        nx = nxVals(nxi);
        [A,B,C,D,K] = subid(allDataCat, [], horizon, nx, [], [], 1);  % Subspace sys id
        [~, yHat] = kalmanFilter(allDataCat, A, C, K);      % Kalman filter to get states
        AIC(nxi) = calcAIC(allDataCat, yHat, (2*size(allDataCat, 2) + 1)*nx);
    end
    [~, selectedNxInd] = min(AIC);
    nx = nxVals(selectedNxInd);
    fprintf('LSSM state dimensions: %d\n', nx);
    % Extract latent states
    [A,B,C,D,K] = subid(allDataCat, [], horizon, nx, [], [], 1);  % Subspace sys id
    predictors = [];
    for ri = 1:size(allChanPredictors, 1)
        y = permute(allChanPredictors(ri, :, :), [3 2 1]);
        y = (y-mu)./repmat(sigma, [size(y, 1), 1]);  
        xHat = kalmanFilter(y, A, C, K);      % Kalman filter to get states
        predictors = cat(1, predictors, xHat(NOkFeatSamples(ri), :));
    end
else
    predictors = mean(allChanPredictors, 3, 'omitnan');
end

% z-score the predictors
predictors = zscore(predictors);

%%
X = predictors(isNonNaN,:);

settings = struct;
settings.numFolds = [];             % Cross-validation folds (leave empty for leave-one-out)
settings.predictor = 'LDA';         % Can be 'logistic' or 'SVM' or 'ridge' or 'LDA'
if fitSettings.useLSSM
    settings.predictor = 'LDA';  % Comment if you want to use LSSM with other regressors
end
settings.SVMSearchPoints = 50;      % For SVM, points to search to select hyperparameters
settings.titleHead = sprintf('%s %s (N=%d) %s\n', allData.subject, measureName, size(X, 1), settings.predictor);
settings.trainMinDistFromTest = 4; % Set to some number greater than 1, to leave data around the test from training

fNames = fieldnames(fitSettings);
for fi = 1:numel(fNames)
    settings.(fNames{fi}) = fitSettings.(fNames{fi});
end

if ismember(settings.predictor, {'logistic', 'SVM', 'LDA'})
    Y = scoreVals >= median(scoreVals);  
else
    Y = scoreVals;
    Y = zscore(Y);
end
    
performCVedModelFit(X,Y, settings);   % Cross-validated classification
[~, ~, results] = performPermTest(X,Y, settings);       % Find chance level with permutation

end
