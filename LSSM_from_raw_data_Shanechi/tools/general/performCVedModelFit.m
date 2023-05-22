% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Omid Sani, Maryam Shanechi
%  Shanechi Lab, University of Southern California, 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%performCVedModelFit Performs cross-validated model fitting
%   Inputs:
%     - (1) X: predictors
%     - (2) Y: labels
%     - (4) settings: struct with additional settings. For details, see
%   Outputs:
%     - (1) ax
%   Usage example:
%       [foldModels, foldLabels, foldACCs, foldProbs, stats, foldIndexes] = performCVedModelFit(X, Y, settings)

function [foldModels, foldPreds, stats, foldIndexes] = performCVedModelFit(X, Y, settings)

if nargin < 3, settings = struct; end
if ~isfield(settings, 'numFolds'), settings.numFolds = 10; end  % If empty, will do leave-one-out CV
if ~isfield(settings, 'predictor'), settings.predictor = 'SVM'; end
if ~isfield(settings, 'SVMKernelScaleRange'), settings.SVMKernelScaleRange = [1e-3,1e3]; end
if ~isfield(settings, 'SVMBoxConstraintRange'), settings.SVMBoxConstraintRange = [1e-3,1e3]; end
if ~isfield(settings, 'SVMSearchPoints'), settings.SVMSearchPoints = 30; end
if ~isfield(settings, 'SVMShowPlots'), settings.SVMShowPlots = false; end
if ~isfield(settings, 'verbose'), settings.verbose = 0; end
if ~isfield(settings, 'savePath'), settings.savePath = ''; end
if ~isfield(settings, 'perfMeasure')
    if ismember(settings.predictor, {'SVM','logistic','LDA'})
        settings.perfMeasure = 'AUC';
    else
        settings.perfMeasure = 'R2';
    end
end

if ~isequal(numel(Y), size(X, 1)), error('Expected same number of rows in X and Y!\n'); end

if isempty(settings.numFolds), settings.numFolds = numel(Y); end

% Convert class labels into 0 and 1
if ismember(settings.predictor, {'SVM','logistic'})
    [classNames, ~, classInds] = unique(Y);
    YBU = Y;
    Y = classInds - 1;
else
    classNames = {};
end
if isfield(settings,'CVFoldInds')
    foldIndexes = settings.CVFoldInds;
else
    foldIndexes = genCVFoldIndexes( numel(Y), settings.numFolds, false );
end

for i=1:settings.numFolds
    isTest  =  ismember((1:numel(Y))', foldIndexes{i});
    if settings.trainMinDistFromTest > 1
        minTestDist = arrayfun(@(di)( min(abs(foldIndexes{i}-di)) ), (1:numel(Y))');
        isCloseToTest = isTest | (minTestDist < settings.trainMinDistFromTest);
    else
        isCloseToTest = isTest;
    end
    isTrain = ~isCloseToTest;
    % figure; plot([isTrain, isTest]);
    
    trainData = X(isTrain, :);
    trainTrueOutput = Y(isTrain);
    testData = X(isTest, :);
    testLabel = Y(isTest);
    
    if strcmpi(settings.predictor, 'SVM')
        % Train
        Params = hyperparameters('fitcsvm',trainData,trainTrueOutput);
        Params(1, 1).Range = settings.SVMBoxConstraintRange;
        Params(2, 1).Range = settings.SVMKernelScaleRange;
        foldModels{i} = fitcsvm(trainData,trainTrueOutput,'KernelFunction','rbf','Standardize','on' ...
            ,'OptimizeHyperparameters',Params ...
             ,'HyperparameterOptimizationOptions',struct('ShowPlots', settings.SVMShowPlots, 'MaxObjectiveEvaluations', settings.SVMSearchPoints, 'Verbose', settings.verbose));
        % Test
        [~,score{i}] = predict(foldModels{i},testData);
        foldPreds{i} = sigmoid(score{i}(:, end));
    elseif strcmpi(settings.predictor, 'logistic')
        % Train
        linkFunc = 'logit';
        [foldModels{i},dev,stats] = glmfit(trainData,trainTrueOutput,'binomial','link', linkFunc,'constant','on'); % Logistic regression
        % Test
        foldPreds{i} = glmval(foldModels{i}, testData,linkFunc); % Logistic regression
    elseif strcmpi(settings.predictor, 'LDA')
        [~,shrinkage] = cov1para(trainData);
        foldModels{i} = fitcdiscr(trainData, trainTrueOutput, 'DiscrimType', 'linear', 'prior', 'uniform', 'Gamma', shrinkage);
        [label, score, cost] = predict( foldModels{i}, testData );
        posClass = 1;
        if label(1) == posClass
            [~, posClassCol] = max(score(1, :));
        else
            [~, posClassCol] = min(score(1, :));
        end
        foldPreds{i} = score(:, posClassCol);
    elseif strcmpi(settings.predictor, 'ridge')
        lambda = selectRidgeHyperparamWithInnerCV(trainData, trainTrueOutput);
        foldModels{i} = ridge(trainTrueOutput, trainData, lambda, 0);
        foldPreds{i} = [ones(size(testData, 1), 1), testData] * foldModels{i};
    end
end

% Combine predictions from all folds to evaluate CVed prediction
trueVals = Y;
preds = cell2mat(foldPreds(:));

if strcmpi(settings.perfMeasure, 'AUC')
    [~,~,~,perf] = perfcurve(trueVals, preds, 1);    
else
    err = preds - trueVals;
    perf = ones(1, size(trueVals, 2)) - mean(err.^2, 1) ./ mean((trueVals-repmat(mean(trueVals, 1), [size(trueVals, 1) 1, 1])).^2);
end

stats = struct;
stats.perfMeasure = settings.perfMeasure;
stats.perf = perf;
stats.classNames = classNames;
stats.allTrueVals = trueVals;
stats.allPreds = preds;

if nargout < 1 || ~isempty(settings.savePath)
    if strcmpi(settings.perfMeasure, 'AUC')
        ax = plotAUC(trueVals, preds, 1, settings);
        figH = ax.Parent;
    else
        figH = figure; ax = axes; 
        scatter(ax, trueVals, preds, 'Marker', 'x', 'MarkerEdgeColor', 'k');
        title(ax, sprintf('%s = %.2g', stats.perfMeasure, stats.perf));
        h = lsline(ax);
    end
    if ~isempty(settings.savePath)
        fprintf('Saving figure as %s\n', settings.savePath);
        saveas(figH, settings.savePath);
    end
end

end

function y = sigmoid(x)

    y = 1 ./ (1+exp(-x));

end

function k = selectRidgeHyperparamWithInnerCV(XTr, YTr, settings)
    if nargin < 3, settings = struct; end
    if ~isfield(settings, 'lambdaVals')
        settings.lambdaVals = logspace(-4,5, 100);
    end

    YTrHatForKs = zeros(length(YTr), length(settings.lambdaVals));
    % Inner cross-validation to find optimum k
    folds = min(10, size(YTr, 1)); 
    if folds < size(YTr, 1)
        CVFoldIndexes = genCVFoldIndexes( size(YTr, 1), folds, true);
    else
        CVFoldIndexes = num2cell(1:size(YTr, 1))';
    end
    for j = 1:folds
        innerTsIndex = CVFoldIndexes{j};
        innerTrIndex = find(~ismember(1:size(YTr, 1), innerTsIndex));
        beta = ridge(YTr(innerTrIndex,:), XTr(innerTrIndex,:), settings.lambdaVals, 0);
        YTrHatForKs(innerTsIndex, :) = [ones(size(XTr(innerTsIndex, :), 1), 1), XTr(innerTsIndex, :)]*beta;
    end
    
    e = (repmat(YTr, 1, length(settings.lambdaVals)) - YTrHatForKs);
    RMSE = sqrt(mean(e.^2, 1));
    STDE = std(e, [], 1)./sqrt(size(e, 1));
    [minRMSE, minIndex] = min(RMSE);
    lVals = settings.lambdaVals;
    lVals( RMSE > (minRMSE+STDE(minIndex)) ) = nan;
    [~, kIndex] = min(lVals); 
%     [~, kIndex] = max(lVals); 
%     kIndex = minIndex;
    
%     figure; plot(RMSE); hold on; plot([RMSE'+STDE', RMSE'-STDE'], '--'); scatter([minIndex; kIndex], RMSE([minIndex; kIndex])); 
    
    k = settings.lambdaVals(kIndex);
end