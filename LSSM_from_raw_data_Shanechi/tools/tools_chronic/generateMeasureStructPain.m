% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Omid Sani, Maryam Shanechi
%   Shanechi Lab, University of Southern California, 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%generateMeasureStructPain Generates the score structure array for the subject
%   Inputs:
%     - (1) content of Pain data file
%     - (2) factor (default: 'painVAS'): pain measure to return
%   Outputs:
%     - (1) measureVals: the structure containing scores and their info
%   Usage example:
%       measureVals = generateMeasureStructPain( PData, 'painVAS' ); % Returns painVAS

function [ measureVals ] = generateMeasureStructPain( PData, factor )

if nargin < 2, factor = 'relativeVAS'; end

if isstruct(PData)
    fD = PData;
else
    rawDataFile = PData;
    fD = load(rawDataFile);
end

fNames = fieldnames(fD);
dataFName = fNames(contains(fNames, 'LFPmeta'));
allD = fD.(dataFName{1});

if ~isfield(allD, 'time') && isfield(allD, 'left')
    allD = allD.left;
elseif ~isfield(allD, 'time') && isfield(allD, 'right')
    allD = allD.right;
end

N = numel(allD.time);

if contains(factor, 'smooth')
    ind = strfind(factor, '_smooth');
    sm = factor(ind:end);
    factor = factor(1:(ind-1));
    A = sscanf(sm, '_smooth%d');
    smoothHist = A;
else
    smoothHist = 1;
end

if isfield(allD, factor)
    mValue = allD.(factor);
    IMSTimeNum = datenum(allD.time(:));
elseif ismember(factor, fieldnames(allD.autopain))
    mValue = allD.autopain.(factor);
    IMSTimeNum = datenum(allD.autopaintime(:));
end

nonNaN = find(~isnan(mValue));
if smoothHist > 1 && any(nonNaN)
    x = [mValue(nonNaN(1))*ones(smoothHist-1, 1); mValue(nonNaN)];
    H = ones(smoothHist, 1)/smoothHist;
    y = conv(x, H);
    yS = y(smoothHist:(end-smoothHist+1));
    % figure; plot([mValue(nonNaN), yS]);
    mValue(nonNaN) = yS;
    fprintf('Using smoothed "%s"\n', factor);
end

mTimeZoneOffset = zeros(size(IMSTimeNum));
% [ offset ] = calcTimeZoneOffsetFromDatenumInSomeOtherOffset( GMTDatenum, 'America/Los_Angeles', 'UTC' );
mTimeVec = datevec(IMSTimeNum); 
mDuration = nan(size(mValue));
mType = repmat({''}, size(mValue));
mDuringStim = false(size(mValue));
mStimInfo = repmat({''}, size(mValue));

if exist('mValue', 'var')&&~isempty(mValue)
    if ~exist('mDuration','var'), mDuration = nan(size(mValue)); end
    clear IMS
    for i = 1:length(mValue)
        measureVals(i, 1) = struct('value', mValue(i, :), 'time', IMSTimeNum(i), 'timevec', mTimeVec(i, :), 'tzone', mTimeZoneOffset(i), ...
                           'duration', mDuration(i), 'type', mType{i}, 'durstim', mDuringStim(i), 'stiminfo', mStimInfo{i});
    end
else
    measureVals = struct('value', {}, 'time', {}, 'timevec', {}, 'tzone', {}, ...
                 'duration', {}, 'type', {}, 'durstim', {}, 'stiminfo', {});
end

end
