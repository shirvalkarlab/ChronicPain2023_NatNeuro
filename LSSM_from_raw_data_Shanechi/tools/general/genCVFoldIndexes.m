% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Omid Sani, Maryam Shanechi
%   Shanechi Lab, University of Southern California, 2016
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%genCVFoldIndexes Generates some indexs for n-fold cross-validation
%   Inputs:
%     - (1) numOfSamples: number of data points
%     - (2) folds: number of folds
%     - (3) rndOrder: If true, will randomize the order of data points 
%                     (default: false)
%   Outputs:
%     - (1) foldIndexes: cell array of index arrays. 
%   Usage example:
%       foldIndexes = genCVFoldIndexes( 100, 5, true );

function foldIndexes = genCVFoldIndexes( numOfSamples, folds, rndOrder )
if nargin < 3, rndOrder = false; end

if (folds > numOfSamples), error('More folds than samples!\n'); end

allFoldSamples = zeros(folds, 1);
lastUsedIndex = 0;
for i = 1:folds
    foldSamples = round( (numOfSamples - lastUsedIndex) / (folds - i + 1) );
    if (foldSamples < 1), keyboard; foldSamples = 1; end
    allFoldSamples(i) = foldSamples;
    lastUsedIndex = lastUsedIndex + foldSamples;
end

allFoldSamples = sort(allFoldSamples, 'descend');
foldIndexes = cell(folds, 1);

lastUsedIndex = 0;
for i = 1:folds
    foldIndexes{i} = (lastUsedIndex)+(1:allFoldSamples(i));
    lastUsedIndex = foldIndexes{i}(end);
end

if (rndOrder)
    sampleOrder = randperm(numOfSamples);
    for i = 1:folds
        foldIndexes{i} = sampleOrder(foldIndexes{i});
    end
end


end
