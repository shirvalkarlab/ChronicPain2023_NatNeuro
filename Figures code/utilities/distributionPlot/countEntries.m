function [uniqueEntries,numberOfOccurences,whereIdx] = countEntries(m,isRow, keepNaN)
%COUNTENTRIES returns all unique entries (sorted) in the array m and how many times the respective entries occured
%
%SYNOPSIS [uniqueEntries,numberOfOccurences,whereIdx] = countEntries(m,isRow)
%
%INPUT  m          : any matrix (not cells or structs)
%       isRow(opt) : should rows be counted or not [1/{0}]
%                       (if it's cols, transpose m before calling the function!)
%       keepNaN (opt) : count NaN as entry? [{1}/0] If 0, NaNs (or
%                       NaN-containing rows) are removed after sorting, so
%                       that whereIdx still refers to the original position
%                       of the uniqueEntries in the input array.
%
%OUTPUT uniqueEntries : unique(m)
%                       if only one output argument is requested,
%                       countEntries returns [uniqueEntries,#ofOcc]
%       numberOfOccurences : how many times the unique entries appear in m
%       whereIdx      : where in m do the entries appear? (m = uniqueEntries(whereIdx,:))
%
%
%c: 11/03, jonas
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%---test input
if iscell(m) || isstruct(m)
    error('cells and structs are not supportet as input');
end

if nargin < 2 || isempty(isRow)
    doRow = 0;
else
    if isRow == 1;
        doRow = 1;
    elseif isRow == 0
        doRow = 0;
    else
        error('input argument isRow has to be 1 or 0!')
    end
end
if nargin < 3 || isempty(keepNaN)
    keepNaN = true;
end
%---end test input



if ~doRow %do the fast method
    
    %make m into a vector
    m = m(:);
    
    % new approach: remove NaNs, find uniques, use hist for counting, then
    % use ismember and find(isnan()) for whereIdx if requested. The old
    % approach won't work with the changes to unique :(
    mFinite = m(isfinite(m));
    uniqueEntries = unique(mFinite);
    % use histc instead of hist, since the latter adds zeros if the entries
    % are all equal and >1 (the bins are interpreted as number of bins, not
    % a vector of bins!)
    numberOfOccurences = histc(mFinite,uniqueEntries);
    

    
    % check for inf
    infIdx = m==inf;
    if any(infIdx)
        uniqueEntries(end+1) = inf;
        numberOfOccurences(end+1) = sum(infIdx);
    end
    
    needWhere = nargout > 2;
    nanIdx = isnan(m);
    [~,whereIdx] = ismember(m(~nanIdx),uniqueEntries);
    
    % add back NaNs at the end
    if keepNaN
        
        if any(nanIdx)
            uniqueEntries(end+1) = NaN;
            numberOfOccurences(end+1) = sum(nanIdx);
            if needWhere
                tmp = whereIdx;
                whereIdx = ones(size(m))*length(uniqueEntries);
                whereIdx(~nanIdx) = tmp;
            end
        end
    end
    
   
    
    
else %do it the complicated way
    
    %we do not care about the ordering of the matrix here: if the user
    %specified rows, he/she wanted a columnVector as output (or should read the help)
    [uniqueEntries, dummy, uniqueIdx] = unique(m,'rows');
    
    %rember output
    whereIdx = uniqueIdx;
    
    if ~keepNaN
        % remove NaN, inf
        badIdx = find(any(~isfinite(uniqueEntries),2));
        uniqueEntries(badIdx,:) = [];
        whereIdx(ismember(whereIdx,badIdx)) = [];
        uniqueIdx = whereIdx;
    end
    
    %uniqueIdx returns the indexList where uniqueEntriy #x occurs.
    %We will now sort this list and take a diff to find where this index
    %changes.
    %adding zero and length(uniqueIndex) to the vector, we can now via
    %another diff see how many entries there are (see example)
    
    %example m: [11,11,22,33,33,22,22,22,44,11]
    %corresponding uniqueEntries, uniqueIdx: [11,22,33,44] / [1 1 2 3 3 2 2 2 4 1]
    
    %sort: [1     1     1     2     2     2     2     3     3     4]
    sortedIdx = sort(uniqueIdx);
    
    %diff: [0     0     1     0     0     0     1     0     1]
    sortedIdxDiff = diff(sortedIdx);
    
    %find and add entries: [0     3     7     9    10]
    changeValueIdx = [0;find(sortedIdxDiff);length(uniqueIdx)];
    
    %diff again for the numberOfOccurences: [3     4     2     1]
    numberOfOccurences = diff(changeValueIdx);
end

if nargout < 2
    uniqueEntries = [uniqueEntries,numberOfOccurences];
end