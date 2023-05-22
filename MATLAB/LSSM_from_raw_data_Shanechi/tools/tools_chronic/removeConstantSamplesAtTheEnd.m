% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Omid Sani, Shanechi Lab, University of Southern California, 2016
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%removeConstantSamplesAtTheEnd Removes last few samples of data if all
%channels have constant signal. Useful for trimming end of recordings.
%   Inputs:
%     - (1) x: signal, each column is considered as a separate signal
%     - (2) t (optional): time vector (will be trimmed similar to x)
%     - (3) refChannel (optional): will only consider these channels
%     - (4) tol (optional): tolerance for equality
%   Outputs:
%     - (1) xNew: x with constant samples at the end removed
%     - (2) tNew (optional): t trimmed similar to x
%     - (3) removedIndices (optional): index of removed samples
%   Usage example:
%       x = removeConstantSamplesAtTheEnd( x, t )

function [x, t, toRemoveInd] = removeConstantSamplesAtTheEnd( x, t, refChannels, tol )

if nargin < 2, t = []; end
if nargin < 3 || isempty(refChannels), refChannels = 1:size(x, 2); end
if nargin < 4, tol = 0; end

toRemoveInd = [];

if ~isempty(t)&&(numel(t)~=size(x, 1)), error('Size mismatch of x and t!\n'); end

% Make the special case of no const sample at the end fast
if size(x, 1)>1 && all(abs(diff(x((end-1):end, refChannels), 1, 1)) <= tol, 2)
    
    lastValStartsAt = 1 + find(abs(diff(x(:, refChannels(1), 1), 1, 1)) > tol, 1, 'last');

    if isempty(lastValStartsAt), lastValStartsAt = 1; end

    lastValInAllChansAt = find(any(abs(diff(x(lastValStartsAt:end, refChannels), 1, 1)) > tol, 2), 1, 'last');
    if ~isempty(lastValInAllChansAt)
        lastValStartsAt = lastValStartsAt + lastValInAllChansAt;
    end
    if (lastValStartsAt < size(x, 1))
        % fprintf('Removing %d constant samples from the end.\n', size(x, 1) - lastValStartsAt);
        % Chop off the last few samples that are constant
        toKeepInd = 1:lastValStartsAt;
        toRemoveInd = find(~ismember(1:size(x, 1), toKeepInd));
        x = x(toKeepInd, :);
        if ~isempty(t)
            t = t(toKeepInd);
        end        
    end
end

end
