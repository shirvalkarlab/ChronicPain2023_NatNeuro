% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Omid Sani, Maryam Shanechi
%  Shanechi Lab, University of Southern California, 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%plotAUC Plots AUC
%   Inputs:
%     - (1) labels: true labels of data
%     - (2) scores: scores from classifier
%     - (3) posclass: label of positive class (score should be associated with this class)
%     - (4) settings: struct with additional settings. For details, see
%   Outputs:
%     - (1) ax
%   Usage example:
%       plotAUC(labels, scores, posclass);

function ax = plotAUC(labels, scores, posclass, settings)

if nargin < 4, settings = struct; end
if ~isfield(settings, 'ax'), settings.ax = []; end
if ~isfield(settings, 'titleHead'), settings.titleHead = ''; end

[X,Y,T,AUC] = perfcurve(labels, scores, posclass);

if isempty(settings.ax)
    figure; ax = axes; 
else
    ax = settings.ax;
end
hold(ax, 'on'); axis(ax, 'square');
plot(ax, X, Y, 'LineWidth', 2, 'color', 'k'); 
line(ax, [0 1], [0 1], 'LineStyle', '--', 'Color', 'k');
xlabel(ax, 'False positive rate (FPR)');
ylabel(ax, 'True positive rate (TPR)');

titleStr = sprintf('AUC: %.3f', AUC);
title(ax, sprintf('%s%s', settings.titleHead, titleStr));

thr = 0.5;
% [~, ind] = min(abs(T-thr));
% thr = T(ind);

predLabel = (scores >= thr);
trueLabels = (labels == posclass);
ACC = sum( predLabel == trueLabels ) / numel(labels);
FPR = sum( predLabel == 1 & trueLabels == 0 ) / sum(trueLabels == 0);
TPR = sum( predLabel == 1 & trueLabels == 1 ) / sum(trueLabels == 1);

h = scatter(ax, FPR, TPR, 'filled', 'ro');
legStr = sprintf('TPR: %.3g\nFPR: %.3g\nACC: %3g', TPR, FPR, ACC);
legend(ax, h, legStr, 'Location', 'SE', 'AutoUpdate', false);


end
