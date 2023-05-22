% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Omid Sani, Maryam Shanechi
%   Shanechi Lab, University of Southern California, 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script runs the cross-validated model fitting with some demo 
% MATLAB datasets to confirm its validity

addpath(genpath('./'));

%%
% Demo dataset 1
load fisheriris
inds = ~strcmp(species,'setosa');
X = meas(inds,3:4);
Y = species(inds);
% [labelStrs, ~, Y] = unique(YStr);

settings = struct;
settings.predictor = 'logistic';

performCVedModelFit(X,Y, settings);
performPermTest(X,Y, settings);

%%
% Demo dataset 2
load ionosphere
NumFolds=4;
NumGrid=30;

settings = struct;
settings.predictor = 'logistic';

performCVedModelFit(X,Y, settings);
performPermTest(X,Y, settings);

%%
% Demo dataset 3
N = 100;  % Samples
Ny = 1; % Outputs
X = randn(N, 5);
b = randn(5, Ny);
n = randn(N, Ny)/10;
Y = X * b + n;

settings = struct;
settings.predictor = 'ridge';

performCVedModelFit(X,Y, settings);
performPermTest(X,Y, settings);

