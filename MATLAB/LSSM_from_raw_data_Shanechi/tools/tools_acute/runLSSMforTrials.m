% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parima Ahmadipour, Maryam Shanechi
% Shanechi Lab, University of Southern California, 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% runLSSMforTrials extracts LSSM states for each trial by running the Kalman filter starting from the previous
% trial
%   Inputs:
%     - (1) X: predictors
%     - (2) Y: labels
%     - (3) TrialInds: cell array containing indices for each extracted trial   
%     - (4) nx: model order (dimension of the latent state for fitting
%     LSSM)
%   Outputs:
%     - (1) xPred_no1stFold : estimated states with the Kalman filter for all trials, excluding first
%     trial (fold) 
%     - (2) Y_no1stFold: labels, excluding labels of the first trial(fold)
%     - (3) TrialInds_no1stFold: indices for all trials, excluding the
%     1st trial(fold)

function [xPred_no1stFold,Y_no1stFold,TrialInds_no1stFold]=runLSSMforTrials(X,Y,TrialInds,nx)
X.data=(X.data-nanmean(X.data,1));%subtracting the mean
horizon = 10;
[A,~,C,~,K] = subid(X.data, [], horizon, nx, [], [], 1);  % Subspace sys id
xPred = nan(size(X.data, 1),nx);
start_2ndFold=TrialInds{2}(1,1);
TrialInds_no1stFold=cell(length(TrialInds)-1,1);
for fi = 2:numel(TrialInds)
    testFeatInds = [TrialInds{fi-1},TrialInds{fi}];
    xPred_temp = kalmanFilter(X.data(testFeatInds,:), A, C, K); % Kalman filter to estimate the states 
    xPred(TrialInds{fi}, :)=xPred_temp(length(TrialInds{fi-1})+1:end,:); % just keeping the estimated states of the "fi"th trial (excluding the states from the trial before it.)
    TrialInds_no1stFold{fi-1}=TrialInds{fi}-start_2ndFold+1;
    
end


xPred_no1stFold.data=xPred(start_2ndFold:end,:);
xPred_no1stFold.time=X.time(start_2ndFold:end);
if ~isempty(Y)
Y_no1stFold.data=Y.data(start_2ndFold:end);
Y_no1stFold.time=Y.time(start_2ndFold:end);
else
Y_no1stFold=[];
end
end