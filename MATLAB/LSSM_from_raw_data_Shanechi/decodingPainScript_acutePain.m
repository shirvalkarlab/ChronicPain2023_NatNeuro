% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parima Ahmadipour, Maryam Shanechi
%   Shanechi Lab, University of Southern California, 2020
%  edited by PShirvalkar 2023 UCSF  -  for file directory friendliness
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script runs the cross-validated model fitting with the acute pain
% data for decoding pain
% Change settings.predictor, useLSSM, avgSamplesAroundPeak, and subID variables to try different decoding methods/features for different subjects.

addpath(genpath('./'));

clear
subID=4; % subject ID  % CHANGE THIS to run a particular participant
fileDir='G:/2023_Shirvalkar_NatNeuro_ChronicPainBiomarkers/MATLAB/'; %Path to data files (change this based on your local machine)


%% Load Data and extract power features from neural features
[pDataStruct,noisy_samples_time,LFPmetaexport] = load_extractPower_acutePainData(fileDir,subID);  
%% create Pain data struct
events=[LFPmetaexport.behaviortable.events]; 
ps=find(strcmp(events, 'peakstart'));
times_allevents= [LFPmetaexport.behaviortable.time];pain_allevents=[LFPmetaexport.behaviortable.pain];
ps_times=times_allevents(ps);% times of temperature peak start for all trials
ps_pain=pain_allevents(ps);% pain scores for all trials
painDataStruct=struct('data',ps_pain,'time',ps_times);
%% Make noisy power samples NaN
pDataStructC=MakeNanBadPowers(pDataStruct,noisy_samples_time);%C at the end of the struct name mean clean
%% Finding the trial indices and cut the correct part of data based on trials
[pDataTrialBased,~,TrialInds]=trialExtraction(pDataStructC,[],LFPmetaexport,subID);
%% Classification with "logistic" or exact prediction using "ridge", using "power" or "LSSM states" as features.
settings = struct;
settings.predictor='logistic'; % could be 'ridge' or 'logistic'
useLSSM=false;%  Set to true to fit LSSM model and use its states for decoding, otherwise power is used for decoding.
avgSamplesAroundPeak=[0,1,2];% feature indices around the start time of temperature peak to take their average and use as a new feature for decoding
if ~useLSSM
    pFeature_painAlg=featuresAlignedWithPainReport(painDataStruct,pDataTrialBased,avgSamplesAroundPeak,true); % power features aligned with pain report
    X_noZscore=pFeature_painAlg.data;
    Y=painDataStruct.data;
    settings.trainMinDistFromTest=1;
else
    nx=1; % The optimal model order (state dimension) is found in the original code
    [xPred_no1stFold]=runLSSMforTrials(pDataTrialBased,[],TrialInds,nx);% For each trial, LSSM states are extracted by running a Kalman filter starting from the previous trial.
    painDataStruct.data=painDataStruct.data(2:end,:);
    painDataStruct.time=painDataStruct.time(2:end,1);% exluding the 1st trial
    Y=painDataStruct.data;%1st trial excluded
    xFeature_painAlg=featuresAlignedWithPainReport(painDataStruct,xPred_no1stFold,avgSamplesAroundPeak,false);% state features aligned with pain report
    X_noZscore=xFeature_painAlg.data;
    settings.trainMinDistFromTest=2;% Leaving a trial gap between train and test (because states are extracted using previous trials data)
end
settings.numFolds=[];% leave-one-out cross-validation
settings.numPerms=100;% number of random permutations for statistical test
X = (X_noZscore - nanmean(X_noZscore,1))./repmat(nanstd(X_noZscore,0,1),size(X_noZscore,1),1); % z scoring data
% X = X_noZscore;
if strcmp(settings.predictor,'logistic') % Convering labels to 0 and 1 for classification
Y_med=nanmedian(Y);
    for y_index=1:length(Y)
        if ~isnan(Y(y_index))
            if Y(y_index)>Y_med
                Y(y_index)=true;
            else
            	Y(y_index)=false;
            end
        end
    end
end
isNonNaN=~isnan(Y);
X=X(isNonNaN,:);Y=Y(isNonNaN,1);
performCVedModelFit(X,Y, settings);   % Cross-validated classification
performPermTest(X,Y, settings);       % Find chance level with permutation