% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parima Ahmadipour, Maryam Shanechi
%   Shanechi Lab, University of Southern California, 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script runs the cross-validated model fitting with the acute pain
% data for decoding temperature
% Change settings.predictor, useLSSM, and subID variables to try different decoding methods for different subjects.

addpath(genpath('./'));

clear all;
subID=4; % subject ID % change this to run a partiular participant (CP1-4)
fileDir= 'G:/2023_Shirvalkar_NatNeuro_ChronicPainBiomarkers/MATLAB/'; %Path to AcuteQST data files (change this based on your local machine)



%% Load Data and extract power features from neural features
[pDataStruct,noisy_samples_time,LFPmetaexport] = load_extractPower_acutePainData(fileDir,subID);  
%% Align temperature and nerual power features
tTime = LFPmetaexport.temptime;% tTime:temprature time
if size(tTime,1)==1
    tTime=tTime';% column vector
end
tFs = 1/median(diff(tTime));% frequency of tempreture
[b, a] = designFilterLPIIRButter(tFs, 0.5, 4);% Doing a lowpass filter on temperature before low-pass filtering to avoid aliasing
aaFilt = dFilter(b, a, false, tFs);
tData = aaFilt.apply(LFPmetaexport.temperature);
tDataStruct = struct('data', tData, 'time', tTime, 'Fs', tFs);% tDataStruct: structure containing temperature info
[pDataStructS, tDataStructS] = syncDataStructs(pDataStruct, tDataStruct);%"S" at the end of the struct name means synchronyzed with temperature/power
%% Make noisy power samples NaN
pDataStructSC=MakeNanBadPowers(pDataStructS,noisy_samples_time);%C at the end of the struct name mean clean
%% Finding the trial indices and cut the correct part of data based on trials
[pDataTrialBased,tDataTrialBased,TrialInds]=trialExtraction(pDataStructSC,tDataStructS,LFPmetaexport,subID);
%% Classification with "logistic" or exact prediction using "ridge", using "power" or "LSSM states" as features.
settings = struct
settings.predictor='logistic'; % could be 'ridge' or 'logistic'
useLSSM=false;% if true LSSM states are used for decoding, otherwise power features are used.    
if ~useLSSM
    X_noZscore=pDataTrialBased.data;
    Y=tDataTrialBased.data;
    settings.CVFoldInds=TrialInds;
    settings.trainMinDistFromTest=1;
else
    nx=1; % The optimal model order (state dimension) is found in the original code
    [xPred_no1stFold,Y_no1stFold,CVFoldInds_no1stFold]=runLSSMforTrials(pDataTrialBased,tDataTrialBased,TrialInds,nx);% For each trial, LSSM states are extracted by running a Kalman filter starting from the previous trial.
    X_noZscore=xPred_no1stFold.data;
    Y=Y_no1stFold.data;
    settings.CVFoldInds=CVFoldInds_no1stFold;
    for i=1:length(TrialInds)-1
        trial_lengths(i)=length(TrialInds{i}); 
    end
    settings.trainMinDistFromTest=max(trial_lengths);% Leaving a trial gap between train and test (because states are extracted using previous trials data)
end
settings.numFolds=length(settings.CVFoldInds);
settings.numPerms=100;
X = (X_noZscore - nanmean(X_noZscore,1))./repmat(nanstd(X_noZscore,0,1),size(X_noZscore,1),1);
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
performCVedModelFit(X,Y, settings);   % Cross-validated classification
performPermTest(X,Y, settings);       % Find chance level with permutation