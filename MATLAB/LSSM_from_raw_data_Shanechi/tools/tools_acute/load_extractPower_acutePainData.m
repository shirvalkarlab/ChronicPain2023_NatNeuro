% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parima Ahmadipour, Omid Sani, Maryam Shanechi
%   Shanechi Lab, University of Southern California, 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function loads data and extracts power features from neural
%  data, and finds noisy samples of it.
%   Inputs:
%     - (1) fileDir: Path to data files
%     - (2) subID: subjects ID
%   
%   Outputs:
%     - (1) pDataStruct: data structure containing extracted power features and corresponding times
%     - (2) noisy samples_time: times corresponding to detected noisy samples of neural data based
%     on a user defined threshold on variance of neural data
%     - (3) LFPmetaexport: structure containing information about pain score and
%     temperature 

function [pDataStruct,noisy_samples_time,LFPmetaexport] = loadPainData_extracPower_acutePain(fileDir,subID)
%% Load data here:
load(sprintf('%sQSTipsi_CP%d_bilateral.mat',fileDir,subID),'LFPexport','LFPmetaexport');
%% loading neural data
rawdata=struct;

if subID==1
rawdata.data=[LFPexport.rightacc,LFPexport.rightofc];
else
rawdata.data=[LFPexport.leftacc,LFPexport.leftofc,LFPexport.rightacc,LFPexport.rightofc];
end     
rawdata.time=LFPexport.time;

%% Finding the noisy epoch of raw neural data. Later, we replace corresponding samples in power with NaN! 
th_var_subject={[100,500],[1000,500,1000,500],[500,500,500,500],[500,500,500,500]}; % I found this threshold for each subject based on visual inspection
th_var=th_var_subject{1,subID};
[noisy_samples_time]=findNoisySamples(rawdata,th_var);% find noisy segments of neural data
%% Do preprocessing and extracting power
% Frequency bands to extract power from 
bands = [1 4; 4 8; 8 12; 12 20; 20 30; 30 55; 65 100];
Fs=1/median(diff(rawdata.time(:,1)));
% Preprocessing settings
settings = struct( ...
    'desiredFs', Fs, ...
    'causal', true, ...
    'doHPFiltering', false, ...
    'doBPFiltering', false, ...
    'doLNFiltering', false, ...
    'DFilterSpecs', [struct('shape', 'BS', 'Fc', [ 66  70], 'type', 'IIR', 'order', 4); ...
                     struct('shape', 'BS', 'Fc', [ 81  85], 'type', 'IIR', 'order', 4); ...
                     struct('shape', 'BS', 'Fc', [103 107], 'type', 'IIR', 'order', 4); ...
                     struct('shape', 'BS', 'Fc', [103 107], 'type', 'IIR', 'order', 4); ...
                     struct('shape', 'BS', 'Fc', [126 130], 'type', 'IIR', 'order', 4); ...
                     struct('shape', 'HP', 'Fc', [1], 'type', 'IIR', 'order', 4); ...
                    ] ...
);

% Preprocessing
[ dataPrepro, time, source, settings, filters ] = doPreprocessing(rawdata.data, [], settings, Fs, [], [] );
% Extract power features
[powers,windowEndSamples] = getBandPowers(bands, dataPrepro, Fs, struct('windowSize', 1)); % Compute powers
powers_perm=permute(powers,[3,2,1]);
powers_concat=reshape(powers_perm,[size(powers_perm,1),size(powers_perm,2)*size(powers_perm,3)]); % concatentates powers from all channels
time_power=((windowEndSamples-1)/Fs+rawdata.time(1))'; % computing the time corresponding to each power sample
pDataStruct=struct('data', powers_concat, 'time', time_power, 'Fs', 1);

end