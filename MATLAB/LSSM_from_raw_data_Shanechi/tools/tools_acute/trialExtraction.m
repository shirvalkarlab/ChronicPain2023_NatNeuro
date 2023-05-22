% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parima Ahmadipour, Maryam Shanechi
%   Shanechi Lab, University of Southern California, 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  This function extracts trials of the experiment. Based on trial end
%  labels and their corresponding times in LFPmetaExport file. Each extracted trial spans the period from
%  one trial end to the next trial end.
%   Inputs:
%     - (1) pDataStruct: structure containing extracted power (or LSSM
%     states)and correspomding times
%     - (2) tDataStruct: structure containing temperature and correponding
%     times
%     - (3) LFPmetaexport: structure containing pain scores and temperature
%     information
%     - (4) subID: subject id
%   
%   Outputs:
%     - (1) X: pDataStruct cut between start of first extracted trial and end of last extracted trial
%     - (2) Y: tDataStruct cut between start of first extracted trial and end of last extracted trial
%     - (3) foldIndices_shifted: a cell array, containing data indices
%     corresponding to each extracted trial
function [X,Y,foldIndices_shifted]=trialExtraction(pDataStruct,tDataStruct,LFPmetaexport,subID);
delay_sub=[3,10,10,10];% this is based on how much earlier than the first trial start we have neural activity for each of the 4 subjects

events=[LFPmetaexport.behaviortable.events];
ts=find(strcmp(events, 'trialstart'));te=find(strcmp(events, 'trialend'));
times_events= [LFPmetaexport.behaviortable.time];
ts_times=times_events(ts);te_times=times_events(te);



isNonNaN = ~isnan(pDataStruct.data(:,1));
pDataStruct.data=pDataStruct.data(isNonNaN ,:);pDataStruct.time=pDataStruct.time(isNonNaN ,1);
if ~isempty(tDataStruct)
tDataStruct.data=tDataStruct.data(isNonNaN ,1);tDataStruct.time=tDataStruct.time(isNonNaN ,1);
end

foldIndices=cell(length(ts_times),1);
foldIndices_shifted=cell(length(ts_times),1);

for i=1:length(ts_times)
    if i==1
  
        t_start_fold=ts_times(i)-delay_sub(subID);         
    else 
        t_start_fold=te_times(i-1);
    end
        t_end_fold=te_times(i);
        [~,fold_index_start]=min(abs(pDataStruct.time-t_start_fold));
        if i==1
            shift_index=fold_index_start-1;
        end
        [~,fold_index_endp1]=min(abs(pDataStruct.time-t_end_fold));
        foldIndices{i,1}=fold_index_start:fold_index_endp1-1;
        foldIndices_shifted{i,1}=foldIndices{i,1}-shift_index;   
end
start_index=shift_index+1;
end_index=fold_index_endp1-1;
X.data=pDataStruct.data(start_index:end_index,:);
X.time=pDataStruct.time(start_index:end_index);
if ~isempty(tDataStruct)
    Y.data=tDataStruct.data(start_index:end_index);
    Y.time=tDataStruct.time(start_index:end_index);
else
    Y=[];
end

end