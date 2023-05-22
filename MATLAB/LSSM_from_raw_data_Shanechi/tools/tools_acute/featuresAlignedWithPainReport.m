% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parima Ahmadipour, Maryam Shanechi
%   Shanechi Lab, University of Southern California, 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes features at time correpondng to time of pain score in
% PainDataStruct 
%    Inputs:
%     - (1) painDataStruct: structure containing pain scores and
%     corresponding times
%     - (2) feat: feature matrix with dimension:(time length of features)*(number of features)
%     - (3) avgSamplesAroundPeak: samples around the painDataStruct.time (temperature peak start
%     time) that we want to average over and consider as a feature for
%     decoding
%     - (4) AvgPowerBands: if true, the features are also averaged over
%     channels and frequency bands.
%   
%   Outputs:
%     - (1) featPainAlg: extracted features for decoding which are aligned with time of pain
%     scores in painDataStruct


function featPainAlg=featuresAlignedWithPainReport(painDataStruct,feat,avgSamplesAroundPeak,AvgPowerBands)
 
for i=1:length(painDataStruct.data)
   [~,index(i)]= min(abs(painDataStruct.time(i)-feat.time));
   avg_time_indices=index(i)+avgSamplesAroundPeak;
   if AvgPowerBands
      featPainAlg.data(i,:)=nanmean(nanmean(feat.data(avg_time_indices,:),1)); % average over both time and channels/frequencies around the times in painDataStruct
   else
      featPainAlg.data(i,:)=nanmean(feat.data(avg_time_indices,:),1); % average over a few time samples around the temperature peak start times
   end
   featPainAlg.time(i,1)=feat.time(index(i));

end
       
            
end