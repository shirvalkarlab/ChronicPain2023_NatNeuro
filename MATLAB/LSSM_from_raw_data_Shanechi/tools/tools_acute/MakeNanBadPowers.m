% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parima Ahmadipour, Maryam Shanechi
%   Shanechi Lab, University of Southern California, 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function inserts a NaN where the data is noisy (specified by
% bad_samples_time)
%     Inputs:
%     - (1) Feats: features
%     - (2) folds: bad_samples_time: time of noisy samples
% 
%   Outputs:
%     - (1) FeatsNan: featur matrices with NaN, where they are noisy. 
function FeatsNan=MakeNanBadPowers(Feats,bad_samples_time)
FeatsNan=Feats;

for i=1:length(bad_samples_time)
    ind_nan=min(find((Feats.time-bad_samples_time(i))>=0));
    FeatsNan.data(ind_nan,:)=nan;
end

end