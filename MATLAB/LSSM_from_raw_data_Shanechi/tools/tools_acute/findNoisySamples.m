% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Parima Ahmadipour, Maryam Shanechi
%   Shanechi Lab, University of Southern California, 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function finds the noisy samples of data based on comparison with a user
% defined threshold times the variance of data
%   Inputs:
%     - (1) rawdata: raw neural data
%     - (2) th_var: user defined threshold (if a data sample is larger than var_th*var(data) it is considered as noisy) 
%   Outputs:
%     - (1) noisy_samples_time: an array containing time of noisy samples
%     of data
function [noisy_samples_time]=findNoisySamples(rawdata,th_var)
noisy_samples=zeros(size(rawdata.data));
if length(th_var)==1
    th_var=ones(1,size(rawdata.data,2))*th_var;
end
for i=1:size(rawdata.data,2)
    data_ch=rawdata.data(:,i);
    var_ch=var(rawdata.data(:,i));
    noisy_samples(:,i)=(abs(data_ch)>th_var(i)*var_ch);
end
nosiy_samples_all=sum(noisy_samples,2)>0;
noisy_samples_time=rawdata.time(nosiy_samples_all);
end
    