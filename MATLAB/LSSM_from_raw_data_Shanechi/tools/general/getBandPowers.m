% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Omid Sani, Maryam Shanechi
%   Shanechi Lab, University of Southern California, 2019
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%getBandPowers Computes signal power in frequency bands using pwelch
%   Inputs:
%     - (1) bands: frequency bands
%     - (2) data: input time domain data
%     - (3) Fs sampling rate of input data in Hz
%     - (4) settings: struct with additional settings. For details, see
%   Outputs:
%     - (1) power (channel, band, window)
%     - (2) windowEndSamples: sample indices corresponding to the end of each window used
%     for computation of power
%   Usage example:
%      power  = getBandPowers(bands, data, Fs)


function [ power,windowEndSamples] = getBandPowers(bands, data, Fs, settings)

if nargin < 4, settings = struct; end
if ~isfield(settings, 'windowSize'), settings.windowSize = []; end % Window size in [s]
if isempty(settings.windowSize)
    settings.windowSize = size(data, 1)/Fs;
end

if ~isfield(settings, 'method'), settings.method = 'multitaper'; end  % Can be 'multitaper' or 'pwelch'
% pwelch settings
if ~isfield(settings, 'pwWindow'), settings.pwWindow = []; end
if ~isfield(settings, 'pwOverlap'), settings.pwOverlap = []; end
% Multitaper settings
if ~isfield(settings, 'mtNW'), settings.mtNW = max(settings.windowSize, 1.25); end % The time-halfbandwidth product
% Shared settings
if ~isfield(settings, 'NFFT'), settings.NFFT = 2^nextpow2(Fs); end


winSamples = round(settings.windowSize * Fs);

windowEndSamples = winSamples:winSamples:size(data, 1);


for wi = 1:numel(windowEndSamples)
    winInds = (windowEndSamples(wi)-winSamples+1):windowEndSamples(wi);
    if strcmpi(settings.method, 'pwelch')
        [pxx,w]=pwelch(data(winInds, :), settings.pwWindow, settings.pwOverlap, settings.NFFT, 'twosided');
    elseif strcmpi(settings.method, 'multitaper')
        [pxx,w]=pmtm(data(winInds, :), settings.mtNW, settings.NFFT, 'twosided');
    else
        error('Not recognized!');
    end
    pxx = pxx * (2*pi);
    f = w / (2*pi) * Fs;
    for i=1:size(bands,1)
        bandInds=find(f>=bands(i,1)&f<bands(i,2));
        power(:, i, wi)=pow2db(mean(pxx(bandInds,:),1) + eps);
    end
end

end
