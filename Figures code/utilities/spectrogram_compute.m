function LFPspectra = spectrogram_compute(LFP,Fs,specgram_method,time_limits,plot_var,varargin)
% function LFPspectra = spectrogram_compute(LFP,420.11,'mt',[time lims],plot var (1 to plot),[Behavior{y}],[1=trialavg],'z')
% 
% USE THIS FOR QST OR ACUTE TASKS !!!
% 
% 
% % This will compute spectrograms using 1 of 4 possible methods: 
%   'mt' multitaper spectrogram (Chronux) 
%   'wavelet' Morelet Wavelet
%   'welch' periodogram (as good as mt usually)
%   'hilbert' hilbert transform complex amp and phase
% 
%  INPUTS
%   - LFP is the LFP structure containing LFPdata usually LFP.acc or
%   LFP.ofc
%            if multiple trials, they should be organized in columns (samples,trials)
%   - Fs= sampling rate Hz
%   - method to use for spectrogram (string variable, either 'mt','wavelet'
%   etc as above)
%   - time_limits as [start end] 2 element vector to use as start and end
%   times in seconds- IF EMPTY, it will plot all data ***ONLY USEFUL FOR
%   PLOTTING ONE TRIAL, IF COMPUTING MANY TRIALS, TURN THIS OFF**
%   - plot_var -> 0=no, 1 = yes
% 
%    VARARGIN
%  1 - can be either a 'params' structure OR B
% %                     Behavioral Annotation data variable (i.e. Behavior with .events) if plot_var =1
%  2 -second variable input for trialavg 0=no, 1 =yes
%  3 - if 'z' then will zscore spectrogram
% 
%       IF a params structure is specified, it should have these fields:
%             params.winstep= [1 0.05]; %window len (s), window step (s)
%             params.tapers = [3 5]; % [W T(Duration) p] or   [Time x BW,  #slepian tapers] (for timestep [2 .1] use [12 22]
%             params.Fs=420.1;
%             params.pad=1;
%             params.fpass=([0 100]);
%  
% 
% 
% 
% 
%   OUTPUTS
%        - LFPspectra is a structure containing the z-scored spectra across
%        frequency slices (to use relative zscoring, change zscore
%        orientation with zscore(x,0,2)
% 
% GENERAL ------
% Bandpass is set at 1-100 Hz
% For mtspectrum, window is 0.5 sec, with 0.1 sec step
% For welch, time resolution is set to 5 sec
% 
% EXAMPLE
% LFPspectra = spectrogram_compute(LFP,LFPmeta,'welch',[],1,Behavior)
% 
% 
% prasad shirvalkar mdphd sept 2018 
% updated jan 2019 
% updated apr 2020 - added log colorlim to non 'z' scored data and
% alternate 'params' input to adjust window and step size
%   updated oct 2021 - added bandpower timeseries to output
% *** IF YOU MAKE CHANGES to this file PLEASE NOTE THEM HERE**

if nargin>5 && isfield(varargin{1},'tapers') && strcmp(specgram_method,'mt')
    paraminput = varargin{1};
winstep = paraminput.winstep; %window len (s), window step (s)
params.tapers = paraminput.tapers; % [W T(Duration) p] or   [Time x BW,  #slepian tapers] (for timestep [2 .1] use [12 22]
params.Fs = 420.1;
params.pad=paraminput.pad;
params.fpass=paraminput.fpass;
else 
winstep = [1 0.1];
params.tapers = [3 5];
params.pad=1;
params.fpass=([0 100]);
fs=420.1; %actual sampling ratewinstep= [0.1 0.05]; %window len (s), window step (s)
end




 if size(varargin,2)>1 && varargin{2}==1 %? trial averaging
     params.trialave=1;
     disp('Averaging all trials...')
 else
     params.trialave=0;
 end
 
 
if ~isempty(time_limits) %set time limits for processing
Si=time_limits(1);
Fi=time_limits(2);
start_ind=round(Fs*Si+1);
fin_ind = round(Fs*Fi);
else
    
    Si=0;
    Fi=round(length(LFP.acc)/Fs);
    start_ind=1;
    fin_ind =size(LFP.acc,1);
end

t_lfp=linspace(Si,Fi,fin_ind-start_ind+1); %LFP time vector
 

 if winstep(1)>= t_lfp(end)
     error(['Window length of ' num2str(winstep(1)) 'sec is too long for the short bit of LFP being computed. Please adjust...'])
 end
 

 switch specgram_method
     case 'mt'
% % ================== mtspecgram =================

disp('running mt spectrograms')
fprintf('Time window/ step is  %0.1f and %0.2f.  TW and Tapers are %d and %d \n', winstep(1), winstep(2), params.tapers(1), params.tapers(2));

[SaL,tl,fq1]= mtspecgramc(LFP.acc(start_ind:fin_ind,:),winstep,params);
[SoL,~,~]= mtspecgramc(LFP.ofc(start_ind:fin_ind,:),winstep,params);
    
   if size(varargin,2)>2 &&  strcmp(varargin{3},'z')
         SaccL = zscore(log10(SaL)); SofcL = zscore(log10(SoL)); %zscore output
   else
       SaccL = (log10(SaL)); SofcL = (log10(SoL));
%        lim_color=10; %color limits for plotting
   end
   



            plotSaccL = SaccL;
            plotSofcL  = SofcL;

%     Relative power
% SaccL = zscore(log10(SaL),0,2); SofcL = zscore(log10(SoL),0,2);

% % or 
% SaccL = log10(SaL)./sum(log10(SaL),2);
% SofcL = log10(SoL)./sum(log10(SoL),2);
%     

tl=tl+Si;

% =========================================

     case 'wavelet'
% ==========  Wavelet ============
disp('running wavelet')
tic

Fs= params.Fs;
F=1:1:100;  
fq1=F;

parfor x = 1:size(LFP.acc,2)
[SaL{x}] = cwt(LFP.acc(start_ind:fin_ind,x),centfrq('cmor1 -1 ')*Fs./F,'cmor1 -1 ');
% ,'NumOctaves',8,'VoicesPerOctave',48);  % will give min Hz about 1 (if want 0.5 need to use 9 octaves.  voices per octave gives the resolution.  this is the max allowed
[SoL{x}] = cwt(LFP.ofc(start_ind:fin_ind,x),centfrq('cmor1 -1 ')*Fs./F,'cmor1 -1 ');
% ,'NumOctaves',8,'VoicesPerOctave',48);

end

for y=1:size(LFP.acc,2)
% Collect outputs for each trial
SaccL(:,:,y) = SaL{y}'; % to make the dims = time, fq, trials
SofcL(:,:,y) = SoL{y}';
end

tl=linspace(Si,Fi,size(SaccL,2));


            plotSaccL = SaccL;
            plotSofcL  = SofcL;
toc
% YTicks = 2.^(round(log2(min(freqs))):round(log2(max(freqs))));
% fq1=log2(freqs);
% SaccL=SaccL';SofcL=SofcL';
% =======================================

     case 'welch'
% ================== Pwelch power spectrum =================
tic
Timeres=1; %sec
for x=1:size(LFP.acc,2)
[Sa,fq1,tl]= pspectrum(LFP.acc(start_ind:fin_ind,x),params.Fs,'spectrogram','FrequencyLimits',[0 100],'TimeResolution',Timeres);
[So,~,~]= pspectrum(LFP.ofc(start_ind:fin_ind,x),params.Fs,'spectrogram','FrequencyLimits',[0 100],'TimeResolution',Timeres);

SaL(:,:,x)=Sa';
SoL(:,:,x)=So';
end

    
SaccL = (log10(SaL)); SofcL = (log10(SoL)); %remove ' to see zscore within time slice, to see TENS artifact well, and high beta/ low gamma
% 
% % Relative scaling (% freq bandpower)
% SaccL = log10(SaL)./sum(log10(SaL),2);
% SofcL = log10(SoL)./sum(log10(SoL),2);
% SaccL=SaccL';SofcL=SofcL';

            plotSaccL = SaccL;
            plotSofcL  = SofcL;

tl=tl+Si;
toc
% =========================================

    case 'hilbert'
% ================== Hilbert Analytic Amplitude =================
fprintf('\n')
tic
% Calculate just the bands of interest for burst detection
bands= {'delta','theta','alpha','beta','Lgamma','Hgamma'};
q.delta =  [1 4];
q.theta = [4 8];
q.alpha = [8 12];
q.beta = [12 30];
q.Lgamma =[30 70];
q.Hgamma = [70 100];


a=LFP.acc(start_ind:fin_ind,:);
b=LFP.ofc(start_ind:fin_ind,:);
            for bb=1: numel(bands)

            lo= q.(bands{bb})(1);
            hi= q.(bands{bb})(2);

            a2=eegfilt_ps(a',Fs,lo,hi);
            b2=eegfilt_ps(b',Fs,lo,hi);
            ha=hilbert(a2');hb=hilbert(b2');
            lfpfiltacc(:,bb,:) = a2;
            lfpfiltofc(:,bb,:) = b2;

            SaccL(:,bb,:)=ha;
            SofcL(:,bb,:)=hb;

            plotSaccL =  abs(SaccL);
            plotSofcL  = abs(SofcL);
            end
tl=linspace(Si,Fi,size(SaccL,1));
freqs=[4,8,12,30,70,100];
% YTicks = 2.^(round(log2(min(freqs))):round(log2(max(freqs))));
fq1=(freqs);
toc
% =========================================
 end %switch case
 
 %%  DECLARE OUTPUT VARIABLES

LFPspectra.acc=SaccL;
LFPspectra.ofc=SofcL;
LFPspectra.method=specgram_method;

        %  Adjust the spectrogram times if LFP has been adjusted
        if isfield(LFP,'time')
            lfp1=LFP.time(1);
            lfpend=LFP.time(end);
            tl=linspace(lfp1,lfpend,length(tl));
            t_lfp=LFP.time;
        end
        
LFPspectra.times=tl;
LFPspectra.lfptimes=t_lfp;
LFPspectra.fq=fq1;

        if isfield(LFP,'trialtimes')
            LFPspectra.trialtimes=linspace(LFP.trialtimes(1),LFP.trialtimes(end),numel(tl));
        else
            LFPspectra.trialtimes=linspace(t_lfp(1),t_lfp(end),numel(tl));
        end



%% Add the bandpwr time series 
bands= {'delta','theta','alpha','beta','Lgamma','Hgamma'};
LFPspectra.bands=bands;
fq=fq1;

q.delta = (fq>=1 & fq<=4);
q.theta = (fq>4 & fq<=8);
q.alpha =(fq>8 & fq<=12);
q.beta = (fq>12 & fq<=30);
q.Lgamma =(fq>30 & fq<=70);
q.Hgamma = (fq>70 & fq<=100);  
   

switch specgram_method
    case 'mt'
LFPspectra.accdelta = squeeze(mean(LFPspectra.acc(:,q.delta,:),2));
LFPspectra.acctheta = squeeze(mean(LFPspectra.acc(:,q.theta,:),2));
LFPspectra.accalpha = squeeze(mean(LFPspectra.acc(:,q.alpha,:),2));
LFPspectra.accbeta = squeeze(mean(LFPspectra.acc(:,q.beta,:),2));
LFPspectra.accLgamma = squeeze(mean(LFPspectra.acc(:,q.Lgamma,:),2));
LFPspectra.accHgamma = squeeze(mean(LFPspectra.acc(:,q.Hgamma,:),2));

LFPspectra.ofcdelta = squeeze(mean(LFPspectra.ofc(:,q.delta,:),2));
LFPspectra.ofctheta = squeeze(mean(LFPspectra.ofc(:,q.theta,:),2));
LFPspectra.ofcalpha = squeeze(mean(LFPspectra.ofc(:,q.alpha,:),2));
LFPspectra.ofcbeta = squeeze(mean(LFPspectra.ofc(:,q.beta,:),2));
LFPspectra.ofcLgamma = squeeze(mean(LFPspectra.ofc(:,q.Lgamma,:),2));
LFPspectra.ofcHgamma = squeeze(mean(LFPspectra.ofc(:,q.Hgamma,:),2));

    case 'wavelet'
LFPspectra.accdelta = squeeze(mean(abs(LFPspectra.acc(:,q.delta,:)),2));
LFPspectra.acctheta = squeeze(mean(abs(LFPspectra.acc(:,q.theta,:)),2));
LFPspectra.accalpha = squeeze(mean(abs(LFPspectra.acc(:,q.alpha,:)),2));
LFPspectra.accbeta = squeeze(mean(abs(LFPspectra.acc(:,q.beta,:)),2));
LFPspectra.accLgamma = squeeze(mean(abs(LFPspectra.acc(:,q.Lgamma,:)),2));
LFPspectra.accHgamma = squeeze(mean(abs(LFPspectra.acc(:,q.Hgamma,:)),2));

LFPspectra.ofcdelta = squeeze(mean(abs(LFPspectra.ofc(:,q.delta,:)),2));
LFPspectra.ofctheta = squeeze(mean(abs(LFPspectra.ofc(:,q.theta,:)),2));
LFPspectra.ofcalpha = squeeze(mean(abs(LFPspectra.ofc(:,q.alpha,:)),2));
LFPspectra.ofcbeta = squeeze(mean(abs(LFPspectra.ofc(:,q.beta,:)),2));
LFPspectra.ofcLgamma = squeeze(mean(abs(LFPspectra.ofc(:,q.Lgamma,:)),2));
LFPspectra.ofcHgamma = squeeze(mean(abs(LFPspectra.ofc(:,q.Hgamma,:)),2));

    case 'hilbert'
LFPspectra.accdelta = squeeze(LFPspectra.acc(:,1,:));
LFPspectra.acctheta = squeeze(LFPspectra.acc(:,2,:));
LFPspectra.accalpha = squeeze(LFPspectra.acc(:,3,:));
LFPspectra.accbeta = squeeze(LFPspectra.acc(:,4,:));
LFPspectra.accLgamma = squeeze(LFPspectra.acc(:,5,:));
LFPspectra.accHgamma = squeeze(LFPspectra.acc(:,6,:));

LFPspectra.ofcdelta = squeeze(LFPspectra.ofc(:,1,:));
LFPspectra.ofctheta = squeeze(LFPspectra.ofc(:,2,:));
LFPspectra.ofcalpha = squeeze(LFPspectra.ofc(:,3,:));
LFPspectra.ofcbeta = squeeze(LFPspectra.ofc(:,4,:));
LFPspectra.ofcLgamma = squeeze(LFPspectra.ofc(:,5,:));
LFPspectra.ofcHgamma = squeeze(LFPspectra.ofc(:,6,:));
 
% for filtered LFP output
LFPspectra.lfpfiltacc = lfpfiltacc;
LFPspectra.lfpfiltofc = lfpfiltofc;


end

%%  For plotting annotated events
 if nargin>5 && isfield(varargin{1},'events')
     Behavior=varargin{1}; 
     Bh=1;
 
% Create linemarkers at Behavioral events for visualization.  
% Plotting variables
X = repmat(Behavior.time',2,1);  
xx=(Behavior.time);
Y= repmat([0 105]',1,length(X));

 else 
     Bh=0;
 end


if plot_var == 1
disp('only plotting first trial if trial avg is off')


figure; 
ax(1)=axes('Position',[.05 .87 .83 .1]); %for LFP

%PLOT ACC LFP =====================================
             plot(LFP.acc(start_ind:fin_ind,1)); ylim([-0.3 0.3])
%              xlim([0 t_lfp(end)])
% set(a1,'XTickLabelMode','manual')
% set(a1,'XTickLabel',{[]})

            
%   PLOT SPECTROGRAM    =====================================     
% Set color limits same for ACC and oFC
  %set color lim based on percentile 
    clims = prctile(plotSaccL(:),[1 99]);
   lim_color1 = clims(1);
   lim_color2 =  clims(2);



         ax(2) = axes('Position',[.05 .58 .9 .3]); %for imagesc              
    imagesc(tl,fq1,(plotSaccL(:,:,1))'); colorbar; colormap('jet'); set(gca,'ydir','normal'); 
    max(tl)
            caxis([lim_color1  lim_color2]);
title('ACC'); 
                            
                             if any(strcmp(specgram_method,{'wavelet'}))
                                AX = gca;
                                AX.YTick = log2(YTicks);
                                AX.YTickLabels = YTicks;
                                ylim([0 max(fq1)]);
                            
                            end
                                    if Bh==1
                                    % set line and texts of events
                                    line(X,Y,'color','k'); text(xx,linspace(1,100,length(xx)),Behavior.events,'color','k');
                                    end

            
ax(3)=axes('Position',[.05 .40 .83 .1]);
 %plot the OFC LFP   =====================================
             plot(LFP.ofc(start_ind:fin_ind,1)); ylim([-0.3 0.3])
%              xlim([0 t_lfp(end)])
%              set(a2,'XTickLabelMode','manual')
% set(a2,'XTickLabel',{[]})s
             title('OFC'); 
%    PLOT SPECTROGRAM          =====================================
             ax(4) =  axes('Position',[.05 0.1 .9 .3]); %for imagesc  
    imagesc(tl,(fq1),(plotSofcL(:,:,1))'); colorbar; colormap('jet'); set(gca,'ydir','normal');   
            caxis([lim_color1  lim_color2]);

           
                     if Bh==1
                    % set line and texts of events
                    line(X,Y,'color','k'); text(xx,linspace(1,100,length(xx)),Behavior.events,'color','k');
                     end
                
                            if any(strcmp(specgram_method,{'wavelet'}))
                                AX = gca;
                                AX.YTick = log2(YTicks);
                                AX.YTickLabels = YTicks;
                                ylim([0 max(fq1)]);
                            end
                            
                            
                            set(gcf,'Position',[440 150 727 655])
                             linkaxes(ax,'x');
                             
%                             make widths the same 
%             ax(1).Position(3) = ax(2).Position(3);
%              ax(3).Position(3) = ax(4).Position(3);
                            
else
    disp(['Done computing spectrograms with ' specgram_method])
       
end







end %function



             