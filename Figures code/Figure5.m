% For Nature Neuroscience article: 
% Shirvalkar et al, 2023
% Figure 5
%

% 
% To answer Q2 from Reviewer 2 re time dynamics
% Underlying neurophysiology. What is missing is more sophisticated neurophysiological
% interpretation of why these signals can be decoded. What aspect of the signal (which frequency
% bands, increases/decreases) was used by decoders in a given area? Do these power changes appear in
% bursts or transiently or is there a persistent increase/decrease over the 30s? Are the authors
% proposing that if a feature is important (Fig. 4B), power in this band is chronically elevated
% over such a long period of time? Given there are only 4 subjects, something like Fig 2C could be
% plotted for each subject, averaging across the features in a weighted form according to which were
% used. In the absence of an interpretation of the features that the decoders rely on I find it hard
% to interpret the decoding result. As an example, in Fig S2 it is noted that “summarizes the
% information in neural features into a latent state, which is then provided…“. But what neural
% features end up being used?
%
%
%
% Look at the time dynamics  over 30 seconds 
%
%
% March 2023 
% Prasad Shirvalkar MD, PHD 
% UCSF
% 
%%  LOAD data input signals, filter them based on BP time, and Look at the power time series 
% QST spectograms  are precomputed,  here we only recompute chronic Specgrams

%go to main code directory and add subfolders to path
painfield = {'mayoNRS','mayoNRS','mayoNRS','painVAS','painVAS','painVAS','painVAS'};
activityname = 'home';


% COMPARE TO THE JSON Bandpower tables (from LDA Analysis) and take only the timepoints with Righttimes = json times
[BP, featurenames]= load_jsons_to_timetables(folderPath);


params.winstep = [1 0.1]; %same as for QST data (total samples should = 10* QST samples)
% WINSTEP should be [0.5 0.05];
params.tapers = [3 5];  % same as QST data
params.pad=1;
params.fpass=([0 100]);
fs=420.1; %actual sampling rate

for x = 1:4

    if x==1
        sidenames= {'right'};
    else
        sidenames = {'left','right'};
    end

    PATIENTID = ['CP' num2str(x)];
    load([PATIENTID '_home_LFPbilateral.mat']);
    disp(['Data loaded for ' PATIENTID])
    disp(['Most recent date included is: ' char(LFPmetabl.(LFPmetabl.refside).time(end))])

    bptimes = BP{x}.time;

    for s=1:numel(sidenames)
        %    Define the raw LFPs to use before filtering
        inputLFP{x,s}.acc = LFPbl.([sidenames{s} 'acc']);
        inputLFP{x,s}.ofc = LFPbl.([sidenames{s} 'ofc']);

        % =================================================x
        %         FILTER THE LFPs to KEEP
        %         index recordings to keep if they match the  time in bptimes
        [idxspec,~]= ismembertol(datenum(LFPmetabl.right.time'),datenum(bptimes),1,'DataScale',1/24/30);

        if  ~(sum(idxspec) == numel(bptimes))
             [idxspec,~]= ismembertol(datenum(LFPmetabl.left.time'),datenum(bptimes),1,'DataScale',1/24/30);
        end

        if x==1 %for CP1 take the only first 30 seconds to match other patients since some clips are slightly longer 
            maxdur  = 12604; %(fs *30s)
            %         LFPs to keep
            unsortedlfp(x).([sidenames{s} 'acc'])  = inputLFP{x,s}.acc(1:maxdur,idxspec);
            unsortedlfp(x).([sidenames{s} 'ofc'])  = inputLFP{x,s}.ofc(1:maxdur,idxspec);
        else
            unsortedlfp(x).([sidenames{s} 'acc'])  = inputLFP{x,s}.acc(:,idxspec);
            unsortedlfp(x).([sidenames{s} 'ofc'])  = inputLFP{x,s}.ofc(:,idxspec);
        end
        disp('Filtering done...')
        % =================================================
        % REAL SIGNAL
        holdLFP{x,s}.acc = unsortedlfp(x).([sidenames{s} 'acc']); %for computing spectrogram need to combine acc and ofc vars
        holdLFP{x,s}.ofc = unsortedlfp(x).([sidenames{s} 'ofc']);

        % Calculate the spectrogram
        LFPspecgram = spectrogram_compute(holdLFP{x,s},fs,'mt',[],0,params);
        specgrams.(PATIENTID).(sidenames{s}) = LFPspecgram;
    end
end

disp('done w all spectrograms')


%% Sort specgrams by fieldnames 
% to keep and organize them similar to the bandbwrs, then plot the images with Z scored mean label next to it

clear hold* new*
clc

for x = 1:4
    PATIENTID = ['CP' num2str(x)];
    if x==1
        sidenames= {'right'};
    else
        sidenames = {'left','right'};
    end

    for y= 1: numel(sidenames)
        disp(['CP' num2str(x) ' ' sidenames{y}  ])
        bands = specgrams.(PATIENTID).(sidenames{y}).bands;

        % =================================================
        %         SORT THE SPECGRAMS by side/ feature
        % =================================================
        for b = 1:numel(bands)
            holdspec(x).([sidenames{y} 'acc' bands{b}]) = specgrams.(PATIENTID).(sidenames{y}).(['acc' bands{b}]);
            holdspec(x).([sidenames{y} 'ofc' bands{b}]) = specgrams.(PATIENTID).(sidenames{y}).(['ofc' bands{b}]);
            % ============ ACC  ====================================
            newspec(x).([sidenames{y} 'acc' bands{b}])  = holdspec(x).([sidenames{y} 'acc' bands{b}]);
        end

        % **Keep  ACC and OFC for loops separate to keep the featurenames in the correct order
        for b= 1:numel(bands)
            % ============ OFC  ===================================
            newspec(x).([sidenames{y} 'ofc' bands{b}])  = holdspec(x).([sidenames{y} 'ofc' bands{b}]);
        end

        ft(x).fq = specgrams.(PATIENTID).right.fq;
        ft(x).times = specgrams.(PATIENTID).right.trialtimes(specgrams.(PATIENTID).right.trialtimes<=30);
    end

end
clear hold*
disp('done sorting')


%%  PLOT THE TOP 5 features time series for each patient 

% load preprocessed data
load powertimeseries.mat
load qstpowertimeseries.mat
load featweights.mat 
feattorun = 'top5';


% define sem fxn.
sem = @(n) nanstd(n)./ sqrt(size(n,1));
fs= size(newspec(2).leftaccalpha,1)/30; %291 samples in 30 s 
qstfs =  size(QSTspecmt(1).ipsi.rightaccalpha,1)/3; %60 samples in 3 seconds 

featnum=[12,24,24,24];
group = {'chronic','qstipsi'};






for g=1:2
%     1:numel(group)
    figure
t=tiledlayout('flow','TileSpacing','tight','Padding','none');
for  x=1:4 

    
   % for QST pain scores s
pkind  = strcmp(QSTbehavior(x).ipsi.events,'peakstart');
qstnrs = QSTbehavior(x).ipsi.pain(pkind);

    hipain = BP{x}.nrs>= median(BP{x}.nrs);
    lopain = ~hipain;
 
    qsthipain = qstnrs>= nanmedian(qstnrs);
    qstlopain =  qstnrs< nanmedian(qstnrs);  %to avoid some nans
 
    if strcmp(feattorun,'all')
        numfeats = featnum(x);
    else strcmp(feattorun,'top5')
    numfeats = numel(top5{x}); 
    end

   


  for f=1:numfeats

%       FOR CHRONIC PAIN DATA 
if g==1 
        feature = top5{x}{f};
        TOPweight = top5weights;
hilines = newspec(x).(feature)(:,hipain)';
lolines = newspec(x).(feature)(:,lopain)';
time1=linspace(0,30,size(hilines,2));
fs= size(newspec(2).leftaccalpha,1)/30; %291 samples in 30 s 


        % FOR ACUTE PAIN affected
elseif g==2
        feature = top5Ai{x}{f};
        TOPweight =top5weightsAi;
hilines = QSTspecmt(x).ipsi.(feature)(1:end-5,qsthipain)';     % to subtract last 5 samples, make end-5 , due to edge artifact
lolines = QSTspecmt(x).ipsi.(feature)(1:end-5,qstlopain)';
time1= linspace(0,3,size(hilines,2));
fs =  size(QSTspecmt(1).ipsi.rightaccalpha,1)/3; %60 samples in 3 seconds
end

%  GROUP AVERAGED CURVES
% Now calculate when the higher avg line crosses above the other and for how long 
if TOPweight{x}(f)>0 %the hipain power > lowpain power 
newcurve = nanmean(hilines) - nanmean(lolines); %Do this and calculate $zero crossings, avg duration of zero crossings
trialcurves =  (hilines) - nanmean(lolines); % for each trial compared to mean of other pain class   (always in terms of hilines # of trials)
elseif TOPweight{x}(f) <0
newcurve = nanmean(lolines) - nanmean(hilines); 
trialcurves =  -( (hilines) - nanmean(lolines) ); % for each trial compared to mean of other pain class (always in terms of hilines # of trials)
end

% CALCULATE STATS ON THE Averaged CROSSINGS
[idx_over0] =  find(newcurve > 0);
%     # Initialize values in identified period
bursting = zeros(size(newcurve));
bursting(idx_over0) = 1;



% TRIAL WISE CURVES
[idx_over0trial] =  find(trialcurves > 0);
%     # Initialize values in identified period
burstingtrial = zeros(size(trialcurves));
burstingtrial(idx_over0trial) = 1;

    change = diff(bursting);
    [idcs] = find(~(change==0));
    idcs = idcs+1;%   # Get indices following the change.
    idcs= idcs';

        clear  starts ends durations

        if bursting(1)
            %         # If the first sample is part of a burst, prepend a 1.
            idcs = [1; idcs];
        end
        if bursting(end)
            %         # If the last sample is part of a burst, append an index corresponding
            %         # to the length of signal.
            idcs = [idcs; length(bursting)];
        end     

% **** TRIAL CALCS! this takes the  high pain trials, but calculation above is based on feature
% weight (+ or -ve) 
numtrials = size(burstingtrial,1);

for t=1:numtrials
%     organized with trials in rows, features in cols
changetrials = diff(burstingtrial(t,:));
    [idcstrial] = find(~(changetrials==0));
  idcstrial = idcstrial+1; %   # Get indices following the change.
idcstrial=idcstrial';
   
           if burstingtrial(t,1)>0
            %         # If the first sample is part of a burst, prepend a 1.
            idcstrial = [1; idcstrial];
     
           end
        if burstingtrial(t,end)
            %         # If the last sample is part of a burst, append an index corresponding
            %         # to the length of signal.
            idcstrial = [idcstrial; length(bursting)];
        end     

        starttrial = idcstrial(1:2:end);
        endtrial = idcstrial(2:2:end);
        durationtrial{g,x,f}{t} = (endtrial-starttrial) ./fs;
        
%        For individual trials, compared to avg of the 'other' group (high vs low pain)
    trial.n{g,x}(f,t) = numel(durationtrial{g,x,f}{t});
    trial.meandur{g,x}(f,t) = nanmean(durationtrial{g,x,f}{t});
    trial.percent{g,x}(f,t) = 100 * sum(burstingtrial(t,:))/ length(burstingtrial(t,:));
    trial.tottime{g,x}(f,t) = sum(burstingtrial(t,:)) / fs;
end


    starts = idcs(1:2:end);
    ends = idcs(2:2:end);
    durations = (ends - starts) ./ fs; 
      
%     For trial avg in high/low pain groups
    n_bursts{g}(x,f) = numel(durations);
    duration_mean{g}(x,f) = nanmean(durations);
    percent_burst{g}(x,f) = 100 * sum(bursting) / length(bursting);
    time_w_bursts{g}(x,f)  = sum(bursting) / fs;
    burstsamp{g}(x,f) =   sum(bursting);
    lensamp{g}(x,f) =   numel(bursting);

% PLOTS
nt = nexttile;
hold on 

shadedErrorBar(time1,hilines,{@nanmean,sem},'lineprops','r');
shadedErrorBar(time1,lolines,{@nanmean,sem},'lineprops','b');
T= title(feature);
T.FontSize = 14;
ymin  = min(nt.YLim);
ymax  = max(nt.YLim);

txtstr = sprintf('bouts above = %d \n mean dur above = %0.1f sec \n total dur above = %0.1f sec',n_bursts{g}(x,f), duration_mean{g}(x,f),time_w_bursts{g}(x,f));
text(0,0.1,txtstr,'Units','normalized')
  end

end
sgtitle(group{g})
set(gcf,'Position',[1000*g 672 1365 825])
end

%% get TRIAL AVGed stats for above metrics
clc
x=1:4
a=1:2
% for x=1:4
% fprintf('\n \n CP%d \n',x)

fprintf('Avg Total time \n')
% Avg total time 
disp('CHRONIC')
chronicpt  = percent_burst{1}(x,:);
percenttottime =  nanmean(chronicpt(:))
stdtottime = nanstd(chronicpt(:))
fprintf('\n')

disp('ACUTE')
acute2pt =percent_burst{2}(a,:);
percenttottime =  nanmean(acute2pt(:))
stdtottime = nanstd(acute2pt(:))

% Avg mean time above
fprintf('\n')
fprintf('Avg mean time \n')
disp('CHRONIC')
chronicpt  = duration_mean{1}(x,:);
meantime = nanmean(chronicpt(:))
stdmeantime = nanstd(chronicpt(:))
fprintf('\n')

disp('ACUTE')
acute2pt =duration_mean{2}(a,:);
meantime = nanmean(acute2pt(:))
stdmeantime = nanstd(acute2pt(:))
% end


% total  dur
fprintf('\n \n Total Dur\n')
disp('CHRONIC')
chronicpt  = time_w_bursts{1}(x,:);
meantime = nanmean(chronicpt(:))
stdmeantime = nanstd(chronicpt(:))
fprintf('\n')

disp('ACUTE')
acute2pt =time_w_bursts{2}(a,:);
meantime = nanmean(acute2pt(:))
stdmeantime = nanstd(acute2pt(:))
% end

%mean bouts

% total  dur
fprintf('\n \n num bouts per sec\n')
disp('CHRONIC')
chronicpt  = n_bursts{1}(x,:);
meanfq= nanmean(chronicpt(:))./30
stdfq = nanstd(chronicpt(:)./30)
fprintf('\n')

disp('ACUTE')
acute2pt =n_bursts{2}(a,:);
meanfq = nanmean(acute2pt(:))./3
stdfq = nanstd(acute2pt(:)./3)





%% Chi2 test on raw data for proportions of samples >0 

% clc
for x=1:4
fprintf('\n \n CP%d \n',x)
        c1 = sum(burstsamp{1}(x,:));    C1 =  sum(lensamp{1}(x,:));
       a1 = sum(burstsamp{2}(x,:)); A1 = sum(lensamp{2}(x,:));
       x1 = [repmat('c',C1,1); repmat('a',A1,1)];
       x2 = [repmat(1,c1,1); repmat(2,C1-c1,1); repmat(1,a1,1); repmat(2,A1-a1,1)];
       [tbl,chi2stat,pval(x),labels] = crosstab(x1,x2)
  

end

pcor = fdr(pval);


%% Compare the  TRIALWISE  % above etc metrics (ttests) from acute to chronic for CP1-2 across trials
% Trialwise compared to the mean of the other condition


clc

clear p pallcor pcor


% (divide chronicdata by 30  and acutedata by 3  to get # of increases or decreases per second)
usemetric{1}  = trial.percent;
usemetric{2}  = trial.meandur;
usemetric{3}  = trial.n;

    
titlestr{1} = '% time';
titlestr{2} = 'mean duration';
titlestr{3} = 'num bouts';

max1{1} = 100;
max1{2} = 10;
max1{3} = 3;



blues  =brewermap(9, 'Blues'); 
metric  = {'percent','meandur'};
bands ={'delta','theta','alpha','beta','Lgamma','Hgamma'};
regions ={'rightacc','rightofc','leftacc','leftofc'};
markertype = {'o','s','+','^'};
regionind = {1:6,7:12,13:18,19:24}; 
xind =  {[1 1.8],[3 3.8],[5 5.8],[7 7.8]};

figure
for u=1:3
    nexttile
    for  x=1:4

        if strcmp(feattorun,'all')
            regions = (allfeats{x});
            pcor = nan(24,4);
            featnum=[12,24,24,24];
        elseif strcmp(feattorun,'top5')
            regions = top5{x};
            pcor = nan(5,4);
            featnum=[5,5,5,5];
        end


        hold  on


        if u==3 %normalize by seconds

            chronicdata = usemetric{u}{1,x}'./30;
            acutedata = usemetric{u}{2,x}'./3;

        else

            chronicdata = usemetric{u}{1,x}';
            acutedata = usemetric{u}{2,x}';

        end



%         Collect the number of samples for reporting 
numsamples_chronic(u,x) = numel(chronicdata(:));
numsamples_acute(u,x) = numel(acutedata(:));


        X1= xind{x}(1);
        Y1=chronicdata(:);
        X2 = xind{x}(2);
        Y2= acutedata(:);


        % Chronic
       
        H=notBoxPlot(Y1,X1,'jitter',0.4);
        set([H.data],...
            'MarkerFaceColor',[1,1,1]*0.25,...
            'markerEdgeColor',[1,1,1]*0.25,...
            'MarkerSize',5,...
            'Marker','.')
        set([H.mu],'color','w')
        set([H.semPtch],...
            'FaceColor', 'none',...
            'EdgeColor','none')
        set([H.sdPtch],...
            'FaceColor','none',...
            'EdgeColor','none')
        hold on
         boxplot(Y1,'positions',X1,'whisker',5,'BoxStyle','filled')
                 h=findobj(gca,'tag','Outliers'); % Get handles for outlier lines.
       set(h(1),'MarkerEdgeColor','w'); % Change color for one group

        if x<3 %only for patients with sig acute decoding  (CP1 and 2)
            % Acute
                    
            H2=notBoxPlot(Y2,X2,'jitter',0.4);
            set([H2.data],...
                'MarkerFaceColor',[1 0 0],...
                'markerEdgeColor',[1 0 0 ],...
                'MarkerSize',5,...
                'Marker','.')
                 set([H2.mu],'color','w')
            set([H2.semPtch],...
                'FaceColor', 'none',...
                'EdgeColor','none')
            set([H2.sdPtch],...
                'FaceColor','none',...
                'EdgeColor','none')
        hold on
        boxplot(Y2,'positions',X2,'whisker',5,'BoxStyle','filled')
        h=findobj(gca,'tag','Outliers'); % Get handles for outlier lines.
       set(h(1),'MarkerEdgeColor','w'); % Change color for one group
 
                    hold on
       
       


        end

        ylim([0 max1{u}])


        % Across all features
        [p(x),~,STATS{x}] = ranksum(chronicdata(:),acutedata(:));

        title(titlestr{u})

    end
    [~,pcor,padj] = fdr_calc(p)
    pfdr{u} = pcor;
    statsout{u} = STATS;

end

set(gcf,'Position',[1640 1014 900 683])


