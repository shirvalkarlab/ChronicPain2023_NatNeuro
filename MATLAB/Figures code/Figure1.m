% For Nature Neuroscience article:
% Shirvalkar et al, 2023,
% 
% 
% Figure 1C,D,E
% 1C- Histograms of pain NRS
% 1D- NRS / VAS correlations
% 1E - partial autocorrelations
% 
% 
% Prasad Shirvalkar MD, PhD
% UCSF
% March 2023



PATIENTID= {'CP1','CP2','CP3','CP4'};
load painsurveys.mat


%% Plot the histograms for all participants' NRS pain distributions
% plot all histograms in one figure

close all
EDGES = [0:0.6:10];
plotopt = {'mayoNRS','mayoNRS','mayoNRS','mayoNRS'};
k=[1,1,1,1];
bincounts = zeros(4,numel(EDGES)-1);
painvals = [];

%  process each patient
for x=1:4
    holdpainvals = painsurveytable.(PATIENTID{x}).(plotopt{x});
    holdpainvals = holdpainvals(~isnan(holdpainvals));

    painvals1 = [holdpainvals(:),repmat(x,numel(holdpainvals),1)];  %one column pain vals,  second ptID
    painvals = [painvals;painvals1];

    [bincounts(x,:),~] =  histcounts(holdpainvals,EDGES,'Normalization','pdf');

    % SHADED histogram (1)
    hold on
    plot_histogram_shaded(holdpainvals,'edges',EDGES(2:end));
end


set(gca,'XTick',[-0.5:1:9.5])
set(gca,'XTickLabel',[0:1:10])
xlabel('Pain NRS')
ylabel('Probability')
title('Fig 1 C- Histograms of all patients'' pain NRS')


%% Plot matrix of NRS vs VAS scores (see if there is a high Correlation) 

figure

paintype_VAS= {'painVAS','unpleasantVAS'};
paintype_NRS = {'mayoNRS','unpleasantNRS'};
numplots = numel(paintype_VAS);
filledstyles={'filled','filled'};

allvas = cell(numplots,1);
allnrs = cell(numplots,1);

for x=1:4
    %    collect all subjects' data
    for p=1:numplots
        allvas{p} = [allvas{p}; painsurveytable.(PATIENTID{x}).(paintype_VAS{p})(:)];
        allnrs{p} = [allnrs{p}; painsurveytable.(PATIENTID{x}).(paintype_NRS{p})(:)];
    end
end


for p =1:numplots
    exclidx = (allvas{p} == 0 | allvas{p} == 100 | allnrs{p} <=1 | isnan(allvas{p}) );  % exclude error values
    allvas{p}(exclidx)=[];
    allnrs{p}(exclidx)=[];

    scatter(allvas{p},allnrs{p}-(0.05*(p-1)),50,filledstyles{p})
    hold on
end

legend({'pain intensity','pain unpleasantness'})
line([0 100], [0 10],'Color','k','LineStyle','--')
xlim([0 100])
ylim([0 10])
xlabel('Visual Analog Score')
ylabel('Numerical Rating Score')

title('Fig 1 D- VAS vs NRS correlations')

% calculate the Pearson's R for  VAS vs NRS
fprintf(' PEARSON''S R = \n')
for p = 1:2
    nanBS= isnan(allnrs{p});
    allvas{p}(nanBS)=[];
    allnrs{p}(nanBS)=[];

    disp(paintype_VAS{p})
    corrcoef(allvas{p},allnrs{p})
end



%% Plot the partial autocorrelation (in hours) of each patients pain scores, after resampling
paintype = 'mayoNRS';
% plot 72 hours worth of data
figure

for x=1:4
    tt=painsurveytable.(PATIENTID{x}).times; %pain times
    pp=painsurveytable.(PATIENTID{x}).(paintype); %pain score

    % Resample to uniform time using pchip
    [Y,Ty] = resample(pp,tt,'pchip');

    %get the time lag in hours
    t_int = mean(diff(Ty));
    t_intnum = datevec(t_int);
    Hint = t_intnum(4)+(t_intnum(5)/60);
    %round off to nearest 0.1
    Hint = round(Hint*10)/10;
    numlags = round(72/Hint);


    % PLOT ALL PATIENTS ON ONE FIGURE
    [pacf1{x},lags{x},bounds(:,x)] =parcorr(Y,'NumLags',numlags);
    [PARCOR,sig1,cil,ciu]= pacf(Y',numlags);
    idx_pac = find(abs(pacf1{x})>bounds(1,x));
    idx_pac = idx_pac(2:end); %ignore pacf of 1 (zero lag)

    ax = gca;
    pac_xvals = round(0:Hint*10:Hint*numlags*10)/10;
    allstem(x) = stem(pac_xvals,pacf1{x});
    hold on

    % make the significant values BOLD
    idx_bold = abs(pacf1{x})>bounds(1,x);
    stem(pac_xvals(idx_bold),pacf1{x}(idx_bold),'LineWidth',2,'Color',allstem(x).Color);

end
set(ax,'XTick',0:8:72)
legend({'CP1','','CP2','','CP3','','CP4'});
xlabel('Hours')
title('Fig 1 E- partial autocorr of pain NRS')

