% For Nature Neuroscience article: 
% Shirvalkar et al, 2023, 
% Figure 4
%
%
% This script will compare feature importances / weights after LDA for each
% participant across chronic and acute pain models to determine WHICH REGION the best features come from

%Prasad Shirvalkar MD PhD UCSF
%use to compare ACC vs OFC among Chronic vs Acute pain
%March 2023    
 
% Define folder path for python outputs  (t extract 
Pypath = ['/Users/pshirvalkar/Dropbox (UCSF Department of Neurological Surgery)/SUBNETS Dropbox/Chronic Pain - Activa and Summit 2.0/DATA ANALYSIS/' ...
    'PCS_Paindecoding_Nature_Neuroscience_Revisions/classification_02_02_2022/pcs_exports'];



%% MAIN FIGURE 4 OFC-ACC diffs  - Get the difference in feature weights  between  ACC and OFC , compare to shuffled feature weight diffs
PATIENTID={'CP1','CP2','CP3','CP4'};

%%%%%%
edges=0:.01:1;
normtype = 'probability';
bands = {'delta','theta','alpha','beta','Lgamma','Hgamma'};
band_labels = {'delta','theta','alpha','beta','Lgamma','Hgamma'};
sides = {'right','left'};
leftnm = cellfun(@(x) sprintf('left%s',x),bands,'UniformOutput',false);
rightnm = cellfun(@(x) sprintf('right%s',x),bands,'UniformOutput',false);
regions = [rightnm, leftnm];
%%%%%%


allbandsCacc = [];
allbandsCofc = [];
allbandsAacc = [];
allbandsAofc = [];
Caccgroup = {};
Cofcgroup={};
Aaccgroup={};
Aofcgroup={};

numbands=numel(band_labels);
featdiff = cell(1,4);
featdiff_surr = cell(1,4);
featdiffA = cell(1,4);
featdiff_surrA = cell(1,4);

for x=1:4
    if x==1
        numsides = 1;
    else
        numsides = 2;
    end


    % import LDA coeffs for only significanto models (Acuteaffected and Chronic) 
    Atableaffected = readtable(fullfile(Pypath,[PATIENTID{x} '_acute_ipsi_pain_LDA_coefs.csv']));
    Ctable  = readtable(fullfile(Pypath,[PATIENTID{x} '_chronic_nrs_LDA_coefs.csv']));


    for s=1:numsides

        clear Ca* Co* Cm* Aa* Ao* Am*
        for b=1:numel(band_labels)
            holddiff=[];
            Aholddiff = [];


%           Take all bands from each brain region (e.g. LeftACC)
            Cacc = contains(Ctable.Properties.VariableNames,'acc') & contains(Ctable.Properties.VariableNames,band_labels{b}) & contains(Ctable.Properties.VariableNames,sides{s});
            Cofc = contains(Ctable.Properties.VariableNames,'ofc') & contains(Ctable.Properties.VariableNames,band_labels{b}) & contains(Ctable.Properties.VariableNames,sides{s});

            Aacc = contains(Atableaffected.Properties.VariableNames,'acc') & contains(Atableaffected.Properties.VariableNames,band_labels{b}) & contains(Atableaffected.Properties.VariableNames,sides{s});
            Aofc = contains(Atableaffected.Properties.VariableNames,'ofc') & contains(Atableaffected.Properties.VariableNames,band_labels{b}) & contains(Atableaffected.Properties.VariableNames,sides{s});

%           Take magnitude and  Normalize to make max = 1
% Chronic
            Caccmat = abs(table2array(Ctable(:,Cacc))); % contains chronic bands in one side ACC, all trials
            Cofcmat = abs(table2array(Ctable(:,Cofc)));  %Same for OFC 
                maxbothC  = max([Caccmat;Cofcmat]);
            Caccmat = Caccmat ./maxbothC; % NORMALIZING TO THE MAX ACROSS BOTH REGIONS! (to maintain relative scaling) so max = 1
            Cofcmat = Cofcmat./maxbothC; %** ? normalize?
            Caccmat = Caccmat(:);
            Cofcmat = Cofcmat(:);
            
% acute affected
            Aaccmat = abs(table2array(Atableaffected(:,Aacc))); % contains chronic bands in one side ACC, all trials
            Aofcmat = abs(table2array(Atableaffected(:,Aofc)));  %Same for OFC
                maxbothA  = max([Aaccmat;Aofcmat]);
            Aaccmat = Aaccmat ./ maxbothA; % NORMALIZING TO THE MAX so max = 1 ?? (What if we dont?)
            Aofcmat = Aofcmat./ maxbothA; %** ? normalize?
            Aaccmat = Aaccmat(:);
            Aofcmat = Aofcmat(:);



            %  Get trialwise differences in OFC- ACC feature weights for each brain region
            holddiff = Cofcmat - Caccmat;
            featdiff{x} = [featdiff{x},holddiff];

            Aholddiff = Aofcmat - Aaccmat;
            featdiffA{x} = [featdiffA{x},Aholddiff];


            %  Then shuffle the regions (ACC or OFC) and do same difference
            %  put all ofc and acc values into a vector, and then randomly choose some to be ACC and some to be
            %  OFC.

            Cmixmat = [Caccmat;Cofcmat];
            szCmix = size(Cmixmat,1);
            numsamp = size(Caccmat,1);

            Amixmat = [Aaccmat;Aofcmat];
            szAmix = size(Amixmat,1);
            Anumsamp = size(Aaccmat,1);

            holdsurr=[]; holdsurr2=[]; holdsurrA= []; holdsurrA2=[];

            for z=1:1000   % 1000 shuffles
                %      Chronic
                perm1 = randperm(szCmix,numsamp);
                perm2 = randperm(szCmix,numsamp);
                mixind = zeros(szCmix,1);
                mixind2 = zeros(szCmix,1);
                mixind(perm1) = 1;
                mixind2(perm2)=1;
                mixind=  logical(mixind);
                mixind2 = logical(mixind2);
                surrogateCofc = Cmixmat(mixind,:);
                surrogateCacc = Cmixmat(mixind2,:);
                holdsurr = surrogateCofc - surrogateCacc;
                holdsurr2 = [holdsurr2;holdsurr];

                % Acute affected
                Aperm1 = randperm(szAmix,Anumsamp);
                Aperm2 = randperm(szAmix,Anumsamp);
                Amixind = zeros(szCmix,1);
                Amixind2 = zeros(szCmix,1);
                Amixind(Aperm1) = 1;
                Amixind2(Aperm2)=1;
                Amixind=  logical(Amixind);
                Amixind2 = logical(Amixind2);
                surrogateAofc = Amixmat(Amixind,:);
                surrogateAacc = Amixmat(Amixind2,:);
                holdsurrA = surrogateAofc - surrogateAacc;
                holdsurrA2 = [holdsurrA2;holdsurrA];

            end

            featdiff_surr{x}= [featdiff_surr{x},holdsurr2];
            featdiff_surrA{x}= [featdiff_surrA{x},holdsurrA2];

        end
    end

    % PLOT the OFC-ACC difference vs shuffled
    edges = -1:0.05:1;
    h1 = figure(101);
    %CHRONIC
    n1(x) = nexttile;

    H = distributionPlot(featdiff{x},'showMM',3,'histOri','left','distWidth',0.9,'widthDiv',[2 1],'xNames',regions);
    H{2}.Color = 'b';
    H = distributionPlot(featdiff_surr{x},'showMM',3,'histOri','right','color',[0.5 0.5 0.5],'distWidth',0.9,'widthDiv',[2 2],'xNames',regions);   % SHUFFLED DATA PLOTTED IN GREY ON RIGHT!
        H{2}.Color = 'b';
    title(['CP' num2str(x)])
    ylim([-1 1])
    ylabel('OFC - ACC feature weight')


    % ACUTE
    h2 = figure(102);
    n2(x) = nexttile;
    H = distributionPlot(featdiffA{x},'showMM',3,'histOri','left','distWidth',0.9,'widthDiv',[2 1],'xNames',regions);
        H{2}.Color = 'b';
    H = distributionPlot(featdiff_surrA{x},'showMM',3,'histOri','right','color',[0.5 0.5 0.5],'distWidth',0.9,'widthDiv',[2 2],'xNames',regions);   % SHUFFLED DATA PLOTTED IN GREY ON RIGHT!
        H{2}.Color = 'b';
    title(['CP' num2str(x)])
    ylim([-1 1])
    ylabel('OFC - ACC feature weight')

    for y = 1:size(featdiff{x},2)
        [pval(y),~,stats1] = ranksum(featdiff{x}(:,y),featdiff_surr{x}(:,y),'tail','both');
        [pvalA(y),~,stats1a] = ranksum(featdiffA{x}(:,y),featdiff_surrA{x}(:,y),'tail','both');

    Zval(x,y) = stats1.zval;
    ZvalA(x,y) = stats1a.zval;
    

    end
 
%   Multiple comparisons correction with FDR
    fdrcorr = fdr(pval');
    fdrmatC(1:numel(fdrcorr),x) = fdrcorr;

    fdrcorrA = fdr(pvalA');
    fdrmatA(1:numel(fdrcorr),x) = fdrcorrA;



end %x

sgtitle(h1,'Chronic');
sgtitle(h2,'Acute affected')


fnames = Ctable.Properties.VariableNames(2:end);
fdrmask = zeros(size(fdrmatC));
fdrmask(fdrmatC>=0.05)=1;

% put asterisks for * sig
for x=1:4
    for f=1:12

        %     chronic sig*
        if fdrmatC(f,x)<0.05 && fdrmatC(f,x)>0 % to exclude empty vals for CP1
            text(n1(x),f-0.2,0,'*','FontSize',30,'Color','r')
        end



        fnames = Atableaffected.Properties.VariableNames(2:end);
        fdrmaskA = zeros(size(fdrmatA));
        % ACUTE SIG
        % put asterisks for * sig
        if fdrmatA(f,x)<0.05 && fdrmatA(f,x)>0 % to exclude empty vals for CP1
            text(n2(x),f-0.2,0,'*','FontSize',30,'Color','r')
        end

    end
end


