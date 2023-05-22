function [BP, featurenames,LDAcoef]  = load_jsons_to_timetables(folderPath)

PATIENTID={'CP1','CP2','CP3','CP4'};

dataframe_path = [folderPath(1:end-6) 'PYTHON'];
addpath(genpath(dataframe_path));

for x = 1:4
fileName = [PATIENTID{x} '_bndpwrbl.json']; % filename in JSON extension


% for Bandpwrs
str = fileread(fullfile(dataframe_path,fileName)); % dedicated for reading files as text
tmpdata = jsondecode(str); % Using the jsondecode function to parse JSON from string
  bphold = struct2table(tmpdata);

% for LDA feature coefs 
featFN =  [PATIENTID{x} '_chronic_nrs_LDA_coefs.csv'];
tmpcoef = readtable(featFN); 
featcoef = mean(table2array(tmpcoef(:,2:end)));
LDAcoef{x} = array2table(featcoef,'VariableNames',tmpcoef.Properties.VariableNames(2:end));


tt=datetime(bphold.righttimes); %pain times

% find the variables to keep
varnames  = bphold.Properties.VariableNames;
% separate the brain regions to keep order consistent
idxBPracc = find(contains(varnames,'rightacc'));
idxBProfc = find(contains(varnames,'rightofc'));
idxBPlacc = find(contains(varnames,'leftacc'));
idxBPlofc = find(contains(varnames,'leftofc'));
idxNRS = find(strcmp(varnames,'nrs')); 
% Sort the BP names and place after NRS

featurenames{x}= varnames([idxNRS,idxBPracc,idxBProfc,idxBPlacc,idxBPlofc]);
pp=bphold(:,[idxNRS,idxBPracc,idxBProfc,idxBPlacc,idxBPlofc]); %all features including nrs
pp.time =tt;

if iscell(pp.nrs) %kill all  empty rows
    idxkill = cellfun('isempty',pp.nrs);
    pp(idxkill,:) = [];
    pp.nrs = cell2mat(pp.nrs);

end

BP{x} = sortrows(table2timetable(pp));

end

