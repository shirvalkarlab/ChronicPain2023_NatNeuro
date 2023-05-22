% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Omid Sani, Parima Ahmadipour, Maryam Shanechi 
%   Shanechi Lab, University of Southern California, 2020
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script runs the cross-validated model fitting with the chronic pain data
% Change measureName, settings.predictor, and useLSSM variables for
% different decoding methods.

addpath(genpath('./tools'));

fileDir= 'G:/2023_Shirvalkar_NatNeuro_ChronicPainBiomarkers/MATLAB/' ;  % Path to MAT data files (change this based on your local machine)
   fDir = 'G:/2023_Shirvalkar_NatNeuro_ChronicPainBiomarkers/PYTHON/'; % Path to JSON files from UCSF

subjects = {'CP1', 'CP2', 'CP3', 'CP4'};
subjectsIndsToRun = 1:numel(subjects);
% subjectsIndsToRun = [3]; % Uncomment and modify to run only some subjects

measureNames = {'autoNRS', 'painVAS', 'MPQsum', 'unpleasantNRS', 'unpleasantVAS'};   % Pain score, can be "painVAS", "relativeVAS", "autoNRS", etc
measureNameIndsToRun = 1:numel(measureNames);
% measureNameIndsToRun = [3, 5]; % Uncomment and modify to run only some pain scores

for mi = measureNameIndsToRun
    measureName = measureNames{mi};
    
    for si = subjectsIndsToRun
        subject = subjects{si};
        fileName = sprintf('%s_home_LFPbilateral.mat', subject);
        
        % If want to remove all data not included in the json file from Prasad
     
        refFilePath = fullfile(fDir, sprintf('%s_bndpwrbl.json', sprintf('%s', fileName(1:3)))); % Path to json file from UCSF, to exactly keep the same recordings
        % refFilePath = '';   % Uncomment to keep all data points
        
        saveDir = './resultsF';
        if ~exist(saveDir, 'dir'), mkdir(saveDir); end
        savePath = fullfile(saveDir, sprintf('%s_%s.fig', subject, measureName));
        
        fitSettings = struct('savePath', savePath, 'numPerms', 10000);
        
        results = decodePain_chronicPain(fileDir, fileName, refFilePath, measureName, fitSettings);
        
        if ~isempty(results)
            savePath = fullfile(saveDir, sprintf('%s_%s.mat', subject, measureName));
            save(savePath, '-struct', 'results');
            fprintf('Saved results for %s %s to %s\n', subject, measureName, savePath);
            
            sheetFile = fullfile(saveDir, sprintf('summary.xlsx'));
            colLetter = char( 65 + find(contains(measureNames, measureName)) ); 
            rowNumber = find(contains(subjects, subject));
            % AUC table
            writecell({measureName}, sheetFile, 'Range', sprintf('%s%d', colLetter, 1));
            writecell({subject}, sheetFile, 'Range', sprintf('A%d', rowNumber+1));
            writecell({results.trueStats.perf}, sheetFile, 'Range', sprintf('%s%d', colLetter, rowNumber+1));
            % p-value table
            writecell({measureName}, sheetFile, 'Range', sprintf('%s%d', colLetter, 10));
            writecell({subject}, sheetFile, 'Range', sprintf('A%d', rowNumber+10));
            writecell({results.pValue}, sheetFile, 'Range', sprintf('%s%d', colLetter, rowNumber+10));
            fprintf('Added results for %s %s to %s\n', subject, measureName, sheetFile);
        end
    end
end
