% Run this script or individual cells to generate figures from the paper that rely on MATLAB scripts
% (except LSSM which is separate)
%
% You can also open scripts for individual Figures run them separately 
% (located in the 'Figures code' folder)
% 
% 
%  Below we add paths and call scripts for the individual figures
% 
% Prasad Shirvalkar MD PhD 
% UCSF
% May 2023
% 


%% add paths
% get this script's path and add the article files' to the working search path
clear
close all
clc

filePath = matlab.desktop.editor.getActiveFilename; 
[folderPath,~,~] = fileparts(filePath);
addpath(genpath(folderPath));
cd(folderPath);
%%  FIGURE 1

Figure1
 
%% FIGURES 2 and 3 
% for Figures 2 and 3, and most of the Supplementary Figures, please use the Python code 

%%  FIGURE 4

Figure4

%% FIGURE 5 (including related supplementary figures)

Figure5