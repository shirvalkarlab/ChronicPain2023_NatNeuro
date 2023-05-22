This Repository contains original data and code to reproduce results and figures from the article: 

Shirvalkar, P., Prosky, J., Chin, G. et al. First-in-human prediction of chronic pain state using intracranial neural biomarkers. Nat Neurosci (2023). https://doi.org/10.1038/s41593-023-01338-z 
Open Access Link: https://rdcu.be/dcICQ
=====================================================================================
NOTE!!! Due to github's petulance, there are 2 *.mat files missing from the database on github. 
TO DOWNLOAD ALL CODE AND DATA FILES PLEASE USE THE FOLLOWING LINK: https://www.dropbox.com/sh/q0nyyoipudk28wu/AAAgUolShLa7_ddVdgL3n9kMa?dl=0
=====================================================================================
 
1. MATLAB based code and data is in the MATLAB folder including: 
- raw LFP and pain survey data + preprocessing and analysis for Figures 1, 4 and 5 and much of the Supplementary Results 
- code for the Linear State Space Modeling (LSSM) results 

* All code can be successfully run on MATLAB 2021b or 2022b with the following toolboxes:
'Signal Processing Toolbox'			'8.7'		
'Statistics and Machine Learning Toolbox'	'12.2'		
'Econometrics Toolbox'				'5.7'
'Image Processing Toolbox'			'11.4'	
'Statistics and Machine Learning Toolbox'	'12.2'	
'Curve Fitting Toolbox'				'3.6' 

* Other dependent files are provided in utilities folder
* 		TO START MAIN analysis : Run 'Main_script_all_analysis.m'
*		TO START LSSM analysis : See LSSM_README.txt in that folder






2. PYTHON code and data is in the PYTHON folder including: 
- all LDA based results in main and supplementary figures
- all LASSO based results in main and supplementary figures

* Note, to generate the difference based (Pain fluctuation) analysis, please set the variable in cell # 3 of the Jupyter notebook as follows: DIFFERENCE = True


==========================
Prasad Shirvalkar MD, PhD
UCSF 
May 20, 2023
