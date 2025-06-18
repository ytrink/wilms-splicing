% finding archetypes using ParTI - all samples
% 

%% 
clear; clc


% make sure you are in preprocessing folder


cd 'G:\Yaron\TARGET_paper\11_03_2024_scripts_for_revision\preprocessing'



%% 
addpath('../ParTI_uri_alon_lab','../ParTI_uri_alon_lab/PCHA','../ParTI_uri_alon_lab/ADVMM_and_SDVMM_codes',...
   '../ParTI_uri_alon_lab/sisal_demo','../ParTI_uri_alon_lab/SeDuMi_1_3')

counts  = readtable('../data/tcga_log_normalized_counts.csv','ReadRowNames',true);
filePath = '../data/tcga_log_normalized_counts.csv';

%%

expressionArray=table2array(counts)';

[arc,arcOrig,pc,errs,pval]=ParTI(expressionArray);






%% three archetypes

% from alon website: We recommend log-transforming gene expression data 
%(e.g. log-mRNA abundance or log-fold change) prior to analysis with ParTI.


% select three archetypes when prompted
% if you want to run using the matrix you created then change the filepath

[arcs,data] = find_Arcs(filePath); 

sampleNames = [counts.Properties.VariableNames {'Arc1','Arc2','Arc3'}];

data.Properties.VariableNames = sampleNames;

writetable(data,'../data/tcga_counts_all_threeArcs.csv','WriteRowNames',true);


