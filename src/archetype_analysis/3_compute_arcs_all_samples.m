% finding archetypes using ParTI - all samples
%%
clear; clc

% --- FIX: use this script's directory as root instead of hard-coded cd ---
thisScript = mfilename('fullpath');
scriptDir  = fileparts(thisScript);
cd(scriptDir);  % now we're in .../preprocessing

% --- FIX: robust addpath using fullfile relative to scriptDir ---
addpath(fullfile(scriptDir, '..', 'ParTI_uri_alon_lab'), ...
        fullfile(scriptDir, '..', 'ParTI_uri_alon_lab', 'PCHA'), ...
        fullfile(scriptDir, '..', 'ParTI_uri_alon_lab', 'ADVMM_and_SDVMM_codes'), ...
        fullfile(scriptDir, '..', 'ParTI_uri_alon_lab', 'sisal_demo'), ...
        fullfile(scriptDir, '..', 'ParTI_uri_alon_lab', 'SeDuMi_1_3'));

% input / output paths (relative)
countsPath = fullfile(scriptDir, '..', 'data', 'tcga_log_normalized_counts.csv');
outPath    = fullfile(scriptDir, '..', 'data', 'tcga_counts_all_threeArcs.csv');

% sanity check
if ~isfile(countsPath)
    error('Counts file not found: %s', countsPath);
end

% read counts (genes as row names, samples as columns)
counts = readtable(countsPath, 'ReadRowNames', true);


[arcs, data] = find_Arcs(countsPath);  


k = width(data) - width(counts);  

% write output
writetable(data, outPath, 'WriteRowNames', true);

fprintf('Wrote table with archetypes to: %s\n', outPath);