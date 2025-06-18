function [arcs,geneExpressionArcs] = find_Arcs(filename)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% [file_ge,path_ge] = uigetfile('*.csv','Select gene expression matrix');
% %selpath = uigetdir('./','Pick folder to save file into');  % probably best
% %not to save in function but after
% %filenameOutput = inputdlg('Please provide a name for the output file','s');
% 
% filenameInput = [path_ge,file_ge];
expressionTable = readtable(filename,'ReadRowNames',true);

expressionArray=table2array(expressionTable)';

%run PARTIlite on DKFZ expression set, select 3 archetypes
[arc,arcOrig,pc]=ParTI_lite(expressionArray);
expressionArray = [expressionArray; arcOrig]; %concatenate archetypes to DKFZ array


arcTable=array2table(arcOrig');

%table of gene expression with archetypes                            
geneExpressionArcs=[expressionTable arcTable];                            
arcs = arcOrig;
%filename = [selpath filenameOutput]; filename = strjoin(filename);
%writetable(geneExpressionArcs,filename,'WriteRowNames',true);

end

