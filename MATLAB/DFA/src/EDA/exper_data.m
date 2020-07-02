%{
Experimental Data
Author: Surabhi Joshi

Table of Contents
Experimental Data
Summary
Create tables for averages and standard deviations.  
Normalize the Data
Separate the Data: GOT1 vs NT Data
Create Numeric Arrays of all Tables
Create the coefficient of variation, CV 
DFA Filter 1: CV thresholds + NaN: 0.1
DFA Filter 1: CV thresholds + NaN: 0.3
DFA Filter 1: CV thresholds: 0.6
DFA Filter 2: NaN Filter
%}

%% Summary
%This script reads in some experimental data and performs exploratory initial statistical analysis to understand how each condition affects the metabolomics before running DFA. 
% download the experimental data and create two tables, one for each
% cancer cell line TU8902 and HCT116.
experimental_data_TU8902 = readtable('Barb Nelson Pos Table-20180410 - HL.xlsx', ...
'Sheet', 'TU8902');
experimental_data_HCT116 = readtable('Barb Nelson Pos Table-20180410 - HL.xlsx', ... 
 'Sheet',"HCT116");
% delete row 2 and column 2 from both tables
experimental_data_TU8902.transition = [];
experimental_data_HCT116.transition = [];

%delete the second histidine row
experimental_data_HCT116(108,:) = [];
experimental_data_TU8902(108, :) = [];
experimental_data_HCT116(192,:) = [];
experimental_data_TU8902(192, :) = [];

%clear away all the extra clutter in the metabolites column of the above
%tables
experimental_data_TU8902.metabolite = strrep(experimental_data_TU8902.metabolite, 'RP_Pos_', '');
experimental_data_TU8902.metabolite = strrep(experimental_data_TU8902.metabolite, '_', '');

experimental_data_HCT116.metabolite = strrep(experimental_data_HCT116.metabolite, 'RP_Pos_', '');
experimental_data_HCT116.metabolite = strrep(experimental_data_HCT116.metabolite, '_', '');

%download metabolomics and clusterring data
metabolomics_TU8902 = readtable('Barb Nelson Pos Table-20180410 - HL.xlsx', ... 
    'Sheet', 'Clustering');
metabolomics_TU8902(1,:) = [];
metabolomics_TU8902(108,:) = [];
metabolomics_TU8902(192,:) = [];
metabolomics_TU8902 = table2array(metabolomics_TU8902);
metabolomics_TU8902 = string(metabolomics_TU8902);
metabolomics_TU8902(:, 1) = strrep(metabolomics_TU8902(:, 1), 'RP_Pos_', '');
metabolomics_TU8902(:, 1) = strrep(metabolomics_TU8902(:, 1), '_', '');

metabolomics_HCT116 = metabolomics_TU8902;


%have a categorical array of all the conditions present
all_conditions = ["shNT no dox"; "shNT 2d dox"; "shNT 4d dox"; "shGOT1 no dox"; "shGOT1 2d dox"; "shGOT1 4d dox"; "shNT no dox media"; "shNT 2d dox media"; "shNT 4d dox media"; "shGOT1 no dox media"; "shGOT1 2d dox media"; "shGOT1 4d dox media"];
got1_conditions = ["shGOT1 no dox"; "shGOT1 2d dox"; "shGOT1 4d dox"; "shGOT1 no dox media"; "shGOT1 2d dox media"; "shGOT1 4d dox media"];
nt_conditions = ["shNT no dox"; "shNT 2d dox"; "shNT 4d dox"; "shNT no dox media"; "shNT 2d dox media"; "shNT 4d dox media"];

%% Create tables for averages and standard deviations.  
%This here is going to take the averages and standard deviations from the original experimental data and fill two new tables with their values. Any NaN will be turned to 0 prior to finding the mean and standard deviations. 
Tabl_average_TU8902 = computeReplicateMetrics(experimental_data_TU8902, 'mean');
Tabl_average_TU8902(225,:) = [];
Tabl_std_TU8902 = computeReplicateMetrics(experimental_data_TU8902, 'stddev');
Tabl_std_TU8902(225,:) = [];
Tabl_average_HCT116 = computeReplicateMetrics(experimental_data_HCT116, 'mean');
Tabl_std_HCT116 = computeReplicateMetrics(experimental_data_HCT116, 'stddev');

%% Normalize the Data
%Use matlab's normalize function to scale the datasets. Scale the data using range such that the data is between 0 and 1. 
% Scale the data
Tabl_average_HCT116(:,2:end) = normalize(Tabl_average_HCT116(:, 2:end), 'range'); 
Tabl_std_HCT116(:,2:end) = normalize(Tabl_std_HCT116(:, 2:end), 'range');

Tabl_average_TU8902(:,2:end) = normalize(Tabl_average_TU8902(:, 2:end), 'range');
Tabl_std_TU8902(:, 2:end) = normalize(Tabl_std_TU8902(:, 2:end), 'range');

%% Separate the Data: GOT1 vs NT Data
% CELL LINE: HCT116
Tabl_avg_HCT116_GOT1 = Tabl_average_HCT116(:, {'metabolite', 'HCT116_shGOT1_no_dox_1','HCT116_shGOT1_2d_dox_1', 'HCT116_shGOT1_4d_dox_1', 'HCT116_shGOT1_no_dox_media_1', 'HCT116_shGOT1_2d_dox_media_1', 'HCT116_shGOT1_4d_dox_media_1'});
Tabl_avg_HCT116_NT = Tabl_average_HCT116(:, {'metabolite', 'HCT116_shNT_no_dox_1','HCT116_shNT_2d_dox_1', 'HCT116_shNT_4d_dox_1', 'HCT116_shNT_no_dox_media_1', 'HCT116_shNT_2d_dox_media_1', 'HCT116_shNT_4d_dox_media_1'});

Tabl_std_HCT116_GOT1 = Tabl_std_HCT116(:, {'metabolite', 'HCT116_shGOT1_no_dox_1','HCT116_shGOT1_2d_dox_1', 'HCT116_shGOT1_4d_dox_1', 'HCT116_shGOT1_no_dox_media_1', 'HCT116_shGOT1_2d_dox_media_1', 'HCT116_shGOT1_4d_dox_media_1'});
Tabl_std_HCT116_NT = Tabl_std_HCT116(:, {'metabolite', 'HCT116_shNT_no_dox_1','HCT116_shNT_2d_dox_1', 'HCT116_shNT_4d_dox_1', 'HCT116_shNT_no_dox_media_1', 'HCT116_shNT_2d_dox_media_1', 'HCT116_shNT_4d_dox_media_1'});

% CELL LINE: TU9802
Tabl_avg_TU8902_GOT1 = Tabl_average_TU8902(:, {'metabolite', 'TU8902_shGOT1_no_dox_1','TU8902_shGOT1_2d_dox_1', 'TU8902_shGOT1_4d_dox_1', 'TU8902_shGOT1_no_dox_media_1', 'TU8902_shGOT1_2d_dox_media_1', 'TU8902_shGOT1_4d_dox_media_1'});
Tabl_avg_TU8902_NT = Tabl_average_TU8902(:, {'metabolite', 'TU8902_shNT_no_dox_1','TU8902_shNT_2d_dox_1', 'TU8902_shNT_4d_dox_1', 'TU8902_shNT_no_dox_media_1', 'TU8902_shNT_2d_dox_media_1', 'TU8902_shNT_4d_dox_media_1'});

Tabl_std_TU8902_GOT1 = Tabl_std_TU8902(:, {'metabolite', 'TU8902_shGOT1_no_dox_1','TU8902_shGOT1_2d_dox_1', 'TU8902_shGOT1_4d_dox_1', 'TU8902_shGOT1_no_dox_media_1', 'TU8902_shGOT1_2d_dox_media_1', 'TU8902_shGOT1_4d_dox_media_1'});
Tabl_std_TU8902_NT = Tabl_std_TU8902(:, {'metabolite', 'TU8902_shNT_no_dox_1','TU8902_shNT_2d_dox_1', 'TU8902_shNT_4d_dox_1', 'TU8902_shNT_no_dox_media_1', 'TU8902_shNT_2d_dox_media_1', 'TU8902_shNT_4d_dox_media_1'});

%% Create Numeric Arrays of all Tables
% Cell Line: TU8902
%all conditions
mean_metabolomic_array_TU8902 = table2array(Tabl_average_TU8902(:,2:end));
std_metabolomic_array_TU8902 = table2array(Tabl_std_TU8902(:,2:end));

%GOT1 and NT
avg_TU8902_GOT1 = table2array(Tabl_avg_TU8902_GOT1(:, 2:end));
avg_TU8902_NT = table2array(Tabl_avg_TU8902_NT(:, 2:end));

std_TU8902_GOT1 = table2array(Tabl_std_TU8902_GOT1(:, 2:end));
std_TU8902_NT = table2array(Tabl_std_TU8902_NT(:, 2:end));

% Cell Line: HCT116
%all conditions
mean_metabolomic_array_HCT116 = table2array(Tabl_average_HCT116(:,2:end));
std_metabolomic_array_HCT116 = table2array(Tabl_std_HCT116(:,2:end));

%GOT1 and NT
avg_HCT116_GOT1 = table2array(Tabl_avg_HCT116_GOT1(:, 2:end));
avg_HCT116_NT = table2array(Tabl_avg_HCT116_NT(:, 2:end));

std_HCT116_GOT1 = table2array(Tabl_std_HCT116_GOT1(:, 2:end));
std_HCT116_NT = table2array(Tabl_std_HCT116_NT(:, 2:end));

%% Create the coefficient of variation, CV 
%CV matrix for all combinations

CV_TU8902 = std_metabolomic_array_TU8902 ./ mean_metabolomic_array_TU8902;
CV_TU8902_GOT1 = std_TU8902_GOT1 ./ avg_TU8902_GOT1;
CV_TU8902_NT = std_TU8902_NT ./ avg_TU8902_NT;

CV_HCT116 = std_metabolomic_array_HCT116 ./ mean_metabolomic_array_HCT116;
CV_HCT116_GOT1 = std_HCT116_GOT1 ./ avg_HCT116_GOT1;
CV_HCT116_NT = std_HCT116_NT ./ avg_HCT116_NT;

%% DFA Filter 1: CV thresholds + NaN: 0.1
[filtered_table_TU8902_01, metabolite_delete_TU8902_01,metabolomics_TU8902_01] = all_CV_Filter(CV_TU8902,metabolomics_TU8902,Tabl_average_TU8902,0.1, all_conditions);
[filtered_table_HCT116_01, metabolite_delete_HCT116_01,metabolomics_HCT116_01] = all_CV_Filter(CV_HCT116,metabolomics_HCT116,Tabl_average_HCT116,0.1, all_conditions);


bool = isnan(table2array(filtered_table_HCT116_01(:, 2:end)));
% If there are more than 6 empty values row-wise, then we'll remove the metabolite. Of course, if you find that there are still a lot NaNs, you can make this filter more stringent.
metabolite_deleteT = metabolomics_HCT116_01(sum(double(bool),2) > 5);
filtered_table_HCT116_01(sum(double(bool), 2) > 5, :) = [];
metabolomics_HCT116_01(sum(double(bool), 2) > 5, :) = [];

%filename = 'CVnanFilter0.6CancerHCT116.xlsx';
%writetable(filtered_table_HCT116,filename,'WriteRowNames', true, "FileType","spreadsheet");

bool = isnan(table2array(filtered_table_TU8902_01(:, 2:end)));

% If there are more than 6 empty values row-wise, then we'll remove the metabolite. Of course, if you find that there are still a lot NaNs, you can make this filter more stringent.
metabolite_deleteT = metabolomics_TU8902_01(sum(double(bool),2) > 5);
filtered_table_TU8902_01(sum(double(bool), 2) > 5, :) = [];
metabolomics_TU8902_01(sum(double(bool),2) > 5,:) = [];
%filename = 'CVnanFilter0.6CancerTU8902.xlsx';
%writetable(filtered_table_TU8902,filename,'WriteRowNames', true, "FileType","spreadsheet");

%% DFA Filter 1: CV thresholds + NaN: 0.3
[filtered_table_TU8902_03, metabolite_delete_TU8902_03,metabolomics_TU8902_03] = all_CV_Filter(CV_TU8902,metabolomics_TU8902,Tabl_average_TU8902,0.3, all_conditions);
[filtered_table_HCT116_03, metabolite_delete_HCT116_03,metabolomics_HCT116_03] = all_CV_Filter(CV_HCT116,metabolomics_HCT116,Tabl_average_HCT116,0.3, all_conditions); 

bool = isnan(table2array(filtered_table_HCT116_03(:, 2:end)));
% If there are more than 6 empty values row-wise, then we'll remove the metabolite. Of course, if you find that there are still a lot NaNs, you can make this filter more stringent.
metabolite_deleteT = metabolomics_HCT116_03(sum(double(bool),2) > 5);
filtered_table_HCT116_03(sum(double(bool), 2) > 5, :) = [];
metabolomics_HCT116_03(sum(double(bool), 2) > 5, :) = [];

%filename = 'CVnanFilter0.6CancerHCT116.xlsx';
%writetable(filtered_table_HCT116,filename,'WriteRowNames', true, "FileType","spreadsheet");

bool = isnan(table2array(filtered_table_TU8902_03(:, 2:end)));

% If there are more than 6 empty values row-wise, then we'll remove the metabolite. Of course, if you find that there are still a lot NaNs, you can make this filter more stringent.
metabolite_deleteT = metabolomics_TU8902_03(sum(double(bool),2) > 5);
filtered_table_TU8902_03(sum(double(bool), 2) > 5, :) = [];
metabolomics_TU8902_03(sum(double(bool),2) > 5,:) = [];
%filename = 'CVnanFilter0.6CancerTU8902.xlsx';
%writetable(filtered_table_TU8902,filename,'WriteRowNames', true, "FileType","spreadsheet");

%% DFA Filter 1: CV thresholds: 0.6
[filtered_table_TU8902_06, metabolite_delete_TU8902_06,metabolomics_TU8902_06] = all_CV_Filter(CV_TU8902,metabolomics_TU8902,Tabl_average_TU8902,0.6, all_conditions);
[filtered_table_HCT116_06, metabolite_delete_HCT116_06,metabolomics_HCT116_06] = all_CV_Filter(CV_HCT116,metabolomics_HCT116,Tabl_average_HCT116,0.6, all_conditions);


bool = isnan(table2array(filtered_table_HCT116_06(:, 2:end)));
% If there are more than 6 empty values row-wise, then we'll remove the metabolite. Of course, if you find that there are still a lot NaNs, you can make this filter more stringent.
metabolite_deleteT = metabolomics_HCT116_06(sum(double(bool),2) > 5);
filtered_table_HCT116_06(sum(double(bool), 2) > 5, :) = [];
metabolomics_HCT116_06(sum(double(bool), 2) > 5, :) = [];

%filename = 'CVnanFilter0.6CancerHCT116.xlsx';
%writetable(filtered_table_HCT116,filename,'WriteRowNames', true, "FileType","spreadsheet");

bool = isnan(table2array(filtered_table_TU8902_06(:, 2:end)));

% If there are more than 6 empty values row-wise, then we'll remove the metabolite. Of course, if you find that there are still a lot NaNs, you can make this filter more stringent.
metabolite_deleteT = metabolomics_TU8902_06(sum(double(bool),2) > 5);
filtered_table_TU8902_06(sum(double(bool), 2) > 5, :) = [];
metabolomics_TU8902_06(sum(double(bool),2) > 5,:) = [];
%filename = 'CVnanFilter0.6CancerTU8902.xlsx';
%writetable(filtered_table_TU8902,filename,'WriteRowNames', true, "FileType","spreadsheet");

%filename = 'cvFilter0.1CancerHCT116.xlsx';
%writetable(filtered_table_HCT116,filename,'WriteRowNames', true, "FileType","spreadsheet");

%filename = 'cvFilter0.1CancerTU8902.xlsx';
%writetable(filtered_table_TU8902,filename2, 'WriteRowNames', true ,"FileType","spreadsheet");

%% DFA Filter 2: NaN Filter
Now, we will use the Filter 2 to create an excel sheet for DFA. 
We will first locate all the locations that have only NaN's in tables and then delete the NaN rows. The HCT116 cell line should delete 7 rows while the TU8902 cell line should delete 6 rows. 
NaN_Filtered_table_HCT116 = Tabl_average_HCT116;
% Construct an N x M boolean matrix that finds NaNs for each element
bool = isnan(table2array(Tabl_average_HCT116(:, 2:end)));

% If there are more than 6 empty values row-wise, then we'll remove the metabolite. Of course, if you find that there are still a lot NaNs, you can make this filter more stringent.
metabolite_deleteH = metabolomics_HCT116(sum(double(bool),2) >= 6);
NaN_Filtered_table_HCT116(sum(double(bool), 2) >= 6, :) = [];
%filename = 'nanFilterCancerHCT116.xlsx';
%writetable(Tabl_average_HCT116,filename,'WriteRowNames', true, "FileType","spreadsheet");
Repeat for the TU9802 Cell line
NaN_Filtered_table_TU8902 = Tabl_average_TU8902;
% Construct an N x M boolean matrix that finds NaNs for each element
bool = isnan(table2array(Tabl_average_TU8902(:, 2:end)));

% If there are more than 6 empty values row-wise, then we'll remove the metabolite. Of course, if you find that there are still a lot NaNs, you can make this filter more stringent.
metabolite_deleteT = metabolomics_TU8902(sum(double(bool),2) >= 6);
NaN_Filtered_table_TU8902(sum(double(bool), 2) >= 6, :) = [];
%filename = 'nanFilterCancerTU8902.xlsx';
%writetable(Tabl_average_TU8902,filename,'WriteRowNames', true, "FileType","spreadsheet");
