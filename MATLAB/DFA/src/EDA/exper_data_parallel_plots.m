%{ 
Experimental Data - Parallel Coordinate Plots
Author: Surabhi Joshi

Table of Contents
Experimental Data - Parallel Coordinate Plots
Load Exper_data.mlx 
Parallel Coordinate Plots - Mean, Unfiltered
Parallel Coordinate Plots - CV, Unfiltered
Parallel Coordinate Plots - CV, Filtered CV = 0.1 + NaN
Parallel Coordinate Plots - CV, Filtered CV = 0.3 + NaN
Parallel Coordinate Plots - CV, Filtered CV = 0.6 + NaN
%}

%% Load Exper_data.mlx 
%Start by loading exper_data.mlx to have all of the data tables for the parallel coordinate plots
exper_data;

%% Parallel Coordinate Plots - Mean, Unfiltered
% TU8902 
figure (1);
aminoAcidIndex = findClusterType(metabolomics_TU8902,'amino acid');

subplot(3,1,1);
parallelcoords(mean_metabolomic_array_TU8902(aminoAcidIndex,:), 'group', metabolomics_TU8902(aminoAcidIndex,1), 'Labels', all_conditions);
xtickangle(45);
title('Mean Values - TU9802 Amino Acids')
legend('location', 'bestoutside');
ylabel('Averages')

subplot(3,1,2);
parallelcoords(avg_TU8902_GOT1(aminoAcidIndex,:), 'group', metabolomics_TU8902(aminoAcidIndex,1), 'Labels', got1_conditions);
xtickangle(45);
title('Mean Values - TU9802 Amino Acids GOT1')
legend('location', 'bestoutside');
ylabel('Averages')

subplot(3,1,3);
parallelcoords(avg_TU8902_NT(aminoAcidIndex,:), 'group', metabolomics_TU8902(aminoAcidIndex,1), 'Labels', nt_conditions);
xtickangle(45);
title('Mean Values - TU9802 Amino Acids NT')
legend('location', 'bestoutside');
ylabel('Averages')

figure (2);
hold on
subplot(3,1,1);
parallelcoords(mean_metabolomic_array_TU8902, 'group', metabolomics_TU8902(:,2), 'labels', all_conditions);
xtickangle(45);
title('Mean of the metabolomics: TU8902');
legend('location', 'bestoutside');
ylabel('Averages')

subplot(3,1,2);
parallelcoords(avg_TU8902_GOT1, 'group', metabolomics_TU8902(:,2), 'labels', got1_conditions);
xtickangle(45);
title('Mean of the metabolomics: TU8902 GOT1');
legend('location', 'bestoutside');
ylabel('Averages')

subplot(3,1,3);
parallelcoords(avg_TU8902_NT, 'group', metabolomics_TU8902(:,2), 'labels', nt_conditions);
xtickangle(45);
title('Mean of the metabolomics: TU8902 NT');
legend('location', 'bestoutside');
ylabel('Averages')
hold off


%HCT116 
figure (3);
aminoAcidIndex = findClusterType(metabolomics_HCT116,'amino acid');

subplot(3,1,1);
parallelcoords(mean_metabolomic_array_HCT116(aminoAcidIndex,:), 'group', metabolomics_HCT116(aminoAcidIndex,1), 'Labels', all_conditions);
xtickangle(45);
title('Mean Values - HCT116 Amino Acids')
legend('location', 'bestoutside');
ylabel('Averages')

subplot(3,1,2);
parallelcoords(avg_HCT116_GOT1(aminoAcidIndex,:), 'group', metabolomics_HCT116(aminoAcidIndex,1), 'Labels', got1_conditions);
xtickangle(45);
title('Mean Values - HCT116 Amino Acids GOT1')
legend('location', 'bestoutside');
ylabel('Averages')

subplot(3,1,3);
parallelcoords(avg_HCT116_NT(aminoAcidIndex,:), 'group', metabolomics_HCT116(aminoAcidIndex,1), 'Labels', nt_conditions);
xtickangle(45);
title('Mean Values - HCT116 Amino Acids NT')
legend('location', 'bestoutside');
ylabel('Averages')

figure (4)
subplot(3,1,1);
parallelcoords(mean_metabolomic_array_HCT116, 'group', metabolomics_TU8902(:,2), 'Labels', all_conditions);
xtickangle(45);
title('Mean metabolomics Data: HCT116');
legend('location', 'bestoutside');
ylabel('Averages')

subplot(3,1,2);
parallelcoords(avg_HCT116_GOT1, 'group', metabolomics_TU8902(:,2), 'Labels', got1_conditions);
xtickangle(45);
title('Mean metabolomics Data: HCT116 GOT1');
legend('location', 'bestoutside');
ylabel('Averages')

subplot(3,1,3);
parallelcoords(avg_HCT116_NT, 'group', metabolomics_TU8902(:,2), 'Labels', nt_conditions);
xtickangle(45);
title('Mean metabolomics Data: HCT116 NT');
legend('location', 'bestoutside');
ylabel('Averages')
hold off

%% Parallel Coordinate Plots - CV, Unfiltered
figure (5);
subplot(2,1,1);
parallelcoords(CV_TU8902, 'group', metabolomics_TU8902(:,2), 'Labels', all_conditions);
xtickangle(45);
title('Coefficient of variation: TU8902');
legend('location', 'bestoutside');
ylabel('Coefficient Of Variation')

subplot(2,1,2);
parallelcoords(CV_HCT116, 'group', metabolomics_TU8902(:,2), 'Labels', all_conditions);
xtickangle(45);
title('Coefficient of Variation: HCT116');
legend('location', 'bestoutside');
ylabel('Coefficient Of Variation')

figure (6);
subplot(2,2,1);
parallelcoords(CV_TU8902_GOT1, 'group', metabolomics_TU8902(:,2), 'Labels', got1_conditions);
xtickangle(45);
title('Coefficient of variation: TU8902 GOT1');
legend('location', 'bestoutside');
ylabel('Coefficient Of Variation')

subplot(2,2,2);
parallelcoords(CV_TU8902_NT, 'group', metabolomics_TU8902(:,2), 'Labels', nt_conditions);
xtickangle(45);
title('Coefficient of variation: TU8902 NT');
legend('location', 'bestoutside');
ylabel('Coefficient Of Variation')

subplot(2,2,3);
parallelcoords(CV_HCT116_GOT1, 'group', metabolomics_TU8902(:,2), 'Labels', got1_conditions);
xtickangle(45);
title('Coefficient of Variation: HCT116 GOT1');
legend('location', 'bestoutside');
ylabel('Coefficient Of Variation')

subplot(2,2,4);
parallelcoords(CV_HCT116_NT, 'group', metabolomics_TU8902(:,2), 'Labels', nt_conditions);
xtickangle(45);
title('Coefficient of Variation: HCT116 NT');
legend('location', 'bestoutside');
ylabel('Coefficient Of Variation')
Now we will look specifically at the amino acids.
figure (7);
hold on
subplot(2,1,1);
parallelcoords(CV_TU8902(aminoAcidIndex,:), 'group', metabolomics_TU8902(aminoAcidIndex,1), 'Labels', all_conditions);
xtickangle(45);
title('Coefficient of Variation Values - TU9802 Amino Acids')
legend('location', 'bestoutside');
ylabel('Coefficient Of Variation')

subplot(2,1,2);
parallelcoords(CV_HCT116(aminoAcidIndex,:), 'group', metabolomics_TU8902(aminoAcidIndex,1), 'Labels', all_conditions);
xtickangle(45);
title('Coefficient of Variation Values - HCT116 Amino Acids')
legend('location', 'bestoutside');
ylabel('Coefficient Of Variation')
hold off
% Now we will look at nucleosides
nucleosideIndex = findClusterType(metabolomics_TU8902, 'nucleoside');

figure (8);
hold on
subplot(2,1,1);
parallelcoords(CV_TU8902(nucleosideIndex,:), 'group', metabolomics_TU8902(nucleosideIndex,1), 'Labels', all_conditions);
xtickangle(45);
title('Coefficient of Variation Values - TU9802 nucleoside')
legend('location', 'bestoutside');
ylabel('Coefficient Of Variation')

subplot(2,1,2);
parallelcoords(CV_HCT116(nucleosideIndex,:), 'group', metabolomics_TU8902(nucleosideIndex,1), 'Labels', all_conditions);
xtickangle(45);
title('Coefficient of Variation Values - HCT116 nucleoside')
legend('location', 'bestoutside');
ylabel('Coefficient Of Variation')
hold off
% Now we will look at nucleotides
nucleotideIndex = findClusterType(metabolomics_TU8902, 'nucleotide');

figure (9);
hold on
subplot(2,1,1);
parallelcoords(CV_TU8902(nucleotideIndex,:), 'group', metabolomics_TU8902(nucleotideIndex,1), 'Labels', all_conditions);
xtickangle(45);
title('Coefficient of Variation Values - TU9802 nucleotide')
legend('location', 'bestoutside');
ylabel('Coefficient Of Variation')

subplot(2,1,2);
parallelcoords(CV_HCT116(nucleotideIndex,:), 'group', metabolomics_TU8902(nucleotideIndex,1), 'Labels', all_conditions);
xtickangle(45);
title('Coefficient of Variation Values - HCT116 nucleotide')
legend('location', 'bestoutside');
ylabel('Coefficient Of Variation')
hold off

%% Parallel Coordinate Plots - CV, Filtered CV = 0.1 + NaN
%TU8902
AA_Index = findClusterType(metabolomics_TU8902_01, 'amino acid');
glycolysis_Index = findClusterType(metabolomics_TU8902_01, 'glycolysis');
ppp_Index = findClusterType(metabolomics_TU8902_01, 'PPP');

figure (10);
subplot(3,1,1);
parallelcoords(table2array(filtered_table_TU8902_01(AA_Index,:)), 'group', metabolomics_TU8902_01(AA_Index), 'labels', all_conditions);
title('CV Values - TU8902 Amino Acids');
ylabel('CV');
legend('Location',"bestoutside");
xtickangle(45);

subplot(3,1,2);
parallelcoords(table2array(filtered_table_TU8902_01(glycolysis_Index,:)), 'group', metabolomics_TU8902_01(glycolysis_Index), 'labels', all_conditions);
title('CV Values - TU8902 Glycolysis');
ylabel('CV');
legend('Location',"bestoutside");
xtickangle(45);

subplot(3,1,3);
parallelcoords(table2array(filtered_table_TU8902_01(ppp_Index,:)), 'group', metabolomics_TU8902_01(ppp_Index), 'labels', all_conditions);
title('CV Values - TU8902 PPP');
ylabel('CV');
legend('Location',"bestoutside");
xtickangle(45);



%HCT116
AA_Index = findClusterType(metabolomics_HCT116_01, 'amino acid');
glycolysis_Index = findClusterType(metabolomics_HCT116_01, 'glycolysis');
ppp_Index = findClusterType(metabolomics_HCT116_01, 'PPP');

figure (11);
subplot(3,1,1);
parallelcoords(table2array(filtered_table_HCT116_01(AA_Index,:)), 'group', metabolomics_HCT116_01(AA_Index), 'labels', all_conditions);
title('CV Values - HCT116 Amino Acids');
ylabel('CV');
legend('Location',"bestoutside");
xtickangle(45);

subplot(3,1,2);
parallelcoords(table2array(filtered_table_HCT116_01(glycolysis_Index,:)), 'group', metabolomics_HCT116_01(glycolysis_Index), 'labels', all_conditions);
title('CV Values - HCT116 Glycolysis');
ylabel('CV');
legend('Location',"bestoutside");
xtickangle(45);

subplot(3,1,3);
parallelcoords(table2array(filtered_table_HCT116_01(ppp_Index,:)), 'group', metabolomics_HCT116_01(ppp_Index), 'labels', all_conditions);
title('CV Values - HCT116 PPP');
ylabel('CV');
legend('Location',"bestoutside");
xtickangle(45);

%% Parallel Coordinate Plots - CV, Filtered CV = 0.3 + NaN
AA_Index = findClusterType(metabolomics_TU8902_03, 'amino acid');
glycolysis_Index = findClusterType(metabolomics_TU8902_03, 'glycolysis');
ppp_Index = findClusterType(metabolomics_TU8902_03, 'PPP');

figure (14);
subplot(3,1,1);
parallelcoords(table2array(filtered_table_TU8902_03(AA_Index,:)), 'group', metabolomics_TU8902_03(AA_Index), 'labels', all_conditions);
title('CV Values - TU8902 Amino Acids');
ylabel('CV');
legend('Location',"bestoutside");
xtickangle(45);

subplot(3,1,2);
parallelcoords(table2array(filtered_table_TU8902_03(glycolysis_Index,:)), 'group', metabolomics_TU8902_03(glycolysis_Index), 'labels', all_conditions);
title('CV Values - TU8902 Glycolysis');
ylabel('CV');
legend('Location',"bestoutside");
xtickangle(45);

subplot(3,1,3);
parallelcoords(table2array(filtered_table_TU8902_03(ppp_Index,:)), 'group', metabolomics_TU8902_03(ppp_Index), 'labels', all_conditions);
title('CV Values - TU8902 PPP');
ylabel('CV');
legend('Location',"bestoutside");
xtickangle(45);



%HCT116
AA_Index = findClusterType(metabolomics_HCT116_03, 'amino acid');
glycolysis_Index = findClusterType(metabolomics_HCT116_03, 'glycolysis');
ppp_Index = findClusterType(metabolomics_HCT116_03, 'PPP');

figure (15);
subplot(3,1,1);
parallelcoords(table2array(filtered_table_HCT116_03(AA_Index,:)), 'group', metabolomics_HCT116_03(AA_Index), 'labels', all_conditions);
title('CV Values - HCT116 Amino Acids');
ylabel('CV');
legend('Location',"bestoutside");
xtickangle(45);

subplot(3,1,2);
parallelcoords(table2array(filtered_table_HCT116_03(glycolysis_Index,:)), 'group', metabolomics_HCT116_03(glycolysis_Index), 'labels', all_conditions);
title('CV Values - HCT116 Glycolysis');
ylabel('CV');
legend('Location',"bestoutside");
xtickangle(45);

subplot(3,1,3);
parallelcoords(table2array(filtered_table_HCT116_03(ppp_Index,:)), 'group', metabolomics_HCT116_03(ppp_Index), 'labels', all_conditions);
title('CV Values - HCT116 PPP');
ylabel('CV');
legend('Location',"bestoutside");
xtickangle(45);

%% Parallel Coordinate Plots - CV, Filtered CV = 0.6 + NaN
AA_Index = findClusterType(metabolomics_TU8902_06, 'amino acid');
glycolysis_Index = findClusterType(metabolomics_TU8902_06, 'glycolysis');
ppp_Index = findClusterType(metabolomics_TU8902_06, 'PPP');

figure (18);
subplot(3,1,1);
parallelcoords(table2array(filtered_table_TU8902_06(AA_Index,:)), 'group', metabolomics_TU8902_06(AA_Index), 'labels', all_conditions);
title('CV Values - TU8902 Amino Acids');
ylabel('CV');
legend('Location',"bestoutside");
xtickangle(45);

subplot(3,1,2);
parallelcoords(table2array(filtered_table_TU8902_06(glycolysis_Index,:)), 'group', metabolomics_TU8902_06(glycolysis_Index), 'labels', all_conditions);
title('CV Values - TU8902 Glycolysis');
ylabel('CV');
legend('Location',"bestoutside");
xtickangle(45);

subplot(3,1,3);
parallelcoords(table2array(filtered_table_TU8902_06(ppp_Index,:)), 'group', metabolomics_TU8902_06(ppp_Index), 'labels', all_conditions);
title('CV Values - TU8902 PPP');
ylabel('CV');
legend('Location',"bestoutside");
xtickangle(45);



%HCT116
AA_Index = findClusterType(metabolomics_HCT116_06, 'amino acid');
glycolysis_Index = findClusterType(metabolomics_HCT116_06, 'glycolysis');
ppp_Index = findClusterType(metabolomics_HCT116_06, 'PPP');

figure (19);
subplot(3,1,1);
parallelcoords(table2array(filtered_table_HCT116_06(AA_Index,:)), 'group', metabolomics_HCT116_06(AA_Index), 'labels', all_conditions);
title('CV Values - HCT116 Amino Acids');
ylabel('CV');
legend('Location',"bestoutside");
xtickangle(45);

subplot(3,1,2);
parallelcoords(table2array(filtered_table_HCT116_06(glycolysis_Index,:)), 'group', metabolomics_HCT116_06(glycolysis_Index), 'labels', all_conditions);
title('CV Values - HCT116 Glycolysis');
ylabel('CV');
legend('Location',"bestoutside");
xtickangle(45);

subplot(3,1,3);
parallelcoords(table2array(filtered_table_HCT116_06(ppp_Index,:)), 'group', metabolomics_HCT116_06(ppp_Index), 'labels', all_conditions);
title('CV Values - HCT116 PPP');
ylabel('CV');
legend('Location',"bestoutside");
xtickangle(45);