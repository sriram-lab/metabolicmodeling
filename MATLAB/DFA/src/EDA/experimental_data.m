 exper_data = readtable('Barb Nelson Pos Table-20180410 - HL.xlsx', ...
    'Sheet', 'TU8902');
%delete row 2 and column 2
exper_data.transition = [];

Tabl_aver = exper_data;
Tabl_std = exper_data;

for n = 2:13;
    Tabl_aver(:,[n]) = num2cell(mean(exper_data{:,n:n+2},2));
    Tabl_std(:,[n]) = num2cell(std(exper_data{:,n:n+2},0,2));
    Tabl_aver(:,[n+2]) = [];
    Tabl_std(:,[n+2]) = [];
    Tabl_aver(:,[n+1]) = [];
    Tabl_std(:,[n+1]) = [];
end 

%heatmaps
