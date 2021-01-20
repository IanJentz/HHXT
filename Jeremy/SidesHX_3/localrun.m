results = readtable('input_data_quick.txt');
results = table2array(results);
for i=1:length(results)
    SideHX_CHTC(results(i,1),results(i,2),results(i,3),results(i,4),results(i,5))
end