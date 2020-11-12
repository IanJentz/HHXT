space = [];
NTU = logspace(0,2,11);
alpha = linspace(0.0001,.9999,11);
C_r = linspace(0,1,11);

for i = 1:11
    for j = 1:11
        for k = 1:11
    space = [space;NTU(i) alpha(j) C_r(k)];
        end
    end
end
space = [linspace(1,length(space),length(space))' space];
csvwrite('input_data.txt',space)