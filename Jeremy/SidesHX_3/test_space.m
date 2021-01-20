space = [];
NTU = logspace(0,1.8,11);
% alpha = linspace(0.1,0.5,5); 
alpha = 0.1;
% alpha(1) = 0.0001;
% alpha(end) = 0.9999;
% C_r = linspace(0,1,11);
% C_r(1) = 0.01;
C_r = [0.9 1];
for i = 1:11
    for j = 1:1
        for k = 1:2
            %if ismember(i,[1,5,11]) && ismember(j,[1,5,11]) && ismember(k,[1,5,11])
             %   space = [space;[1 alpha(j) NTU(i) C_r(k)]];
            %else
               space = [space;[0  alpha(j) NTU(i) C_r(k)]];
            %end
    
        end
    end
end
space = [linspace(1,length(space),length(space))' space];
csvwrite('input_data_quick.txt',space)