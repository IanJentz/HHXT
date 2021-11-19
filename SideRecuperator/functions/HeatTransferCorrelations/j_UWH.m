function j = j_UWH(Re,Pr)  
    
    alpha = 80;
    cosa = 1/cos(pi*alpha/360);
%     Re = cosa*Re;
    
%     % mean value, not H or C
%     j_mean = 0.009528045454545;
%     j = 2.691883641449576*Re.^(-0.663806836429960);
%     j(j>j_mean) = j_mean;
    
    a = 2.331077636630432;
    b = -0.633047975298456;

    a = 5.313298567093888;
    b = -0.718546630088044;

    % cold value
    j_mean = 0.009528045454545;
    j = a*Re.^(b);
    j(j>j_mean) = j_mean;

   
    
end