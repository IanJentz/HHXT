function j = j_UWC(Re,Pr)  
    
    alpha = 80;
    cosa = 1/cos(pi*alpha/360);
%     Re = cosa*Re;
    
%     % mean value, not H or C
%     j_mean = 0.009528045454545;
%     j = 2.691883641449576*Re.^(-0.663806836429960);
%     j(j>j_mean) = j_mean;
    
    a = 3.945987730432754;
    b = -0.708312671462686;
    
    % cold value
    j_mean = 0.009528045454545;
    j = a*Re.^(b);
    j(j>j_mean) = j_mean;
   
    
end