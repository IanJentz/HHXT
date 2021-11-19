function Nu = func_Nu_Nikitin(pos,time,Re,Pr)
    cosa = 1.305407289332279;
    Re = cosa*Re;
    Re(Re<1000) = 1000;
%     Re(Re>12100) = 12100;
%     Re(Re>60000) = 60000;
%     Nu = Nu_UW(Re).*Pr.^0.4;
    Nu = Nu_Nikitin(Re,Pr);
%     Nu = 1.305407289332279*Nu_Nikitin(Re,Pr);
    a = 0.784647712172;
    Nu = a*Nu;
end

