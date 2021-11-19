function Nu = func_Nu_Katz_He(pos,time,Re,Pr)
    cosa = 1.252135658156226;
    Re = cosa*Re;
    Re(Re<500) = 500;
%     Re(Re>2500) = 2500;
    Nu = (0.8072*Re.^(0.4359).*Pr.^(1/3));
%     Nu = cosa*(0.8072*Re.^(0.4359).*Pr.^(1/3));
end