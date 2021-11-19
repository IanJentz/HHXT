function Nu = func_Nu_Katz_CO2(pos,time,Re,Pr)
    cosa = 1.252135658156226;
    Re = cosa*Re;
    Re(Re<2000) = 2000;
%     Re(Re>12000) = 12000;
    Nu = (0.02609*Re.^(0.8765).*Pr.^(1/3));
%     Nu = cosa*(0.02609*Re.^(0.8765).*Pr.^(1/3));
    a = 0.832438990976;
    Nu = a*Nu;
end