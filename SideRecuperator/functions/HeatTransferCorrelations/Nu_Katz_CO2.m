function Nu = Nu_Katz_CO2(Re,Pr)
    cosa = 1.252135658156226;
    Re = cosa*Re;
    
%     Re(Re<370) = 370;
%     Re(Re>17056) = 17056;
    Re(Re<1000) = 1000;
    Re(Re>100000) = 100000;

    Nu = (0.02609*Re.^(0.8765).*Pr.^(1/3));
%     Nu = cosa*(0.02609*Re.^(0.8765).*Pr.^(1/3));
%     a = 0.832438990976;
%     Nu = a*Nu;
    
    a = 1.60290485103;
    b = -0.0637885211338;
    Nu = a*Nu.*Re.^b;
    
end