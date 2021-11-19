function Nu = Nu_Nikitin(Re,Pr)  
    %cold side channels from Nikitin
    Re = 1.305407289332279*Re;
    Re(Re<1000) = 1000;
%     Re(Re>12100) = 12100;
    Re(Re>100000) = 100000;
%    Nu = C_Nu*5.49*Re.^0.625;
    im = 2;%C_Nu; 
    %Dhdx = [0.0548,0.0616,0.0683]; % Dh/k [m2-K/W]
    Dhdx = [0.0335,0.0376,0.0417]; % Dh/k [m2-K/W]
    Nu = (Dhdx(im)*5.49*Re.^0.625); % 1
%     a = 0.848361835135;
%     Nu = a*Nu;
    
%     a = 1.54646642583;
%     b = -0.0567293119117;
%     Nu = a*Nu.*Re.^b;
end