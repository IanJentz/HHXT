function f = f_Katz_He(Re)
    cosa = 1.252135658156226;   
    Re = cosa*Re;
    Re(Re<390) = 390;
%     Re(Re>3130) = 3130;
%     f = cosa*4*(0.508*Re.^(-0.276)); %1/cos
%     f = cosa^2*4*(0.508*Re.^(-0.276)); %(1/cos)^2
    f = cosa^3*4*(0.508*Re.^(-0.276)); %(1/cos)^3
end