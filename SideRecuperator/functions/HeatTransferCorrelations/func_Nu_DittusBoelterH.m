function Nu = func_Nu_DittusBoelterH(pos,time,Re,Pr)
    % Dittus-Boelter relation for a cooled fluid
    Nu = 0.023*Re.^0.8.*Pr.^0.3;
end

