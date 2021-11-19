function Nu = func_Nu_DittusBoelterC(pos,time,Re,Pr)
    % Dittus-Boelter relation for a heated fluid
    Nu = 0.023*Re.^0.8.*Pr.^0.4;
end