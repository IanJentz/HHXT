function val = interpinbounds(T,P,y,temp,press)
    % provides 2D interpolation of temperature and pressure
    % with values falling outside of the T and P data being taken as the
    % value at the closest boundary (always returns a valid value)
    val = interp2(T,P,y,temp,press,'linear');
    
    bind = 1*(temp<T(1));
    bind = bind + 2*(temp>T(end));
    bind = bind + 3*(press<P(1));
    bind = bind + 6*(press>P(end));
    
    val(bind==4) = y(end,1);
    val(bind==5) = y(end,end);
    val(bind==7) = y(1,1);
    val(bind==8) = y(1,end);
    val(bind==1) = interp1(P(:,1),y(:,1),press(bind==1));
    val(bind==2) = interp1(P(:,end),y(:,end),press(bind==2));
    val(bind==3) = interp1(T(end,:),y(end,:),temp(bind==3));
    val(bind==6) = interp1(T(1,:),y(1,:),temp(bind==6));
end