function fmatr = func_fmatr_h_Katz_CO2(Re)
    f = f_Katz_CO2( Re(1,:) );
    fmatr = [f;f*1000000];
end