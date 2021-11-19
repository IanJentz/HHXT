function fmatr = func_fmatr_h_Katz_He(Re)
    f = f_Katz_He( Re(1,:) );
    fmatr = [f;f*1000000];
end