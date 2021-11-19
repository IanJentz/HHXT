function fmatr = func_fmatr_h_Nikitin(Re)
    f = f_Nikitin( Re(1,:) );
    fmatr = [f;f*1000000];
end