






% f_D = 0.4507; % we will keep friction factor constant
f_D_lo = 0.01; f_D_hi = 10; % range of Darcy Friction factor
j_C_lo = 0.0001; j_C_hi = .1; % range of Colburn heat transfer coefficients
N_f = 11; % number of friction factors to run
N_j = 11; % number of heat transfer coeffs to run
N = N_f*N_j*N_j; %total number of points
f_D = logspace(log10(f_D_lo),log10(f_D_hi),N_f); % distribute f_D logarithmically
j_C_C = logspace(log10(j_C_lo),log10(j_C_hi),N_j); % distribute j_C_C logarithmically
j_C_H = logspace(log10(j_C_lo),log10(j_C_hi),N_j); % distribute j_C_H logarithmically

[j_C_C_i,j_C_H_i,f_D_i] = meshgrid(j_C_C,j_C_H,f_D);
f_D_i = f_D_i(:);
j_C_C_i = j_C_C_i(:);
j_C_H_i = j_C_H_i(:);

% preallocate the result variables as column arrays
zcol = zeros(N,1);
rmsd_fiber = zcol;
rmsd_TC = zcol;
rmsd_DP = zcol;
rmsd_Tout = zcol;
rmsd_qdot = zcol;
Re_C = zcol;
Re_H = zcol;
DP_C = zcol;
DP_H = zcol;
T_C_out = zcol;
T_H_out = zcol;
qdot_C = zcol;
qdot_H = zcol;
Re_C_std = zcol;
Re_H_std = zcol;

% run the many instances of ZZvsModel on a distributed system
i = 690;
for i = 1:N
        
[rmsd_fiber(i),rmsd_TC(i),rmsd_DP(i),rmsd_Tout(i),rmsd_qdot(i),Re_C(i),Re_H(i),...
    DP_C(i),DP_H(i),T_C_out(i),T_H_out(i),qdot_C(i),qdot_H(i),...
    Re_C_std(i),Re_H_std(i)] = ...
    ZZvsModel([j_C_C_i(i),j_C_H_i(i)],f_D_i(i));

end