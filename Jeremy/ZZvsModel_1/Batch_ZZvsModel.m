






f_D = 0.4507; % we will keep friction factor constant
j_C_lo = 0.005; j_C_hi = .5; % range of Colburn heat transfer coefficients
N_j = 5; % number of heat transfer coeffs to run
j_C = logspace(log10(j_C_lo),log10(j_C_hi),N_j); % distribute j_C logarithmically

% preallocate the result variables as column arrays
zcol = zeros(N_j,1);
a_rmsd_fiber = zcol;
rmsd_TC = zcol;
rmsd_DP = zcol;
rmsd_Tout = zcol;
rmsd_qdot = zcol;
Re_C = zcol;
Re_H = zcol;

% run the many instances of ZZvsModel on a distributed system

for i = 1:N_j
    
[rmsd_fiber(i),rmsd_TC(i),rmsd_DP(i),rmsd_Tout(i),rmsd_qdot(i),Re_C(i),Re_H(i)] = ZZvsModel(j_C(i),f_D);

end