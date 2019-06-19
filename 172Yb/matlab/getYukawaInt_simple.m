function YukawaInt = getYukawaInt_simple(rwfn,m)
r = rwfn.r.';
P = rwfn.P.';
Q = rwfn.Q.';

rho = P.^2 + Q.^2;

dw = (r(3:end) - r(1:end-2))/2;
r_int = r(2:end-1);
rho_int = rho(2:end-1);

% normalization test
norm_test = dw*rho_int.';
if abs(norm_test - 1) > 0.01
    error('Wrong normalization: %g',norm_test)
end
rho_int = rho_int/norm_test;

[M,R] = meshgrid(m,r_int);
% f_Yukawa = @(r,m) exp(-m.*r)./r;
V_Yukawa = exp(-M.*R)./R;
YukawaInt = (dw.*rho_int)*V_Yukawa;
end