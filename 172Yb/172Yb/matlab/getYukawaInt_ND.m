function Y = getYukawaInt_ND(r,ND,m)
% % normalization test
% norm_test = dw*rho_int.';
% if abs(norm_test - 1) > 0.01
%     error('Wrong normalization: %g',norm_test)
% end
% rho_int = rho_int/norm_test;
r = r(2:end);
ND = ND(2:end);

[M,NNDD] = meshgrid(m,ND);
[M,R] = meshgrid(m,r);
% f_Yukawa = @(r,m) exp(-m.*r)./r;
V_Yukawa = exp(-M.*R)./R;
% V_Yukawa = exp(-M.*r)./r;
% f_Y = @(m) trapz(r,f_Yukawa(r,m).*ND);
% Y = arrayfun(f_Y,m);
Y = trapz(r,V_Yukawa.*NNDD,1);

if iscolumn(m)
    Y = Y.';
end
end