function YukawaInt = getYukawaCorrInt_simple(rwfn1,rwfn2,m)
r1 = rwfn1.r.'; P1 = rwfn1.P.'; Q1 = rwfn1.Q.';
r2 = rwfn2.r.'; P2 = rwfn2.P.'; Q2 = rwfn2.Q.';

rho1 = P1.^2 + Q1.^2;
rho2 = P2.^2 + Q2.^2;

dw1 = (r1(3:end) - r1(1:end-2))/2;
dw2 = (r2(3:end) - r2(1:end-2))/2;

% normalization test
norm_test1 = dw1*rho1(2:end-1).';
norm_test2 = dw2*rho2(2:end-1).';
if abs(norm_test1 - 1) > 0.01, error('Wrong normalization 1: %g',norm_test1); end
if abs(norm_test2 - 1) > 0.01, error('Wrong normalization 2: %g',norm_test2); end
P1 = P1/sqrt(norm_test1); Q1 = Q1/sqrt(norm_test1); 
P2 = P2/sqrt(norm_test1); Q2 = Q2/sqrt(norm_test1); 

N = min(length(r1),length(r2));
r = r1(1:N);
dw = (r(3:end) - r(1:end-2))/2;
r_int = r(2:end-1);
rcorr = P1(2:N-1).*P2(2:N-1) + Q1(2:N-1).*Q2(2:N-1);



[M,R] = meshgrid(m,r_int);
% f_Yukawa = @(r,m) exp(-m.*r)./r;
V_Yukawa = exp(-M.*R)./R;
YukawaInt = (dw.*rcorr)*V_Yukawa;

if iscolumn(m)
    YukawaInt = YukawaInt.';
end
end