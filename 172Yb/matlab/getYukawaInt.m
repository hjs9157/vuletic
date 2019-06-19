function YukawaInt = getYukawaInt(rwfn,gamma,m)
r = rwfn.r;
P = rwfn.P;
Q = rwfn.Q;

rho = P.^2 + Q.^2;
rhor = rho./r.^(2*gamma);
rhor(1) = rhor(2);

% normalization test
dw = (r(3:end) - r(1:end-2))/2;
norm_test = dw.'*rho(2:end-1);
if abs(norm_test - 1) > 0.01;
    error('Wrong normalization: %g',norm_test)
end

% integration for each m
N_int_min = 10000;
YukawaInt = nan*m;
for minx = 1:length(m)
    intrange = 10/m(minx);
    %     if sum(r < intrange) >= N_int_min || sum(r < intrange) == length(r)
    if intrange >= r(end)
        r_int = r;
        rhor_int = rhor;
    else
        r_int = linspace(0,intrange,N_int_min+2).';
        rhor_int = interp1(r,rhor,r_int);
    end
    
    dw = (r_int(3:end) - r_int(1:end-2))/2;
    rho_int = rhor_int.*r_int.^(2*gamma);
    rhodw = dw.*rho_int(2:end-1);
    
    V_Yukawa = exp(-m(minx)*r_int.')./r_int.';
    
    YukawaInt(minx) = V_Yukawa(:,2:end-1)*rhodw;
end

end