stateID = '1S0';
DHForCI = 'DHF';
T_orbital = getOrbitals(stateID,DHForCI);

plotOrbital = {'6s','5p','4d','4f-','4f'};

orbital_fig = figure; hold on;
legendstr = {};
clear h;
for oinx = 1:length(plotOrbital)
    r = T_orbital{plotOrbital{oinx},'rwfn'}{1}.r;
    P = T_orbital{plotOrbital{oinx},'rwfn'}{1}.P;
    Q = T_orbital{plotOrbital{oinx},'rwfn'}{1}.Q;
    gamma = T_orbital{plotOrbital{oinx},'gamma'};
    
    rho = P.^2 + Q.^2;
%     rho = P.^2;
%     rho = Q.^2;

    rhor = rho./r.^(2*gamma);
    
    xplot = r;
    yplot = rhor;
%     yplot = yplot/yplot(2);
    
    h(oinx) = plot(xplot,yplot,'.-');
    ax = h(oinx).Parent;
    ax.XScale = 'log';
    ax.YScale = 'log';
    
    legendstr = [legendstr,plotOrbital(oinx)];
end

a0 = 52.917721067; % pm; 5.2917721067(12)×10−11 m (2014 CODATA)
r_nucleus_fm = 5.2192449328117521; % fm; rms radius set in grasp
r_nucleus = r_nucleus_fm/a0/1e3;

yrange = ax.YLim;
plot([1,1]*r_nucleus,yrange,'k--')

legend(h,legendstr{:})