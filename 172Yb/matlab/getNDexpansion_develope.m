numberDensity_fermions_2S12_2D52;
drawnow;
ND = dND;
% % ND = dND_valence;

% numberDensity_fermions_2S12;
% ND = ND_valence;

% stateID = {'2S12','2D52'};
% confstr_core = ['1s(2)2s(2)2p-(2)2p(4)',...
%     '3s(2)3p-(2)3p(4)3d-(4)3d(6)',...
%     '4s(2)4p-(2)4p(4)4d-(4)4d(6)',...
%     '5s(2)5p-(2)5p(4)4f-(6)4f(8)'];
% confstr_valence1 = '6s(1)';
% confstr_valence2 = '5d(1)';
% stateinfo{1} = {stateID{1},confstr_core,confstr_valence1}; % {stateID, confstr_core, confstr_valence}
% stateinfo{2} = {stateID{2},confstr_core,confstr_valence2};
% % DHForCI = 'DHF';
% DHForCI = 'CI';

% %% calculate
% [R,dND,dND_core,dND_pop,dND_corr] = getNumberDensity_transition(stateinfo{:},DHForCI);
% dND_valence = dND_pop + dND_corr;


%%
r_n_fm = 5.3; %fm
a_0 = 52.917721067; % pm; Bohr radius = 5.2917721067(12)×10−11 m (2014 CODATA)
r_n = r_n_fm/a_0;

order_fit = 3; % order of polinomial fit (2*n, n=1:order_fit)
InxCap_fit = [50:5:150,151:1:160,160:5:220];
% InxCap_fit = [50:5:220];
% InxCap_fit = 150;
InxBottom_fit = 10;

w = ([R(2:end);R(end)]-[R(1);R(1:end-1)])/2; % binsize of R
w_fit = w(InxBottom_fit:end);
R_fit = R(InxBottom_fit:end);
ND_fit = ND(InxBottom_fit:end);

Rcap_fit = R(InxCap_fit).';

% constuct fitting model
n_fit = (1:order_fit).';
ft = fittype(sprintf('poly%u',order_fit-1));
fitopts = fitoptions(sprintf('poly%u',order_fit-1));

% perform fitting
N_fit = length(InxCap_fit);
xi = nan(order_fit,N_fit); % ND = xi{fi}(1)*r^2 + xi{fi}(2)*r^4 + xi{fi}(3)*r^6 + ...
xi_err = nan(order_fit,N_fit);
fitresult = struct;
for fi = 1:N_fit
% for fi = 7
    icap = InxCap_fit(fi);
    R_fit_fi = R_fit(2:icap); % exclude r=0
    ND_fit_fi = ND_fit(2:icap); % exclude r=0
    
    Rsq_fit = R_fit_fi.^2;
    
    x_fit = Rsq_fit;
    y_fit = ND_fit_fi./x_fit;
    
    % scale x_fit
    sigma_x = std(x_fit); sigma_y = std(y_fit);
    x_fit_sc = x_fit/sigma_x; y_fit_sc = y_fit/sigma_y;
    
%     w = ([0;x_fit_sc] + [x_fit_sc;0])/2; % size of bins
    weight = w_fit(2:icap);
    
%     [p,S] = polyfit(x_fit_sc,y_fit_sc,order_fit-1);
    fitopts.Weights = weight;
    [fobj,gof,output] = fit(x_fit_sc,y_fit_sc,ft,fitopts);
    fitresult(fi).fobj = fobj;
    fitresult(fi).gof = gof;
    fitresult(fi).output = output;
    
    xi_sc = fliplr(coeffvalues(fobj)).';
    xi_err_sc = (diff(confint(fobj,erf(1/sqrt(2))))/2).';
    
    xi(:,fi) = xi_sc./sigma_x.^(n_fit-1)*sigma_y;
    xi_err(:,fi) = xi_err_sc./sigma_x.^(n_fit-1)*sigma_y;
    
end

fitfun = @(fi,x) xi(:,fi).'*x.^(n_fit-1);


fx_plot = @(x) x;
fy_plot = @(y) y./R_fit.^2;
xrange = fx_plot(Rcap_fit([1,end]).');
% xscale = 'linear';
xscale = 'log';
yscale = 'linear';
yscale = 'log';

Position_fig = [560 50 700 900];

polyexp_fig = figure;
polyexp_fig.Position = Position_fig;
h = gobjects(0);
clear ax
% wavefunction
subplot(order_fit+3+1,1,[1,2]); hold on;
% h(end+1) = plot(x_fit,y_fit,'k.-'); ax(1) = h(1).Parent;
h(end+1) = plot(fx_plot(R_fit),fy_plot(ND_fit),'k.-'); ax(1) = h(end).Parent;
title({sprintf('Radial number density: n(r) = \\Sigma_{k=1}^{%u}\\xi^{(2k)}r^{2k}',order_fit),[],...
    'Number density'})
ax(end).XScale = xscale;
% ax(end).YScale = yscale;
ax(end).XLim = xrange;
ylabel('n(r)/r^2 (a_0^{-3})')
xline(fx_plot(r_n),'k--');
ax(end).YLimMode = 'manual';
yrange = ax(end).YLim;
text(fx_plot(r_n),yrange(2),...
    'r_n',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

% (weighted) RMSE in nucleus
w_n = w(R <= r_n).';
R_n = R(R <= r_n).';
NDIR2_n = ND(R <= r_n).'./R_n.^2; 
R_n = R_n(2:end); NDIR2_n = NDIR2_n(2:end); w_n = w_n(2:end);
RMSE_n = nan(1,N_fit);
for fi = 1:N_fit
    NDIR2_n_est = fitfun(fi,R_n.^2);
    err_n = NDIR2_n - NDIR2_n_est;
%     RMSE_n(fi) = sum(err_n.^2)/length(R_n);
    RMSE_n(fi) = err_n.^2*w_n.'/sum(w_n);
end
subplot(order_fit+3+1,1,4); hold on;
h(end+1) = plot(fx_plot(Rcap_fit),RMSE_n,'r.-');
title('RMSE in nucleus')
ax(end+1) = h(end).Parent;
ax(end).XScale = xscale;
ax(end).XLim = xrange;
% ax(end).YScale = yscale;
ax(end).YScale = 'log';

% best fit
fi_best = find(RMSE_n == min(RMSE_n));
Rcap_best = Rcap_fit(fi_best);
RMSE_n_best = RMSE_n(fi_best);
xi_best = xi(:,fi_best);
xi_err_best = xi_err(:,fi_best);
subplot(order_fit+3+1,1,4)
xline(r_n,'k--');
xline(fx_plot(Rcap_best),'r-');
yrange = ax(end).YLim;
t = text(fx_plot(Rcap_best),geomean(yrange),...
    sprintf(' r = %g\nRMSE = %g',Rcap_best,RMSE_n_best),...
    'HorizontalAlignment','left',...
    'VerticalAlignment','middle');
t.Color = 'r';
subplot(order_fit+3+1,1,[1,2])
xline(fx_plot(Rcap_best),'r-');
h_fit = plot(fx_plot(R_fit.'),fitfun(fi_best,R_fit.'.^2),'r-');

% residuals
subplot(order_fit+3+1,1,3); hold on;
% h(end+1) = plot(x_fit,y_fit,'k.-'); ax(1) = h(1).Parent;
residual = ND_fit./R_fit.^2 - fitfun(fi_best,R_fit.'.^2).';
Inx_nucleus = R_fit <= 1.1*r_n;
h(end+1) = plot(fx_plot(R_fit(Inx_nucleus)),residual(Inx_nucleus),'k.-'); ax(end+1) = h(end).Parent;
title('Residuals')
ax(end).XScale = xscale;
% ax(end).YScale = yscale;
ax(end).XLim = xrange;
xline(fx_plot(r_n),'k--');
xline(fx_plot(Rcap_best),'r-');
yline(0,'r-');

n_subplot = length(polyexp_fig.Children)+1;
% fit params vs rcap
F = nan(order_fit,1);
for ni = 1:order_fit
    subplot(order_fit+n_subplot,1,ni+n_subplot); hold on;
%     h(ni+n_subplot) = plot(fx_plot(Rcap_fit),xi(ni,:),'.-');
    y_plot = xi(ni,:);
    if strcmp(yscale,'log')
        y_plot(sign(xi_best(ni))*y_plot < 0) = NaN;
    end
%     h(ni+n_subplot) = errorbar(fx_plot(Rcap_fit),abs(y_plot),xi_err(ni,:),'.-');
    h(ni+n_subplot) = plot(fx_plot(Rcap_fit),y_plot,'.-');
%     h(ni+n_subplot).CapSize = 0;
%         title(sprintf('\\xi^{(%u)} = %g \\pm %.3f a_0^{%i}',...
%             2*ni,xi(ni,fi_best),xi_err(ni,fi_best),-2*ni-1))
    
    ax(ni+n_subplot) = h(ni+n_subplot).Parent;
    ax(ni+n_subplot).XScale = xscale;
    ax(ni+n_subplot).XLim = xrange;
    ax(ni+n_subplot).YScale = yscale;
    
    xline(fx_plot(r_n),'k--');
    xline(fx_plot(Rcap_best),'r-');
    ylabel(sprintf('a_0^{%i}',-2*ni-1))
    
    alpha = 0.0072973525664; % 0.0072973525664(17) (2014 CODATA)
    h_plank = 6.62607015e-34; % Plank constant; 6.62607015×10−34 Js (exact; 2019 ISO)
    E_h = 4.359744650e-18; % Hartree energy; 4.359 744 650(54)e−18 (2014 CODATA)
    E_h_GHz = E_h/h_plank/1e9;
    a_0 = 52.917721067; % pm; Bohr radius = 5.2917721067(12)×10−11 m (2014 CODATA)
    Z = 70;
    
    F(ni) = Z*E_h_GHz/(2*ni)/(2*ni+1)*xi(ni,fi_best)/(a_0*1e3)^(2*ni);
    F_err(ni) = Z*E_h_GHz/(2*ni)/(2*ni+1)*xi_err(ni,fi_best)/(a_0*1e3)^(2*ni);
    fprintf('\t - F^(%u) = %g +- %.3e GHz/fm^(%u)\n',2*ni,F(ni),F_err(ni),2*ni)
%     text(xrange(1),geomean(ax(ni+n_subplot).YLim),...
%         sprintf(' F^(%u) = %g GHz/fm^(%u)\n',2*ni,F(ni),2*ni),...
%         'HorizontalAlignment','left',...
%         'VerticalAlignment','middle');

%     title({sprintf('\\xi^{(%u)} = %g \\pm %.3g a_0^{%i}',...
%         2*ni,xi(ni,fi_best),xi_err(ni,fi_best),-2*ni),...
%         sprintf('F^{(%u)} = %g \\pm %.3g GHz/fm^{(%u)}\n',...
%         2*ni,F(ni),F_err(ni),2*ni)})

      title([sprintf('\\xi^{(%u)} = %g \\pm %.3g a_0^{%i}',...
        2*ni,xi(ni,fi_best),xi_err(ni,fi_best),-2*ni-1),'   ',...
        sprintf('F^{(%u)} = %g \\pm %.3g GHz/fm^%u\n',...
        2*ni,F(ni),F_err(ni),2*ni)])

end

xlabel('r (a_0)')
