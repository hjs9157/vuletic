numberDensity_fermions_2S12_2D52;
drawnow;
ND = dND;

%%
% nuclear size
r_n_fm = 5.3; %fm
a_0 = 52.917721067; % pm; Bohr radius = 5.2917721067(12)×10−11 m (2014 CODATA)
r_n = r_n_fm/a_0/1e3;

% Fermi distribution parameters
c_fm = 6.227;
a_fm = 2.18/4/log(3); % for 172Yb
c = c_fm/a_0/1e3; a = a_fm/a_0/1e3; 

k = [2,3,4,5]; % orders of polinomial fit (2*k, n=1:order_fit).
InxCap_fit = find(r > r_n*1e-3 & r < r_n*1e1);
% InxCap_fit = [50:5:150,151:1:160,160:5:220];
% InxCap_fit = [50:5:220];
% InxCap_fit = 150;
InxBottom_fit = 10;

% take binsize as weight of fit
% w = ([R(2:end);R(end)]-[R(1);R(1:end-1)])/2; % binsize of R
w = ([r(2:end);r(end)]-[r(1);r(1:end-1)])/2; % binsize of r
rCap_fit = r(InxCap_fit).';
rBottom_fit = r(InxBottom_fit);

% perform fittingr
N_fit = length(InxCap_fit);
xi = nan(length(k),N_fit);
xi_err = nan(length(k),N_fit);
for fi = 1:N_fit
% for fi = 92
%     icap = 155;
    icap = InxCap_fit(fi);
    
    x_fit = r(InxBottom_fit:icap);
%     x_fit = x_fit;
    y_fit = ND(InxBottom_fit:icap);
    y_fit = y_fit./x_fit.^2;
    k_fit = (k - 2);
    
    % scale x_fit
    sigma_x = std(x_fit); sigma_y = std(y_fit);
    x_fit_sc = x_fit/sigma_x; y_fit_sc = y_fit/sigma_y;
    
    weight = w(InxBottom_fit:icap);
%     weight = ones(size(weight));
    
    % fit
    [p,se] = wpolyfit(x_fit_sc,y_fit_sc,weight,k_fit);
    xi_sc = p.';
    xi_err_sc = se.';
    
    xi(:,fi) = xi_sc./sigma_x.^k_fit.'*sigma_y;
    xi_err(:,fi) = xi_err_sc./sigma_x.^k_fit.'*sigma_y;    
end

% fitfun = @(fi,x) xi(:,fi).'*x.^k_fit;
fitfun = @(fi,x) x.^k_fit*xi(:,fi);

fx_plot = @(x) x;
% xrange = fx_plot(Rcap_fit([1,end]).');
xrange = fx_plot([rBottom_fit,max(rCap_fit)]);
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
subplot(length(k)+3+1,1,[1,2]); hold on;
% h(end+1) = plot(x_fit,y_fit,'k.-'); ax(1) = h(1).Parent;
h(end+1) = plot(fx_plot(r),ND./r.^2,'k.-'); ax(1) = h(end).Parent;
title({sprintf('Radial number density: n(r) = \\Sigma_{k=1}^{%u}\\xi^{(2k)}r^{2k}',length(k)),[],...
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

% nuclear charge distribution
if ~isnan(c) && ~isnan(a)
    yyaxis right
    if a == 0
        f_fermi = @(r) heaviside(c-r);
    else
        f_fermi = @(r) 1./(1 + exp((r-c)/a));
    end
    rho_0 = 1/integral(@(r) 4*pi*r.^2.*f_fermi(r),0,Inf);
    rho_fermi = @(r) rho_0.*f_fermi(r);
    h_nucleus = fplot(rho_fermi,xrange,'b-');
    ax_nucleus = h_nucleus.Parent;
    ax_nucleus.YColor = 'b';
    ylabel('\rho_n(r) (a_0^{-3} Z^{-1})')
    yyaxis left
end

% (weighted) RMSE in nucleus
Inx_nucleus = find(r >= rBottom_fit & r <= max(rCap_fit) & r <= r_n);
w_n = w(Inx_nucleus);
r_nucleus = r(Inx_nucleus);
NDIR2_n = ND(Inx_nucleus)./r_nucleus.^2; 
RMSE_n = nan(1,N_fit);
for fi = 1:N_fit
    NDIR2_n_est = fitfun(fi,r_nucleus);
    err_n = NDIR2_n - NDIR2_n_est;
%     RMSE_n(fi) = sum(err_n.^2)/length(R_n);
    RMSE_n(fi) = err_n.'.^2*w_n/sum(w_n);
end
subplot(length(k)+3+1,1,4); hold on;
h(end+1) = plot(fx_plot(rCap_fit),RMSE_n,'r.-');
title('RMSE in nucleus')
ax(end+1) = h(end).Parent;
ax(end).XScale = xscale;
ax(end).XLim = xrange;
% ax(end).YScale = yscale;
ax(end).YScale = 'log';

% best fit
fi_best = find(RMSE_n == min(RMSE_n));
rCap_best = rCap_fit(fi_best);
RMSE_n_best = RMSE_n(fi_best);
xi_best = xi(:,fi_best);
xi_err_best = xi_err(:,fi_best);
subplot(length(k)+3+1,1,4)
xline(r_n,'k--');
xline(fx_plot(rCap_best),'r-');
yrange = ax(end).YLim;
t = text(fx_plot(rCap_best),geomean(yrange),...
    sprintf(' r = %g\nRMSE = %g',rCap_best,RMSE_n_best),...
    'HorizontalAlignment','left',...
    'VerticalAlignment','middle');
t.Color = 'r';
subplot(length(k)+3+1,1,[1,2])
xline(fx_plot(rCap_best),'r-');
h_fit = plot(fx_plot(r.'),fitfun(fi_best,r),'r-');

% residuals
subplot(length(k)+3+1,1,3); hold on;
% h(end+1) = plot(x_fit,y_fit,'k.-'); ax(1) = h(1).Parent;
residual = ND./r.^2 - fitfun(fi_best,r);
% Inx_nucleus = R_fit <= 1.1*r_n;
h(end+1) = plot(fx_plot(r(Inx_nucleus)),residual(Inx_nucleus),'k.-'); ax(end+1) = h(end).Parent;
title('Residuals')
ax(end).XScale = xscale;
% ax(end).YScale = yscale;
ax(end).XLim = xrange;
xline(fx_plot(r_n),'k--');
xline(fx_plot(rCap_best),'r-');
yline(0,'r-');

n_subplot = length(polyexp_fig.Children)+1;
% fit params vs rcap
F = nan(length(k),1);
for ni = 1:length(k)
    subplot(length(k)+n_subplot,1,ni+n_subplot); hold on;
%     h(ni+n_subplot) = plot(fx_plot(Rcap_fit),xi(ni,:),'.-');
    y_plot = xi(ni,:);
    if strcmp(yscale,'log')
        y_plot(sign(xi_best(ni))*y_plot < 0) = NaN;
    end
%     h(ni+n_subplot) = errorbar(fx_plot(Rcap_fit),abs(y_plot),xi_err(ni,:),'.-');
    h(ni+n_subplot) = plot(fx_plot(rCap_fit),y_plot,'.-');
%     h(ni+n_subplot).CapSize = 0;
%         title(sprintf('\\xi^{(%u)} = %g \\pm %.3f a_0^{%i}',...
%             2*ni,xi(ni,fi_best),xi_err(ni,fi_best),-2*ni-1))
    
    ax(ni+n_subplot) = h(ni+n_subplot).Parent;
    ax(ni+n_subplot).XScale = xscale;
    ax(ni+n_subplot).XLim = xrange;
    ax(ni+n_subplot).YScale = yscale;
    
    xline(fx_plot(r_n),'k--');
    xline(fx_plot(rCap_best),'r-');
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
