numberDensity_fermions_2D52;
% ND = ND_valence;
% ND = ND_core;

%% inputs
% function Y = getYukawaInt_ND_powerexp(r,ND,m,Nmin_int,k,varargin)
hbar = 1.054571800e-34; % 1.054571800(13)×10−34 J s  (2014 CODATA)
c = 299792458; % m/s
a0 = 5.2917721067e-11; % 5.2917721067(12)×10−11 (2014 CODATA)
eV = 1.6021766208e-19; % J; 1eV = 1.6021766208(98)×10−19 J (2014 CODATA)

mc2_a0 = hbar*c/a0; % J
eV2a0 = eV/mc2_a0; % 1 eV ~ 1/3729 mc2_a0

m_eV = logspace(0,20,1000);
m = m_eV*eV2a0;

r_cut = 1e-6;
N_int_cut = 120; % number of point for integration under r_cut
k = [2,4,6,8];
showPowerexp = true;
showYukawaint = true;
int_m_ratio = 10; % (integration range to be accrate)/(1/m)
N_bottom = 10;

fprintf('r_cut = %e\n',r_cut)
fprintf('Nmin_int = %u\n',Nmin_int)

% processing inputs
Nmin_int = sum(r <= r_cut);
r = r(N_bottom:end);
ND = ND(N_bottom:end);

%% power expansion
w = ([r(2:end);r(end)]-[r(1);r(1:end-1)])/2; % binsize of r
Inx_fit = 1:Nmin_int;
r_fit = r(Inx_fit);
ND_fit = ND(Inx_fit);
w_fit = w(Inx_fit);

x_fit = r_fit;
y_fit = ND_fit./r_fit.^2;
k_fit = k - 2;

% [p,se] = wpolyfit(x_fit,y_fit,w_fit,k_fit);

xscale_fit = std(x_fit);
x_fit_sc = x_fit/xscale_fit;
yscale_fit = 1;
y_fit_sc = y_fit/yscale_fit;

[p_sc,se_sc] = wpolyfit(x_fit_sc,y_fit_sc,w_fit,k_fit);
p = yscale_fit*p_sc.*(xscale_fit.^(-k_fit(:)).');
se = yscale_fit*se_sc.*(xscale_fit.^(-k_fit(:)).');

fitfun = @(x) x(:).^(k_fit(:).')*p.';
fND_fit = @(r) r(:).^(k(:).')*p.';

if showPowerexp
    ND_est = fND_fit(r_fit(:));
    
    figure;
    ax = subplot(2,1,1); hold on
    h = plot(r_fit,ND_fit./r_fit.^2,'.-');
    h_fit = plot(r_fit,ND_est./r_fit.^2,'r-');
    ax.XScale = 'log';
    xlim(r_fit([1,end]))
    
    ax_res = subplot(2,1,2); hold on
    yline(0,'r-');
    residual = (ND_fit - ND_est)./r_fit.^2;
    h_res = plot(r_fit,residual);
    ax_res.XScale = 'log';
    xlim(r_fit([1,end]))
end

%% Yukawa integral
f_Yukawa = @(r,m) exp(-m.*r)./r;

Y = nan(size(m));
isPowerexpUsed = false(size(m));
for mi = 1:length(m)
    % if int_m_ratio/m > r_cut then do simple integration
    r_cap = int_m_ratio/m(mi);
    if r_cap >= r_cut
%         isPowerexpUsed(mi) = false;
        
        r_int  = r;
        ND_int = ND;
    else % r_cap < r_cut
        isPowerexpUsed(mi) = true;
        
        r1 = linspace(0,r_cap,N_int_cut).';
        r1 = r1(2:end);
        ND1 = fND_fit(r1);
        
%         Inx2 = r >= r_cut;
%         r2 = r(Inx2);
%         ND2 = ND(Inx2);

        r2 = [];
        ND2 = [];
        
        r_int = [r1;r2];
        ND_int = [ND1;ND2];
    end
    
    V_Yukawa_int = f_Yukawa(r_int,m(mi));
    Y(mi) = trapz(r_int,V_Yukawa_int.*ND_int,1);
end

% plot Y vs m
Y_fig = figure; hold on;
ax = Y_fig.CurrentAxes;
h_Y = plot(m_eV,Y,'.-');
ax.XScale = 'log';
ax.YScale = 'log';
yrange = ax.YLim;
m_eV_powerexpUsed = m_eV(isPowerexpUsed);
temp = 0.1;
y_powerexpUsed = yrange(1)*(1-temp) + yrange(2)*temp;
plot(m_eV_powerexpUsed,y_powerexpUsed*ones(size(m_eV_powerexpUsed)),'g.-')
text(max(m_eV_powerexpUsed),y_powerexpUsed,...
    'power expansion used',...
    'HorizontalAlignment','right',...
    'VerticalAlignment','bottom',...
    'Color','g')


% end