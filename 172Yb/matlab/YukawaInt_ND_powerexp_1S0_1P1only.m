% numberDensity_fermions_2D52;

numberDensity_fermions_1S0_1P1only;
ND = dND;
ND_core = dND_core;
ND_valence = dND_valence;

drawnow;

%%

% % nuclear size
% r_n_fm = 5.3; % fm; nuclear size 
% a_0 = 52.917721067; % pm; Bohr radius = 5.2917721067(12)×10−11 m (2014 CODATA)
% r_n = r_n_fm/a_0/1e3;

% % Fermi distribution parameters
% c_fm = 6.227;
% a_fm = 2.18/4/log(3); % for 172Yb
% c = c_fm/a_0/1e3; a = a_fm/a_0/1e3; 

% order_fit = 3;
% InxCap_fit = find(R > r_n*1e-2 & R < r_n*1e1);
% InxBottom_fit = 10;
% 
% [S,NDexp_fig] = getNDexpansion(R,ND,order_fit,r_n,InxCap_fit,InxBottom_fit);

hbar = 1.054571800e-34; % 1.054571800(13)×10−34 J s  (2014 CODATA)
c = 299792458; % m/s
a0 = 5.2917721067e-11; % 5.2917721067(12)×10−11 (2014 CODATA)
eV = 1.6021766208e-19; % J; 1eV = 1.6021766208(98)×10−19 J (2014 CODATA)
mc2_a0 = hbar*c/a0; % J
eV2a0 = eV/mc2_a0; % 1 eV ~ 1/3729 mc2_a0
m_eV = logspace(0,15,10000);
m = m_eV*eV2a0;

Y = getYukawaInt_ND_powerexp(r,ND,m,[2,4,6,8],1e-6,1000,20,10,false,false);
Y_core = getYukawaInt_ND_powerexp(r,ND_core,m,[2,4,6,8],1e-6,1000,20,10,false,false);
% Y_pop = getYukawaInt_ND_powerexp(r,ND_pop,m,[6,8],1e-6,1000,20,10,false,false);
% Y_corr = getYukawaInt_ND_powerexp(r,ND_corr,m,[6,8],1e-6,1000,20,10,false,false);
Y_valence = getYukawaInt_ND_powerexp(r,ND_valence,m,[2,4,6,8],1e-6,1000,20,10,false,false);

% plot
h = gobjects(0);
legendstr = {};
Y_fig = figure; hold on;
ax = Y_fig.CurrentAxes;
ax.XScale = 'log';
ax.YScale = 'log';
h(end+1) = plot(m_eV,abs(Y),'k-'); legendstr{end+1} = 'Total';
h(end+1) = plot(m_eV,abs(Y_valence),'r--'); legendstr{end+1} = 'valence';
h(end+1) = plot(m_eV,abs(Y_core),'b--'); legendstr{end+1} = 'core';
% yrange = ax.YLim;
% h(end+1) = plot(m_eV,abs(Y_pop),':'); legendstr{end+1} = 'pop';
% h(end+1) = plot(m_eV,abs(Y_corr),':'); legendstr{end+1} = 'corr';
% ax.YLim = yrange;

xline(1/eV2a0,'k:');
r_n_fm = 5.3; % rms nuclear size
xline(r_n_fm/a0/1e3/eV2a0,'k:');


xlim(m_eV([1,end]))

% legendstr = {'Total','core','pop','corr'};
legend(legendstr{:})

xlabel('m_\phi (eV)')
ylabel('\langle e^{-m_\phi r}/r \rangle')
% title(sprintf('Yukawa integral: %s - %s transition, %s\n through number density function',stateID{:},DHForCI))
