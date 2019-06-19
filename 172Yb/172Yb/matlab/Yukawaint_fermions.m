hbar = 1.054571800e-34; % 1.054571800(13)×10−34 J s  (2014 CODATA)
c = 299792458; % m/s
a0 = 5.2917721067e-11; % 5.2917721067(12)×10−11 (2014 CODATA)
eV = 1.6021766208e-19; % J; 1eV = 1.6021766208(98)×10−19 J (2014 CODATA)

mc2_a0 = hbar*c/a0; % J
eV2a0 = eV/mc2_a0; % 1 eV ~ 1/3729 mc2_a0

m_eV = logspace(-10,10,1000);
m = m_eV*eV2a0;

stateID = '1S0';
% stateID = '1P1';
DHForCI = 'DHF';
% DHForCI = 'CI';

% electron configuration of state
confstr_core = ['1s(2)2s(2)2p-(2)2p(4)',...
    '3s(2)3p-(2)3p(4)3d-(4)3d(6)',...
    '4s(2)4p-(2)4p(4)4d-(4)4d(6)',...
    '5s(2)5p-(2)5p(6)4f-(6)4f(8)'];
confstr_valence = '6s(2)';
% confstr_valence = '6s(1)6p(1)';

[Y,Y_core,Y_pop,Y_corr] = ...
    getYukawaInt_fermions(stateID,confstr_core,confstr_valence,DHForCI,m);

Y_fig = figure; hold on;
h = plot(m,Y,'k-');
h_core = plot(m,Y_core,'r--');
h_pop = plot(m,Y_pop,':');
h_corr = plot(m,abs(Y_corr),':');

ax = h.Parent;
ax.XScale = 'log';
ax.YScale = 'log';

legendstr = {'Total','core','pop','corr'};
legend(legendstr{:})