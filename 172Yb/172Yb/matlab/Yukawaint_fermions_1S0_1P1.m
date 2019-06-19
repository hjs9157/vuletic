hbar = 1.054571800e-34; % 1.054571800(13)×10−34 J s  (2014 CODATA)
c = 299792458; % m/s
a0 = 5.2917721067e-11; % 5.2917721067(12)×10−11 (2014 CODATA)
eV = 1.6021766208e-19; % J; 1eV = 1.6021766208(98)×10−19 J (2014 CODATA)

mc2_a0 = hbar*c/a0; % J
eV2a0 = eV/mc2_a0; % 1 eV ~ 1/3729 mc2_a0

m_eV = logspace(0,9,1000);
m = m_eV*eV2a0;

stateID = {'1S0','1P1'};
confstr_core = ['1s(2)2s(2)2p-(2)2p(4)',...
    '3s(2)3p-(2)3p(4)3d-(4)3d(6)',...
    '4s(2)4p-(2)4p(4)4d-(4)4d(6)',...
    '5s(2)5p-(2)5p(4)4f-(6)4f(8)'];
confstr_valence1 = '6s(2)';
confstr_valence2 = '6s(1)6p(1)';
stateinfo{1} = {stateID{1},confstr_core,confstr_valence1}; % {stateID, confstr_core, confstr_valence}
stateinfo{2} = {stateID{2},confstr_core,confstr_valence2};
% DHForCI = 'DHF';
DHForCI = 'CI';

%% calculate
% fprintf('%s:\n',stateID{1})
% [Y1,Y_core1,Y_pop1,Y_corr1] = getYukawaInt_fermions(stateID{1},DHForCI,m);
% fprintf('%s:\n',stateID{2})
% [Y2,Y_core2,Y_pop2,Y_corr2] = getYukawaInt_fermions(stateID{2},DHForCI,m);
% 
% dY = Y1 - Y2;
% dY_core = Y_core1 - Y_core2;
% dY_pop = Y_pop1 - Y_pop2;
% dY_corr = Y_corr1 - Y_corr2;
[dY,dY_core,dY_pop,dY_corr] = getYukawaint_transition(stateinfo{:},DHForCI,m);
dY_valence = dY_pop + dY_corr;

%% plot
legendstr = {};
Y_fig = figure; hold on;
% h = plot(m,dY);
% h_core = plot(m,dY_core,'--');
% h_pop = plot(m,dY_pop,':');
% h_corr = plot(m,abs(dY_corr));
h(1) = plot(m_eV,abs(dY),'k-'); legendstr = [legendstr,{'Total'}];
h(2) = plot(m_eV,abs(dY_valence),'r--'); legendstr = [legendstr,{'valence'}];
h(3) = plot(m_eV,abs(dY_core),'b--'); legendstr = [legendstr,{'core'}];
h(4) = plot(m_eV,abs(dY_pop),':'); legendstr = [legendstr,{'pop'}];
h(5) = plot(m_eV,abs(dY_corr),':'); legendstr = [legendstr,{'corr'}];

ax = h.Parent;
ax.XScale = 'log';
ax.YScale = 'log';

% legendstr = {'Total','core','pop','corr'};
legend(legendstr{:})

xlabel('m_\phi (eV)')
ylabel('\langle e^{-m_\phi r}/r \rangle')
title(sprintf('Yukawa integral: %s - %s transition, %s',stateID{:},DHForCI))

% Julian's data for comparison
data_6s2_Julian_CI = [% [m (keV), Y]
0.001	-0.532212
0.01	-0.527415
0.1	-0.482056
1	-0.207475
10	-0.00682125
100	-0.00490513
1000	-0.00136634
10000	-5.33672E-05
100000	-8.67982E-07
];

data_6s6p_Julian_CI = [
0.001	-0.463106
0.01	-0.458316
0.1	-0.413572
1	-0.163934
10	-0.00470111
100	-0.00324696
1000	-0.000832215
10000	-3.23869E-05
100000	-5.26712E-07
];

data_6s2_Julian_CIMBPT = [
0.001		-0.551729
0.01		-0.546931
0.1		-0.501408
1		-0.218901
10		-0.00627391
100		-0.00640285
1000		-0.00180431
10000		-7.07004E-05
100000		-1.15002E-06
];

data_6s6p_Julian_CIMBPT = [
0.001		-0.480467
0.01		-0.475675
0.1		-0.430739
1		-0.173045
10		-0.00412872
100		-0.00427583
1000		-0.00107745
10000		-4.20315E-05
100000		-6.83678E-07
];


m_julian = data_6s2_Julian_CI(:,1)*1e3;
dY_julian_CI = abs(data_6s2_Julian_CI(:,2) - data_6s6p_Julian_CI(:,2));
h_julian_CI = plot(m_julian,dY_julian_CI,'o','Color',[1,0.5,0]);
dY_julian_CIMBPT = abs(data_6s2_Julian_CIMBPT(:,2) - data_6s6p_Julian_CIMBPT(:,2));
h_julian_CIMBPT = plot(m_julian,dY_julian_CIMBPT,'o','Color',[0,1,0.5]);

legend([h,h_julian_CI,h_julian_CIMBPT],...
    legendstr{:},'Julian\_CI','Julian\_CIMBTP',...
    'Location','southwest')

