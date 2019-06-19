a_0 = 52.917721067; % pm; Bohr radius = 5.2917721067(12)×10−11 m (2014 CODATA)
r_n_fm = 5.298; % fm; rms nuclear radius of 172Yb

r_n = r_n_fm/a_0/1e3;

stateID = {'1S0','3P0only'};
confstr_core1 = ['1s(2)2s(2)2p-(2)2p(4)',...
    '3s(2)3p-(2)3p(4)3d-(4)3d(6)',...
    '4s(2)4p-(2)4p(4)4d-(4)4d(6)',...
    '5s(2)5p-(2)5p(4)4f-(6)4f(8)'];
confstr_valence1 = '6s(2)';
confstr_core2 = ['1s(2)2s(2)2p-(2)2p(4)',...
    '3s(2)3p-(2)3p(4)3d-(4)3d(6)',...
    '4s(2)4p-(2)4p(4)4d-(4)4d(6)',...
    '5s(2)5p-(2)5p(4)'];
confstr_valence2 = '4f-(6)4f(8)6s(1)6p(1)';
stateinfo{1} = {stateID{1},confstr_core1,confstr_valence1}; % {stateID, confstr_core, confstr_valence}
stateinfo{2} = {stateID{2},confstr_core2,confstr_valence2};
% DHForCI = 'DHF';
DHForCI = 'CI';

%% calculate
[r,dND,dND_core,dND_pop,dND_corr] = getNumberDensity_transition(stateinfo{:},DHForCI);
dND_valence = dND_pop + dND_corr;

%% plot
fx_plot = @(x) x;
% fy_plot = @(y) y;
fy_plot = @(y) y./r.^2;

h = gobjects(0);
legendstr = {};
ND_fig = figure; hold on;
h(end+1) = plot(fx_plot(r),fy_plot(dND),'k-'); legendstr{end+1} = 'Total';
h(end+1) = plot(fx_plot(r),fy_plot(dND_valence),'r--'); legendstr{end+1} = 'valence';
h(end+1) = plot(fx_plot(r),fy_plot(dND_core),'b--'); legendstr{end+1} = 'core';
h(end+1) = plot(fx_plot(r),fy_plot(dND_pop),':'); legendstr{end+1} = 'pop';
h(end+1) = plot(fx_plot(r),fy_plot(dND_corr),':'); legendstr{end+1} = 'corr';

ax = h.Parent;
ax.XScale = 'log';
% ax.YScale = 'log';

% legend(legendstr{:},'Location','southwest')

xlabel('r (a_0)')
% ylabel('\Deltan(r)')
ylabel('\Deltan(r)/r^2')
title(sprintf('Radial number density difference:\n %s - %s transition, %s',stateID{:},DHForCI))

xrange = [1e-5,1e2];
xlim(xrange)

% yrange = ax.YLim;
% plot([1,1]*r_n/a_0/1e3,yrange,'k--')
lx = xline(r_n,'k--');
yrange = ax.YLim;
text(r_n,yrange(2),...
    'r_n',...
    'HorizontalAlignment','center',...
    'VerticalAlignment','bottom')

legend(h,legendstr{:},'Location','best')