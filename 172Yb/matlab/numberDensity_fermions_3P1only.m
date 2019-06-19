
% stateID = '1S0';
stateID = '3P1only';
% DHForCI = 'DHF';
DHForCI = 'CI';

% electron configuration of state
confstr_core = ['1s(2)',...
    '2s(2)2p-(2)2p(4)',...
    '3s(2)3p-(2)3p(4)3d-(4)3d(6)',...
    '4s(2)4p-(2)4p(4)4d-(4)4d(6)',...
    '5s(2)5p-(2)5p(4)4f-(6)4f(8)'];
% confstr_valence = '6s(2)';
confstr_valence = '6s(1)6p(1)';

[r,ND,ND_core,ND_pop,ND_corr] = ...
    getNumberDensity_fermions(stateID,confstr_core,confstr_valence,DHForCI);
ND_valence = ND_pop + ND_corr;

%% Plot
legendstr = {};
h = gobjects(0);
ND_fig = figure; hold on;
h(end+1) = plot(r,ND,'k-'); legendstr{end+1} = 'Total';
h(end+1) = plot(r,ND_valence,'r--'); legendstr{end+1} = 'valence';
h(end+1) = plot(r,ND_core,'b--'); legendstr{end+1} = 'core';
h(end+1) = plot(r,ND_pop,':'); legendstr{end+1} = 'pop';
h(end+1) = plot(r,abs(ND_corr),':'); legendstr{end+1} = 'corr';

ax = h.Parent;
ax.XScale = 'log';
xlim([1e-5,1e2])
ax.YScale = 'log';

legend(legendstr{:})