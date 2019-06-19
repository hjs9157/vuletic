% evalc('T_orbital = getOrbitals(''2S12'')');

m_eV = logspace(0,10,100); % eV
a0toeV = 3728.9398377647327968373543999879; % eV for reduced compton wavelength of a_0
m = m_eV/a0toeV;

% orbital = {{'2S12','6s'},{'2D32','5d-'},{'2D52','5d'},{'2F72','6s'},{'2F72','4f-'},{'2F72','4f'}};
orbital = {{'1S0','6s'},{'1S0','5d-'},{'1S0','5d'},...
    {'1S0','5p-'},{'1S0','5p'},{'1S0','7s'}};
DHForCI = 'DHF';

YukawaInt = {};
rwfn_fig = figure; hold on
int_fig = figure; hold on
clear h_int h_rwfn
for oi=1:length(orbital)
    stateID = orbital{oi}{1};
    subshell = orbital{oi}{2};
    
    T_orbital = getOrbitals(stateID,DHForCI);
%     evalc(['T_orbital = getOrbitals(''',stateID,''')']);
    
    rwfn = T_orbital.rwfn{subshell};
    gamma = T_orbital.gamma(subshell);
    
    figure(rwfn_fig)
    h_rwfn(oi) = plot(rwfn.r,rwfn.P.^2+rwfn.Q.^2,'.-');
    
    YukawaInt{oi} = getYukawaInt(rwfn,gamma,m);
    figure(int_fig)
    h_int(oi) = plot(m_eV,YukawaInt{oi},'-');
    
    ax = h_int(oi).Parent;
    ax.XScale = 'log';
    ax.YScale = 'log';
    
end
figure(int_fig); ylim([1e-8,1e1])
figure(rwfn_fig); ylim([1e-8,1e1])

legendstr = cellfun(@(a) [a{1},'\_',a{2}],orbital,'UniformOutput',false);
figure(rwfn_fig); legend(legendstr{:})
figure(int_fig); legend(legendstr{:})