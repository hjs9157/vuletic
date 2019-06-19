stateID = {'2S12','2D52'}; % states in transition

F_state = nan(1,2);
for sti = 1:2
    evalc(['T_orbital = getOrbitals(''',stateID{sti},''')']);
    
    corestr = '1s(2)2s(2)2p(6)3s(2)3p(6)3d(10)4s(2)4p(6)4d(10)5s(2)5p(6)4f(14)';
    valencestr = '6s(1)';
    
    exp_orbital = '[spd-z]';
    
    % expression = ['(\d+)(',exp_orbital,')\((\d+)\)'];
    expression = ['(?<n>\d+)(?<l>',exp_orbital,')\((?<occupn>\d+)\)'];
    % core = regexp(corestr,expression,'tokens');
    occpn_core = regexp(corestr,expression,'names');
    % valence = regexp(valencestr,expression,'tokens');
    occpn_valence = regexp(valencestr,expression,'names');
    
    subshells = T_orbital.Properties.RowNames;
    
    % get F for core
    F_core = 0;
    occpn = occpn_core;
    for occi = 1:length(occpn)
        if strcmp(occpn(occi).l,'s')
            subshell = sprintf('%s%s',occpn(occi).n,occpn(occi).l);
            shInx = find(cellfun(@(c) ~isempty(c),regexp(subshells,[subshell,'-?']))); % index of row in T_orbital
            shInx = shInx.';
            
            for shi = shInx
                %         rho2 = T_orbital.P2(shi)^2 + T_orbital.Q2(shi)^2;
                r2 = T_orbital.rwfn{shi}.r(2);
                P2 = T_orbital.rwfn{shi}.P(2);
                Q2 = T_orbital.rwfn{shi}.Q(2);
                rho2 = P2^2 + Q2^2;
                
                F_subshell = rho2/r2^2;
                F_core = F_core + F_subshell*T_orbital.occupation(shi);
                
            end
        end
    end
    
    % get F for valence
    F_valence = 0;
    occpn = occpn_valence;
    for occi = 1:length(occpn)
        if strcmp(occpn(occi).l,'s')
            subshell = sprintf('%s%s',occpn(occi).n,occpn(occi).l);
            shInx = find(cellfun(@(c) ~isempty(c),regexp(subshells,[subshell,'-?']))); % index of row in T_orbital
            shInx = shInx.';
            
            for shi = shInx
                %         rho2 = T_orbital.P2(shi)^2 + T_orbital.Q2(shi)^2;
                r2 = T_orbital.rwfn{shi}.r(2);
                P2 = T_orbital.rwfn{shi}.P(2);
                Q2 = T_orbital.rwfn{shi}.Q(2);
                rho2 = P2^2 + Q2^2;
                
                F_subshell = rho2/r2^2;
                F_valence = F_valence + F_subshell*T_orbital.occupation(shi);
                
            end
        end
    end
    
    F_state(sti) = F_core + F_valence;
end

F_transition = F_state(1) - F_state(2)